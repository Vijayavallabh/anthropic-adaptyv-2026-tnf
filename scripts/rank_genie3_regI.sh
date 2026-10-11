#!/bin/bash
# Real ranking for Genie3-regI binders (the Anthropic/Overath recipe): Boltz-2 co-fold vs the TWO-PROTOMER
# groove (binder A + human/mouse TNF as chains B AND C), then ipSAE_min (interface-focused, 1.4x better than
# af2ig for TNF) + groove-pose-gate (binder must contact the hotspots) + cross-react = min(human,mouse).
# First run on the 12 CRO-validated binders (runs/genie3_regI_calib.csv) CALIBRATES the ranker (POS vs NEG).
# Usage: rank_genie3_regI.sh <survivors.csv name,sequence> <tag> "<gpus 'a b'>"
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
CSV=$1; TAG=$2; IFS=' ' read -ra GA <<< "${3:-4 5 6 7}"; NG=${#GA[@]}
ESMC="$REPO/envs/uplifting-esmc/bin/python"
IPSAE_PY="$REPO/third_party/LigandMPNN/.venv/bin/python"
HUM=$(tr -d '\n ' < "$REPO/runs/partial_denovo2/funnel/xreact/human_tnf.txt")
MOU=$(tr -d '\n ' < "$REPO/runs/buried2/mouse_tnf_mature.txt")
OUT="$REPO/runs/${TAG}"; mkdir -p "$OUT"
# 1) write YAMLs: A=binder, B=C=TNF (homotrimer groove), msa empty; one dir per species
"$ESMC" - "$CSV" "$OUT" "$HUM" "$MOU" <<'PY'
import sys,csv,os
csvf,out,hum,mou=sys.argv[1:5]
for sp,tgt in (("h",hum),("m",mou)):
    d=os.path.join(out,f"in_{sp}"); os.makedirs(d,exist_ok=True)
    for r in csv.DictReader(open(csvf)):
        nm=r.get("name") or r.get("design"); seq=r["sequence"].split(":")[0]
        open(os.path.join(d,nm+".yaml"),"w").write(
          "version: 1\nsequences:\n"
          f"  - protein: {{id: A, sequence: {seq}, msa: empty}}\n"
          f"  - protein: {{id: B, sequence: {tgt}, msa: empty}}\n"
          f"  - protein: {{id: C, sequence: {tgt}, msa: empty}}\n")
    print("wrote",sum(1 for _ in csv.DictReader(open(csvf))),"yaml ->",d)
PY
# 2) Boltz co-fold, shard species x gpu
run_boltz(){ local IN=$1 O=$2 G=$3
  ( unset VIRTUAL_ENV 2>/dev/null||true; source "$REPO/envs/uplifting-boltz2/bin/activate"
    set -a; source "$REPO/envs/uplifting.env" 2>/dev/null||true; set +a
    cd "$REPO/third_party/uplifting-biomolecular-modeling/boltz2"
    CUDA_VISIBLE_DEVICES=$G bash run.sh pred --config h100 --mode fast --input "$IN" --out_dir "$O" ) ; }
g0=${GA[0]}; g1=${GA[$((NG>1?1:0))]}
run_boltz "$OUT/in_h" "$OUT/boltz_h" "$g0" > "$OUT/boltz_h.log" 2>&1 &
ph=$!
run_boltz "$OUT/in_m" "$OUT/boltz_m" "$g1" > "$OUT/boltz_m.log" 2>&1 &
pm=$!
wait $ph; wait $pm
echo "[rank] boltz done $(date -u +%FT%TZ)"
# 3) ipSAE (A-B, A-C -> min) + groove pose-gate per species; cross-react = min(h,m)
"$ESMC" - "$CSV" "$OUT" "$IPSAE_PY" <<'PY'
import sys,csv,os,glob,subprocess,numpy as np
csvf,out,ipy=sys.argv[1:4]
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
HOT={("B",69),("B",70),("B",92),("C",109),("C",110)}
def find(sp,nm):
    for cif in glob.glob(f"{out}/boltz_{sp}/**/{nm}_model_0.cif",recursive=True):
        pae=os.path.join(os.path.dirname(cif),f"pae_{nm}_model_0.npz")
        if os.path.exists(pae): return cif,pae
    return None,None
def ipsae_pair(cif,pae):
    txt=cif.replace(".cif","_10_10.txt")
    if not os.path.exists(txt):
        subprocess.run([ipy,f"{R}/tools/ipsae.py",pae,cif,"10","10"],capture_output=True,cwd=os.path.dirname(cif))
    vals={}
    if os.path.exists(txt):
        for ln in open(txt):
            f=ln.split()
            # cols: Chn1 Chn2 PAE Dist Type ipSAE ...  -> Type=f[4], ipSAE=f[5]
            if len(f)>6 and f[4]=="asym":
                try: vals[(f[0],f[1])]=float(f[5])
                except: pass
    # binder (A) interfaces to targets B,C -> min; fall back to any A-pair
    cand=[v for (a,b),v in vals.items() if "A" in (a,b) and v is not None]
    return (min(cand) if cand else None), vals
def pose_contacts(cif):
    # CA of binder (A) vs hotspot residues
    A=[];HOTxyz=[]
    for l in open(cif):
        pass
    # cif parse: use gemmi-free simple columns if PDB; cif is mmCIF -> parse _atom_site
    import re
    ca={}
    with open(cif) as fh:
        cols=[];inloop=False
        for ln in fh:
            if ln.startswith("_atom_site."): cols.append(ln.strip()); inloop=True; continue
            if inloop and (ln.startswith("ATOM") or ln.startswith("HETATM")):
                f=ln.split()
                rec=dict(zip([c.split('.')[1] for c in cols],f))
                if rec.get("label_atom_id")=="CA":
                    ch=rec.get("auth_asym_id") or rec.get("label_asym_id")
                    rid=int(rec.get("auth_seq_id") or rec.get("label_seq_id"))
                    xyz=np.array([float(rec["Cartn_x"]),float(rec["Cartn_y"]),float(rec["Cartn_z"])])
                    ca[(ch,rid)]=xyz
            elif inloop and ln.startswith("#"): inloop=False
    A=np.array([v for (c,r),v in ca.items() if c=="A"])
    hot=np.array([ca[k] for k in HOT if k in ca])
    if len(A)==0 or len(hot)==0: return 0,99.0
    d=np.sqrt(((A[:,None,:]-hot[None,:,:])**2).sum(-1))
    return int((d<8).any(0).sum()), float(d.min())
rows=[]
for r in csv.DictReader(open(csvf)):
    nm=r.get("name") or r.get("design")
    res={}
    for sp in ("h","m"):
        cif,pae=find(sp,nm)
        if cif: res[sp]=(ipsae_pair(cif,pae)[0], *pose_contacts(cif))
        else: res[sp]=(None,0,99.0)
    ih=res["h"][0]; im=res["m"][0]
    xr=min(ih,im) if (ih is not None and im is not None) else None
    rows.append((nm,ih,im,xr,res["h"][1],res["h"][2],r.get("kind","")))
rows.sort(key=lambda x:-(x[3] or -1))
print(f"{'name':40}{'ipSAE_h':>8}{'ipSAE_m':>8}{'xr_min':>8}{'hotCov':>7}{'minD':>6}  kind")
for nm,ih,im,xr,hc,md,k in rows:
    f=lambda v:(f'{v:.3f}' if isinstance(v,float) else 'NA')
    print(f"{nm[:40]:40}{f(ih):>8}{f(im):>8}{f(xr):>8}{hc:>7}{md:>6.1f}  {k}")
import csv as _c
_c.writer(open(f"{out}/ranked.csv","w")).writerows(
    [["name","ipsae_h","ipsae_m","xr_min","hot_cov","min_d","kind"]]+
    [[nm,ih,im,xr,hc,round(md,1),k] for nm,ih,im,xr,hc,md,k in rows])
print(f"\n-> {out}/ranked.csv")
PY
echo "RANK_GENIE3_REGI_DONE $(date -u +%FT%TZ)"
