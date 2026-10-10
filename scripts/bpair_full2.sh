#!/bin/bash
# rfd2 downstream funnel (2nd de-novo generation wave for the best-20-on-all-metrics campaign):
# waits for the 8 rfd2 RFd streams, then MPNN (8 GPU) -> full-atom af2in -> af2ig gate (+LCB3 control)
# -> unique passers -> Boltz human+mouse cross-react, into *_full2 dirs (rfd1 results preserved).
# pgrep wait is SAFE here because this is a SCRIPT FILE: its argv is "bash .../bpair_full2.sh", so the
# "rfdiffusion1_opt" pattern inside the file never appears in this process's own cmdline (no self-match).
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
log(){ echo "[bpair-full2 $(date -u +%H:%M:%S)] $*"; }
log "waiting for rfd2 RFd streams..."
sleep 60
while pgrep -f "rfdiffusion1_opt.*bpair/rfd2/des" >/dev/null 2>&1; do sleep 90; done
NB=$(ls runs/bpair/rfd2/des_*.pdb 2>/dev/null | wc -l); log "rfd2 done: $NB backbones"
[ "$NB" -ge 1 ] || { log "NO backbones - abort"; echo "BPAIR_FULL2_DONE $(date -u +%FT%TZ)"; exit 0; }
L="$R/third_party/LigandMPNN"; PY="$L/.venv/bin/python"; CKPT="$L/model_params/proteinmpnn_v_48_020.pt"
W="$R/runs/bpair/mpnn_full2"; AIN="$R/runs/bpair/af2in_full2"; AIG="$R/runs/bpair/af2ig_full2"; rm -rf "$W" "$AIN" "$AIG"; mkdir -p "$W"
mapfile -t ALL < <(ls runs/bpair/rfd2/des_*.pdb)
GL="0 1 2 3 4 5 6 7"; NG=8; declare -A LS; gi=0
for p in "${ALL[@]}"; do g=$(echo $GL|cut -d' ' -f$((gi%NG+1))); LS[$g]+="\"$(readlink -f "$p")\","; gi=$((gi+1)); done
pids=()
for g in $GL; do mkdir -p "$W/gpu$g"; printf '[%s]' "${LS[$g]%,}" > "$W/gpu$g/multi.json"
  ( cd "$L" && CUDA_VISIBLE_DEVICES=$g "$PY" run.py --model_type protein_mpnn --checkpoint_protein_mpnn "$CKPT" \
     --pdb_path_multi "$W/gpu$g/multi.json" --chains_to_design A --omit_AA C \
     --out_folder "$W/gpu$g" --batch_size 6 --number_of_batches 1 --temperature 0.15 ) > "$W/gpu$g.log" 2>&1 &
  pids+=($!); done
for p in "${pids[@]}"; do wait "$p"||true; done
NS=$(find "$W" -path '*/backbones/*.pdb'|wc -l); log "MPNN backbones: $NS"
"$R/.venv/bin/python" -I scripts/build_fullatom_af2in_tnf.py --mpnn "$W/gpu*/backbones/*.pdb" \
  --target data/targets/1TNF_B.pdb --target-chains "B" --binder-chain A --out "$AIN" 2>&1 | tail -1
cp runs/_pc_denovo/in/7JZM_lcb3.pdb "$AIN"/ 2>/dev/null && log "added LCB3 control" || true
NA=$(ls "$AIN"/*.pdb 2>/dev/null|wc -l); log "af2in complexes: $NA"
"$R/.venv/bin/python" -I scripts/run_af2ig.py --pdbdir "$AIN" --out "$AIG" --gpus "0-7" --mode fast > runs/bpair/af2ig_full2_run.log 2>&1
log "=== af2ig gate (pae<10 & plddt_binder>85 & rmsd<2) ==="
cat "$AIG"/gpu*/out.sc "$AIG"/out.sc 2>/dev/null | awk '$1=="SCORE:"&&$4!="nan"{n++;if($4+0<10&&$6+0>85&&$2+0<2)p++}END{printf "scored=%d PASS=%d\n",n,p+0}'
cat "$AIG"/gpu*/out.sc "$AIG"/out.sc 2>/dev/null | awk '$1=="SCORE:"&&$NF~/lcb3/{print "  LCB3 control: pae="$4" plddt="$6" rmsd="$2}' | head -1
cat "$AIG"/gpu*/out.sc "$AIG"/out.sc 2>/dev/null | awk '$1=="SCORE:"&&$4!="nan"&&$4+0<10&&$6+0>85&&$2+0<2&&$NF!~/lcb3/{print $4,$6,$2,$NF}' | sort -n | awk '!s[$NF]++' > runs/bpair/af2ig_full2_passers.txt
log "unique passers: $(wc -l < runs/bpair/af2ig_full2_passers.txt)"
"$R/.venv/bin/python" -I - <<'PY'
import os,glob,csv
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
aa3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def seqA(pdb):
    s={}
    for l in open(pdb):
        if l[:4]=="ATOM" and l[21]=="A" and l[12:16].strip()=="CA": s[int(l[22:26])]=aa3.get(l[17:20].strip(),"X")
    return "".join(s[k] for k in sorted(s))
rows=[]
for ln in open(f"{R}/runs/bpair/af2ig_full2_passers.txt"):
    p=ln.split()
    if len(p)<4: continue
    nm=p[3].replace("_af2pred","")
    c=glob.glob(f"{R}/runs/bpair/af2in_full2/{nm}.pdb")
    if c:
        s=seqA(c[0])
        if s and "X" not in s: rows.append((nm,s))
with open(f"{R}/runs/bpair/xreact_full2.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in rows]
print("xreact_full2 seqs:",len(rows))
PY
NXR=$(($(wc -l < runs/bpair/xreact_full2.csv)-1)); log "folding $NXR passers on human+mouse"
if [ "$NXR" -ge 1 ]; then
  nohup bash scripts/run_boltz_tnf.sh runs/bpair/xreact_full2.csv runs/buried2/human_tnf.txt bpf2_h 0 > runs/bpair/bpf2_h.log 2>&1 &
  nohup bash scripts/run_boltz_tnf.sh runs/bpair/xreact_full2.csv runs/buried2/mouse_tnf_mature.txt bpf2_m 1 > runs/bpair/bpf2_m.log 2>&1 &
  wait
  "$R/.venv/bin/python" -I - <<'PY'
import os,subprocess,glob,csv
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"
def ip(did,tag):
    bs=glob.glob(f"{R}/runs/boltz_{tag}/out/by_seed/{did}/s0")
    if not bs: return None
    b=bs[0]; cif=f"{b}/{did}_model_0.cif"
    if not os.path.exists(cif): return None
    out=cif.replace(".cif","_10_10.txt")
    if not os.path.exists(out):
        npz=glob.glob(f"{b}/pae_{did}_model_0.npz")
        if npz: subprocess.run([PY,f"{R}/tools/ipsae.py",os.path.basename(npz[0]),f"{did}_model_0.cif","10","10"],capture_output=True,cwd=b)
    if not os.path.exists(out): return None
    a=[float(l.split()[5]) for l in open(out) if len(l.split())>6 and l.split()[4]=="asym"]
    return min(a) if a else None
res=[]
for r in csv.DictReader(open(f"{R}/runs/bpair/xreact_full2.csv")):
    d=r["name"]; h=ip(d,"bpf2_h"); m=ip(d,"bpf2_m")
    if h is None or m is None: continue
    res.append((d,h,m,min(h,m)))
res.sort(key=lambda x:-x[3])
with open(f"{R}/runs/bpair/xreact_full2_results.txt","w") as fh:
    fh.write(f"{'design':28}{'h':>7}{'m':>7}{'min':>7}\n")
    for d,h,m,mn in res: fh.write(f"{d[:27]:28}{h:>7.3f}{m:>7.3f}{mn:>7.3f}\n")
print("cross-react scored:",len(res))
for d,h,m,mn in res[:15]: print(f"  {d[:27]:28} h{h:.3f} m{m:.3f} min{mn:.3f}")
PY
fi
echo "BPAIR_FULL2_DONE $(date -u +%FT%TZ)"
