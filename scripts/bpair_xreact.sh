#!/bin/bash
# beta-pairing cross-react stage: wait for the af2ig funnel, then fold the af2ig passers on human+mouse
# TNF (Boltz) and compute cross-react min-ipSAE (worst ortholog). beta-pairing binds the conserved edge
# strand BACKBONE -> expect intrinsic cross-react. Writes runs/bpair/xreact_results.txt.
set -uo pipefail
R=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$R"
log(){ echo "[bpair-xreact $(date -u +%H:%M:%S)] $*"; }
log "waiting for BPAIR_FUNNEL_DONE..."
while ! grep -q BPAIR_FUNNEL_DONE runs/bpair/funnel.log 2>/dev/null; do sleep 60; done
AIN="$R/runs/bpair/af2in"
# 1. build csv (name,sequence) from af2in chain-A of each af2ig passer (robust: seq from the folded backbone)
"$R/.venv/bin/python" -I - <<'PY'
import os,glob
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"
aa3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def chainA_seq(pdb):
    seq={}
    for l in open(pdb):
        if l[:4]=="ATOM" and l[21]=="A" and l[12:16].strip()=="CA":
            seq[int(l[22:26])]=aa3.get(l[17:20].strip(),"X")
    return "".join(seq[k] for k in sorted(seq))
names=[]
pf=f"{R}/runs/bpair/af2ig_passers.txt"
if os.path.exists(pf):
    for ln in open(pf):
        p=ln.split()
        if len(p)>=4: names.append(p[3].replace("_af2pred",""))
import csv
rows=[]
for nm in names:
    cands=glob.glob(f"{R}/runs/bpair/af2in/{nm}.pdb")
    if not cands: continue
    s=chainA_seq(cands[0])
    if s and "X" not in s: rows.append((nm,s))
with open(f"{R}/runs/bpair/xreact.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["name","sequence"]); [w.writerow(r) for r in rows]
print("xreact seqs:",len(rows))
PY
NS=$(($(wc -l < runs/bpair/xreact.csv)-1)); log "folding $NS passers on human+mouse"
[ "$NS" -lt 1 ] && { log "no passers to fold; done"; echo "BPAIR_XREACT_DONE $(date -u +%FT%TZ)"; exit 0; }
# 2. Boltz human (gpu4) + mouse (gpu5)
nohup bash scripts/run_boltz_tnf.sh runs/bpair/xreact.csv runs/buried2/human_tnf.txt bp_h 4 > runs/bpair/bp_h.log 2>&1 &
nohup bash scripts/run_boltz_tnf.sh runs/bpair/xreact.csv runs/buried2/mouse_tnf_mature.txt bp_m 5 > runs/bpair/bp_m.log 2>&1 &
wait
# 3. cross-react min-ipSAE
"$R/.venv/bin/python" -I - <<'PY'
import os,subprocess,glob,csv
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; PY=f"{R}/third_party/LigandMPNN/.venv/bin/python"
def ipsae(did,tag):
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
for r in csv.DictReader(open(f"{R}/runs/bpair/xreact.csv")):
    d=r["name"]; h=ipsae(d,"bp_h"); m=ipsae(d,"bp_m")
    if h is None or m is None: continue
    res.append((d,h,m,min(h,m)))
res.sort(key=lambda x:-x[3])
with open(f"{R}/runs/bpair/xreact_results.txt","w") as fh:
    fh.write(f"{'design':28}{'h':>7}{'m':>7}{'min':>7}\n")
    for d,h,m,mn in res: fh.write(f"{d[:27]:28}{h:>7.3f}{m:>7.3f}{mn:>7.3f}\n")
print(f"cross-react scored: {len(res)}")
print("top 10 by min-ipSAE (de-novo b-pairing; ref des_59_1 min 0.564):")
for d,h,m,mn in res[:10]: print(f"  {d[:27]:28} h{h:.3f} m{m:.3f} min{mn:.3f}")
PY
echo "BPAIR_XREACT_DONE $(date -u +%FT%TZ)"
