#!/usr/bin/env python3
"""Genie3 full-backbone route, TWO-PROTOMER (Region-I groove) variant of build_genie3_fullbb.py.
Genie3 emits binder=chain A (CA-only UNK) + target chains B,C (CA-only, in Genie3's own frame).
For each Genie3 complex:
  1. extract binder CA, UNK->ALA, PULCHRA-reconstruct -> binder N/CA/C/O
  2. JOINT Kabsch-superpose the REAL full-atom 2-protomer target (chains B+C) onto the Genie3 B+C CA
     positions using ALL B+C CA at once (one rigid transform preserves the groove docking geometry)
  3. write a full-backbone complex: reconstructed binder (A) + real full-atom target (B, C)
Output feeds standard full ProteinMPNN (design chain A) -> build_fullatom_af2in_tnf.py --target-chains B,C.
Usage: build_genie3_fullbb_bc.py <pulchra_bin> <backbone_glob> <out_dir> [target_pdb] [chains=B,C]
"""
import sys, os, glob, subprocess, tempfile, re
import numpy as np
REPO="/home/prachh/v/anthropic-adaptyv-2026-tnf"
PULCHRA=sys.argv[1]; GLOB=sys.argv[2]; OUT=sys.argv[3]
TARGET=sys.argv[4] if len(sys.argv)>4 else f"{REPO}/data/targets/1TNF_regI_BC.pdb"
CHAINS=(sys.argv[5] if len(sys.argv)>5 else "B,C").split(",")
os.makedirs(OUT, exist_ok=True)

def atom_lines(path, chain=None):
    return [ln for ln in open(path)
            if ln.startswith(("ATOM","HETATM")) and (chain is None or ln[21]==chain)]

def ca_xyz(lines):
    return np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])]
                     for l in lines if l[12:16].strip()=="CA"])

def set_xyz(ln, xyz):
    return f"{ln[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{ln[54:]}"

def kabsch(P, Q):
    """Return transform mapping P onto Q (least squares): moved = (P-Pc)@R.T + Qc."""
    Pc=P.mean(0); Qc=Q.mean(0)
    H=(P-Pc).T@(Q-Qc); U,S,Vt=np.linalg.svd(H)
    d=np.sign(np.linalg.det(Vt.T@U.T)); D=np.diag([1.0,1.0,d])
    R=Vt.T@D@U.T
    return R, Pc, Qc

# real full-atom target: gather per-chain lines + CA, in CHAINS order, plus all-atom coords
TGT_LINES=[]; TGT_CA=[]
for ch in CHAINS:
    ls=atom_lines(TARGET, ch); TGT_LINES+=ls
    TGT_CA.append(ca_xyz(ls))
TGT_CA=np.vstack(TGT_CA)
TGT_ALL=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in TGT_LINES])

files=[]
for pat in GLOB.split(): files+=glob.glob(pat)
# optional sharding for parallel workers: process files[SHARD_I::SHARD_N]
SI=int(os.environ.get("SHARD_I","0")); SN=int(os.environ.get("SHARD_N","1"))
files=sorted(files)[SI::SN]
n_ok=0; n_fail=0; n_skip=0
for p in files:
    # include the generation shard in the name: shard_g4/.../tnf_regI_0.pdb -> g4_tnf_regI_0
    _m=re.search(r'shard_(g\d+)', p); _pfx=(_m.group(1)+"_" if _m else "")
    name=_pfx+os.path.basename(p).replace(".pdb","")
    if os.path.exists(os.path.join(OUT,f"{name}.pdb")): n_skip+=1; n_ok+=1; continue  # reuse already-built
    binderA=atom_lines(p,"A")
    # genie3 target CA in CHAINS order (must match real target CA count/order)
    g_ca=[]
    ok=True
    for ch in CHAINS:
        c=ca_xyz(atom_lines(p, ch))
        if len(c)==0: ok=False; break
        g_ca.append(c)
    if not ok:
        n_fail+=1; continue
    g_ca=np.vstack(g_ca)
    b_ca=ca_xyz(binderA)
    if len(b_ca)==0 or g_ca.shape!=TGT_CA.shape:
        n_fail+=1; continue
    # 1) reconstruct binder backbone via PULCHRA
    with tempfile.TemporaryDirectory() as td:
        bca=os.path.join(td,"b.pdb")
        with open(bca,"w") as fh:
            for l in binderA:
                if l[12:16].strip()=="CA": fh.write(f"{l[:17]}ALA{l[20:]}")
            fh.write("END\n")
        subprocess.run([PULCHRA,bca],capture_output=True,text=True,cwd=td)
        reb=os.path.join(td,"b.rebuilt.pdb")
        if not os.path.exists(reb): n_fail+=1; continue
        binder_bb=[l for l in atom_lines(reb) if l[12:16].strip() in ("N","CA","C","O")]
    # 2) JOINT superpose real target (B+C) onto Genie3 (B+C) CA frame
    R,Pc,Qc=kabsch(TGT_CA, g_ca)        # map real-target CA (P) onto genie3 B+C CA (Q)
    moved=(TGT_ALL-Pc)@R.T+Qc
    tgt_moved=[set_xyz(TGT_LINES[i], moved[i]) for i in range(len(TGT_LINES))]
    # 3) write complex: binder A + real target (original chain IDs preserved, e.g. B,C)
    with open(os.path.join(OUT,f"{name}.pdb"),"w") as fh:
        for l in binder_bb: fh.write(f"{l[:21]}A{l[22:]}")
        fh.write("TER\n")
        prev=None
        for l in tgt_moved:
            if prev is not None and l[21]!=prev: fh.write("TER\n")
            fh.write(l); prev=l[21]
        fh.write("TER\nEND\n")
    n_ok+=1
print(f"built {n_ok-n_skip} (+{n_skip} reused) two-protomer full-backbone complexes -> {OUT}  (failed {n_fail}; chains={CHAINS}; shard {SI}/{SN})")
