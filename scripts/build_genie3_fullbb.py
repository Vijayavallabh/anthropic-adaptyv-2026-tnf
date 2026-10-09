#!/usr/bin/env python3
"""Genie3 full-backbone route (the untried fix for CA-only cross-react failure).
For each Genie3 complex (binder=chain A CA-only UNK, target=chain B CA-only, positioned):
  1. extract binder CA, rename UNK->ALA, PULCHRA-reconstruct -> binder N/CA/C/O
  2. Kabsch-superpose the REAL full-atom 1TNF_B target onto Genie3 chain B CA positions
  3. write a full-backbone complex: reconstructed binder (A) + real full-atom target (B)
Output complexes feed STANDARD full ProteinMPNN (not CA-only), optionally multi-state tied.
Usage: build_genie3_fullbb.py <pulchra_bin> <backbone_glob> <out_dir> [target_pdb]
"""
import sys, os, glob, subprocess, tempfile
import numpy as np
REPO="/home/prachh/v/anthropic-adaptyv-2026-tnf"
PULCHRA=sys.argv[1]
GLOB=sys.argv[2]
OUT=sys.argv[3]
TARGET=sys.argv[4] if len(sys.argv)>4 else f"{REPO}/data/targets/1TNF_B.pdb"
os.makedirs(OUT, exist_ok=True)

def atom_lines(path, chain=None):
    out=[]
    for ln in open(path):
        if ln.startswith(("ATOM","HETATM")) and (chain is None or ln[21]==chain):
            out.append(ln)
    return out

def ca_xyz(lines):
    return np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])]
                     for l in lines if l[12:16].strip()=="CA"])

def set_xyz(ln, xyz):
    return f"{ln[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{ln[54:]}"

def kabsch(P, Q):
    """Return R,t mapping P onto Q (least squares): Q ~ P@R.T + t."""
    Pc=P.mean(0); Qc=Q.mean(0)
    H=(P-Pc).T@(Q-Qc); U,S,Vt=np.linalg.svd(H)
    d=np.sign(np.linalg.det(Vt.T@U.T)); D=np.diag([1.0,1.0,d])
    R=Vt.T@D@U.T
    return R, Pc, Qc

# real full-atom target, in its own frame
TGT_LINES=atom_lines(TARGET)
TGT_CA=ca_xyz(TGT_LINES)
TGT_ALL=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in TGT_LINES])

files=[]
for pat in GLOB.split():
    files+=glob.glob(pat)
n_ok=0; n_fail=0
for p in sorted(files):
    name=os.path.basename(p).replace(".pdb","")
    binderA=atom_lines(p,"B".replace("B","A"))  # chain A
    targetB=atom_lines(p,"B")
    b_ca=ca_xyz(binderA); t_ca=ca_xyz(targetB)
    if len(b_ca)==0 or len(t_ca)!=len(TGT_CA):
        n_fail+=1; continue
    # 1) reconstruct binder: UNK->ALA CA pdb -> PULCHRA
    with tempfile.TemporaryDirectory() as td:
        bca=os.path.join(td,"b.pdb")
        with open(bca,"w") as fh:
            for l in binderA:
                if l[12:16].strip()=="CA":
                    fh.write(f"{l[:17]}ALA{l[20:]}")
            fh.write("END\n")
        r=subprocess.run([PULCHRA,bca],capture_output=True,text=True,cwd=td)
        reb=os.path.join(td,"b.rebuilt.pdb")
        if not os.path.exists(reb):
            n_fail+=1; continue
        # keep only N/CA/C/O of rebuilt binder
        binder_bb=[l for l in atom_lines(reb) if l[12:16].strip() in ("N","CA","C","O")]
    # 2) superpose real target onto Genie3 chain B CA frame
    R,Pc,Qc=kabsch(TGT_CA, t_ca)       # map real-target CA (P) onto genie3 B CA (Q)
    moved=(TGT_ALL-Pc)@R.T+Qc
    tgt_moved=[set_xyz(TGT_LINES[i], moved[i]) for i in range(len(TGT_LINES))]
    # 3) write complex: binder chain A + target chain B
    with open(os.path.join(OUT,f"{name}.pdb"),"w") as fh:
        for l in binder_bb:
            fh.write(f"{l[:21]}A{l[22:]}")
        fh.write("TER\n")
        for l in tgt_moved:
            fh.write(f"{l[:21]}B{l[22:]}")
        fh.write("TER\nEND\n")
    n_ok+=1

print(f"built {n_ok} full-backbone complexes -> {OUT}  (failed {n_fail})")
