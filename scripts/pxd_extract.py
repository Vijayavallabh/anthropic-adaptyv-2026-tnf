#!/usr/bin/env python3
"""PXDesign backbone outputs -> staged PDBs for SolubleMPNN -> af2ig.

PXDesign (fast mode) is a BACKBONE generator: each predictions/<task>_sample_<i>.cif holds the
TNF trimer (chains A0/B0/C0, full-atom, real seq) + N binder backbones (poly-'xpb', backbone-only
N/CA/C/O) placed around the trimer. Many binders are OFF the intended epitope (diffusion scatters
them), so we KEEP only backbones that actually contact the task's hotspot cluster, then stage each
as binder=chain A (GLY backbone, MPNN will design) + TNF trimer=chains B/C/D (full-atom, the af2ig
target).  Dedup near-identical poses; cap per task.

Epitopes (1TNF numbering, chain A0 unless noted):
  recsite: {32,65,73,75,92,144,145,146}              (kit's proven TNF example patch, one protomer)
  groove : A0{74,75,97} + B0{114,115}                 (inter-subunit receptor cleft; the TMB-1 site)

Usage: pxd_extract.py <pxd_out_dir> <stage_out_dir> [min_contacts=5] [contact_A=8.0] [cap_per_task=16]
"""
import sys, os, glob, itertools
import numpy as np
import gemmi

HOTSPOT = {
    "recsite": {"A0": [32,65,73,75,92,144,145,146]},
    "groove":  {"A0": [74,75,97], "B0": [114,115]},
}
BB = ("N","CA","C","O")

def epitope_for(task):
    for k in HOTSPOT:
        if k in task: return HOTSPOT[k]
    return None

def ca_coords(ch):
    out={}
    for r in ch:
        a=r.find_atom("CA","*")
        if a: out[r.seqid.num]=np.array([a.pos.x,a.pos.y,a.pos.z])
    return out

def binder_bb(ch):
    """list of (resnum, {atom:xyz}) backbone atoms for an xpb chain"""
    res=[]
    for r in ch:
        d={}
        for at in r:
            if at.name in BB: d[at.name]=np.array([at.pos.x,at.pos.y,at.pos.z])
        if "CA" in d: res.append((r.seqid.num,d))
    return res

def main():
    pxd_out, stage = sys.argv[1], sys.argv[2]
    min_contacts = int(sys.argv[3]) if len(sys.argv)>3 else 5
    contact_A    = float(sys.argv[4]) if len(sys.argv)>4 else 8.0
    cap          = int(sys.argv[5]) if len(sys.argv)>5 else 16
    os.makedirs(stage, exist_ok=True)
    manifest=[]
    for taskdir in sorted(glob.glob(f"{pxd_out}/tnf_*")):
        task=os.path.basename(taskdir)
        epi=epitope_for(task)
        if epi is None: continue
        cifs=sorted(glob.glob(f"{taskdir}/seed_*/predictions/*_sample_*.cif"))
        kept=[]  # (ncontact, label, binder_res, tnf_chains_structure)
        for cif in cifs:
            st=gemmi.read_structure(cif); m=st[0]
            # TNF protomers (full-atom) + hotspot CA coords
            tnf={ch.name: ch for ch in m if any(r.name!='xpb' for r in ch) and len(ch)>50}
            hot=[]  # list of hotspot CA xyz
            for cname,resl in epi.items():
                if cname in tnf:
                    cas=ca_coords(tnf[cname])
                    for rn in resl:
                        if rn in cas: hot.append(cas[rn])
            if not hot: continue
            hot=np.array(hot)
            sidx=os.path.basename(cif).split("_sample_")[1].split(".")[0]
            for ch in m:
                if not (len(ch) and ch[0].name=='xpb'): continue
                bb=binder_bb(ch)
                if len(bb)<30: continue
                ca=np.array([d["CA"] for _,d in bb])
                # per-binder-residue min dist to any hotspot CA
                dmin=np.sqrt(((ca[:,None,:]-hot[None,:,:])**2).sum(-1)).min(1)
                ncontact=int((dmin<contact_A).sum())
                if ncontact>=min_contacts:
                    kept.append((ncontact, f"{task}_s{sidx}_{ch.name}", bb, tnf))
        # dedup by pose (CA centroid proximity) keeping highest-contact, cap
        kept.sort(key=lambda x:-x[0])
        chosen=[]
        cents=[]
        for nc,label,bb,tnf in kept:
            cen=np.mean([d["CA"] for _,d in bb],axis=0)
            if any(np.linalg.norm(cen-c)<6.0 for c in cents): continue
            cents.append(cen); chosen.append((nc,label,bb,tnf))
            if len(chosen)>=cap: break
        # write staged PDB: binder=A (GLY bb), TNF=B,C,D (full-atom)
        for nc,label,bb,tnf in chosen:
            out=gemmi.Structure(); mod=gemmi.Model("1")
            # binder chain A, GLY, renumber 1..L
            ca_ch=gemmi.Chain("A")
            for i,(_,atoms) in enumerate(bb,1):
                r=gemmi.Residue(); r.name="GLY"; r.seqid=gemmi.SeqId(str(i))
                for an in BB:
                    if an in atoms:
                        a=gemmi.Atom(); a.name=an; a.element=gemmi.Element(an[0])
                        a.pos=gemmi.Position(*atoms[an]); r.add_atom(a)
                ca_ch.add_residue(r)
            mod.add_chain(ca_ch)
            # TNF protomers as B,C,D (full-atom, keep)
            for newname,cname in zip("BCD", sorted(tnf)):
                src=tnf[cname]; tc=gemmi.Chain(newname)
                for r in src:
                    if r.name=='xpb': continue
                    tc.add_residue(r.clone())
                if len(tc): mod.add_chain(tc)
            out.add_model(mod)
            pdb=f"{stage}/{label}.pdb"
            out.write_pdb(pdb)
            L=len(bb)
            manifest.append(dict(design_id=label, pdb=pdb, chain="A", length=L, ncontact=nc, task=task))
    # write manifest
    import csv
    with open(f"{stage}/manifest.csv","w") as fh:
        w=csv.DictWriter(fh, fieldnames=["design_id","pdb","chain","length","ncontact","task"])
        w.writeheader()
        for r in manifest: w.writerow(r)
    print(f"staged {len(manifest)} on-target binder backbones -> {stage}/manifest.csv")
    import collections
    print("per task:", dict(collections.Counter(r['task'] for r in manifest)))
    if manifest:
        ncs=[r['ncontact'] for r in manifest]
        print(f"contacts: min {min(ncs)} max {max(ncs)} median {int(np.median(ncs))}")

if __name__=="__main__":
    main()
