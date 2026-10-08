#!/usr/bin/env python3
"""Thread SolubleMPNN sequences onto PXDesign binder backbones -> af2ig input PDBs.

For each staged complex (binder=chain A GLY backbone + TNF trimer B/C/D full-atom) and each MPNN
designed sequence, write a PDB with chain A renamed to the designed sequence (backbone kept; af2ig
builds side chains) + the TNF protomers the binder actually contacts (full-atom, the af2ig target).
Binder=A renumbered 1.., target protomers=B,C renumbered 1.. . Target MUST stay full-atom (af2ig
templates it; backbone-only target -> plddt_target~40 -> 0 binders).

Usage: pxd_make_af2ig.py <stage_dir> <mpnn_dir> <out_af2in_dir> [contact_A=10] [max_target_chains=2]
"""
import sys, os, glob, csv
import numpy as np
import gemmi

AA3={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY','H':'HIS',
     'I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR','W':'TRP',
     'Y':'TYR','V':'VAL'}

def designed_seqs(fa):
    """return list of designed binder sequences (first ':'-field of each id=N record)."""
    out=[]
    lines=open(fa).read().splitlines()
    for i,ln in enumerate(lines):
        if ln.startswith(">") and "id=" in ln:
            seq=lines[i+1].split(":")[0].strip()
            out.append(seq)
    return out

def ca(chain):
    out=[]
    for r in chain:
        a=r.find_atom("CA","*")
        if a: out.append([a.pos.x,a.pos.y,a.pos.z])
    return np.array(out) if out else np.zeros((0,3))

def main():
    stage, mpnn, out = sys.argv[1], sys.argv[2], sys.argv[3]
    contact_A = float(sys.argv[4]) if len(sys.argv)>4 else 10.0
    maxt = int(sys.argv[5]) if len(sys.argv)>5 else 2
    os.makedirs(out, exist_ok=True)
    rows=list(csv.DictReader(open(f"{stage}/manifest.csv")))
    n=0
    for r in rows:
        did=r["design_id"]; pdb=r["pdb"]
        fa=f"{mpnn}/{did}/seqs/{did}.fa"
        if not (os.path.exists(pdb) and os.path.exists(fa)): continue
        seqs=designed_seqs(fa)
        if not seqs: continue
        st=gemmi.read_structure(pdb); m=st[0]
        binder=m["A"]; bca=ca(binder)
        # pick contacted TNF protomers
        contacts=[]
        for ch in m:
            if ch.name=="A": continue
            tca=ca(ch)
            if len(tca)==0: continue
            d=np.sqrt(((bca[:,None,:]-tca[None,:,:])**2).sum(-1))
            nc=int((d<contact_A).any(1).sum()) if d.size else 0  # binder residues contacting this protomer
            if nc>=3: contacts.append((nc,ch.name))
        contacts.sort(reverse=True)
        keep=[c for _,c in contacts[:maxt]]
        if not keep: keep=[sorted(ch.name for ch in m if ch.name!="A")[0]]  # fallback: 1 chain
        L=len(bca)
        for si,seq in enumerate(seqs):
            if len(seq)!=L: continue
            out_st=gemmi.Structure(); mod=gemmi.Model("1")
            # binder chain A with designed seq
            bc=gemmi.Chain("A")
            idx=0
            for res in binder:
                if res.find_atom("CA","*") is None: continue
                nr=gemmi.Residue(); nr.name=AA3.get(seq[idx],"ALA"); nr.seqid=gemmi.SeqId(str(idx+1))
                for at in res:
                    if at.name in ("N","CA","C","O"):
                        a=gemmi.Atom(); a.name=at.name; a.element=at.element; a.pos=at.pos; nr.add_atom(a)
                bc.add_residue(nr); idx+=1
            mod.add_chain(bc)
            # target protomers full-atom, relabel B,C..
            for newname,cname in zip("BCDE", keep):
                src=m[cname]; tc=gemmi.Chain(newname); j=0
                for res in src:
                    if res.find_atom("CA","*") is None: continue
                    j+=1; nr=res.clone(); nr.seqid=gemmi.SeqId(str(j)); tc.add_residue(nr)
                mod.add_chain(tc)
            out_st.add_model(mod)
            out_st.write_pdb(f"{out}/{did}_m{si+1}.pdb")
            n+=1
    print(f"wrote {n} af2ig input PDBs -> {out}")

if __name__=="__main__":
    main()
