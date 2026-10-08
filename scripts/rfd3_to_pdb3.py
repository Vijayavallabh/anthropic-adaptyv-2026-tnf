import sys,gzip,glob,os
# RFd3 cif.gz -> PDB keeping binder chain A + BOTH target protomers (B,C). For the Region-I two-protomer
# epitope. af2ig/MPNN-ready. Usage: rfd3_to_pdb3.py "<glob of cif.gz>" <outdir>
outdir=sys.argv[2]; os.makedirs(outdir,exist_ok=True)
KEEP=['A','B','C']
n=0
for fn in sorted(glob.glob(sys.argv[1])):
    rows={c:[] for c in KEEP}
    op=gzip.open if fn.endswith('.gz') else open
    for ln in op(fn,'rt'):
        if not ln.startswith('ATOM'): continue
        f=ln.split()
        atom=f[2]; resn=f[4]; ch=f[5]; seqid=int(f[7]); el=f[1]
        x,y,z=float(f[18]),float(f[19]),float(f[20])
        if ch not in rows: continue
        rows[ch].append((seqid,atom,resn,el,x,y,z))
    if not rows['A'] or not rows['B'] or not rows['C']: continue
    out=[]; serial=1
    for ch in KEEP:
        for seqid,atom,resn,el,x,y,z in rows[ch]:
            out.append(f"ATOM  {serial:>5} {atom:<4} {resn:>3} {ch}{seqid:>4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}\n")
            serial+=1
        out.append("TER\n")
    out.append("END\n")
    base=os.path.basename(os.path.dirname(fn))+'__'+os.path.basename(fn).replace('.cif.gz','').replace('.cif','')
    open(os.path.join(outdir,base+'.pdb'),'w').write("".join(out)); n+=1
print(f"converted {n} -> {outdir} (binder A + target B,C)")
