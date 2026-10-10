#!/usr/bin/env python3
"""ESMFold single-seq MONOMER folder: fold each binder ALONE (no MSA) -> PDB + mean-CA-pLDDT.
Per Garcia/Dixit/Rocklin 2026 (PMC12817478) monomer pLDDT is the best de-novo success predictor
(expressed+soluble+monomeric+folded, AUROC ~0.72). Also yields the single-seq structure for the
novelty proxy. Run with envs/uplifting-esmc/bin/python -I.
Usage: esmfold_monomer.py <csv name,sequence> <outdir>"""
import sys, os, csv, json
import torch
from transformers import AutoTokenizer, EsmForProteinFolding
csvp, outdir = sys.argv[1], sys.argv[2]; os.makedirs(outdir, exist_ok=True)
print("loading ESMFold v1 ...", flush=True)
model=EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True).cuda().eval()
try: model.esm=model.esm.half()
except Exception: pass
torch.backends.cuda.matmul.allow_tf32=True
rows=list(csv.DictReader(open(csvp)))
res={}
for r in rows:
    name=r.get("name") or r.get("design_id"); seq=(r.get("sequence") or "").split(":")[0].strip()
    if not name or not seq: continue
    try:
        with torch.no_grad(): pdb=model.infer_pdb(seq)
    except Exception as e:
        print(f"{name}: FOLD FAIL {str(e)[:80]}", flush=True); continue
    open(f"{outdir}/{name}.pdb","w").write(pdb)
    bf=[float(l[60:66]) for l in pdb.splitlines() if l[:4]=="ATOM" and l[12:16].strip()=="CA"]
    plddt=sum(bf)/len(bf) if bf else 0.0
    res[name]=round(plddt,2); print(f"{name:26} pLDDT={plddt:6.2f} L={len(seq)}", flush=True)
json.dump(res, open(f"{outdir}/plddt.json","w"), indent=1)
print(f"\nwrote {len(res)} monomer folds + plddt.json -> {outdir}", flush=True)
