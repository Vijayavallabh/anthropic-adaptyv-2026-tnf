#!/usr/bin/env python3
"""Build the two-state MD manifest for cb1-epitope His-switch variants that RETAINED binding.
Reads the switch-af2ig scorefile, keeps pae_interaction<14 & plddt_binder>80 & rmsd<2.5, parses the
installed His position from the variant name (..._A<pos><aa>2H_vsK<anchor>), locates the af2pred
complex, and writes a manifest: name<TAB>his_positions<TAB>lineage<TAB>pdb_path  (one His per row).
Usage: cb1ep_md_manifest.py <switch_af2ig_dir> <out_manifest.tsv> [max_designs]
"""
import sys, os, re, glob
scdir, out = sys.argv[1], sys.argv[2]
maxd = int(sys.argv[3]) if len(sys.argv) > 3 else 40
sc = os.path.join(scdir, "out.sc")
rows = []
hdr = None
for ln in open(sc):
    p = ln.split()
    if p[:1] == ["SCORE:"] and "description" in ln: hdr = p[1:]
    elif p[:1] == ["SCORE:"] and hdr:
        d = dict(zip(hdr, p[1:]))
        try:
            pae = float(d["pae_interaction"]); plb = float(d["plddt_binder"]); rmsd = float(d["binder_aligned_rmsd"])
        except (ValueError, KeyError): continue
        if pae != pae or plb != plb or rmsd != rmsd: continue  # nan
        if pae < 14 and plb > 80 and rmsd < 2.5:
            rows.append((pae, d["description"]))
rows.sort()
# dedup by PARENT BACKBONE (distinct de-novo fold) -> keep best-pae switch variant per backbone,
# so the MD-gated set is fold-diverse rather than many variants of one good backbone.
seen = set(); n = 0
with open(out, "w") as fh:
    for pae, desc in rows:
        m = re.search(r"_A(\d+)[A-Z]2H_vsK", desc)   # installed His position on binder chain A
        if not m: continue
        bb = re.search(r"(cb1_epitope_\d+_model_\d+)", desc)
        bbid = bb.group(1) if bb else desc
        if bbid in seen: continue
        pdb = glob.glob(os.path.join(scdir, "**", "pdbs", desc + ".pdb"), recursive=True) or \
              glob.glob(os.path.join(scdir, "pdbs", desc + ".pdb"))
        if not pdb: continue
        seen.add(bbid)
        fh.write(f"{desc}\t{m.group(1)}\tcb1ep\t{pdb[0]}\n")
        n += 1
        if n >= maxd: break
print(f"wrote {n} MD-manifest rows (best switch variant per distinct backbone) -> {out}")
