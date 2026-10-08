#!/usr/bin/env python3
"""Select best-pae af2pred complex per distinct cb1ep2 backbone (pae<THR), copy to a binders dir,
run scan_switch on each, and report the switchable yield (backbones with >=1 His variant).
Usage: cb1ep2_select_scan.py <af2ig_dir> <out_root> [pae_thr]
"""
import sys, os, re, glob, shutil, subprocess
afdir, outroot = sys.argv[1], sys.argv[2]
THR = float(sys.argv[3]) if len(sys.argv) > 3 else 10.0
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
binders = os.path.join(outroot, "binders"); scan = os.path.join(outroot, "switch_scan")
for d in (binders, scan):
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
# best-pae desc per distinct backbone
best = {}  # backbone -> (pae, desc)
for ln in open(os.path.join(afdir, "out.sc")):
    p = ln.split()
    if p[:1] != ["SCORE:"] or len(p) < 4 or p[3] == "pae_interaction": continue
    try: pae = float(p[3])
    except ValueError: continue
    if pae != pae or pae >= THR: continue
    desc = p[-1]
    m = re.search(r"cb1_epitope_\d+_model_\d+", desc)
    if not m: continue
    bb = m.group(0)
    if bb not in best or pae < best[bb][0]: best[bb] = (pae, desc)
print(f"distinct backbones with pae<{THR}: {len(best)}")
# locate pdb, copy, scan
pdbmap = {os.path.basename(f)[:-4]: f for f in glob.glob(os.path.join(afdir, "gpu*/pdbs/*.pdb"))}
switchable, total_var, missing = 0, 0, 0
for bb, (pae, desc) in sorted(best.items(), key=lambda kv: kv[1][0]):
    src = pdbmap.get(desc)
    if not src: missing += 1; continue
    dst = os.path.join(binders, desc + ".pdb"); shutil.copy(src, dst)
    before = len(glob.glob(os.path.join(scan, "*.pdb")))
    subprocess.run(["uv", "run", "python", "scripts/scan_switch.py", dst, scan, "6"],
                   cwd=R, capture_output=True, text=True)
    nv = len(glob.glob(os.path.join(scan, "*.pdb"))) - before
    total_var += nv
    if nv > 0: switchable += 1
print(f"binders copied: {len(os.listdir(binders))}  missing-pdb: {missing}")
print(f"SWITCHABLE backbones (>=1 His var): {switchable}  total His variants: {total_var}")
