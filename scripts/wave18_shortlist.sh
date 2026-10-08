#!/bin/bash
# WAVE 18 efficient takeover (run AFTER results/novelty_w18.csv exists, with the pipeline stopped).
# Instead of the expensive His433 switch install-scan on ALL ~1785 novelty-survivors, narrow the 2380
# af2ig binders to a DIVERSE shortlist (spread across seed lineages + folds; user directive: distinct
# designs, not cousins of a few) and run the switch only on the shortlist. Priority downstream of this:
# switch > cross-react > pKAI > MD, all on the shortlist. The portal is the novelty oracle.
# Usage: wave18_shortlist.sh [N=200] [MAX_PER_LINEAGE=15]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
MD="$REPO/envs/md/bin/python"; PR="$REPO/envs/pyrosetta/bin/python"
N=${1:-200}; MAXPL=${2:-15}
export TMPDIR="$REPO/runs/tmp"; mkdir -p "$TMPDIR"

# 1) pool = af2 binding (pae) joined with novelty (TM/level); drop only CLEAR structural fails (TM>=0.90)
"$MD" - <<'PY'
import csv, os
af = {r["design"]: r for r in csv.DictReader(open("results/rfd_af2_passers_w18.csv"))}
nov = {r["name"]: r for r in csv.DictReader(open("results/novelty_w18.csv"))} if os.path.exists("results/novelty_w18.csv") else {}
rows = []
for name, a in af.items():
    n = nov.get(name)
    if n:
        try: tm = float(n.get("struct_tm", "") or 0)
        except ValueError: tm = 0.0
        if tm >= 0.90:            # clear structural fail -> drop (no known passer sits here)
            continue
        lvl, pss, tmv = n.get("level", ""), n.get("pass", ""), n.get("struct_tm", "")
    else:
        lvl, pss, tmv = "", "", ""   # no novelty score -> keep (portal decides)
    rows.append({"name": name, "pae_interaction": a["pae_interaction"], "plddt_binder": a["plddt_binder"],
                 "struct_tm": tmv, "novelty_level": lvl, "novelty_pass": pss})
with open("results/pool_w18.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["name", "pae_interaction", "plddt_binder", "struct_tm", "novelty_level", "novelty_pass"])
    w.writeheader(); w.writerows(rows)
print(f"pool_w18.csv: {len(rows)} binders kept (dropped clear novelty-fails TM>=0.90)")
PY

# 2) diverse shortlist (spread across lineages + folds), ranked within by binding (low pae)
"$MD" scripts/select_diverse.py --pool results/pool_w18.csv --pdbs runs/novelty_w18_struct \
    --n "$N" --max-per-lineage "$MAXPL" --score-col pae_interaction --lower-better \
    --out results/shortlist_w18.csv --work runs/diverse_sel/w18

# 3) gather the shortlist COMPLEXES (binder+target) + canon-renumber
sl=runs/shortlist_w18_cplx; mkdir -p "$sl"; rm -f "$sl"/*.pdb
"$MD" - <<'PY'
import csv, shutil, os
k = 0
for r in csv.DictReader(open("results/shortlist_w18.csv")):
    src = f"runs/rfd_passers_w18/{r['name']}.pdb"
    if os.path.exists(src):
        shutil.copy(src, "runs/shortlist_w18_cplx/"); k += 1
print("shortlist complexes copied:", k)
PY
"$MD" scripts/canon_renumber.py --indir runs/shortlist_w18_cplx --outdir runs/shortlist_w18_canon 2>&1 | tail -2

# 4) His433 pH-switch (PROPKA install-scan) on the shortlist ONLY -- fast vs the full pool
export PATH="$REPO/envs/pyrosetta/bin:$PATH"
"$PR" scripts/obj1_switch.py --pdbs runs/shortlist_w18_canon/*.pdb \
    --outdir runs/rfd_obj1_w18 --out results/rfd_obj1_switch_w18.csv --workers 48 2>&1 | tail -4

echo "WAVE18_SHORTLIST_SWITCH_DONE $(date -u +%FT%TZ)"
echo "switch rows: $(($(wc -l < results/rfd_obj1_switch_w18.csv 2>/dev/null || echo 1)-1))"
echo "NEXT: cross-react (run_xreact.sh runs/rfd_obj1_w18 w18 0-7) + pKAI on switch-passers (ddg_pH>0), then final diverse-12."
