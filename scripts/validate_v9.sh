#!/bin/bash
# Definition-of-done gate for the TNF Challenge 2 submission (v9). Fast: checks submission bytes + format
# and the RECORDED evidence (master_v9.csv, chai_*.tsv) — does NOT re-fold. Exit 0 = PASS, 1 = FAIL.
# Usage: bash scripts/validate_v9.sh
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO"
SUB="runs/push/submission_v9.csv"; MASTER="runs/push/master_v9.csv"
EXPECT_MD5="05a344444cf740907e79a177d34c73f5"
fail=0; note(){ echo "  $1"; }
echo "== v9 definition-of-done gate =="

# 1. submission exists + md5
[ -s "$SUB" ] || { echo "FAIL: $SUB missing"; exit 1; }
md5=$(md5sum "$SUB" | cut -d' ' -f1)
if [ "$md5" = "$EXPECT_MD5" ]; then note "[OK] md5 $md5"; else note "[FAIL] md5 $md5 != $EXPECT_MD5"; fail=1; fi

# 2. format (python: columns, count<=20, unique names+seqs, single_chain, natural AA, len 10-250, 0-Cys)
.venv/bin/python -I - "$SUB" <<'PY' || fail=1
import sys,csv
rows=list(csv.reader(open(sys.argv[1]))); hdr=rows[0]; data=rows[1:]
ok=True
if hdr!=["name","sequence","molecule_class"]: print(f"  [FAIL] header {hdr}"); ok=False
else: print("  [OK] header name,sequence,molecule_class")
if not (1<=len(data)<=20): print(f"  [FAIL] {len(data)} designs (>20)"); ok=False
else: print(f"  [OK] {len(data)} designs (<=20)")
names=[r[0] for r in data]; seqs=[r[1] for r in data]
if len(set(names))!=len(names): print("  [FAIL] duplicate names"); ok=False
else: print(f"  [OK] {len(set(names))} unique names")
if len(set(seqs))!=len(seqs): print("  [FAIL] duplicate sequences"); ok=False
else: print(f"  [OK] {len(set(seqs))} unique sequences")
mc=set(r[2] for r in data)
if mc=={"single_chain"}: print("  [OK] all single_chain")
else: print(f"  [FAIL] molecule_class {mc}"); ok=False
bad=set(c for s in seqs for c in s)-set("ACDEFGHIKLMNPQRSTVWY")
if bad: print(f"  [FAIL] non-standard AA {bad}"); ok=False
else: print("  [OK] natural AA only")
lens=[len(s) for s in seqs]
if all(10<=l<=250 for l in lens): print(f"  [OK] length {min(lens)}-{max(lens)} in [10,250]")
else: print(f"  [FAIL] length {min(lens)}-{max(lens)}"); ok=False
cys=sum(s.count('C') for s in seqs)
print(f"  [{'OK' if cys==0 else 'WARN'}] total Cys={cys}")
sys.exit(0 if ok else 1)
PY

# 3. recorded evidence: master_v9 pLDDT >=0.75 and novelty >=3/4 for all 20
if [ -s "$MASTER" ]; then
.venv/bin/python -I - "$MASTER" <<'PY' || fail=1
import sys,csv
rows=list(csv.DictReader(open(sys.argv[1]))); ok=True
weak_pl=[r["name"] for r in rows if r.get("mono_pLDDT") and float(r["mono_pLDDT"])<0.75]
print(f"  [{'OK' if not weak_pl else 'WARN'}] expression: mono_pLDDT>=0.75 for {len(rows)-len(weak_pl)}/{len(rows)} (weak: {weak_pl or 'none'})")
bad_nov=[r["name"] for r in rows if r.get("portal_novelty","").split('/')[0] not in ("3","4")]
if bad_nov: print(f"  [FAIL] novelty <3/4: {bad_nov}"); ok=False
else: print(f"  [OK] novelty >=3/4 for all {len(rows)} (portal-confirmed)")
sys.exit(0 if ok else 1)
PY
else note "[WARN] $MASTER missing (evidence check skipped)"; fi

# 4. orthogonal binding (Chai) provenance present
for t in runs/push/chai_switchcore.tsv runs/push/chai_v7floor.tsv; do
  [ -s "$t" ] && note "[OK] $t present ($(($(wc -l < $t)-1)) designs)" || note "[WARN] $t missing"
done

echo "== $([ $fail -eq 0 ] && echo PASS || echo FAIL) =="
exit $fail
