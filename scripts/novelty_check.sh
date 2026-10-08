#!/bin/bash
# Track-3 novelty gate: no design >70% sequence identity to a natural protein (PDB seqres).
# Self-contained: fetches mmseqs2 (static) + pdb_seqres once, then mmseqs easy-search the submission.
# FAIL-LOUD: every fetch uses curl -f, the search exit code is checked, and a POSITIVE CONTROL (the
# EGFR DIII target, which is in the PDB) must hit -- otherwise an empty hits.m8 could be a silent tool
# failure masquerading as "all novel" (a disqualification risk). Usage: novelty_check.sh [submission.csv]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
CSV=${1:-results/submission_track3_v2.csv}
NOV="$REPO/tools/novelty"; DB="$NOV/db"; mkdir -p "$NOV" "$DB" runs/novelty
export PATH="$NOV/mmseqs/bin:$PATH"
die(){ echo "NOVELTY_CHECK_FAILED: $*" >&2; exit 2; }

if ! command -v mmseqs >/dev/null 2>&1; then
  echo "fetching mmseqs2..."; curl -fsSL https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz | tar xz -C "$NOV" \
    || die "mmseqs fetch failed"
fi
command -v mmseqs >/dev/null 2>&1 || die "mmseqs still not on PATH after fetch"
if [ ! -f "$DB/pdb_seqres.fasta" ]; then
  echo "fetching pdb_seqres..."
  curl -fsSL https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz -o "$DB/pdb_seqres.txt.gz" \
    || die "pdb_seqres fetch failed"
  gunzip -kf "$DB/pdb_seqres.txt.gz" || die "pdb_seqres gunzip failed"
  awk '/^>/{p=($0 ~ /mol:protein/)} {if(p)print}' "$DB/pdb_seqres.txt" > "$DB/pdb_seqres.fasta"
  echo "pdb_seqres protein entries: $(grep -c '^>' "$DB/pdb_seqres.fasta")"
fi
[ -s "$DB/pdb_seqres.fasta" ] || die "pdb_seqres.fasta missing/empty"

# query = submitted designs + a POSITIVE CONTROL (EGFR DIII target, must hit the PDB)
envs/md/bin/python -c "
import csv
for r in csv.DictReader(open('$CSV')):
    print('>'+r['name']); print(r['sequence'].split(':')[0])
" > runs/novelty/final20.fasta
NQ=$(grep -c '^>' runs/novelty/final20.fasta)
[ "$NQ" -gt 0 ] || die "no query sequences from $CSV"
cp runs/novelty/final20.fasta runs/novelty/query.fasta
printf ">POSCTRL_egfr_dom3\n%s\n" "$(tr -d '\n ' < data/targets/human_domainIII.seq)" >> runs/novelty/query.fasta

echo "searching $NQ designs (+1 positive control) vs PDB ..."
rm -rf runs/novelty/tmp
mmseqs easy-search runs/novelty/query.fasta "$DB/pdb_seqres.fasta" runs/novelty/hits.m8 runs/novelty/tmp \
  --format-output query,target,pident,alnlen,qlen,tlen,qcov,evalue -s 7.5 --max-seqs 50 \
  1>runs/novelty/mmseqs.log 2>&1 || die "mmseqs easy-search exited non-zero (see runs/novelty/mmseqs.log)"
grep -q "^POSCTRL_egfr_dom3" runs/novelty/hits.m8 \
  || die "positive control (EGFR DIII) produced NO PDB hit -> search is broken, result untrustworthy"
echo "NOVELTY_SEARCH_DONE (positive control hit -> search verified working)"
envs/md/bin/python -c "
best={}
for ln in open('runs/novelty/hits.m8'):
    f=ln.rstrip('\n').split('\t'); q,t,pid=f[0],f[1],float(f[2]); qcov=float(f[6])
    if q=='POSCTRL_egfr_dom3' or qcov<0.5: continue
    if q not in best or pid>best[q][0]: best[q]=(pid,t,qcov)
qs=[l[1:].strip() for l in open('runs/novelty/final20.fasta') if l.startswith('>')]
flagged=0
for q in qs:
    if q in best:
        pid,t,qc=best[q]; mark=' <== >70% NOT NOVEL' if pid>70 else ''
        flagged+= pid>70
        print(f'{q[:34]:34s} max_pident={pid:5.1f}% qcov={qc:.2f} {t}{mark}')
    else:
        print(f'{q[:34]:34s} NO PDB HIT (novel)')
print(f'=== NOVELTY: {flagged}/{len(qs)} flagged >70% identity (0 = all de novo, PASS) ===')
"
