#!/bin/bash
# Sequence-novelty gate: no binder >70% id to a natural protein (PDB seqres). TNF target = fail-loud positive control (MUST hit).
set -uo pipefail; cd "$(dirname "$0")/.."
NOV=runs/novelty; mkdir -p "$NOV"; export PATH="$NOV/mmseqs/bin:$PATH"
die(){ echo "NOVELTY-FAIL: $*"; exit 1; }
command -v mmseqs >/dev/null 2>&1 || { echo "fetching mmseqs..."; curl -fsSL https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz | tar xz -C "$NOV" || die "mmseqs fetch"; }
command -v mmseqs >/dev/null 2>&1 || die "mmseqs not on PATH"
# query = binders + TNF mature (positive control, UniProt P01375 mature)
TNF="VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL"
cp runs/submission/binders.fasta "$NOV/query.fasta" || die "runs/submission/binders.fasta missing (gate would otherwise pass on positive control alone)"
printf ">POSCTRL_TNF_target\n%s\n" "$TNF" >> "$NOV/query.fasta"
DB="$NOV/pdb_seqres.fasta"
[ -s "$DB" ] || { echo "fetching pdb_seqres..."; curl -fsSL https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz | gunzip > "$DB" || curl -fsSL https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz | gunzip > "$DB" || die "pdb_seqres fetch"; }
echo "searching $(grep -c '>' "$NOV/query.fasta") seqs vs PDB ($(grep -c '>' "$DB") entries)..."
mmseqs easy-search "$NOV/query.fasta" "$DB" "$NOV/hits.m8" "$NOV/tmp" --format-output "query,target,pident,alnlen,evalue" -s 7.5 1>"$NOV/mmseqs.log" 2>&1 || die "mmseqs search"
grep -qi POSCTRL "$NOV/hits.m8" || die "positive control (TNF) produced NO PDB hit -> search broken"
echo "=== NOVELTY (max %id to PDB per binder; want <70%, POSCTRL must be ~100%) ==="
awk '{if($3>m[$1])m[$1]=$3} END{for(q in m) printf "  %-22s max_id=%.1f%%\n",q,m[q]}' "$NOV/hits.m8" | sort
echo "(binders with NO line above = no PDB hit at all = fully novel)"
echo "NOVELTY_DONE"
