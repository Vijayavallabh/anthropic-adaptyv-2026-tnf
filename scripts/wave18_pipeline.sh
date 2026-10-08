#!/bin/bash
# WAVE 18 novelty-AWARE downstream. Unlike the generic wave17 pipeline (binding + His433 switch only),
# this encodes the constraints that actually gate this submission, in priority order:
#   0-Cys (soluble-MPNN) -> binding (af2ig) -> NOVELTY (local proxy; drop clear structural fails,
#   annotate + rank the rest; the Adaptyv portal is the final novelty oracle) -> His433 pH-switch
#   (PROPKA) -> [cross-react + pKAI + deep MD run by hand on the shortlist].
# The novelty filter sits AFTER af2ig (needs the folded binder) and BEFORE the expensive switch scan,
# per the recovery plan. Proxy is a RANKER, not a gate: it only drops foldseek-TM>=0.90 / seq>70 (safe:
# no known Adaptyv-passer sits there); everything else is kept for validation and the portal decides.
# Usage: wave18_pipeline.sh <wavelog> <rfd_dirs_glob> <tag> <max_pae>
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-egfr; cd "$REPO"
WAVELOG=${1:-scratchpad/win/waveB_pd.log}; GLOB=${2:-runs/rfd/pdnov_egfr_phsw_*}; TAG=${3:-w18}; MAXPAE=${4:-11}
MD="$REPO/envs/md/bin/python"
echo "[$TAG] waiting for generation (PD_DIVERSE_WAVE2_DONE in $WAVELOG) ... $(date -u +%FT%TZ)"
for i in $(seq 1 1200); do
  grep -q PD_DIVERSE_WAVE2_DONE "$WAVELOG" 2>/dev/null && { echo "[$TAG] generation done after ~$((i*30))s"; break; }
  sleep 30
done

# 1) 0-Cys soluble-MPNN + full-atom-target af2ig (binding + expression gate)
bash scripts/run_betA_downstream.sh "$GLOB" "0-7" 2 "$TAG" A

# 2) collect af2ig binders (pae<MAXPAE)
"$MD" scripts/collect_passers.py --runs "runs/rfd_af2ig_$TAG" --outdir "runs/rfd_passers_$TAG" \
    --manifest "results/rfd_af2_passers_$TAG.csv" --max-pae "$MAXPAE"
npass=$(ls runs/rfd_passers_$TAG/*.pdb 2>/dev/null | wc -l)
echo "[$TAG] af2ig binders: $npass"
[ "$npass" -eq 0 ] && { echo "[$TAG] no binders; stop"; exit 0; }

# 3) NOVELTY: extract binder monomers (chain A) + score with the local proxy
nd="runs/novelty_${TAG}_struct"; mkdir -p "$nd"
"$MD" - "$TAG" "$nd" <<'PY'
import sys, glob, os, csv
tag, nd = sys.argv[1], sys.argv[2]
fa = open(f"runs/novelty_{tag}.fasta", "w")
aa3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
 'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
for p in sorted(glob.glob(f"runs/rfd_passers_{tag}/*.pdb")):
    name = os.path.splitext(os.path.basename(p))[0]
    seq = []; seen = set()
    with open(p) as fh, open(f"{nd}/{name}.pdb", "w") as o:
        for l in fh:
            if l[:4] == "ATOM" and l[21] == "A":
                o.write(l)
                rid = l[22:27]
                if l[12:16].strip() == "CA" and rid not in seen:
                    seen.add(rid); seq.append(aa3to1.get(l[17:20].strip(), "X"))
        o.write("END\n")
    if seq: fa.write(f">{name}\n{''.join(seq)}\n")
fa.close()
print(f"extracted {len(glob.glob(f'{nd}/*.pdb'))} binder monomers")
PY
python3 scripts/novelty_proxy.py --pdbs "$nd" --seqs "runs/novelty_${TAG}.fasta" \
    --out "results/novelty_${TAG}.csv" --work "runs/novelty/proxy_${TAG}" --threads 32 || echo "[novelty] proxy warn"

# 4) keep designs that are NOT clear novelty-fails (foldseek-TM<0.90 and seq<=70), rank by TM
"$MD" - "$TAG" <<'PY'
import sys, csv, os, shutil
tag = sys.argv[1]
keep = []
try:
    for r in csv.DictReader(open(f"results/novelty_{tag}.csv")):
        if float(r["struct_tm"]) < 0.90 and float(r["seqid"]) <= 70:
            keep.append(r["name"])
except FileNotFoundError:
    pass
src = f"runs/rfd_passers_{tag}"; dst = f"runs/rfd_passers_{tag}_novel"; os.makedirs(dst, exist_ok=True)
n = 0
for name in keep:
    p = os.path.join(src, name + ".pdb")
    if os.path.exists(p): shutil.copy(p, dst); n += 1
print(f"[novelty] kept {n} not-clear-fail binders (of {len(os.listdir(src))}) for switch scoring -> {dst}")
PY

# 5) canonicalize + His433 pH-switch (PROPKA) on the novelty-surviving binders
export TMPDIR="$REPO/runs/tmp"; mkdir -p "$TMPDIR"; export PATH="$REPO/envs/pyrosetta/bin:$PATH"
PR="$REPO/envs/pyrosetta/bin/python"
"$MD" scripts/canon_renumber.py --indir "runs/rfd_passers_${TAG}_novel" --outdir "runs/rfd_passers_${TAG}_novel_canon" 2>&1 | tail -2
ncan=$(ls runs/rfd_passers_${TAG}_novel_canon/*.pdb 2>/dev/null | wc -l)
[ "$ncan" -gt 0 ] && "$PR" scripts/obj1_switch.py --pdbs runs/rfd_passers_${TAG}_novel_canon/*.pdb \
    --outdir "runs/rfd_obj1_$TAG" --out "results/rfd_obj1_switch_$TAG.csv" --workers 48 2>&1 | tail -3

echo "WAVE18_PIPELINE_DONE $(date -u +%FT%TZ)"
echo "[$TAG] binders=$npass | novelty-surviving=$ncan | switch rows=$(($(wc -l < results/rfd_obj1_switch_$TAG.csv 2>/dev/null || echo 1)-1))"
echo "[$TAG] NEXT (by hand on the shortlist): pKAI + deep MD quad + cross-react mechanism + rank by priority; the portal is the novelty oracle."
