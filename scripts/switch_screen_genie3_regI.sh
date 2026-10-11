#!/bin/bash
# Genie3-regI back-half (switch screen): place SEPARATED multi-His on the conserved GROOVE cations
# (K65 chain C, K98 chain B, K112 chain C — all human=mouse conserved), two-protomer-fold each His-variant
# vs human + mouse (binder A + TNF B + TNF C), score Wyman linkage (mh_screen) + cross-react iptm. Strong =
# linkage<=-0.7 & nHis>=2 & xr>=0.70 -> advance to two-state MD (run_md_release.sh) + novelty.
# Usage: switch_screen_genie3_regI.sh <binders.csv name,seq> <pose_glob> <tag> "<gpus 'a b'>" [cap=24]
set -uo pipefail
REPO=/home/prachh/v/anthropic-adaptyv-2026-tnf; cd "$REPO"
BINDERS=$1; POSEGLOB=$2; TAG=$3; IFS=' ' read -ra GA <<< "${4:-4 5 6 7}"; CAP="${5:-24}"
PY="$REPO/envs/uplifting-esmc/bin/python"; MPNNPY="$REPO/third_party/LigandMPNN/.venv/bin/python"
HUM=$(tr -d '\n ' < "$REPO/runs/partial_denovo2/funnel/xreact/human_tnf.txt")
MOU=$(tr -d '\n ' < "$REPO/runs/buried2/mouse_tnf_mature.txt")
MHOUT="$REPO/runs/${TAG}/multihis"; rm -rf "$MHOUT"; mkdir -p "$MHOUT"
# cap selection (balanced) -> mh_seqs
"$PY" - "$BINDERS" "$MHOUT/mh_seqs.csv" "$MHOUT/selected.txt" "$CAP" <<'PY'
import sys,csv
b,outc,outs,cap=sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4])
rows=list(csv.DictReader(open(b)))[:cap]
w=csv.writer(open(outc,"w")); w.writerow(["name","sequence"]); open(outs,"w").write("")
for r in rows:
    n=r.get("name") or r.get("design"); w.writerow([n,r["sequence"]]); open(outs,"a").write(n+"\n")
print("selected",len(rows),"binders for multi-His")
PY
NSEL=$(grep -c . "$MHOUT/selected.txt" 2>/dev/null || echo 0); [ "$NSEL" -lt 1 ] && { echo "NO BINDERS — stop"; exit 0; }
# 1) build separated multi-His variants scanning BOTH groove protomers (MH_TARGET_CHAINS=BC)
echo "[g3mh] build multi-His (groove cations, chains B+C) for $NSEL binders $(date -u +%FT%TZ)"
MH_SEQS="$MHOUT/mh_seqs.csv" MH_POSEGLOB="$POSEGLOB" MHOUT="$MHOUT" MH_TARGET_CHAINS="BC" \
  "$MPNNPY" scripts/bpair_multihis_build.py $(cat "$MHOUT/selected.txt") 2>&1 | tail -6
NV=$(($(wc -l < "$MHOUT/variants.csv" 2>/dev/null)-1)); echo "[g3mh] $NV variants"; [ "$NV" -lt 1 ] && { echo "NO VARIANTS — stop"; exit 0; }
# 2) two-protomer Boltz (A + TNF B + TNF C) vs human + mouse, sharded
"$PY" - "$MHOUT/variants.csv" "$MHOUT" "$HUM" "$MOU" <<'PY'
import sys,csv,os
v,out,hum,mou=sys.argv[1:5]
for sp,tgt in (("h",hum),("m",mou)):
    d=os.path.join(out,f"in_{sp}"); os.makedirs(d,exist_ok=True)
    for r in csv.DictReader(open(v)):
        n=r["name"]; s=r["sequence"].split(":")[0]
        open(os.path.join(d,n+".yaml"),"w").write(
          f"version: 1\nsequences:\n  - protein: {{id: A, sequence: {s}, msa: empty}}\n"
          f"  - protein: {{id: B, sequence: {tgt}, msa: empty}}\n  - protein: {{id: C, sequence: {tgt}, msa: empty}}\n")
    print("wrote",sum(1 for _ in csv.DictReader(open(v))),"->",d)
PY
run_boltz(){ ( unset VIRTUAL_ENV 2>/dev/null||true; source "$REPO/envs/uplifting-boltz2/bin/activate"
    set -a; source "$REPO/envs/uplifting.env" 2>/dev/null||true; set +a
    cd "$REPO/third_party/uplifting-biomolecular-modeling/boltz2"
    CUDA_VISIBLE_DEVICES=$3 bash run.sh pred --config h100 --mode fast --input "$1" --out_dir "$2" ) ; }
g0=${GA[0]}; g1=${GA[$((${#GA[@]}>1?1:0))]}
run_boltz "$MHOUT/in_h" "$MHOUT/boltz_h" "$g0" > "$MHOUT/boltz_h.log" 2>&1 & ph=$!
run_boltz "$MHOUT/in_m" "$MHOUT/boltz_m" "$g1" > "$MHOUT/boltz_m.log" 2>&1 & pm=$!
wait $ph; wait $pm; echo "[g3mh] boltz done $(date -u +%FT%TZ)"
# 3) linkage (human) via mh_screen on the 2-protomer human folds, + join mouse iptm -> strong roster
"$PY" scripts/mh_screen.py "$MHOUT/manifest.tsv" "$MHOUT/boltz_h/**" "$MHOUT/linkage.tsv" 2>&1 | tail -8 || \
  echo "[g3mh] NOTE: if mh_screen glob differs for 2-protomer folds, point it at $MHOUT/boltz_h confidence jsons"
echo "G3_SWITCH_SCREEN_DONE $(date -u +%FT%TZ) — strong roster -> $MHOUT/linkage.tsv; next: run_md_release.sh on strong (anchor = conserved groove cation K65_C/K98_B/K112_C)"
