#!/usr/bin/env python3
"""Stage the validation-completion pool for the final-20: canonicalize each under-validated
PD/RFd pick's af2pred complex and split into the subsets that still need each validation.

  xreact/boltz set = final picks missing cross_reactive/boltz2
  MD set           = final picks missing md_docc

Outputs:
  runs/cand_final2_canon/<design>.pdb   (all under-validated picks, canon binder=A/target=B)
  runs/cand_xr2/<design>.pdb            (cross-react+boltz subset)
  runs/md_final2_list.txt               (MD subset: design \t canon_pdb, one per line)
  runs/boltz_final2.csv                 (name,sequence for the xr/boltz subset)
"""
import csv, os, shutil, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CANON = REPO/"scripts"/"canon_renumber.py"
PY = sys.executable
D3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H",
      "ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
      "TYR":"Y","VAL":"V","MSE":"M"}


def chainA_seq(pdb):
    return "".join(D3.get(l[17:20].strip(), "X") for l in open(pdb)
                   if l[:4] == "ATOM" and l[21] == "A" and l[12:16].strip() == "CA")


def main():
    rows = [r for r in csv.DictReader(open(REPO/"results"/"submission_manifest_v2.csv"))
            if r["source"] != "baseline" and int(r["rank"]) >= 9]
    canon_dir = REPO/"runs"/"cand_final2_canon"; canon_dir.mkdir(parents=True, exist_ok=True)
    xr_dir = REPO/"runs"/"cand_xr2"; xr_dir.mkdir(parents=True, exist_ok=True)
    md_list, boltz_rows = [], []
    PASSER_DIRS = ["rfd_passers_pd", "rfd_passers_pd2", "rfd_passers_pd3",
                   "rfd_passers_w1fa", "rfd_passers_w2fa"]
    for r in rows:
        d = r["design"]
        raw = next((REPO/"runs"/pd/f"{d}.pdb" for pd in PASSER_DIRS
                    if (REPO/"runs"/pd/f"{d}.pdb").exists()), None)
        if raw is None:
            print(f"WARN: no complex found for {d}; skipping"); continue
        canon_pdb = canon_dir/f"{d}.pdb"
        subprocess.run([PY, str(CANON), str(raw), str(canon_pdb)], check=True, capture_output=True)
        # cross-react + Boltz: stage ALL non-baseline picks and re-run them together, so the shared
        # af2ig_xreact_xr2 / boltz_fin2 dirs hold the CURRENT top-20's results (staging only the
        # "missing" ones would overwrite and drop previously-validated picks that are still selected).
        shutil.copy(canon_pdb, xr_dir/f"{d}.pdb")
        boltz_rows.append((d, chainA_seq(canon_pdb)))
        if not str(r["md_docc"]).strip():   # MD is additive (campaign skips already-scored designs)
            md_list.append((d, str(canon_pdb)))
    (REPO/"runs"/"md_final2_list.txt").write_text(
        "".join(f"{d}\t{p}\n" for d, p in md_list))
    with open(REPO/"runs"/"boltz_final2.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["name", "sequence"])
        for d, s in boltz_rows: w.writerow([d, s])
    print(f"canon pool: {len(rows)} | xreact/boltz: {len(boltz_rows)} | MD: {len(md_list)}")
    print("MD designs:", [d for d, _ in md_list])
    print("xr/boltz designs:", [d for d, _ in boltz_rows])


if __name__ == "__main__":
    main()
