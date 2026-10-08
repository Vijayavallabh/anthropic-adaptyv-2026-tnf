#!/usr/bin/env python3
"""WAVE 19 co-fold: build the chai/boltz input for the NEW MD-selective hits, and parse the outputs
into results/cofold_consensus_w19.csv (a superset of the w18 consensus). The wave19 finalize score uses
bind_consensus = min(chai_iptm_human, boltz2_iptm), so a new hit with no co-fold would be unfairly
penalized against the existing 11; this gives the new hits their third-predictor binding confirmation.

Modes:
  --build-input : read results/md_selective_w19.csv, drop the 11 already in md_selective_w18, pull each
                  new hit's sequence from results/metrics_w18_full.csv, write results/cofold_in_w19.csv.
  --parse       : read chai (runs/chai_w19_human, optional runs/chai_w19_mouse) + boltz (runs/boltz_w19)
                  outputs, merge with results/cofold_consensus_w18.csv, write cofold_consensus_w19.csv.

Between the two, run (GPUs free, after the MD screen):
  scripts/run_chai.sh results/cofold_in_w19.csv w19_human data/targets/human_domainIII.seq 0-7
  scripts/run_boltz.sh results/cofold_in_w19.csv w19 <gpu>
(chai_mouse is optional; it is only reported, not scored, since cross-reactivity is gated by the af2ig
mouse-vs-human delta upstream.)
"""
import argparse, csv, glob, json, os, re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(REPO)


def build_input():
    existing = {r["name"] for r in csv.DictReader(open("results/md_selective_w18.csv"))}
    seqs = {r["name"]: (r.get("sequence") or "") for r in csv.DictReader(open("results/metrics_w18_full.csv"))}
    rows, missing = [], []
    for r in csv.DictReader(open("results/md_selective_w19.csv")):
        n = r["name"]
        if n in existing:
            continue
        s = seqs.get(n, "")
        if s:
            rows.append({"name": n, "sequence": s, "molecule_class": "single_chain"})
        else:
            missing.append(n)
    with open("results/cofold_in_w19.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "sequence", "molecule_class"]); w.writeheader(); w.writerows(rows)
    print(f"wrote results/cofold_in_w19.csv with {len(rows)} new hits for co-fold"
          + (f" | NO SEQUENCE for {len(missing)}: {missing}" if missing else ""))


def chai_iptm(tag):
    """design -> best iptm across chai models. Path: runs/chai_<tag>/out/fast/<name>/seed_*/scores.model_idx_*.npz"""
    try:
        import numpy as np
    except Exception:
        print(f"WARN: numpy unavailable, cannot parse chai_{tag}"); return {}
    out = {}
    files = glob.glob(f"runs/chai_{tag}/out/fast/*/seed_*/scores.model_idx_*.npz")
    for npz in files:
        try:
            name = npz.split("/out/fast/")[1].split("/")[0]
            v = float(np.asarray(np.load(npz)["iptm"]).ravel()[0])   # chai iptm is a (1,) array, not a scalar
            out[name] = max(out.get(name, 0.0), round(v, 3))
        except Exception as e:
            print(f"WARN: chai parse failed for {npz}: {e!r}")
    if files and not out:
        print(f"WARN: chai_{tag}: {len(files)} npz found but 0 parsed -- check the iptm key/shape")
    return out


def boltz_iptm(tag):
    """design -> best iptm. Filename: confidence_<name>_model_<k>.json (robust to the output tree shape)."""
    out = {}
    for j in glob.glob(f"runs/boltz_{tag}/out/**/confidence_*.json", recursive=True):
        try:
            m = re.match(r"confidence_(.+)_model_\d+\.json", os.path.basename(j))
            if not m:
                continue
            name = m.group(1)
            v = float(json.load(open(j))["iptm"])
            out[name] = max(out.get(name, 0.0), round(v, 3))
        except Exception:
            pass
    return out


def parse():
    ch, cm, bz = chai_iptm("w19_human"), chai_iptm("w19_mouse"), boltz_iptm("w19")
    cols = ["name", "chai_iptm_human", "chai_mouse_iptm", "boltz2_iptm"]
    merged = {}
    if os.path.exists("results/cofold_consensus_w18.csv"):
        for r in csv.DictReader(open("results/cofold_consensus_w18.csv")):
            merged[r["name"]] = {k: r.get(k, "") for k in cols}
    new = set(ch) | set(cm) | set(bz)
    for n in new:
        merged[n] = {"name": n, "chai_iptm_human": ch.get(n, ""),
                     "chai_mouse_iptm": cm.get(n, ""), "boltz2_iptm": bz.get(n, "")}
    with open("results/cofold_consensus_w19.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(merged.values())
    print(f"wrote results/cofold_consensus_w19.csv ({len(merged)} designs; new co-folded: {len(new)}; "
          f"chai_human={len(ch)} chai_mouse={len(cm)} boltz={len(bz)})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-input", action="store_true")
    ap.add_argument("--parse", action="store_true")
    a = ap.parse_args()
    if a.build_input:
        build_input()
    if a.parse:
        parse()
    if not (a.build_input or a.parse):
        ap.error("pass --build-input and/or --parse")
