#!/usr/bin/env python
"""
WAVE 16 - Proton-PottsMPNN protonation-state corroboration of the His433 pH switch (objective 1).

A 5th, method-orthogonal pH-selectivity signal for the de novo conditional EGFR binders, using the
protonation-aware Potts inverse-folding model *Proton-PottsMPNN* (Jacobsen et al., 2026;
bioRxiv 10.64898/2026.09.30.755438; code github.com/christian-creator/ProtonPottsMPNN, MIT). The
authors benchmark the model's protonation energy against experimental pKa data (PKAD), so reading a
titratable site's preferred state on a structure the model did NOT design is a supported use.

For each binder:DIII complex we read the model's per-residue protonation-state preference at the
TARGET His433 (the switch) via the SAME conditional-minus-self decode the design engine uses for its
selective field:
    gap(His433) = logP(HIS-P) - logP(HIS-S)        (+ = protonated / binding-ON state favoured)
The ABSOLUTE gap is only a relative preference (the authors caution the state tokens are relative,
not absolute; and no structure-based predictor beats a pKa-6.5 null for a bare histidine), so the
SIGNAL we trust is the binder's CAUSAL contribution, isolated against the apo target:
    delta(His433) = gap_complex(His433) - gap_apo(His433)
delta > 0  <=>  the binder's interface shifts His433 toward the protonated (pH-ON) state -- i.e.
raises its effective pKa toward the tumour window -- an independent, learned-model corroboration of
the engineered switch, orthogonal to PROPKA / pKAI / MD (sign-consensus, not magnitude).

CONTROL: the SAME delta is computed for every other target histidine (5 per DIII crop) that the
binder does NOT engage; His433's shift is ranked against them, so a design that specifically engages
His433 shows delta(His433) above its own non-engaged-His background.

Run with the ProtonPottsMPNN venv (CPU; zero GPU) via scripts/run_potts_phsel.sh.
Outputs: results/wave16_potts_phsel.csv  +  results/wave16_potts_phsel_perhis.json
"""
from __future__ import annotations
import os, sys, json, glob
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")   # CPU-only; zero GPU use
os.environ.pop("DEBUG", None)
import numpy as np
import pandas as pd
from biotite.structure.io.pdb import PDBFile

from mpnn.potts_inference import load_model, prepare_potts_input, compute_log_probs

REPO = os.environ.get("ADAPTYV_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CKPT = os.environ.get("PPM_CKPT",
    f"{REPO}/third_party/ProtonPottsMPNN/checkpoints/potts_v6_afdb_edge_his0.3_acid0.06/epoch-0125.ckpt")
IN_GLOB = os.environ.get("PPM_INGLOB", f"{REPO}/runs/md_bind_gplus/in/egfr_phsw_*.pdb")
SCRATCH = os.environ.get("PPM_SCRATCH", os.path.join(REPO, "runs", "ppm_scratch"))
OUT_CSV = f"{REPO}/results/wave16_potts_phsel.csv"
OUT_JSON = f"{REPO}/results/wave16_potts_phsel_perhis.json"
HIS433_ORD = 100
# v6 token indices (verified from mpnn.transforms.extended_vocab.get_vocab('v6'))
HISP, HISS = 21, 22
ASPP, ASPD, GLUP, GLUD = 24, 25, 27, 28


def target_chain_of(chain: np.ndarray) -> str:
    chs = sorted(set(chain.tolist()))
    return max(chs, key=lambda c: int((chain == c).sum()))   # DIII target (191 res) >> binder


def score_structure(model, structure):
    """Return {his_ordinal: gap(HIS-P - HIS-S)} for every target-chain His, plus binder acid gaps."""
    batch = prepare_potts_input(structure, extended_vocab="v6")
    aa = batch["atom_array"]
    ca = aa[aa.atom_name == "CA"]
    chain = np.array(ca.chain_id); resid = np.array(ca.res_id); rn = np.array(ca.res_name)
    lp = compute_log_probs(model, batch, causality="conditional_minus_self").cpu().numpy()  # [L, V]
    tc = target_chain_of(chain)
    tmask = chain == tc
    boff = int(resid[tmask].min()) - 1
    his = {int(resid[i]) - boff: float(lp[i, HISP] - lp[i, HISS])
           for i in np.where(tmask & (rn == "HIS"))[0]}
    acids = []
    for i in np.where((~tmask) & np.isin(rn, ["ASP", "GLU"]))[0]:
        acids.append(float(lp[i, ASPD] - lp[i, ASPP]) if rn[i] == "ASP"
                     else float(lp[i, GLUD] - lp[i, GLUP]))
    return his, boff, tc, acids


def main():
    os.makedirs(f"{REPO}/results", exist_ok=True)
    os.makedirs(SCRATCH, exist_ok=True)
    hb = os.environ.get("HBPLUS_PATH")
    if not hb or not os.path.exists(hb):
        sys.exit("[FATAL] HBPLUS_PATH must point at a built hbplus binary (needed to featurize).")
    model = load_model(CKPT, extended_vocab="v6")
    print(f"model loaded (v6); scoring {IN_GLOB}")
    rows, perhis = [], {}
    for pdb in sorted(glob.glob(IN_GLOB)):
        name = os.path.splitext(os.path.basename(pdb))[0]
        # complex
        his_c, boff, tc, acids = score_structure(model, pdb)
        # apo target: write the target chain alone to a temp PDB (identical code path as the complex)
        aa_raw = PDBFile.read(pdb).get_structure(model=1)
        tc_raw = target_chain_of(np.array(aa_raw[aa_raw.atom_name == "CA"].chain_id))
        apo_arr = aa_raw[aa_raw.chain_id == tc_raw]
        apo_path = os.path.join(SCRATCH, f"{name}_apo.pdb")
        pf = PDBFile(); pf.set_structure(apo_arr); pf.write(apo_path)
        his_a, _, _, _ = score_structure(model, apo_path)

        if HIS433_ORD not in his_c or HIS433_ORD not in his_a:
            sys.exit(f"[FATAL] {name}: His433 ordinal {HIS433_ORD} missing "
                     f"(complex={HIS433_ORD in his_c}, apo={HIS433_ORD in his_a}).")
        deltas = {o: his_c[o] - his_a[o] for o in his_c if o in his_a}
        d433 = deltas[HIS433_ORD]
        ctrl = np.array([d for o, d in deltas.items() if o != HIS433_ORD], float)
        ranked = sorted(deltas.items(), key=lambda kv: -kv[1])        # most-protonation-shifted first
        rank = [o for o, _ in ranked].index(HIS433_ORD) + 1
        row = {
            "name": name,
            "his433_ord": HIS433_ORD,
            "his433_resid": boff + HIS433_ORD,
            "gap_complex_his433": round(his_c[HIS433_ORD], 3),
            "gap_apo_his433": round(his_a[HIS433_ORD], 3),
            "delta_his433": round(d433, 3),
            "n_target_his": len(deltas),
            "ctrl_delta_mean": round(float(ctrl.mean()), 3) if ctrl.size else float("nan"),
            "ctrl_delta_max": round(float(ctrl.max()), 3) if ctrl.size else float("nan"),
            "his433_minus_ctrlmean": round(float(d433 - ctrl.mean()), 3) if ctrl.size else float("nan"),
            "his433_rank_by_delta": rank,
            "binder_mean_acid_DmP": round(float(np.mean(acids)), 2) if acids else float("nan"),
            "binder_min_acid_DmP": round(float(np.min(acids)), 2) if acids else float("nan"),
        }
        rows.append(row)
        perhis[name] = {"complex": his_c, "apo": his_a, "delta": deltas}
        print(f"{name}: dHis433={row['delta_his433']:+.3f} "
              f"(cplx {row['gap_complex_his433']:+.2f} - apo {row['gap_apo_his433']:+.2f}) "
              f"ctrlDmean={row['ctrl_delta_mean']:+.3f} rank {rank}/{row['n_target_his']} "
              f"his433-ctrl={row['his433_minus_ctrlmean']:+.3f}")

    df = pd.DataFrame(rows).sort_values("name").reset_index(drop=True)
    df.to_csv(OUT_CSV, index=False)
    with open(OUT_JSON, "w") as fh:
        json.dump(perhis, fh, indent=2)
    n = len(df)
    print("\n==== WAVE16 Proton-PottsMPNN His433 pH-switch corroboration ====")
    print(f"n designs: {n}")
    print(f"delta_his433 > 0 (binder shifts His433 toward protonated/ON): "
          f"{int((df.delta_his433 > 0).sum())}/{n}")
    print(f"his433_minus_ctrlmean > 0 (His433 shift exceeds non-engaged-His background): "
          f"{int((df.his433_minus_ctrlmean > 0).sum())}/{n}")
    print(f"His433 is the MOST protonation-shifted target His (rank 1): "
          f"{int((df.his433_rank_by_delta == 1).sum())}/{n}")
    print(f"mean delta_his433 = {df.delta_his433.mean():+.3f}  median = {df.delta_his433.median():+.3f}")
    print(f"mean his433_minus_ctrlmean = {df.his433_minus_ctrlmean.mean():+.3f}")
    print(f"wrote {OUT_CSV}  and  {OUT_JSON}")


if __name__ == "__main__":
    main()
