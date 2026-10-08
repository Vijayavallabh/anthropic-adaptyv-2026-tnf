#!/usr/bin/env python
"""
WAVE 16 (Potts-energy backend) - two-state pH-dependent BINDING energy from Proton-PottsMPNN.

Companion to obj1_potts_phsel.py (decoder/log-prob backend). Here we use the model's Potts
*Hamiltonian* (the paper's headline energy, which predicted experimental recovery) through its
binding-energy thermodynamic cycle:

    g_bind(state) = E_complex(state) - [ E_binder + E_target(state) ]

computed with the TARGET His433 forced into a given protonation state (HIS-P = cationic, the pH-6.5
binding-ON state; HIS-S = neutral, the pH-7.4 OFF state) in BOTH the complex and the isolated-target
pass (E_binder is His433-independent and cancels). The pH switch on the BINDING observable is:

    ddg_bind_pH = g_bind(HIS-P) - g_bind(HIS-S)        (more NEGATIVE = binds better protonated = tumour-ON)

Because the Potts encoder is backbone-only (etab/E_idx do not depend on the sequence), forcing a token
changes only the energy sum, so this is exact and cheap (one featurisation per system). This is the
LEARNED-MODEL analogue of the WAVE-14 MD MM-GBSA two-state binding dG leg -- an orthogonal cross-check
of the same quantity. Signal is sign/consensus, not absolute magnitude (Potts a.u.).

Run with the ProtonPottsMPNN venv (CPU; zero GPU) via scripts/run_potts_phsel.sh --gbind.
Output: results/wave16_potts_gbind.csv
"""
from __future__ import annotations
import os, sys, glob
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")   # CPU-only; zero GPU
os.environ.pop("DEBUG", None)
import numpy as np
import pandas as pd
import torch
from biotite.structure.io.pdb import PDBFile
from atomworks.ml.utils.token import get_token_starts

from mpnn.model.pottsmpnn import PottsMPNN
from mpnn.potts_inference import load_model, prepare_potts_input, _resolve_device

REPO = os.environ.get("ADAPTYV_REPO", "/home/prachh/v/anthropic-adaptyv-2026-egfr")
CKPT = os.environ.get("PPM_CKPT",
    f"{REPO}/third_party/ProtonPottsMPNN/checkpoints/potts_v6_afdb_edge_his0.3_acid0.06/epoch-0125.ckpt")
IN_GLOB = os.environ.get("PPM_INGLOB", f"{REPO}/runs/md_bind_gplus/in/egfr_phsw_*.pdb")
OUT_CSV = f"{REPO}/results/wave16_potts_gbind.csv"
HIS433_ORD = 100
HISP, HISS = 21, 22
DEV = _resolve_device()


def _to_dev(x):
    if torch.is_tensor(x): return x.to(DEV)
    if isinstance(x, dict): return {k: _to_dev(v) for k, v in x.items()}
    return x


def encode(model, structure):
    b = prepare_potts_input(structure, extended_vocab="v6")
    aa = b["atom_array"]
    inp = _to_dev(b["network_input"]["input_features"])
    taa = aa[get_token_starts(aa)]            # one row per token (residue), token order
    etab, E_idx = model.run_potts_encoder(inp)
    return aa, inp, taa, etab, E_idx


def energy(etab, E_idx, S):
    """Native-token Potts energy E(S) = Σ_i cond[i, S_i] (double-counted pair terms cancel in Δ)."""
    cond = PottsMPNN.potts_candidate_energies(etab, E_idx, S.unsqueeze(0))[0]   # [L, V]
    L = int(S.shape[0])
    return float(cond[torch.arange(L, device=cond.device), S].sum().cpu())


def target_chain_of(chain):
    chs = sorted(set([str(c) for c in chain]))
    return max(chs, key=lambda c: int((np.array([str(x) for x in chain]) == c).sum()))


def his433_pos(taa, tc):
    ch = np.array([str(c) for c in taa.chain_id]); rid = np.array([int(r) for r in taa.res_id])
    rn = np.array([str(n) for n in taa.res_name])
    tmask = ch == tc
    off = int(rid[tmask].min()) - 1
    want = off + HIS433_ORD
    hit = np.where(tmask & (rid == want) & (rn == "HIS"))[0]
    if hit.size != 1:
        sys.exit(f"[FATAL] His433 (ordinal {HIS433_ORD}, resid {want}) not uniquely found on chain {tc} (got {hit.size}).")
    return int(hit[0]), want


def main():
    os.makedirs(f"{REPO}/results", exist_ok=True)
    hb = os.environ.get("HBPLUS_PATH")
    if not hb or not os.path.exists(hb):
        sys.exit("[FATAL] HBPLUS_PATH must point at a built hbplus binary.")
    model = load_model(CKPT, extended_vocab="v6")
    print(f"model loaded (v6) on {DEV}; scoring {IN_GLOB}")
    rows = []
    for pdb in sorted(glob.glob(IN_GLOB)):
        name = os.path.splitext(os.path.basename(pdb))[0]
        aa_c, inp_c, taa_c, etab_c, E_idx_c = encode(model, pdb)
        S_c = inp_c["S"][0].clone()
        ch_c = np.array([str(c) for c in taa_c.chain_id])
        rid_c = np.array([int(r) for r in taa_c.res_id])
        tc = target_chain_of(taa_c.chain_id)
        hc, resid433 = his433_pos(taa_c, tc)
        key_S = {(str(c), int(r)): int(s) for c, r, s in zip(taa_c.chain_id, taa_c.res_id, S_c.cpu())}

        # isolated target and binder (re-encoded alone), microstates fixed to the complex assignment
        aa_all = np.array([str(c) for c in aa_c.chain_id])
        aa_t = aa_c[aa_all == tc]; aa_b = aa_c[aa_all != tc]
        _, _, taa_t, etab_t, E_idx_t = encode(model, aa_t)
        _, _, taa_b, etab_b, E_idx_b = encode(model, aa_b)
        ht, _ = his433_pos(taa_t, tc)

        def fix_to_complex(taa):
            S = torch.zeros(len(taa), dtype=S_c.dtype, device=S_c.device)
            for p, (c, r) in enumerate(zip(taa.chain_id, taa.res_id)):
                S[p] = key_S[(str(c), int(r))]
            return S
        S_t = fix_to_complex(taa_t); S_b = fix_to_complex(taa_b)
        E_binder = energy(etab_b, E_idx_b, S_b)           # His433-independent → cancels, but report

        gb = {}
        for state in (HISP, HISS):
            S_c[hc] = state; S_t[ht] = state
            E_complex = energy(etab_c, E_idx_c, S_c)
            E_target = energy(etab_t, E_idx_t, S_t)
            gb[state] = E_complex - (E_binder + E_target)
        ddg = gb[HISP] - gb[HISS]
        rows.append({
            "name": name, "his433_resid": resid433,
            "gbind_hisP": round(gb[HISP], 3), "gbind_hisS": round(gb[HISS], 3),
            "ddg_bind_pH": round(ddg, 3),               # <0 = binds better protonated (tumour-ON)
            "binder_propensity_hisP": round(-gb[HISP], 3),
        })
        print(f"{name}: g_bind(P)={gb[HISP]:+.2f} g_bind(S)={gb[HISS]:+.2f}  ddg_bind_pH={ddg:+.3f} "
              f"{'(ON: binds better protonated)' if ddg < 0 else '(OFF-leaning)'}")

    df = pd.DataFrame(rows).sort_values("name").reset_index(drop=True)
    df.to_csv(OUT_CSV, index=False)
    n = len(df)
    print("\n==== WAVE16 Potts two-state binding pH-switch ====")
    print(f"n designs: {n}")
    print(f"ddg_bind_pH < 0 (binds better protonated = tumour-ON): {int((df.ddg_bind_pH < 0).sum())}/{n}")
    print(f"mean ddg_bind_pH = {df.ddg_bind_pH.mean():+.3f}  median = {df.ddg_bind_pH.median():+.3f}")
    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
