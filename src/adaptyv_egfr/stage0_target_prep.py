"""Stage 0 - target & epitope preparation.

Builds the epitope "contract" every downstream stage targets:
  * load human (P00533) and mouse (Q01279) EGFR; align; compute ECD identity
  * verify + save the exact challenge target region (residues 25-645)
  * check the conserved cross-reactive anchors + cetuximab-divergent positions
  * extract the structural epitope from 6ARU (EGFR chain A contacts to the cetuximab Fab)
    and list candidate Asp/Glu partners for protonated-His salt bridges (the pH switch)
  * write data/epitope/epitope_map.json and data/epitope/structure_epitope.json

Residue numbering: FULL P00533 precursor numbering (signal peptide = residues 1-24).

Run:  uv run python -m adaptyv_egfr.stage0_target_prep
"""
from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TARGETS = REPO / "data" / "targets"
EPITOPE = REPO / "data" / "epitope"

HUMAN_FASTA = TARGETS / "EGFR_human_P00533.fasta"
MOUSE_FASTA = TARGETS / "EGFR_mouse_Q01279.fasta"
TARGET_PDB = TARGETS / "6ARU.pdb"
TARGET_FASTA = TARGETS / "EGFR_human_target_25-645.fasta"

# Human EGFR extracellular target region (full P00533 numbering).
ECD_RANGE = (25, 645)
# Exact target sequence quoted on the challenge page (used as a correctness check).
TARGET_SEQ_EXPECTED = (
    "LEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQ"
    "IIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGS"
    "CQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPL"
    "MLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDS"
    "LSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGR"
    "TKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQVCHALC"
    "SPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHC"
    "VKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPS"
)

# Domain III (L2) approximate range in full numbering - refined by the structural epitope below.
DOMAIN_III_APPROX = (313, 514)
# Conserved (human == mouse) anchors on the Domain III ligand-blocking face to TARGET.
CONSERVED_ANCHORS = [408, 432, 433, 436, 464, 489]
# Cetuximab contacts that DIFFER human vs mouse (why cetuximab misses mouse) - AVOID.
DIVERGENT_AVOID = [377, 442, 467, 491, 492, 497]
# His433 (full numbering) = built-in, mouse-conserved pH sensor for objective 1.
PH_SENSOR_HIS = 433

# 6ARU: chain A = EGFR ECD, chains B/C = cetuximab Fab (light/heavy).
EGFR_CHAIN = "A"
FAB_CHAINS = ("B", "C")
EPITOPE_CUTOFF = 5.0  # angstroms (any-atom contact)


def read_fasta(path: Path) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))


def _aligner(end_gap_free: bool = False):
    from Bio.Align import PairwiseAligner, substitution_matrices

    a = PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11
    a.extend_gap_score = -1
    if end_gap_free:  # place a fragment into a full sequence without end-gap penalties
        try:
            a.end_gap_score = 0.0  # convenience setter for both end gaps (current Biopython)
        except Exception:  # pragma: no cover - older Biopython attribute names
            a.target_end_gap_score = 0.0
            a.query_end_gap_score = 0.0
    return a


def align(human: str, mouse: str):
    return _aligner().align(human, mouse)[0]


def human_to_mouse_map(alignment) -> dict[int, tuple[str, str]]:
    """Map human 1-based position -> (human_aa, mouse_aa or '-')."""
    idx = alignment.indices  # (2, L); -1 marks a gap
    human_seq, mouse_seq = alignment.sequences[0], alignment.sequences[1]
    out: dict[int, tuple[str, str]] = {}
    for col in range(idx.shape[1]):
        hi, mi = int(idx[0, col]), int(idx[1, col])
        if hi < 0:
            continue
        out[hi + 1] = (human_seq[hi], mouse_seq[mi] if mi >= 0 else "-")
    return out


def build_epitope_map(h2m: dict[int, tuple[str, str]]) -> dict:
    lo, hi = ECD_RANGE
    pairs = [h2m[p] for p in range(lo, hi + 1) if p in h2m]
    identity = round(100.0 * sum(a == b for a, b in pairs) / max(1, len(pairs)), 2)

    def annotate(positions):
        return [
            {"pos": p, "human": h2m.get(p, ("?", "?"))[0], "mouse": h2m.get(p, ("?", "?"))[1],
             "conserved": h2m.get(p, ("?", "?"))[0] == h2m.get(p, ("?", "?"))[1]}
            for p in positions
        ]

    his = h2m.get(PH_SENSOR_HIS, ("?", "?"))
    return {
        "target": "human EGFR extracellular region (UniProt P00533-1)",
        "numbering": "full P00533 precursor numbering (signal peptide 1-24)",
        "reference_structure": "PDB 6ARU chain A (EGFR + cetuximab Fab)",
        "ecd_range": list(ECD_RANGE),
        "domain_III_approx": list(DOMAIN_III_APPROX),
        "ecd_identity_human_vs_mouse_pct": identity,
        "conserved_anchors": annotate(CONSERVED_ANCHORS),
        "divergent_avoid": annotate(DIVERGENT_AVOID),
        "ph_sensor": {
            "pos": PH_SENSOR_HIS, "human": his[0], "mouse": his[1],
            "note": "protonated His imidazole (pKa ~6.3) => ON @ 6.5 / OFF @ 7.4; conserved => serves objectives 1 & 2",
        },
        "objective_priority": ["ph_selective", "cross_reactive", "human_affinity"],
    }


def slice_target(human: str) -> tuple[str, bool]:
    """Write the exact 25-645 target region and verify it matches the challenge page."""
    lo, hi = ECD_RANGE
    s = human[lo - 1:hi]
    ok = s == TARGET_SEQ_EXPECTED
    body = "\n".join(s[i:i + 60] for i in range(0, len(s), 60))
    TARGET_FASTA.write_text(
        f">EGFR_human_P00533_25-645 target extracellular region ({len(s)} aa)\n{body}\n"
    )
    return s, ok


def _three_to_one() -> dict[str, str]:
    from Bio.Data.IUPACData import protein_letters_3to1_extended

    return {k.upper(): v for k, v in protein_letters_3to1_extended.items()}


def build_structure_epitope(h2m: dict[int, tuple[str, str]]) -> dict:
    """EGFR (chain A) residues contacting the cetuximab Fab in 6ARU, mapped to P00533."""
    from Bio.PDB import NeighborSearch, PDBParser

    t1 = _three_to_one()
    model = PDBParser(QUIET=True).get_structure("6aru", str(TARGET_PDB))[0]

    # observed chain A sequence + residue ids (standard amino acids only)
    obs_seq, obs_resids = [], []
    for res in model[EGFR_CHAIN]:
        if res.id[0] != " ":
            continue
        aa = t1.get(res.resname.upper())
        if aa and len(aa) == 1:
            obs_seq.append(aa)
            obs_resids.append(res.id[1])
    obs_seq = "".join(obs_seq)

    # map 6ARU chain A numbering -> P00533 full numbering (end-gap-free fragment alignment)
    human = read_fasta(HUMAN_FASTA)
    aln = _aligner(end_gap_free=True).align(obs_seq, human)[0]
    idx = aln.indices
    resid_to_p00533: dict[int, int] = {}
    for col in range(idx.shape[1]):
        oi, hi = int(idx[0, col]), int(idx[1, col])
        if oi >= 0 and hi >= 0:
            resid_to_p00533[obs_resids[oi]] = hi + 1

    # contacts to the Fab
    fab_atoms = [a for ch in FAB_CHAINS if ch in model for a in model[ch].get_atoms()]
    ns = NeighborSearch(fab_atoms)
    epitope = []
    for res in model[EGFR_CHAIN]:
        if res.id[0] != " ":
            continue
        aa = t1.get(res.resname.upper())
        if not aa or len(aa) != 1:
            continue
        if any(ns.search(atom.coord, EPITOPE_CUTOFF) for atom in res):
            rid = res.id[1]
            pos = resid_to_p00533.get(rid)
            mouse_aa = h2m.get(pos, ("?", "?"))[1] if pos else "?"
            epitope.append({
                "resid_6aru": rid, "aa": aa, "p00533_pos": pos,
                "mouse_aa": mouse_aa, "conserved": (pos is not None and aa == mouse_aa),
            })

    epi_pos = [e["p00533_pos"] for e in epitope if e["p00533_pos"]]
    conserved = [e for e in epitope if e["conserved"]]
    # Asp/Glu in the epitope = candidate protonated-His salt-bridge partners for the pH switch.
    acidic_partners = [e for e in epitope if e["aa"] in ("D", "E")]
    return {
        "reference_structure": "6ARU",
        "egfr_chain": EGFR_CHAIN,
        "fab_chains": list(FAB_CHAINS),
        "contact_cutoff_angstrom": EPITOPE_CUTOFF,
        "numbering": "p00533_pos = full P00533; resid_6aru = 6ARU chain A numbering",
        "n_epitope_residues": len(epitope),
        "epitope_p00533_positions": sorted(p for p in epi_pos),
        "epitope": epitope,
        "conserved_epitope_p00533": sorted(e["p00533_pos"] for e in conserved if e["p00533_pos"]),
        "acidic_his_partner_candidates": [
            {"p00533_pos": e["p00533_pos"], "aa": e["aa"], "conserved": e["conserved"]}
            for e in acidic_partners
        ],
    }


def main() -> None:
    EPITOPE.mkdir(parents=True, exist_ok=True)
    human, mouse = read_fasta(HUMAN_FASTA), read_fasta(MOUSE_FASTA)
    h2m = human_to_mouse_map(align(human, mouse))

    m = build_epitope_map(h2m)
    (EPITOPE / "epitope_map.json").write_text(json.dumps(m, indent=2))
    print(f"wrote {EPITOPE / 'epitope_map.json'}")
    print(f"ECD identity human vs mouse (res {ECD_RANGE[0]}-{ECD_RANGE[1]}): "
          f"{m['ecd_identity_human_vs_mouse_pct']}%")
    print("conserved anchors:")
    for a in m["conserved_anchors"]:
        print(f"  {a['human']}{a['pos']}: human={a['human']} mouse={a['mouse']}  "
              f"[{'OK' if a['conserved'] else 'DIVERGENT!'}]")
    print("cetuximab-divergent (avoid):")
    for a in m["divergent_avoid"]:
        print(f"  pos {a['pos']}: human={a['human']} mouse={a['mouse']}  conserved={a['conserved']}")
    hs = m["ph_sensor"]
    print(f"pH sensor His{hs['pos']}: human={hs['human']} mouse={hs['mouse']}")

    seq, ok = slice_target(human)
    print(f"\ntarget region {ECD_RANGE[0]}-{ECD_RANGE[1]}: {len(seq)} aa; "
          f"matches challenge page: {'YES' if ok else 'NO (!!)'}  -> {TARGET_FASTA.name}")

    try:
        se = build_structure_epitope(h2m)
        (EPITOPE / "structure_epitope.json").write_text(json.dumps(se, indent=2))
        print(f"\nstructural epitope (6ARU, EGFR within {EPITOPE_CUTOFF} A of cetuximab Fab): "
              f"{se['n_epitope_residues']} residues")
        print(f"  epitope P00533 positions: {se['epitope_p00533_positions']}")
        print(f"  conserved epitope positions (human==mouse): {se['conserved_epitope_p00533']}")
        print("  acidic His-partner candidates (Asp/Glu in epitope):")
        for a in se["acidic_his_partner_candidates"]:
            print(f"    {a['aa']}{a['p00533_pos']}  conserved={a['conserved']}")
    except Exception as e:  # best-effort; do not fail the verification gate
        print(f"\nWARN: structural epitope extraction skipped: {e}")


if __name__ == "__main__":
    main()
