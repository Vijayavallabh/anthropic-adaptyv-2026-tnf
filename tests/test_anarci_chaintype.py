"""Antibody-format tests: ANARCI numbering + kappa/lambda isotype match.

These exercise the hard filter for nanobody/scFv/Fab classes via anarcism.
Skipped automatically if anarcism is unavailable.
"""
import pytest

from adaptyv_egfr.submission import anarci_available, validate_rows

pytestmark = pytest.mark.skipif(not anarci_available(), reason="anarcism not installed")

VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
      "ADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
VL_KAPPA = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDF"
            "TLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
VL_LAMBDA = ("QSVLTQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNSNRPSGVPDRFSGSKSGT"
             "SASLAITGLQAEDEADYYCQSYDSSLSGSVFGGGTKLTVL")
LINKER = "GGGGSGGGGSGGGGS"


def _one(name, seq, mclass):
    return validate_rows([{"name": name, "sequence": seq, "molecule_class": mclass}], track=3)


def test_fab_kappa_ok():
    rep = _one("f", f"{VH}:{VL_KAPPA}", "fab_kappa")
    assert rep.ok, [e for r in rep.rows for e in r.errors]


def test_fab_lambda_ok():
    rep = _one("f", f"{VH}:{VL_LAMBDA}", "fab_lambda")
    assert rep.ok, [e for r in rep.rows for e in r.errors]


def test_fab_kappa_with_lambda_vl_fails():
    rep = _one("f", f"{VH}:{VL_LAMBDA}", "fab_kappa")
    assert not rep.ok
    assert any("isotype" in e for r in rep.rows for e in r.errors)


def test_nanobody_single_heavy_ok():
    rep = _one("n", VH, "nanobody")
    assert rep.ok, [e for r in rep.rows for e in r.errors]


def test_scfv_ok():
    rep = _one("s", VH + LINKER + VL_KAPPA, "scfv")
    assert rep.ok, [e for r in rep.rows for e in r.errors]


def test_nanobody_rejects_paired_chains():
    # A scFv-shaped sequence declared as nanobody must fail (two domains found).
    rep = _one("n", VH + LINKER + VL_KAPPA, "nanobody")
    assert not rep.ok
