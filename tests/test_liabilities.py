"""Tests for the CDR/sequence liability scanner (stage 5)."""
from adaptyv_tnf.stage5_validate import liability_count, scan_liabilities


def test_n_glycosylation_sequon():
    assert scan_liabilities("AANGTAA")["n_glycosylation"] == [3]   # N-G-T matches N[^P][ST]
    assert scan_liabilities("AANPTAA")["n_glycosylation"] == []    # N-P-T excluded (X != P)


def test_unpaired_cysteine():
    assert scan_liabilities("ACA")["unpaired_cys"] is True
    assert scan_liabilities("ACAC")["unpaired_cys"] is False


def test_deamidation_and_isomerization():
    assert 2 in scan_liabilities("ANGA")["deamidation_NG_NS"]      # N-G at pos 2
    assert 2 in scan_liabilities("ADGA")["isomerization_DG_DS_DT"]  # D-G at pos 2


def test_liability_count_nonzero():
    assert liability_count("ANGTCAA") >= 1
