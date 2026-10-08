"""Length limits and scoring sub-category tests."""
from adaptyv_tnf.submission import length_category, validate_rows


def test_length_categories():
    assert length_category(39) == "protein_microbinder"
    assert length_category(40) == "protein_minibinder"
    assert length_category(100) == "protein_minibinder"
    assert length_category(101) == "large_protein"


def test_too_short_fails():
    rep = validate_rows([{"name": "d", "sequence": "AAAAA", "molecule_class": "protein"}])  # 5 aa
    assert not rep.ok


def test_too_long_fails():
    rep = validate_rows([{"name": "d", "sequence": "A" * 251, "molecule_class": "protein"}])
    assert not rep.ok


def test_length_boundaries_pass():
    for n in (10, 250):
        rep = validate_rows([{"name": "d", "sequence": "A" * n, "molecule_class": "protein"}], track=3)
        assert rep.ok, (n, [e for r in rep.rows for e in r.errors])
