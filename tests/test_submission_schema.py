"""Schema/contract tests for the submission validator."""
from adaptyv_egfr.submission import validate_csv, validate_rows, write_submission

GOOD_PROTEIN = "M" + "A" * 49  # 50 aa, all standard residues


def test_valid_protein_passes():
    rep = validate_rows(
        [{"name": "d1", "sequence": GOOD_PROTEIN, "molecule_class": "protein"}], track=3
    )
    assert rep.ok, rep.errors + [e for r in rep.rows for e in r.errors]


def test_missing_column_fails(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("name,sequence\nd1," + GOOD_PROTEIN + "\n")
    rep = validate_csv(p, track=3)
    assert not rep.ok
    assert any("missing required column" in e for e in rep.errors)


def test_bad_class_fails():
    rep = validate_rows([{"name": "d1", "sequence": GOOD_PROTEIN, "molecule_class": "peptide"}])
    assert not rep.ok


def test_nonstandard_residue_fails():
    rep = validate_rows([{"name": "d1", "sequence": "AAAXAA" + "A" * 10, "molecule_class": "protein"}])
    assert not rep.ok


def test_colon_in_non_fab_fails():
    rep = validate_rows([{"name": "d1", "sequence": "AAAAA:AAAAA", "molecule_class": "protein"}])
    assert not rep.ok


def test_duplicate_name_fails():
    rep = validate_rows(
        [
            {"name": "d", "sequence": GOOD_PROTEIN, "molecule_class": "protein"},
            {"name": "d", "sequence": "C" + GOOD_PROTEIN, "molecule_class": "protein"},
        ]
    )
    assert not rep.ok
    assert any("duplicate name" in e for e in rep.errors)


def test_duplicate_sequence_fails():
    rep = validate_rows(
        [
            {"name": "d1", "sequence": GOOD_PROTEIN, "molecule_class": "protein"},
            {"name": "d2", "sequence": GOOD_PROTEIN, "molecule_class": "protein"},
        ]
    )
    assert not rep.ok
    assert any("duplicate sequence" in e for e in rep.errors)


def test_track_cap_enforced():
    rows = [{"name": f"d{i}", "sequence": "A" * (10 + i), "molecule_class": "protein"} for i in range(21)]
    rep = validate_rows(rows, track=3)
    assert not rep.ok
    assert any("exceeds Track 3 cap" in e for e in rep.errors)


def test_write_submission_roundtrip(tmp_path):
    rows = [{"name": "d1", "sequence": GOOD_PROTEIN, "molecule_class": "protein", "rank": 1}]
    p = write_submission(rows, tmp_path / "s.csv", extra_columns=["rank"])
    rep = validate_csv(p, track=3)
    assert rep.ok, rep.errors + [e for r in rep.rows for e in r.errors]
