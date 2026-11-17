"""Day 05: pytest tests for the business logic in primer_analysis.py.

These tests add to test_primer_analysis.py (written with unittest on day 03).
They focus on edge cases, boundary values and the command-line program.

Run from the day04 folder:
    python -m pytest -v
"""

import csv
import sys

import pytest

from primer_analysis import (
    GC_RANGE,
    TM_RANGE,
    analyze,
    gc_content,
    main,
    melting_temperature,
    molecular_weight,
    parse_fasta,
    pcr_warnings,
    read_fasta,
    reverse_complement,
    validate,
)

ALL_FUNCTIONS = [gc_content, melting_temperature, molecular_weight,
                 reverse_complement, pcr_warnings]


# ---------- validate ----------

@pytest.mark.parametrize("seq", ["A", "ACGT", "GATTACA" * 10])
def test_validate_accepts_dna(seq):
    validate(seq)  # must not raise


@pytest.mark.parametrize("seq", [
    "ACGU",      # RNA base
    "ACGN",      # ambiguous base
    "acgt",      # lower case (parse_fasta upper-cases before validating)
    "AC GT",     # space inside
    "ACGT\n",    # newline
    "123",
])
def test_validate_rejects_other_characters(seq):
    with pytest.raises(ValueError, match="Invalid bases"):
        validate(seq)


@pytest.mark.parametrize("func", ALL_FUNCTIONS)
def test_empty_sequence_raises_value_error(func):
    # Before the fix gc_content("") crashed with ZeroDivisionError and
    # molecular_weight("") returned a negative weight.
    with pytest.raises(ValueError, match="Empty sequence"):
        func("")


# ---------- gc_content ----------

@pytest.mark.parametrize("seq, expected", [
    ("G", 100),
    ("A", 0),
    ("GCGCAT", 200 / 3),
    ("ATATATATGC", 20),
])
def test_gc_content(seq, expected):
    assert gc_content(seq) == pytest.approx(expected)


@pytest.mark.parametrize("seq", ["ACGTTG", "GGGAAATTTCCC", "GTAAAACGACGGCCAGT"])
def test_gc_content_same_for_reverse_complement(seq):
    # G pairs with C, so the reverse complement has the same GC content.
    assert gc_content(reverse_complement(seq)) == pytest.approx(gc_content(seq))


# ---------- melting_temperature ----------

def test_wallace_rule_used_up_to_13_bases():
    assert melting_temperature("A" * 13) == 26          # 2 * 13
    assert melting_temperature("G" * 13) == 52          # 4 * 13


def test_gc_formula_used_from_14_bases():
    # 64.9 + 41 * (0 - 16.4) / 14
    assert melting_temperature("A" * 14) == pytest.approx(16.871, abs=0.001)
    # 64.9 + 41 * (14 - 16.4) / 14
    assert melting_temperature("G" * 14) == pytest.approx(57.871, abs=0.001)


def test_more_gc_means_higher_tm():
    low = melting_temperature("ATATATATATATATATATAT")
    high = melting_temperature("GCGCGCGCGCGCGCGCGCGC")
    assert high > low


# ---------- molecular_weight ----------

@pytest.mark.parametrize("base, mass", [
    ("A", 313.21), ("T", 304.20), ("G", 329.21), ("C", 289.18),
])
def test_molecular_weight_single_bases(base, mass):
    assert molecular_weight(base) == pytest.approx(mass - 61.96)


def test_molecular_weight_grows_with_length():
    assert molecular_weight("ACGTACGT") > molecular_weight("ACGT")


def test_molecular_weight_does_not_depend_on_order():
    assert molecular_weight("AACCGGTT") == pytest.approx(molecular_weight("TGCATGCA"))


# ---------- reverse_complement ----------

@pytest.mark.parametrize("seq, expected", [
    ("A", "T"),
    ("AAAC", "GTTT"),
    ("GAATTC", "GAATTC"),   # EcoRI site is a palindrome
])
def test_reverse_complement(seq, expected):
    assert reverse_complement(seq) == expected


def test_reverse_complement_keeps_length():
    seq = "TAATACGACTCACTATAGGG"
    assert len(reverse_complement(seq)) == len(seq)


# ---------- pcr_warnings ----------

def test_gc_on_range_edges_is_ok():
    # 10 G/C out of 20 = 50%, 8 out of 20 = 40% (lower edge), 12 out of 20 = 60%
    for gc_count in (8, 10, 12):
        seq = "G" * gc_count + "A" * (20 - gc_count)
        assert not any(w.startswith("GC") for w in pcr_warnings(seq)), seq


def test_gc_just_outside_range_warns():
    low = "G" * 7 + "A" * 13      # 35%
    high = "G" * 13 + "A" * 7     # 65%
    assert pcr_warnings(low)[0].startswith("GC 35.0%")
    assert pcr_warnings(high)[0].startswith("GC 65.0%")


def test_warning_text_uses_the_ranges():
    warnings = pcr_warnings("A" * 20)
    assert warnings == [
        f"GC 0.0% outside {GC_RANGE[0]}-{GC_RANGE[1]}%",
        f"Tm 31.3C outside {TM_RANGE[0]}-{TM_RANGE[1]}C",
    ]


# ---------- parse_fasta / read_fasta ----------

def test_parse_fasta_multi_line_and_lower_case():
    lines = [">p1 some description", "acgt", "ACGT", "", ">p2", "gg"]
    assert parse_fasta(lines) == [("p1 some description", "ACGTACGT"), ("p2", "GG")]


def test_parse_fasta_removes_spaces_inside_sequence():
    # Bug found while writing these tests: "AC GT" used to stay "AC GT".
    assert parse_fasta([">p", "AC GT", "\tAA  "]) == [("p", "ACGTAA")]


def test_parse_fasta_header_without_sequence():
    assert parse_fasta([">empty", ">p", "AC"]) == [("empty", ""), ("p", "AC")]


def test_parse_fasta_empty_input():
    assert parse_fasta([]) == []


def test_read_fasta_example_file():
    records = read_fasta("data/primers.fasta")
    assert len(records) == 10
    assert records[0] == ("M13_forward_-20", "GTAAAACGACGGCCAGT")


def test_read_fasta_missing_file():
    with pytest.raises(FileNotFoundError):
        read_fasta("data/no_such_file.fasta")


# ---------- analyze ----------

def test_analyze_has_all_columns():
    row = analyze("t", "ACGT")
    assert list(row) == ["name", "sequence", "length", "gc_percent", "tm_celsius",
                         "mw_g_per_mol", "reverse_complement", "warnings"]


def test_analyze_good_primer_has_empty_warnings():
    assert analyze("T7_terminator", "GCTAGTTATTGCTCAGCGG")["warnings"] == ""


def test_analyze_rounds_to_one_decimal():
    row = analyze("t", "GCGCAT")
    assert row["gc_percent"] == 66.7


# ---------- the command-line program (main) ----------

def run_main(monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", ["primer_analysis.py", *args])
    main()


def test_main_writes_csv(monkeypatch, tmp_path, capsys):
    out = tmp_path / "out.csv"
    run_main(monkeypatch, "data/primers.fasta", str(out))

    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 10
    assert rows[0]["name"] == "M13_forward_-20"
    assert rows[0]["mw_g_per_mol"] == "5228.5"
    assert "Saved 10 results" in capsys.readouterr().out


def test_main_matches_saved_results(monkeypatch, tmp_path):
    out = tmp_path / "out.csv"
    run_main(monkeypatch, "data/primers.fasta", str(out))
    assert out.read_text() == open("data/results.csv").read()


def test_main_wrong_number_of_arguments(monkeypatch):
    with pytest.raises(SystemExit) as exc:
        run_main(monkeypatch, "only_one_argument")
    assert "Usage" in str(exc.value)


def test_main_skips_invalid_primers(monkeypatch, tmp_path, capsys):
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">good\nGCTAGTTATTGCTCAGCGG\n>bad\nACGXN\n>empty\n")
    out = tmp_path / "out.csv"
    run_main(monkeypatch, str(fasta), str(out))

    with open(out, newline="") as fh:
        names = [row["name"] for row in csv.DictReader(fh)]
    assert names == ["good"]
    err = capsys.readouterr().err
    assert "Skipping bad: Invalid bases in sequence: NX" in err
    assert "Skipping empty: Empty sequence" in err


def test_main_no_valid_primers(monkeypatch, tmp_path):
    fasta = tmp_path / "in.fasta"
    fasta.write_text("")
    with pytest.raises(SystemExit) as exc:
        run_main(monkeypatch, str(fasta), str(tmp_path / "out.csv"))
    assert "No valid primers" in str(exc.value)
