"""Day 06: pytest tests for the primer pair checker (primer_pairs.py) and its web page.

Run from the day06 folder:
    python -m pytest -v
"""

import sys

import pytest

from app import app
from primer_pairs import (
    check_pair,
    cross_dimer,
    find_product,
    hairpin_stem,
    main,
    self_dimer,
    three_prime_dimer,
)

M13_FORWARD = "GTAAAACGACGGCCAGT"
M13_REVERSE = "CAGGAAACAGCTATGAC"
TEMPLATE = ("GTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGG"
            "CATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGC")


def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


# ---------- primer-dimers ----------

def test_palindrome_is_a_full_self_dimer():
    # GAATTC is its own reverse complement, so the whole primer pairs with itself.
    assert self_dimer("GAATTCGAATTC") == 12


def test_no_dimer_when_3_prime_base_cannot_pair():
    # The last base A would need a T in the other primer, and there is none.
    assert three_prime_dimer("CCCCA", "GGGG") == 0


def test_dimer_counts_only_3_prime_end():
    # The 3' end "GGCC" can pair (GGCC is a palindrome), the 5' A's cannot.
    assert three_prime_dimer("AAAAGGCC", "TTTTTTTT") == 0
    assert three_prime_dimer("AAAAGGCC", "GGCC") == 4


def test_cross_dimer_checks_both_directions():
    a, b = "TTTTTTAC", "CCCCCCGT"
    # a's 3' end "AC" pairs with b's "GT", and b's 3' end "GT" pairs with a's "AC".
    assert cross_dimer(a, b) == cross_dimer(b, a) == 2


def test_m13_primers_have_small_dimers():
    assert self_dimer(M13_FORWARD) < 4
    assert self_dimer(M13_REVERSE) < 4
    assert cross_dimer(M13_FORWARD, M13_REVERSE) < 4


# ---------- hairpins ----------

def test_hairpin_found():
    # GGGG ... CCCC with a 4-base loop folds into a 4 bp stem.
    assert hairpin_stem("GGGGAAAACCCC") == 4


def test_hairpin_needs_a_loop():
    # Without room for a 3-base loop, the stem is shorter.
    assert hairpin_stem("GGGGCCCC") == 2
    assert hairpin_stem("GGGGCCCC", min_loop=0) == 4


def test_no_hairpin_in_poly_a():
    assert hairpin_stem("AAAAAAAAAA") == 0


def test_hairpin_of_very_short_sequence():
    assert hairpin_stem("A") == 0


# ---------- in-silico PCR ----------

def test_product_on_plus_strand():
    product = find_product(M13_FORWARD, M13_REVERSE, TEMPLATE)
    assert product["strand"] == "+"
    assert product["start"] == 1
    assert product["length"] == 103
    assert product["product"].startswith(M13_FORWARD)
    assert product["product"].endswith(reverse_complement(M13_REVERSE))


def test_product_on_minus_strand():
    # Same primers on the other strand of the template give the same product.
    product = find_product(M13_FORWARD, M13_REVERSE, reverse_complement(TEMPLATE))
    assert product["strand"] == "-"
    assert product["length"] == 103


def test_product_in_the_middle_of_template():
    product = find_product("CCCC", "TTTT", "GGGCCCCAATTAAAAGG")
    assert product == {"strand": "+", "start": 4, "end": 15,
                       "length": 12, "product": "CCCCAATTAAAA"}


def test_no_product_when_reverse_primer_is_missing():
    assert find_product(M13_FORWARD, "ACACACACACACACAC", TEMPLATE) is None


def test_no_product_when_primers_face_away():
    # The reverse primer site is before the forward primer.
    assert find_product("GGGG", "AAAA", "TTTTCCCCGGGG") is None


def test_invalid_template_raises():
    with pytest.raises(ValueError):
        find_product(M13_FORWARD, M13_REVERSE, "ACGTNNACGT")


# ---------- check_pair ----------

def test_good_pair_has_no_warnings():
    result = check_pair(M13_FORWARD, M13_REVERSE, TEMPLATE)
    assert result["warnings"] == []
    assert result["product"]["length"] == 103
    assert result["tm_difference"] == 2.5


def test_bad_pair_gets_warnings():
    result = check_pair("GAATTCGAATTC", "GGGGGGGGGGGGGGGGGGGG", "ACGTACGT")
    warnings = "\n".join(result["warnings"])
    assert "Tm differs" in warnings
    assert "forward primer forms a 3' self-dimer (12 bp)" in warnings
    assert "forward primer forms a hairpin" in warnings
    assert "no PCR product" in warnings


def test_without_template_no_product_check():
    result = check_pair(M13_FORWARD, M13_REVERSE)
    assert result["product"] is None
    assert result["template_checked"] is False
    assert not any("PCR product" in w for w in result["warnings"])


@pytest.mark.parametrize("forward, reverse", [("", "ACGT"), ("ACGT", "ACGX")])
def test_check_pair_rejects_bad_primers(forward, reverse):
    with pytest.raises(ValueError):
        check_pair(forward, reverse)


# ---------- command-line program ----------

def test_main_with_template(tmp_path, monkeypatch, capsys):
    fasta = tmp_path / "template.fasta"
    fasta.write_text(">pUC19\n" + TEMPLATE + "\n")
    monkeypatch.setattr(sys, "argv", ["primer_pairs.py", M13_FORWARD, M13_REVERSE, str(fasta)])
    main()
    out = capsys.readouterr().out
    assert "PCR product:         103 bp" in out
    assert "The primer pair looks good." in out


def test_main_wrong_arguments(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["primer_pairs.py", "ACGT"])
    with pytest.raises(SystemExit) as exc:
        main()
    assert "Usage" in str(exc.value)


def test_main_invalid_primer(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["primer_pairs.py", "ACGX", "ACGT"])
    with pytest.raises(SystemExit) as exc:
        main()
    assert "Invalid bases" in str(exc.value)


# ---------- web page ----------

@pytest.fixture
def client():
    app.config["TESTING"] = True
    return app.test_client()


def test_pair_page_loads(client):
    response = client.get("/pair")
    assert response.status_code == 200
    assert b"Primer Pair Check" in response.data


def test_pair_example(client):
    html = client.get("/pair/example").data.decode()
    assert "PCR product: 103 bp" in html
    assert "The primer pair looks good." in html


def test_pair_form_with_plain_template(client):
    html = client.post("/pair", data={
        "forward": "gtaaaacgac ggccagt", "reverse": M13_REVERSE, "template": TEMPLATE,
    }).data.decode()
    assert "PCR product: 103 bp" in html


def test_pair_form_missing_primer(client):
    html = client.post("/pair", data={"forward": M13_FORWARD, "reverse": ""}).data.decode()
    assert "Please enter both a forward and a reverse primer." in html


def test_pair_form_invalid_primer(client):
    html = client.post("/pair", data={"forward": "ACGTX", "reverse": "ACGT"}).data.decode()
    assert "Invalid bases in sequence: X" in html


def test_pair_form_shows_warnings(client):
    html = client.post("/pair", data={
        "forward": "GAATTCGAATTC", "reverse": M13_REVERSE,
    }).data.decode()
    assert "self-dimer (12 bp)" in html


def test_home_page_links_to_pair_page(client):
    assert b'href="/pair"' in client.get("/").data
