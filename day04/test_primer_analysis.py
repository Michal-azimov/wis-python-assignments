"""Tests for primer_analysis.py.

Run from the day03 folder:
    python -m unittest -v
"""

import os
import tempfile
import unittest

from primer_analysis import (
    analyze,
    gc_content,
    melting_temperature,
    molecular_weight,
    pcr_warnings,
    read_fasta,
    reverse_complement,
    validate,
)

M13_FORWARD = "GTAAAACGACGGCCAGT"
INVALID = "ATGXN"


class TestReadFasta(unittest.TestCase):
    def test_reads_records(self):
        text = ">first\nacgt\nAC\n\n>second\nGGCC\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as fh:
            fh.write(text)
            path = fh.name
        try:
            records = read_fasta(path)
        finally:
            os.remove(path)
        # Lines of one record are joined and upper-cased; blank lines are skipped.
        self.assertEqual(records, [("first", "ACGTAC"), ("second", "GGCC")])


class TestValidate(unittest.TestCase):
    def test_valid_sequence(self):
        validate("ACGT")

    def test_invalid_letters(self):
        with self.assertRaises(ValueError) as ctx:
            validate(INVALID)
        self.assertIn("NX", str(ctx.exception))


class TestGcContent(unittest.TestCase):
    def test_values(self):
        self.assertEqual(gc_content("GGCC"), 100)
        self.assertEqual(gc_content("AATT"), 0)
        self.assertEqual(gc_content("ATGC"), 50)

    def test_invalid_letters(self):
        with self.assertRaises(ValueError):
            gc_content(INVALID)


class TestMeltingTemperature(unittest.TestCase):
    def test_short_primer_uses_wallace_rule(self):
        # 5 A/T and 5 G/C: 2*5 + 4*5 = 30
        self.assertEqual(melting_temperature("ACGTACGTAC"), 30)

    def test_long_primer_uses_gc_formula(self):
        # 9 G/C out of 17: 64.9 + 41 * (9 - 16.4) / 17
        self.assertAlmostEqual(melting_temperature(M13_FORWARD), 47.05, places=2)

    def test_invalid_letters(self):
        with self.assertRaises(ValueError):
            melting_temperature(INVALID)


class TestMolecularWeight(unittest.TestCase):
    def test_single_base(self):
        self.assertAlmostEqual(molecular_weight("A"), 313.21 - 61.96)

    def test_m13_forward_matches_published_value(self):
        self.assertAlmostEqual(molecular_weight(M13_FORWARD), 5228.5, places=1)

    def test_invalid_letters(self):
        with self.assertRaises(ValueError):
            molecular_weight(INVALID)


class TestReverseComplement(unittest.TestCase):
    def test_values(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement(M13_FORWARD), "ACTGGCCGTCGTTTTAC")

    def test_twice_gives_original(self):
        self.assertEqual(reverse_complement(reverse_complement(M13_FORWARD)), M13_FORWARD)

    def test_invalid_letters(self):
        with self.assertRaises(ValueError):
            reverse_complement(INVALID)


class TestPcrWarnings(unittest.TestCase):
    def test_good_primer_has_no_warnings(self):
        self.assertEqual(pcr_warnings("GCTAGTTATTGCTCAGCGG"), [])  # T7 terminator

    def test_low_gc_and_low_tm(self):
        warnings = pcr_warnings("ATTTAGGTGACACTATAG")  # SP6 promoter
        self.assertEqual(len(warnings), 2)
        self.assertTrue(warnings[0].startswith("GC"))
        self.assertTrue(warnings[1].startswith("Tm"))

    def test_high_gc(self):
        self.assertEqual(pcr_warnings("CGCAAATGGGCGGTAGGCGTG"),  # CMV forward
                         ["GC 66.7% outside 40-60%"])

    def test_invalid_letters(self):
        with self.assertRaises(ValueError):
            pcr_warnings(INVALID)


class TestAnalyze(unittest.TestCase):
    def test_result_row(self):
        row = analyze("M13_forward_-20", M13_FORWARD)
        self.assertEqual(row["name"], "M13_forward_-20")
        self.assertEqual(row["length"], 17)
        self.assertEqual(row["gc_percent"], 52.9)
        self.assertEqual(row["tm_celsius"], 47.1)
        self.assertEqual(row["mw_g_per_mol"], 5228.5)
        self.assertEqual(row["reverse_complement"], "ACTGGCCGTCGTTTTAC")
        self.assertEqual(row["warnings"], "Tm 47.1C outside 50-65C")

    def test_invalid_letters(self):
        with self.assertRaises(ValueError):
            analyze("bad", INVALID)


if __name__ == "__main__":
    unittest.main()
