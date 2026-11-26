"""Tests for the Flask web interface (app.py).

Run from the day04 folder:
    python -m unittest -v
"""

import unittest

from app import app, records_from_text


class TestRecordsFromText(unittest.TestCase):
    def test_fasta_text(self):
        text = ">a\nacgt\n>b\nGGCC\n"
        self.assertEqual(records_from_text(text), [("a", "ACGT"), ("b", "GGCC")])

    def test_plain_sequences_get_names(self):
        text = "acgt\n\nGG CC\n"
        self.assertEqual(records_from_text(text),
                         [("primer_1", "ACGT"), ("primer_2", "GGCC")])


class TestWebPages(unittest.TestCase):
    def setUp(self):
        app.config["TESTING"] = True
        self.client = app.test_client()

    def test_home_page(self):
        response = self.client.get("/")
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"DNA Primer Analysis", response.data)

    def test_analyze_one_primer(self):
        response = self.client.post("/", data={"sequences": ">M13\nGTAAAACGACGGCCAGT"})
        html = response.data.decode()
        self.assertIn("ACTGGCCGTCGTTTTAC", html)  # reverse complement
        self.assertIn("5228.5", html)            # molecular weight
        self.assertIn("Tm 47.1C outside 50-65C", html)

    def test_invalid_sequence_shows_error(self):
        response = self.client.post("/", data={"sequences": ">bad\nATGXN"})
        html = response.data.decode()
        self.assertIn("bad: Invalid bases in sequence: NX", html)

    def test_empty_input_shows_error(self):
        response = self.client.post("/", data={"sequences": ""})
        self.assertIn(b"Please enter at least one primer", response.data)

    def test_upload_fasta_file(self):
        from io import BytesIO
        data = {"fasta_file": (BytesIO(b">T7\nTAATACGACTCACTATAGGG\n"), "p.fasta")}
        response = self.client.post("/", data=data, content_type="multipart/form-data")
        self.assertIn(b"CCCTATAGTGAGTCGTATTA", response.data)

    def test_example_page_has_all_primers(self):
        response = self.client.get("/example")
        self.assertIn(b"Results (10 primers)", response.data)

    def test_download_csv(self):
        response = self.client.post("/download", data={"sequences": "GGCC"})
        self.assertEqual(response.mimetype, "text/csv")
        lines = response.data.decode().splitlines()
        self.assertTrue(lines[0].startswith("name,sequence,length"))
        self.assertTrue(lines[1].startswith("primer_1,GGCC,4,100.0"))


if __name__ == "__main__":
    unittest.main()
