"""Web interface for the primer analysis program.

Usage:
    python app.py
Then open http://127.0.0.1:5000 in a browser.
"""

import csv
import io
import os

from flask import Flask, Response, render_template, request

from primer_analysis import GC_RANGE, TM_RANGE, analyze, parse_fasta

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE_FASTA = os.path.join(HERE, "data", "primers.fasta")

app = Flask(__name__)


def records_from_text(text):
    """Turn the text from the form into (name, sequence) pairs.

    The text can be FASTA, or plain sequences with one primer per line
    (those get the names primer_1, primer_2, ...).
    """
    lines = text.splitlines()
    if any(line.strip().startswith(">") for line in lines):
        return parse_fasta(lines)
    sequences = [line.strip().upper().replace(" ", "") for line in lines]
    return [(f"primer_{i}", seq) for i, seq in enumerate(filter(None, sequences), 1)]


def analyze_text(text):
    """Analyze every primer in the text. Return (results, errors)."""
    results, errors = [], []
    for name, seq in records_from_text(text):
        try:
            if not seq:
                raise ValueError("Empty sequence")
            results.append(analyze(name, seq))
        except ValueError as err:
            errors.append(f"{name}: {err}")
    return results, errors


def get_input_text():
    """Read the primers from the uploaded file, or else from the text box."""
    upload = request.files.get("fasta_file")
    if upload and upload.filename:
        return upload.read().decode("utf-8", errors="replace")
    return request.form.get("sequences", "")


@app.route("/", methods=["GET", "POST"])
def index():
    text = ""
    results, errors = [], []
    if request.method == "POST":
        text = get_input_text()
        results, errors = analyze_text(text)
        if not results and not errors:
            errors.append("Please enter at least one primer sequence.")
    return render_template(
        "index.html",
        text=text,
        results=results,
        errors=errors,
        gc_range=GC_RANGE,
        tm_range=TM_RANGE,
    )


@app.route("/example")
def example():
    with open(EXAMPLE_FASTA) as fh:
        text = fh.read()
    results, errors = analyze_text(text)
    return render_template(
        "index.html",
        text=text,
        results=results,
        errors=errors,
        gc_range=GC_RANGE,
        tm_range=TM_RANGE,
    )


@app.route("/download", methods=["POST"])
def download():
    results, _ = analyze_text(request.form.get("sequences", ""))
    out = io.StringIO()
    if results:
        writer = csv.DictWriter(out, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    return Response(
        out.getvalue(),
        mimetype="text/csv",
        headers={"Content-Disposition": "attachment; filename=primer_results.csv"},
    )


if __name__ == "__main__":
    app.run(debug=True)
