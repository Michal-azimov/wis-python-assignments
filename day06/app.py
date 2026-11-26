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
from primer_pairs import check_pair

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE_FASTA = os.path.join(HERE, "data", "primers.fasta")
EXAMPLE_TEMPLATE = os.path.join(HERE, "data", "template.fasta")
# The M13 forward and reverse primers, which amplify the pUC19 cloning site.
EXAMPLE_PAIR = ("GTAAAACGACGGCCAGT", "CAGGAAACAGCTATGAC")

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


def clean(text):
    """Remove spaces and newlines from a sequence and make it upper case."""
    return "".join(text.split()).upper()


def template_from_text(text):
    """Return the template sequence from the form (FASTA or plain), or None."""
    lines = text.splitlines()
    if any(line.strip().startswith(">") for line in lines):
        records = parse_fasta(lines)
        return records[0][1] if records else None
    return clean(text) or None


def render_pair(forward="", reverse="", template_text=""):
    """Check the primer pair (if one was given) and show the pair page."""
    result, error = None, None
    if forward or reverse:
        try:
            result = check_pair(clean(forward), clean(reverse),
                                template_from_text(template_text))
        except ValueError as err:
            error = str(err)
    return render_template("pair.html", forward=forward, reverse=reverse,
                           template_text=template_text, result=result, error=error)


@app.route("/pair", methods=["GET", "POST"])
def pair():
    if request.method == "GET":
        return render_pair()
    forward = request.form.get("forward", "")
    reverse = request.form.get("reverse", "")
    if not clean(forward) or not clean(reverse):
        return render_template("pair.html", forward=forward, reverse=reverse,
                               template_text=request.form.get("template", ""),
                               result=None,
                               error="Please enter both a forward and a reverse primer.")
    return render_pair(forward, reverse, request.form.get("template", ""))


@app.route("/pair/example")
def pair_example():
    with open(EXAMPLE_TEMPLATE) as fh:
        template_text = fh.read()
    return render_pair(*EXAMPLE_PAIR, template_text)


if __name__ == "__main__":
    app.run(debug=True)
