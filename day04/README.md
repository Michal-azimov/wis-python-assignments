# Day 04 – Web Interface for DNA Primer Analysis

This folder starts as a copy of [day03](../day03/) (which stays unchanged) and
adds a web application, written with [Flask](https://flask.palletsprojects.com/),
that uses the same "business logic" from
[`primer_analysis.py`](primer_analysis.py).

## Files

| File | What it is |
|---|---|
| [`primer_analysis.py`](primer_analysis.py) | The business logic and command-line program from day 03 |
| [`app.py`](app.py) | The Flask web application |
| [`templates/index.html`](templates/index.html) | The HTML page (Jinja2 template) |
| [`static/style.css`](static/style.css) | Styling for the page |
| [`test_primer_analysis.py`](test_primer_analysis.py) | Tests for the business logic (from day 03) |
| [`test_app.py`](test_app.py) | Tests for the web application |
| [`data/primers.fasta`](data/primers.fasta) | Example primers |
| [`requirements.txt`](requirements.txt) | Python packages needed (Flask) |

## What the web app does

On the page you can:

- paste primers in FASTA format, or plain sequences with one primer per line
  (they are named `primer_1`, `primer_2`, ...),
- or upload a FASTA file,
- or click **Load example primers** to analyze [`data/primers.fasta`](data/primers.fasta).

For each primer the page shows a table with length, GC content, melting
temperature, molecular weight, reverse complement and PCR warnings — the
same values the command-line program computes. Primers with invalid letters
are listed in an error box instead of crashing the page. The results can be
downloaded as a CSV file, in the same format as the command-line output.

The web app does not duplicate any calculation: it imports `analyze()` and
`parse_fasta()` from `primer_analysis.py`.

### Change to the business logic

In day 03, `read_fasta()` could only read a file from disk. The web app gets
its text from a form, so I split it into two functions:

- `parse_fasta(lines)` — parses FASTA lines (new),
- `read_fasta(path)` — opens the file and calls `parse_fasta()`.

The command-line program works exactly as before.

## How to run

Install Flask once:

```
pip install -r requirements.txt
```

### The web application

```
cd day04
python app.py
```

Then open <http://127.0.0.1:5000> in a browser.

### The command-line program (still works)

```
cd day04
python primer_analysis.py data/primers.fasta data/results.csv
```

The output is the same as in day 03.

## Checking that both systems work

```
cd day04
python -m unittest -v
```

```
Ran 29 tests in 0.085s

OK
```

The 20 day-03 tests check the business logic, and 9 new tests in
`test_app.py` use Flask's test client to check the web app: the home page,
analyzing a primer, an invalid sequence, empty input, uploading a FASTA file,
the example page and the CSV download.

I also checked by hand:

- the original `day03` program and its tests still run and pass,
- the `day04` command-line program gives the same output as `day03`,
- in a browser, the example primers, a mix of a valid and an invalid
  sequence, and the CSV download all work.

## Prompts I used

I worked with an AI coding assistant (Claude). My prompts were:

1. *"Do this assignment pls"*, with a screenshot of the Day 04 assignment.

   From that prompt the assistant copied day03 to day04, wrote the Flask app,
   the HTML template, the CSS and the tests, ran both the old and the new
   programs, checked the page in a browser, and wrote this README.
2. *"Yes, commit and push"*

## Git steps

1. Created a branch `day04` (`git switch -c day04`).
2. Added the `day04` folder and the updated root `README.md`, and committed them.
3. Merged `day04` into `main` and pushed to GitHub.
