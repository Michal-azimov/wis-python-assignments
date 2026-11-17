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
| [`test_business_logic.py`](test_business_logic.py) | More pytest tests for the business logic (day 05 assignment) |
| [`data/primers.fasta`](data/primers.fasta) | Example primers |
| [`requirements.txt`](requirements.txt) | Python packages needed (Flask, pytest) |

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

## Day 05 – More tests for the business logic

The Day 05 assignment was to write tests for the "business logic" of a
previous application. I added [`test_business_logic.py`](test_business_logic.py),
written with [pytest](https://docs.pytest.org/) like the lecturer's
`test_bmi_pytest.py` example. It has 51 tests that cover what the day 03
tests did not:

- **Edge cases:** empty sequences, RNA and ambiguous bases (`U`, `N`), lower
  case, spaces and newlines inside a sequence, an empty FASTA file, a header
  with no sequence, a missing file.
- **Boundary values:** the Tm formula switches from the Wallace rule at 13
  bases to the GC formula at 14 bases; GC content exactly at 40% and 60% is
  OK, but 35% and 65% get a warning.
- **Properties that must always hold:** the reverse complement has the same
  length and GC content as the primer, a palindrome (the EcoRI site `GAATTC`)
  is its own reverse complement, the molecular weight does not depend on the
  order of the bases, and more G/C gives a higher Tm.
- **The command-line program (`main()`):** it writes the expected CSV (the
  same as [`data/results.csv`](data/results.csv)), shows the usage message
  for wrong arguments, skips invalid primers, and stops with a message when
  no primer is valid. These tests use pytest's `tmp_path`, `monkeypatch` and
  `capsys` fixtures.
- `@pytest.mark.parametrize` runs one test with many inputs, so every
  function is checked with an empty sequence.

### Bugs the tests found

I wrote the tests first and 8 of them failed. They showed 4 bugs, which I
then fixed in [`primer_analysis.py`](primer_analysis.py):

| Bug | Before | After |
|---|---|---|
| Empty sequence | `gc_content("")` crashed with `ZeroDivisionError`; `molecular_weight("")` returned −61.96 | `validate()` raises `ValueError("Empty sequence")` |
| Space inside a FASTA sequence line (`AC GT`) | The space stayed, so the primer was rejected | `parse_fasta()` removes all whitespace |
| One invalid primer in the FASTA file | The command-line program crashed with a traceback | It prints `Skipping <name>: <reason>` and analyzes the others |
| FASTA file with no valid primers | Crashed with `IndexError` | Stops with `No valid primers found in <file>` |

The web app already handled empty sequences itself; that check is now done
by the business logic, so I removed it from `app.py`.

### Running all the tests

```
cd day04
pip install -r requirements.txt
python -m pytest -v
```

```
80 passed in 0.22s
```

pytest also runs the older `unittest` tests, so this one command runs all 80
tests: 51 new ones, 20 from day 03 and 9 for the web app.
`python -m unittest -v` still runs the 29 older tests.

I also checked that the command-line program still writes the same
`data/results.csv`, that the web page still shows the example primers, and
that the `day03` folder and its tests are unchanged.

### Prompts I used (day 05)

1. *"I have another assignment. Can you do it for me?"*, with a screenshot of
   the Day 05 assignment.

   The assistant read the course repository, saw that the lecturer used
   pytest, wrote `test_business_logic.py`, ran it to find the bugs above,
   fixed them, ran all the tests and both programs, and updated this README.
2. *"Yes, all three"* (commit, push and open the Day 05 issue).

### Git steps (day 05)

1. Created a branch `day05` (`git switch -c day05`).
2. Committed the new tests, the bug fixes and the README updates.
3. Merged `day05` into `main` and pushed to GitHub.
4. Opened a "Day 05 by Michal Azimov" issue on the course repository with a
   link to this folder.
