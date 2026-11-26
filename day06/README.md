# Day 06 – Primer Pair Checker and In-silico PCR

This folder starts as a copy of my day 05 solution (in [day04](../day04/),
which stays unchanged) and adds a new feature: a **primer pair checker**.

A PCR needs two primers, a forward and a reverse one. Each primer can look
fine on its own (the day 03–05 program checks that) but the pair can still
fail. The new feature checks the two primers *together*:

| Check | Why it matters | Warning when |
|---|---|---|
| **3' self-dimer** | The 3' end of a primer sticks to another copy of the same primer, and the polymerase makes a short "primer-dimer" instead of the real product | ≥ 4 bp |
| **3' cross-dimer** | Same, but between the forward and the reverse primer | ≥ 4 bp |
| **Hairpin** | The primer folds back on itself (needs a loop of at least 3 bases) and cannot bind the template | stem ≥ 4 bp |
| **Tm difference** | Both primers must bind at the same annealing temperature | > 5 °C |
| **In-silico PCR** | Given a template DNA, find where the primers bind (on both strands) and the size and sequence of the PCR product | no product |

## Files

| File | What it is |
|---|---|
| [`primer_pairs.py`](primer_pairs.py) | **New:** primer pair checks and in-silico PCR (business logic + command-line program) |
| [`templates/pair.html`](templates/pair.html) | **New:** web page for the primer pair checker |
| [`test_primer_pairs.py`](test_primer_pairs.py) | **New:** 30 pytest tests for the new feature (logic, command line and web page) |
| [`data/template.fasta`](data/template.fasta) | **New:** example template: the pUC19 cloning site between the M13 primers |
| [`primer_analysis.py`](primer_analysis.py) | Single primer analysis (from day 03, unchanged) |
| [`app.py`](app.py) | The Flask web application, now with the `/pair` page |
| [`templates/index.html`](templates/index.html) | The single-primer page (now with a link to the new page) |
| [`static/style.css`](static/style.css) | Styling |
| [`test_primer_analysis.py`](test_primer_analysis.py), [`test_business_logic.py`](test_business_logic.py), [`test_app.py`](test_app.py) | Older tests (days 03–05) |
| [`data/primers.fasta`](data/primers.fasta) | Example primers |
| [`requirements.txt`](requirements.txt) | Python packages needed (Flask, pytest) |

The new code reuses `analyze()`, `validate()` and `reverse_complement()` from
`primer_analysis.py`, so nothing is computed twice.

## How to run

Install the packages once:

```
pip install -r requirements.txt
```

### Web application

```
cd day06
python app.py
```

Open <http://127.0.0.1:5000/pair> (or click **Primer pair & PCR** at the top of
the home page). Click **Load example** to check the M13 forward and reverse
primers on the pUC19 cloning site: the page shows a 103 bp product with the
primer binding sites highlighted, and says the pair looks good.

### Command line

```
cd day06
python primer_pairs.py GTAAAACGACGGCCAGT CAGGAAACAGCTATGAC data/template.fasta
```

```
forward  GTAAAACGACGGCCAGT  len 17  GC 52.9%  Tm 47.1C
reverse  CAGGAAACAGCTATGAC  len 17  GC 47.1%  Tm 44.6C
Tm difference:       2.5C
3' self-dimers:      forward 2 bp, reverse 1 bp
3' cross-dimer:      2 bp
Hairpin stems:       forward 2 bp, reverse 2 bp
PCR product:         103 bp, positions 1-103 on the + strand
                     GTAAAACGACGGCCAGTGAATTCGAGCTC...GGTCATAGCTGTTTCCTG

The primer pair looks good.
```

The template is optional. An example of a bad pair:

```
python primer_pairs.py GAATTCGAATTC CCCGGGAAATTT
```

```
Warnings:
  - forward primer forms a 3' self-dimer (12 bp)
  - forward primer forms a hairpin (4 bp stem)
  - reverse primer forms a 3' self-dimer (6 bp)
```

The single-primer program still works as before:

```
python primer_analysis.py data/primers.fasta data/results.csv
```

## Tests

```
cd day06
python -m pytest -v
```

```
110 passed in 0.35s
```

That is the 80 tests from days 03–05 and 30 new ones in
`test_primer_pairs.py`. The new tests check:

- dimers: a palindrome (`GAATTCGAATTC`) is a full 12 bp self-dimer, only the
  3' end counts, and the cross-dimer is checked in both directions,
- hairpins: a stem is found, it needs a loop of at least 3 bases, and poly-A
  has none,
- in-silico PCR: a product on the + strand, on the − strand, in the middle of
  a template, no product when a primer is missing or the primers face away
  from each other, and an error for an invalid template,
- the warnings for a good and a bad pair,
- the command-line program (with `tmp_path`, `monkeypatch`, `capsys`),
- the web page: the example, a form with a plain-text template, a missing
  primer, an invalid primer and a bad pair.

I also checked the page in a browser.

### Note about the checks

These are simple checks for teaching, not a replacement for tools like
Primer3: the dimer check only looks for perfect matches at the 3' end, the
hairpin check counts base pairs without computing free energy (ΔG), and the
in-silico PCR only finds exact matches.

## Other changes

- Added a [`.gitignore`](../.gitignore) to the repository so that
  `__pycache__/` and `.pytest_cache/` folders are never committed.

## Reviewing other students' repositories

Part of the assignment was to run two other students' day 05 projects and
open issues on their repositories. Issues I opened:

- [sapirg2079/Assignments#1](https://github.com/sapirg2079/Assignments/issues/1):
  how to install the packages and run the day04 web app.
- [kseniabrt/2026_programming_in_python_using_AI_course_assignments#1](https://github.com/kseniabrt/2026_programming_in_python_using_AI_course_assignments/issues/1):
  which tests belong to Day05, and an idea to avoid the Windows path-length limit.

## Prompts I used

I worked with an AI coding assistant (Claude). My prompts were:

1. *"Do this assignment for me pls"*, with a screenshot of the Day 06
   assignment.

   The assistant copied my day 05 solution to `day06`, suggested and wrote the
   primer pair checker (`primer_pairs.py`), the new web page and the tests,
   ran all the tests and both command-line programs, checked the page in a
   browser, and wrote this README. It also picked two students whose
   repositories had no issues yet, read their day 05 code and drafted the issues.
2. *"Write only 2 points to Sapir, not 4. The instructions said to write
   comments about the app, not about the issues in GitHub."*
3. *"Yes, post it"* (post the issues, push, and open the day06 issue).

## Git steps

1. Created a branch `day06` (`git switch -c day06`).
2. Copied `day04` to `day06` and added the new feature, tests and README.
3. Committed, merged `day06` into `main` and pushed to GitHub.
4. Opened a "day06 Michal Azimov" issue on the course repository with links to
   this repository and to the issues I opened.
