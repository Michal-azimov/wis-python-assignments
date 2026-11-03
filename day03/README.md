# Day 03 – DNA Primer Analysis

## What the program does

PCR and sequencing need short DNA strands called *primers*. Before you order a
primer you check a few of its properties. [`primer_analysis.py`](primer_analysis.py)
reads primers from a FASTA file and calculates, for each primer:

| Property | Function | How it is computed |
|---|---|---|
| GC content (%) | `gc_content()` | (G + C) / length × 100 |
| Melting temperature (°C) | `melting_temperature()` | Wallace rule `2(A+T) + 4(G+C)` for primers shorter than 14 nt; otherwise `64.9 + 41 × (G+C − 16.4) / N` |
| Molecular weight (g/mol) | `molecular_weight()` | Sum of the nucleotide masses − 61.96 (single-stranded, no 5' phosphate) |
| Reverse complement | `reverse_complement()` | Complement each base and reverse the order |
| PCR warnings | `pcr_warnings()` | Flags GC content outside 40–60% or Tm outside 50–65 °C, the ranges usually recommended for PCR primers |

The program rejects sequences that contain characters other than A, C, G, T.

## Data

[`data/primers.fasta`](data/primers.fasta) contains common lab primers (M13,
T7, SP6, BGH, CMV, pGEX, 16S rRNA 1492R) and one short test probe that
checks the Wallace-rule branch. The sequences are the standard published primer
sequences.

The program writes its results to [`data/results.csv`](data/results.csv).

## How to run

```
cd day03
python primer_analysis.py data/primers.fasta data/results.csv
```

Output:

```
name                  len    GC%   Tm(C)        MW
M13_forward_-20        17   52.9    47.1    5228.5  Tm 47.1C outside 50-65C
M13_reverse            17   47.1    44.6    5212.5  Tm 44.6C outside 50-65C
T7_promoter            20   40.0    47.7    6125.1  Tm 47.7C outside 50-65C
T7_terminator          19   52.6    51.1    5834.8
SP6_promoter           18   33.3    41.2    5537.7  GC 33.3% outside 40-60%; Tm 41.2C outside 50-65C
BGH_reverse            18   55.6    50.3    5597.7
CMV_forward            21   66.7    60.2    6552.3  GC 66.7% outside 40-60%
pGEX_5prime            23   65.2    62.4    7136.7  GC 65.2% outside 40-60%
16S_1492R              19   42.1    46.8    5784.8  Tm 46.8C outside 50-65C
short_test_probe       10   50.0      30    3012.0  Tm 30.0C outside 50-65C

Saved 10 results to data/results.csv
```

The M13 forward (−20) molecular weight of 5228.5 g/mol matches the value
that oligo suppliers publish.

Many of the standard primers get a Tm warning. The reason is that the simple
GC-based formula usually gives lower values than the nearest-neighbor method
that oligo suppliers use, so the warnings are a rough first check.

## Prompts I used

I worked with an AI coding assistant (Claude). My prompts were:

1. *"Do this assignment for me"*, with a screenshot of the Day 03 assignment.

   From that prompt the assistant chose a biology computation (primer
   properties), wrote the program and the FASTA data file, ran the program to
   check the output, and wrote this README.
2. *"Add a check that warns when a primer's GC content or melting temperature
   is outside the range usually recommended for PCR."*
   The assistant added `pcr_warnings()` and a `warnings` column in the CSV.

## Git steps

1. Created a branch `day03` (`git switch -c day03`).
2. Added `day03/primer_analysis.py`, `day03/data/primers.fasta`,
   `day03/data/results.csv` and this README.
3. Committed and pushed to GitHub.
4. Opened an issue on the course repository with a link to this folder.
