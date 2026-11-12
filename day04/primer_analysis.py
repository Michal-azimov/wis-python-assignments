"""Analyze DNA primers: GC content, melting temperature and molecular weight.

Usage:
    python primer_analysis.py data/primers.fasta results.csv
"""

import csv
import sys

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

# Average mass (g/mol) of each nucleotide inside a single-stranded DNA chain.
NUCLEOTIDE_MASS = {"A": 313.21, "T": 304.20, "G": 329.21, "C": 289.18}

# Ranges usually recommended for PCR primers.
GC_RANGE = (40, 60)
TM_RANGE = (50, 65)


def parse_fasta(lines):
    """Return a list of (name, sequence) pairs from lines of FASTA text."""
    records = []
    name = None
    parts = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            name = line[1:].strip()
            parts = []
        else:
            parts.append(line.upper())
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def read_fasta(path):
    """Return a list of (name, sequence) pairs from a FASTA file."""
    with open(path) as fh:
        return parse_fasta(fh)


def validate(seq):
    bad = set(seq) - set(COMPLEMENT)
    if bad:
        raise ValueError(f"Invalid bases in sequence: {''.join(sorted(bad))}")


def gc_content(seq):
    """Percent of G and C bases in the sequence."""
    validate(seq)
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)


def melting_temperature(seq):
    """Estimated melting temperature (Tm) in degrees Celsius.

    Short primers (< 14 nt) use the Wallace rule: Tm = 2(A+T) + 4(G+C).
    Longer primers use: Tm = 64.9 + 41 * (G+C - 16.4) / N.
    """
    validate(seq)
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    if len(seq) < 14:
        return 2 * at + 4 * gc
    return 64.9 + 41 * (gc - 16.4) / len(seq)


def molecular_weight(seq):
    """Molecular weight (g/mol) of an unphosphorylated single-stranded DNA oligo."""
    validate(seq)
    return sum(NUCLEOTIDE_MASS[base] for base in seq) - 61.96


def reverse_complement(seq):
    validate(seq)
    return "".join(COMPLEMENT[base] for base in reversed(seq))


def pcr_warnings(seq):
    """Return a list of warnings for values outside the recommended PCR ranges."""
    warnings = []
    gc = gc_content(seq)
    tm = melting_temperature(seq)
    if not GC_RANGE[0] <= gc <= GC_RANGE[1]:
        warnings.append(f"GC {gc:.1f}% outside {GC_RANGE[0]}-{GC_RANGE[1]}%")
    if not TM_RANGE[0] <= tm <= TM_RANGE[1]:
        warnings.append(f"Tm {tm:.1f}C outside {TM_RANGE[0]}-{TM_RANGE[1]}C")
    return warnings


def analyze(name, seq):
    """Compute all properties of one primer and return them as a dict."""
    return {
        "name": name,
        "sequence": seq,
        "length": len(seq),
        "gc_percent": round(gc_content(seq), 1),
        "tm_celsius": round(melting_temperature(seq), 1),
        "mw_g_per_mol": round(molecular_weight(seq), 1),
        "reverse_complement": reverse_complement(seq),
        "warnings": "; ".join(pcr_warnings(seq)),
    }


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    fasta_path, out_path = sys.argv[1], sys.argv[2]

    results = [analyze(name, seq) for name, seq in read_fasta(fasta_path)]

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"{'name':<20}{'len':>5}{'GC%':>7}{'Tm(C)':>8}{'MW':>10}")
    for r in results:
        print(f"{r['name']:<20}{r['length']:>5}{r['gc_percent']:>7}"
              f"{r['tm_celsius']:>8}{r['mw_g_per_mol']:>10}  {r['warnings']}")
    print(f"\nSaved {len(results)} results to {out_path}")


if __name__ == "__main__":
    main()
