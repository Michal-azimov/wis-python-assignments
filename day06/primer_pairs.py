"""Check whether two primers work together as a PCR primer pair.

For a forward and a reverse primer this module checks:
- primer-dimers: can the 3' end of a primer bind to itself or to the other primer?
- hairpins: can a primer fold back on itself?
- Tm difference: the two primers should melt at similar temperatures.
- in-silico PCR: given a template DNA, where do the primers bind and what
  product (amplicon) do they make?

Usage:
    python primer_pairs.py FORWARD REVERSE [template.fasta]
"""

import sys

from primer_analysis import analyze, read_fasta, reverse_complement, validate

# Warn when a 3' dimer or a hairpin stem has at least this many base pairs.
MAX_DIMER = 4
MAX_HAIRPIN = 4
# Smallest loop a hairpin can have.
MIN_LOOP = 3
# Warn when the two primers' Tm differ by more than this (degrees C).
MAX_TM_DIFFERENCE = 5

PAIRS = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}


def three_prime_dimer(primer1, primer2):
    """Number of bases at the 3' end of primer1 that can pair with primer2.

    A 3' dimer is the dangerous kind: the polymerase can extend it and make
    a short "primer-dimer" product instead of the real one.
    """
    validate(primer1)
    validate(primer2)
    target = reverse_complement(primer2)
    best = 0
    for k in range(1, len(primer1) + 1):
        if primer1[-k:] not in target:
            break
        best = k
    return best


def self_dimer(primer):
    """3' dimer of a primer with another copy of itself."""
    return three_prime_dimer(primer, primer)


def cross_dimer(forward, reverse):
    """Worst 3' dimer between the two primers, in either direction."""
    return max(three_prime_dimer(forward, reverse), three_prime_dimer(reverse, forward))


def hairpin_stem(seq, min_loop=MIN_LOOP):
    """Length of the longest hairpin stem (in base pairs) the primer can form.

    A stem pairs seq[i], seq[i+1], ... with seq[j], seq[j-1], ... and must
    leave a loop of at least `min_loop` unpaired bases between the two arms.
    """
    validate(seq)
    best = 0
    n = len(seq)
    for i in range(n):
        for j in range(i + 1, n):
            stem = 0
            while ((j - stem) - (i + stem) - 1 >= min_loop
                   and (seq[i + stem], seq[j - stem]) in PAIRS):
                stem += 1
            best = max(best, stem)
    return best


def find_product(forward, reverse, template):
    """Simulate PCR: find where the primers bind on the template.

    The forward primer must appear in the template, and the reverse
    complement of the reverse primer must appear after it. Both strands of
    the template are searched. Returns a dict with the product, or None.
    Positions are 1-based and refer to the strand given in "strand".
    """
    for seq in (forward, reverse, template):
        validate(seq)
    reverse_site = reverse_complement(reverse)
    for strand, dna in (("+", template), ("-", reverse_complement(template))):
        start = dna.find(forward)
        if start == -1:
            continue
        end = dna.find(reverse_site, start)
        if end == -1:
            continue
        product = dna[start:end + len(reverse_site)]
        return {
            "strand": strand,
            "start": start + 1,
            "end": end + len(reverse_site),
            "length": len(product),
            "product": product,
        }
    return None


def check_pair(forward, reverse, template=None):
    """Analyze a primer pair and return all the results as a dict."""
    fwd = analyze("forward", forward)
    rev = analyze("reverse", reverse)
    result = {
        "forward": fwd,
        "reverse": rev,
        "tm_difference": round(abs(fwd["tm_celsius"] - rev["tm_celsius"]), 1),
        "forward_self_dimer": self_dimer(forward),
        "reverse_self_dimer": self_dimer(reverse),
        "cross_dimer": cross_dimer(forward, reverse),
        "forward_hairpin": hairpin_stem(forward),
        "reverse_hairpin": hairpin_stem(reverse),
        "product": None,
        "template_checked": template is not None,
    }

    warnings = []
    if result["tm_difference"] > MAX_TM_DIFFERENCE:
        warnings.append(f"Tm differs by {result['tm_difference']}C "
                        f"(more than {MAX_TM_DIFFERENCE}C)")
    for name in ("forward", "reverse"):
        if result[f"{name}_self_dimer"] >= MAX_DIMER:
            warnings.append(f"{name} primer forms a 3' self-dimer "
                            f"({result[f'{name}_self_dimer']} bp)")
        if result[f"{name}_hairpin"] >= MAX_HAIRPIN:
            warnings.append(f"{name} primer forms a hairpin "
                            f"({result[f'{name}_hairpin']} bp stem)")
    if result["cross_dimer"] >= MAX_DIMER:
        warnings.append(f"primers form a 3' cross-dimer ({result['cross_dimer']} bp)")
    if template is not None:
        result["product"] = find_product(forward, reverse, template)
        if result["product"] is None:
            warnings.append("no PCR product: the primers do not both bind the template")
    result["warnings"] = warnings
    return result


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    forward = "".join(sys.argv[1].split()).upper()
    reverse = "".join(sys.argv[2].split()).upper()
    template = None
    if len(sys.argv) == 4:
        records = read_fasta(sys.argv[3])
        if not records:
            sys.exit(f"No sequence found in {sys.argv[3]}")
        template = records[0][1]

    try:
        result = check_pair(forward, reverse, template)
    except ValueError as err:
        sys.exit(f"Error: {err}")

    for name in ("forward", "reverse"):
        p = result[name]
        print(f"{name:<8} {p['sequence']}  len {p['length']}  GC {p['gc_percent']}%  "
              f"Tm {p['tm_celsius']}C")
    print(f"Tm difference:       {result['tm_difference']}C")
    print(f"3' self-dimers:      forward {result['forward_self_dimer']} bp, "
          f"reverse {result['reverse_self_dimer']} bp")
    print(f"3' cross-dimer:      {result['cross_dimer']} bp")
    print(f"Hairpin stems:       forward {result['forward_hairpin']} bp, "
          f"reverse {result['reverse_hairpin']} bp")
    product = result["product"]
    if product:
        print(f"PCR product:         {product['length']} bp, positions "
              f"{product['start']}-{product['end']} on the {product['strand']} strand")
        print(f"                     {product['product']}")
    print()
    if result["warnings"]:
        print("Warnings:")
        for warning in result["warnings"]:
            print(f"  - {warning}")
    else:
        print("The primer pair looks good.")


if __name__ == "__main__":
    main()
