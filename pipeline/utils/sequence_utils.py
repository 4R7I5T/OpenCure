"""
DNA sequence manipulation utilities.
"""

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    upper = seq.upper()
    return (upper.count("G") + upper.count("C")) / len(upper)


def clean_sequence(seq: str) -> str:
    """Keep only ACGT characters, uppercased."""
    return "".join(c for c in seq.upper() if c in "ACGT")


def apply_variants_to_sequence(ref_seq: str, variants: list, region_start: int) -> str:
    """
    Apply a sorted list of homozygous-alt variants to a reference sequence
    to produce the patient's actual sequence.

    *variants* must have .pos, .ref, .alt, .zygosity attributes.
    *region_start* is the 1-based genomic coordinate of ref_seq[0].
    """
    seq = list(ref_seq)
    offset = 0

    for var in sorted(variants, key=lambda v: v.pos):
        if var.zygosity not in ("HOM_ALT",):
            continue
        rel = var.pos - region_start + offset
        if rel < 0 or rel >= len(seq):
            continue
        ref_len = len(var.ref)
        alt_allele = var.alt[0] if var.alt else var.ref
        seq[rel : rel + ref_len] = list(alt_allele)
        offset += len(alt_allele) - ref_len

    return "".join(seq)
