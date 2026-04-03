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


def apply_variants_to_sequence(ref_seq: str, variants: list,
                               region_start: int,
                               allele: str = "alt") -> str:
    """
    Apply patient variants to a reference sequence to produce an
    allele-specific patient sequence.

    *variants* must have .pos, .ref, .alt, .zygosity attributes.
    *region_start* is the 1-based genomic coordinate of ref_seq[0].
    *allele*: which allele to build:
        - "alt": apply homozygous-alt variants (original behaviour)
        - "alt_all": apply ALL variants (HOM_ALT + HET) using alt allele
          — used for allele-specific guide design targeting the mutant allele
        - "ref": return reference with no variants applied
          — used for designing guides against the wild-type allele
    """
    if allele == "ref":
        return ref_seq

    seq = list(ref_seq)
    offset = 0

    accepted_zygosities = {"HOM_ALT"}
    if allele == "alt_all":
        accepted_zygosities = {"HOM_ALT", "HET"}

    for var in sorted(variants, key=lambda v: v.pos):
        if var.zygosity not in accepted_zygosities:
            continue
        rel = var.pos - region_start + offset
        if rel < 0 or rel >= len(seq):
            continue
        ref_len = len(var.ref)
        alt_allele_seq = var.alt[0] if var.alt else var.ref
        seq[rel : rel + ref_len] = list(alt_allele_seq)
        offset += len(alt_allele_seq) - ref_len

    return "".join(seq)
