"""
CRISPR guide-RNA design and scoring utilities.

Scoring incorporates Doench 2014/2016-inspired position-weight features
and hard manufacturing filters matching Synthego sgRNA synthesis constraints.

Synthego manufactures chemically modified sgRNAs with 2'-O-methyl and
phosphorothioate modifications at the 3 terminal positions of each end.
Guides must pass manufacturing QC to be orderable.
"""

from dataclasses import dataclass, field
from .sequence_utils import reverse_complement, gc_content


@dataclass
class GuideRNA:
    sequence: str
    pam: str
    strand: str  # "+" or "-"
    position: int  # position within the input sequence
    genomic_start: int = 0  # absolute genomic coordinate
    gc: float = 0.0
    score: float = 0.0
    synthego_compatible: bool = True
    synthego_rejection_reason: str = ""
    offtarget_status: str = "PENDING"
    offtarget_hits: int = 0
    seed_gc: float = 0.0  # GC content of PAM-proximal seed region
    spans_variant: bool = False
    variant_note: str = ""
    gene: str = ""
    target_name: str = ""

    @property
    def full_target(self) -> str:
        return self.sequence + self.pam

    @property
    def synthego_order_sequence(self) -> str:
        """20nt DNA sequence in 5'→3' orientation for Synthego ordering."""
        return self.sequence.upper()


# ── Doench 2014-inspired position-weight scoring ────────────────────────────
#
# Simplified version of position-specific nucleotide preferences from
# Doench et al. (Nat Biotechnol 2014) and Doench et al. (Nat Biotechnol 2016).
# Positions are 1-indexed from PAM-distal (pos 1) to PAM-proximal (pos 20).
#
# Positive = favoured; negative = disfavoured.  Derived from on-target
# activity data across thousands of guides.

# Position-specific single-nucleotide preferences (pos → {base: weight})
# Only the strongest effects are encoded to keep the heuristic interpretable.
_POS_WEIGHTS: dict[int, dict[str, float]] = {
    # PAM-distal positions (1-5): modest effects
    1:  {"G": 0.5, "C": -0.3},
    2:  {"A": 0.3},
    3:  {"C": 0.3},
    # Mid-guide (6-15): weaker effects
    # PAM-proximal seed region (16-20): strong effects
    16: {"A": 0.3, "T": -0.5},
    17: {"A": 0.3, "G": 0.2},
    18: {"G": 0.5, "T": -0.3},
    19: {"G": 0.5, "T": -0.5},
    20: {"G": 1.5, "A": 0.5, "C": -1.0, "T": -0.8},  # strongest position
}


def _max_homopolymer(seq: str) -> int:
    """Return the length of the longest single-nucleotide run in seq."""
    if not seq:
        return 0
    max_run = 1
    current = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
            max_run = max(max_run, current)
        else:
            current = 1
    return max_run


def check_synthego_compatibility(guide: GuideRNA,
                                  min_gc: float = 0.40,
                                  max_gc: float = 0.70,
                                  max_homopolymer: int = 4,
                                  reject_poly_t: bool = True) -> None:
    """
    Apply Synthego sgRNA manufacturing hard-filters.

    Rejects guides that cannot be reliably synthesised or will perform
    poorly.  Sets guide.synthego_compatible and guide.synthego_rejection_reason.
    """
    seq = guide.sequence.upper()
    reasons = []

    # 1. GC content hard limits
    if guide.gc < min_gc:
        reasons.append(f"GC={guide.gc:.0%} below Synthego minimum {min_gc:.0%}")
    if guide.gc > max_gc:
        reasons.append(f"GC={guide.gc:.0%} above Synthego maximum {max_gc:.0%}")

    # 2. Poly-T ≥4 is a Pol III terminator signal -- absolute reject
    if reject_poly_t and "TTTT" in seq:
        reasons.append("Contains TTTT (Pol III terminator signal)")

    # 3. Any homopolymer ≥ max_homopolymer
    hr = _max_homopolymer(seq)
    if hr >= max_homopolymer:
        reasons.append(f"Homopolymer run of {hr} (max allowed: {max_homopolymer - 1})")

    # 4. Low-complexity sequence (≤2 distinct bases in 20mer)
    if len(set(seq)) <= 2:
        reasons.append("Low-complexity sequence (≤2 distinct bases)")

    # 5. Contains BsmBI/BsaI restriction site (common Golden Gate sites
    #    used in cloning -- can cause issues if guides are cloned)
    for site_name, site_seq in [("BsmBI", "CGTCTC"), ("BsaI", "GGTCTC")]:
        if site_seq in seq or reverse_complement(site_seq) in seq:
            reasons.append(f"Contains {site_name} restriction site ({site_seq})")

    if reasons:
        guide.synthego_compatible = False
        guide.synthego_rejection_reason = "; ".join(reasons)
    else:
        guide.synthego_compatible = True
        guide.synthego_rejection_reason = ""


def find_pam_sites(sequence: str, pam: str = "NGG", guide_length: int = 20,
                   region_start: int = 0) -> list[GuideRNA]:
    """
    Find all SpCas9 PAM sites (NGG) in *sequence* and return candidate
    GuideRNA objects for both strands.
    """
    guides: list[GuideRNA] = []
    seq_upper = sequence.upper()

    # Forward strand: guide is upstream of NGG
    for i in range(len(seq_upper) - guide_length - 3 + 1):
        potential_pam = seq_upper[i + guide_length : i + guide_length + 3]
        if len(potential_pam) < 3:
            continue
        if potential_pam[1:] == "GG":
            guide_seq = seq_upper[i : i + guide_length]
            if "N" in guide_seq:
                continue  # skip guides spanning uncovered bases
            guides.append(GuideRNA(
                sequence=guide_seq,
                pam=potential_pam,
                strand="+",
                position=i,
                genomic_start=region_start + i,
                gc=gc_content(guide_seq),
                seed_gc=gc_content(guide_seq[-12:]),
            ))

    # Reverse strand
    rc = reverse_complement(seq_upper)
    seq_len = len(seq_upper)
    for i in range(len(rc) - guide_length - 3 + 1):
        potential_pam = rc[i + guide_length : i + guide_length + 3]
        if len(potential_pam) < 3:
            continue
        if potential_pam[1:] == "GG":
            guide_seq = rc[i : i + guide_length]
            if "N" in guide_seq:
                continue  # skip guides spanning uncovered bases
            # Map position back to forward-strand coordinate
            fwd_pos = seq_len - (i + guide_length)
            guides.append(GuideRNA(
                sequence=guide_seq,
                pam=potential_pam,
                strand="-",
                position=fwd_pos,
                genomic_start=region_start + fwd_pos,
                gc=gc_content(guide_seq),
                seed_gc=gc_content(guide_seq[-12:]),
            ))

    return guides


def score_guide(guide: GuideRNA) -> float:
    """
    Score guide quality using Doench 2014/2016-inspired position-weight
    features plus structural penalties.

    Returns 0-100.  Guides failing Synthego hard-filters are capped at 0.
    """
    if not guide.synthego_compatible:
        return 0.0

    score = 50.0
    seq = guide.sequence.upper()

    # ── 1. GC content (10 points) ──
    # Optimal: 40-60% (sweet spot); acceptable: 60-70%; penalise extremes
    if 0.40 <= guide.gc <= 0.60:
        score += 10
    elif 0.60 < guide.gc <= 0.70:
        score += 5
    elif guide.gc < 0.35 or guide.gc > 0.75:
        score -= 15

    # ── 2. Position-specific nucleotide preferences (up to ~15 points) ──
    for pos_1idx, weights in _POS_WEIGHTS.items():
        if pos_1idx <= len(seq):
            base = seq[pos_1idx - 1]
            score += weights.get(base, 0.0) * 3  # scale factor

    # ── 3. Poly-nucleotide penalties ──
    # TTT is penalised (but TTTT already rejected by Synthego filter)
    if "TTT" in seq:
        score -= 10
    if "GGG" in seq:
        score -= 5
    if "CCC" in seq:
        score -= 3

    # ── 4. Seed region (PAM-proximal 12nt) quality ──
    seed = seq[-12:]
    seed_gc = gc_content(seed)
    # Seed GC 40-75% is optimal for specificity + activity
    if seed_gc < 0.30:
        score -= 10
    elif seed_gc > 0.80:
        score -= 8
    # Penalise low-complexity seed
    if len(set(seed)) <= 2:
        score -= 20

    # ── 5. Self-complementarity / hairpin check ──
    for k in range(4, 8):
        for j in range(len(seq) - k):
            frag = seq[j : j + k]
            rc_frag = reverse_complement(frag)
            if rc_frag in seq[j + k :]:
                score -= 5 * (k - 3)
                break

    # ── 6. Dinucleotide at PAM junction (position 20-PAM) ──
    # G-G at position 20 into PAM (G-NGG) is highly favoured
    if seq[-1] == "G":
        score += 5

    return max(0.0, min(100.0, score))


def score_guides(guides: list[GuideRNA],
                 apply_synthego_filters: bool = True) -> list[GuideRNA]:
    """
    Apply Synthego manufacturing filters, score all guides, and sort
    descending by score.  Synthego-incompatible guides are placed last.
    """
    from ..config import (SYNTHEGO_MIN_GC, SYNTHEGO_MAX_GC,
                          SYNTHEGO_MAX_HOMOPOLYMER, SYNTHEGO_REJECT_POLY_T)

    for g in guides:
        if apply_synthego_filters:
            check_synthego_compatibility(
                g,
                min_gc=SYNTHEGO_MIN_GC,
                max_gc=SYNTHEGO_MAX_GC,
                max_homopolymer=SYNTHEGO_MAX_HOMOPOLYMER,
                reject_poly_t=SYNTHEGO_REJECT_POLY_T,
            )
        g.score = score_guide(g)

    # Sort: Synthego-compatible first (desc by score), then incompatible
    return sorted(guides,
                  key=lambda g: (g.synthego_compatible, g.score),
                  reverse=True)
