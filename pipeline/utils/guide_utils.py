"""
CRISPR guide-RNA design and scoring utilities.
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
    offtarget_status: str = "PENDING"
    offtarget_hits: int = 0
    spans_variant: bool = False
    variant_note: str = ""
    gene: str = ""
    target_name: str = ""

    @property
    def full_target(self) -> str:
        return self.sequence + self.pam


def find_pam_sites(sequence: str, pam: str = "NGG", guide_length: int = 20,
                   region_start: int = 0) -> list[GuideRNA]:
    """
    Find all SpCas9 PAM sites (NGG) in *sequence* and return candidate
    GuideRNA objects for both strands.
    """
    guides: list[GuideRNA] = []

    # Forward strand: guide is upstream of NGG
    for i in range(len(sequence) - guide_length - 3 + 1):
        potential_pam = sequence[i + guide_length : i + guide_length + 3]
        if len(potential_pam) < 3:
            continue
        if potential_pam[1:] == "GG":
            guide_seq = sequence[i : i + guide_length]
            guides.append(GuideRNA(
                sequence=guide_seq,
                pam=potential_pam,
                strand="+",
                position=i,
                genomic_start=region_start + i,
                gc=gc_content(guide_seq),
            ))

    # Reverse strand
    rc = reverse_complement(sequence)
    seq_len = len(sequence)
    for i in range(len(rc) - guide_length - 3 + 1):
        potential_pam = rc[i + guide_length : i + guide_length + 3]
        if len(potential_pam) < 3:
            continue
        if potential_pam[1:] == "GG":
            guide_seq = rc[i : i + guide_length]
            # Map position back to forward-strand coordinate
            fwd_pos = seq_len - (i + guide_length)
            guides.append(GuideRNA(
                sequence=guide_seq,
                pam=potential_pam,
                strand="-",
                position=fwd_pos,
                genomic_start=region_start + fwd_pos,
                gc=gc_content(guide_seq),
            ))

    return guides


def score_guide(guide: GuideRNA) -> float:
    """
    Heuristic scoring for guide quality.
    In production, swap for Doench 2016 / Rule Set 2 / DeepCRISPR.
    Returns 0-100.
    """
    score = 50.0

    # GC content: optimal 40-70%
    if 0.40 <= guide.gc <= 0.70:
        score += 20
    elif guide.gc < 0.25 or guide.gc > 0.80:
        score -= 25
    elif guide.gc < 0.35 or guide.gc > 0.75:
        score -= 10

    seq = guide.sequence

    # Penalise poly-T runs (Pol III terminator signal)
    if "TTTT" in seq:
        score -= 30
    elif "TTT" in seq:
        score -= 10

    # Penalise poly-G (synthesis and secondary-structure issues)
    if "GGGG" in seq:
        score -= 15

    # Reward G at position 20 (PAM-proximal)
    if seq[-1] == "G":
        score += 10

    # Penalise C at position 20
    if seq[-1] == "C":
        score -= 5

    # Reward purines at position 20
    if seq[-1] in "AG":
        score += 5

    # Self-complementarity check (simple hairpin heuristic)
    for k in range(4, 8):
        for j in range(len(seq) - k):
            frag = seq[j : j + k]
            rc_frag = reverse_complement(frag)
            if rc_frag in seq[j + k :]:
                score -= 5 * (k - 3)
                break

    return max(0.0, min(100.0, score))


def score_guides(guides: list[GuideRNA]) -> list[GuideRNA]:
    """Score a list of guides and sort descending by score."""
    for g in guides:
        g.score = score_guide(g)
    return sorted(guides, key=lambda g: g.score, reverse=True)
