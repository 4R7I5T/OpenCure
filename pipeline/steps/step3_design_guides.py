"""
Step 3 — Design personalized CRISPR guide RNAs.

For each target gene:
  1. Fetch or reconstruct the patient's sequence (reference + VCF variants)
  2. Find all PAM sites
  3. Score candidate guides
  4. Flag guides that span patient-specific variants
"""

import json
import logging
import os
from pathlib import Path

from ..config import GUIDE_LENGTH, REFERENCE_GENOME, FLANKING_BASES
from ..utils.guide_utils import find_pam_sites, score_guides, GuideRNA
from ..utils.sequence_utils import apply_variants_to_sequence

log = logging.getLogger(__name__)

TOP_GUIDES_PER_TARGET = 10


def _fetch_reference_sequence(chrom: str, start: int, end: int,
                              ref_path: str,
                              bam_path: str | None = None,
                              bai_path: str | None = None) -> str | None:
    """
    Try multiple strategies to obtain the reference sequence:
      1. Local reference FASTA (via pysam)
      2. BAM consensus (pileup from aligned reads)
      3. NCBI Entrez remote fetch
    """
    # Normalize chrom name variants to try
    chroms = [chrom, chrom.replace("chr", ""), "chr" + chrom.replace("chr", "")]

    # 1. Local reference FASTA
    if os.path.isfile(ref_path):
        try:
            import pysam
            fa = pysam.FastaFile(ref_path)
            for c in chroms:
                try:
                    seq = fa.fetch(c, start, end)
                    fa.close()
                    return seq
                except (KeyError, ValueError):
                    continue
            fa.close()
        except Exception as exc:
            log.debug("Failed to read reference FASTA %s: %s", ref_path, exc)

    # 2. BAM consensus — reconstruct from aligned reads
    if bam_path and os.path.isfile(bam_path):
        seq = _consensus_from_bam(bam_path, bai_path, chroms, start, end)
        if seq:
            return seq

    # 3. NCBI Entrez remote fetch
    try:
        from Bio import Entrez, SeqIO
        ncbi_email = os.getenv("NCBI_EMAIL")
        ncbi_key = os.getenv("NCBI_API_KEY")
        if not ncbi_email:
            log.warning("No local ref, no BAM coverage, no NCBI_EMAIL — "
                        "set NCBI_EMAIL env var to enable remote sequence fetch")
            return None
        Entrez.email = ncbi_email
        if ncbi_key:
            Entrez.api_key = ncbi_key
        # Use GRCh38 primary assembly accessions
        chrom_accessions = {
            "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
            "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
            "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
            "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
            "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
            "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
            "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
            "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
        }
        bare = chrom.replace("chr", "")
        acc = chrom_accessions.get(bare, bare)
        handle = Entrez.efetch(
            db="nucleotide", id=acc,
            seq_start=start, seq_stop=end,
            rettype="fasta", retmode="text",
        )
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)
    except Exception as exc:
        log.warning("Failed to fetch sequence for %s:%d-%d: %s",
                    chrom, start, end, exc)
        return None


def _find_bai(bam_path: str, bai_path: str | None) -> str | None:
    """Locate the BAI index file for a BAM."""
    candidates = []
    if bai_path:
        candidates.append(bai_path)
    candidates.extend([
        bam_path + ".bai",
        bam_path.replace(".bam", ".bai"),
        bam_path.replace(".bam", ".bam.bai"),
    ])
    # Also check for BAI in same directory with different prefix
    bam_dir = os.path.dirname(bam_path)
    if bam_dir:
        import glob
        candidates.extend(glob.glob(os.path.join(bam_dir, "*.bai")))

    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def _consensus_from_bam(bam_path: str, bai_path: str | None,
                        chroms: list[str],
                        start: int, end: int) -> str | None:
    """Build a consensus sequence from BAM pileup at a region."""
    try:
        import pysam

        idx = _find_bai(bam_path, bai_path)
        if not idx:
            log.debug("No BAI index found for %s", bam_path)
            return None

        bam = pysam.AlignmentFile(bam_path, "rb", index_filename=idx)

        for c in chroms:
            try:
                # Check if there are reads in this region first
                read_count = bam.count(c, start, end)
                if read_count == 0:
                    continue

                counts = bam.count_coverage(c, start, end,
                                            quality_threshold=20)
                bases = "ACGT"
                seq = []
                for i in range(len(counts[0])):
                    base_counts = [counts[j][i] for j in range(4)]
                    total = sum(base_counts)
                    if total == 0:
                        seq.append("N")
                    else:
                        best = max(range(4), key=lambda x: base_counts[x])
                        seq.append(bases[best])
                result = "".join(seq)
                n_count = result.count("N")
                # Accept if at least 50% covered
                if n_count < len(result) * 0.5:
                    log.info("  Extracted %d bp consensus from BAM for %s:%d-%d "
                             "(%d uncovered, %d reads)",
                             len(result), c, start, end, n_count, read_count)
                    bam.close()
                    return result
                else:
                    log.debug("  BAM coverage too sparse at %s:%d-%d "
                              "(%d/%d uncovered)", c, start, end,
                              n_count, len(result))
            except (ValueError, KeyError):
                continue

        bam.close()
    except Exception as exc:
        log.debug("BAM consensus extraction failed: %s", exc)
    return None


def design_for_gene(gene: str, coords: dict, profile_entry: dict,
                    ref_path: str, bam_path: str | None = None,
                    bai_path: str | None = None) -> dict:
    """Design guides for a single gene target."""
    chrom = coords["chrom"]
    # Target the promoter region for CRISPRi/CRISPRa (upstream of gene start)
    target_start = max(1, coords["start"] - 1500)
    target_end = coords["start"] + 500

    # Get reference sequence
    ref_seq = _fetch_reference_sequence(chrom, target_start, target_end,
                                        ref_path, bam_path=bam_path,
                                        bai_path=bai_path)
    if not ref_seq:
        return {
            "gene": gene,
            "error": f"Could not obtain sequence for {chrom}:{target_start}-{target_end}",
            "top_guides": [],
        }

    # Trim sequence to expected region size (NCBI may return more)
    expected_len = target_end - target_start
    if len(ref_seq) > expected_len * 2:
        log.debug("  Trimming fetched sequence from %d to %d bp",
                  len(ref_seq), expected_len)
        ref_seq = ref_seq[:expected_len]

    # Build patient-specific sequence by applying variants
    from ..utils.vcf_utils import Variant
    variants = []
    for vd in (profile_entry.get("promoter_variants", []) +
               profile_entry.get("coding_variants", [])):
        variants.append(Variant(
            chrom=chrom, pos=vd["pos"], ref=vd["ref"], alt=vd["alt"],
            qual=vd["qual"], filter_status=vd["filter"],
            var_type=vd["type"], zygosity=vd["zygosity"],
        ))

    # Build both allele sequences for allele-specific guide design
    # Alt-all: applies both HOM_ALT and HET variants (targets mutant allele)
    patient_seq_alt = apply_variants_to_sequence(
        ref_seq, variants, target_start, allele="alt_all")
    # Ref: unmodified reference (targets wild-type allele)
    patient_seq_ref = apply_variants_to_sequence(
        ref_seq, variants, target_start, allele="ref")

    # Design guides on mutant allele sequence (primary)
    raw_guides = find_pam_sites(patient_seq_alt, guide_length=GUIDE_LENGTH,
                                region_start=target_start)

    # Also find reference-allele guides for comparison
    ref_guides = find_pam_sites(patient_seq_ref, guide_length=GUIDE_LENGTH,
                                region_start=target_start)
    ref_guide_seqs = {g.sequence for g in ref_guides}

    # Score and rank (applies Synthego manufacturing filters)
    scored = score_guides(raw_guides)
    top = scored[:TOP_GUIDES_PER_TARGET]

    # Annotate guides
    for guide in top:
        guide.gene = gene
        guide.target_name = f"{gene}_promoter"
        # Check if guide spans a patient variant (allele-specific)
        for vd in variants:
            rel = vd.pos - target_start
            if guide.position <= rel <= guide.position + GUIDE_LENGTH:
                guide.spans_variant = True
                guide.variant_note = (
                    f"Spans patient variant at {vd.pos} "
                    f"({vd.ref}>{vd.alt}, {vd.zygosity}). Guide designed "
                    f"against patient allele — allele-specific."
                )

    synthego_ok = sum(1 for g in top if g.synthego_compatible)

    return {
        "gene": gene,
        "region": f"{chrom}:{target_start}-{target_end}",
        "purpose": coords.get("strategy", coords.get("role", "")),
        "patient_variants_in_region": len(variants),
        "total_pam_sites": len(raw_guides),
        "synthego_compatible_count": synthego_ok,
        "top_guides": [_guide_to_dict(g) for g in top],
    }


def _guide_to_dict(g: GuideRNA) -> dict:
    return {
        "sequence": g.sequence,
        "synthego_order_seq": g.synthego_order_sequence,
        "pam": g.pam,
        "strand": g.strand,
        "position": g.position,
        "genomic_start": g.genomic_start,
        "gc_content": round(g.gc, 3),
        "seed_gc_content": round(g.seed_gc, 3),
        "score": round(g.score, 1),
        "full_target": g.full_target,
        "synthego_compatible": g.synthego_compatible,
        "synthego_rejection_reason": g.synthego_rejection_reason,
        "spans_variant": g.spans_variant,
        "variant_note": g.variant_note,
        "offtarget_status": g.offtarget_status,
    }


def run(targets: dict, profile: dict, output_dir: Path,
        patient_id: str, bam_path: str | None = None,
        bai_path: str | None = None) -> dict:
    """Design guides for all target genes and write output."""
    log.info("[Step 3] Designing personalized guide RNAs")

    ref_path = REFERENCE_GENOME
    all_designs: dict = {}

    for gene, coords in targets.items():
        profile_entry = profile.get(gene, {})
        design = design_for_gene(gene, coords, profile_entry, ref_path,
                                 bam_path=bam_path, bai_path=bai_path)
        all_designs[gene] = design

        n_guides = len(design.get("top_guides", []))
        top_score = (design["top_guides"][0]["score"]
                     if design.get("top_guides") else 0)
        log.info("  %s: %d PAM sites, top %d guides (best score %.1f)",
                 gene, design.get("total_pam_sites", 0), n_guides, top_score)

    out_path = output_dir / f"{patient_id}_guide_designs.json"
    with open(out_path, "w") as f:
        json.dump(all_designs, f, indent=2)
    log.info("  Wrote %s", out_path)

    return all_designs
