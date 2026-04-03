"""
Step 2 — Extract patient-specific variants at target gene loci.

Reads the patient VCF and pulls all variants overlapping each target gene,
including promoter regions (2 kb upstream of gene start). Produces a
per-gene patient profile and a targeting strategy.
"""

import json
import logging
from pathlib import Path

from ..utils.vcf_utils import extract_variants_multi_region

log = logging.getLogger(__name__)

PROMOTER_UPSTREAM = 2000  # bp upstream of gene start


def _passes_30x_quality(v, min_dp: int, min_ad: int,
                        min_qual: float, min_gq: int) -> bool:
    """
    Check if a variant passes quality thresholds appropriate for 30x WGS.

    At 30x mean coverage, heterozygous variants require sufficient depth
    and alt-allele support to distinguish from sequencing errors.
    """
    if v.qual < min_qual:
        return False
    # If depth fields are populated (non-zero), apply filters
    if v.depth > 0 and v.depth < min_dp:
        return False
    if v.alt_depth > 0 and v.alt_depth < min_ad:
        return False
    if v.genotype_qual > 0 and v.genotype_qual < min_gq:
        return False
    return True


def extract_profile(vcf_path: str, targets: dict) -> dict:
    """
    Build a patient variant profile for every target gene.
    Uses multi-region extraction to avoid scanning the VCF multiple times.

    Applies 30x WGS quality filters to ensure variant calls used for
    guide design are reliable at typical clinical WGS depths.

    Returns {gene_name: {coding_variants, promoter_variants, ...}}.
    """
    from ..config import (MIN_VARIANT_DEPTH, MIN_VARIANT_ALT_DEPTH,
                          MIN_VARIANT_QUAL, MIN_GENOTYPE_QUAL)

    # Build all regions to query in one batch
    regions = {}
    for gene, coords in targets.items():
        chrom = coords["chrom"]
        start = coords["start"]
        end = coords["end"]
        prom_start = max(1, start - PROMOTER_UPSTREAM)
        regions[f"{gene}_coding"] = (chrom, start, end)
        regions[f"{gene}_promoter"] = (chrom, prom_start, start - 1)

    # Single extraction pass for all regions
    all_variants = extract_variants_multi_region(vcf_path, regions,
                                                 pass_only=True)

    profile: dict = {}
    for gene, coords in targets.items():
        raw_coding = all_variants.get(f"{gene}_coding", [])
        raw_promoter = all_variants.get(f"{gene}_promoter", [])

        # Apply 30x quality filters
        coding = [v for v in raw_coding
                  if _passes_30x_quality(v, MIN_VARIANT_DEPTH,
                                         MIN_VARIANT_ALT_DEPTH,
                                         MIN_VARIANT_QUAL,
                                         MIN_GENOTYPE_QUAL)]
        promoter = [v for v in raw_promoter
                    if _passes_30x_quality(v, MIN_VARIANT_DEPTH,
                                           MIN_VARIANT_ALT_DEPTH,
                                           MIN_VARIANT_QUAL,
                                           MIN_GENOTYPE_QUAL)]

        filtered_coding = len(raw_coding) - len(coding)
        filtered_promoter = len(raw_promoter) - len(promoter)

        profile[gene] = {
            "chrom": coords["chrom"],
            "start": coords["start"],
            "end": coords["end"],
            "role": coords.get("role", ""),
            "coding_variants": [_variant_to_dict(v) for v in coding],
            "promoter_variants": [_variant_to_dict(v) for v in promoter],
            "total_coding": len(coding),
            "total_promoter": len(promoter),
            "filtered_low_quality_coding": filtered_coding,
            "filtered_low_quality_promoter": filtered_promoter,
        }

    return profile


def _variant_to_dict(v) -> dict:
    return {
        "pos": v.pos,
        "ref": v.ref,
        "alt": v.alt,
        "qual": v.qual,
        "zygosity": v.zygosity,
        "type": v.var_type,
        "filter": v.filter_status,
        "depth": v.depth,
        "alt_depth": v.alt_depth,
        "genotype_qual": v.genotype_qual,
    }


def assess_strategy(profile: dict, targets: dict) -> dict:
    """
    Determine targeting approach per gene based on patient genotype.
    """
    strategies: dict = {}

    for gene, data in profile.items():
        target_info = targets.get(gene, {})
        strategy_text = target_info.get("strategy", "Standard guide design")

        approach = "REFERENCE_MATCHED"
        notes = strategy_text

        if data["total_promoter"] > 0:
            approach = "ALLELE_SPECIFIC"
            notes += (
                f" | Patient has {data['total_promoter']} promoter variant(s) "
                f"— allele-specific guide design required."
            )

        func_vars = [
            v for v in data["coding_variants"]
            if v["type"] in ("snp", "mnp", "missense", "missense_hotspot")
            and v["qual"] > 30
        ]
        if func_vars:
            notes += (
                f" | {len(func_vars)} coding variant(s) may already affect "
                f"protein function."
            )

        strategies[gene] = {
            "gene": gene,
            "role": data["role"],
            "approach": approach,
            "strategy": strategy_text,
            "notes": notes,
            "coding_variant_count": data["total_coding"],
            "promoter_variant_count": data["total_promoter"],
        }

    return strategies


def run(vcf_path: str, targets: dict, output_dir: Path,
        patient_id: str) -> tuple[dict, dict]:
    """Execute step 2 and write outputs."""
    log.info("[Step 2] Extracting patient variants at target loci")

    profile = extract_profile(vcf_path, targets)
    strategy = assess_strategy(profile, targets)

    profile_path = output_dir / f"{patient_id}_gene_profile.json"
    strategy_path = output_dir / f"{patient_id}_targeting_strategy.json"

    with open(profile_path, "w") as f:
        json.dump(profile, f, indent=2, default=str)
    with open(strategy_path, "w") as f:
        json.dump(strategy, f, indent=2, default=str)

    log.info("  Wrote %s", profile_path)
    log.info("  Wrote %s", strategy_path)

    return profile, strategy
