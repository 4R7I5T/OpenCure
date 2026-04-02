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


def extract_profile(vcf_path: str, targets: dict) -> dict:
    """
    Build a patient variant profile for every target gene.
    Uses multi-region extraction to avoid scanning the VCF multiple times.

    Returns {gene_name: {coding_variants, promoter_variants, ...}}.
    """
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
        coding = all_variants.get(f"{gene}_coding", [])
        promoter = all_variants.get(f"{gene}_promoter", [])

        profile[gene] = {
            "chrom": coords["chrom"],
            "start": coords["start"],
            "end": coords["end"],
            "role": coords.get("role", ""),
            "coding_variants": [_variant_to_dict(v) for v in coding],
            "promoter_variants": [_variant_to_dict(v) for v in promoter],
            "total_coding": len(coding),
            "total_promoter": len(promoter),
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
