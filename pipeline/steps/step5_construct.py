"""
Step 5 — Design therapeutic delivery constructs.

Generates AAV vector designs for each target gene, incorporating the
patient-specific guide RNAs selected in previous steps.

Cancer constructs:
  - CRISPRi (dCas9-KRAB) for oncogene silencing
  - CRISPRa (dCas9-VPR) for tumor-suppressor reactivation
  - Cas9 + HDR template for precise mutation correction

Mental-health constructs:
  - CRISPRi/CRISPRa for gene expression modulation
  - Base-editing for pharmacogenomic variant correction
"""

import json
import logging
from pathlib import Path

from ..config import AAV_PACKAGING_LIMIT_BP

log = logging.getLogger(__name__)

# Standard sgRNA scaffold for SpCas9
SGRNA_SCAFFOLD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAAC"
    "TTGAAAAAGTGGCACCGAGTCGGTGCTTTT"
)

# Promoter sizes (bp) for common neuron/tissue-specific promoters
PROMOTER_SIZES = {
    "CamKII": 1300,
    "hSyn1": 470,
    "CMV": 600,
    "EF1a": 1200,
    "CAG": 1700,
    "U6": 250,
    "TPH2": 1800,
    "OXT": 900,
    "TH": 1100,
    "GAD67": 1400,
    "ChAT": 1200,
}

# Effector sizes (bp)
EFFECTOR_SIZES = {
    "dCas9-KRAB": 4200,
    "dCas9-VPR": 5400,
    "SpCas9": 4100,
    "ABE8e": 4800,    # adenine base editor
    "CBE4max": 5200,  # cytosine base editor
}

# DREADD receptor sizes (bp)
DREADD_SIZES = {
    "hM4Di": 1400,    # inhibitory DREADD
    "hM3Dq": 1400,    # excitatory DREADD
}

# Reporter sizes (bp)
REPORTER_SIZES = {
    "IRES-mCherry": 1700,
    "IRES-GFP": 1600,
    "T2A-mCherry": 1000,
    "T2A-GFP": 900,
}


def _infer_strategy(target_info: dict) -> str:
    """
    Infer therapeutic strategy from all available database fields.

    Checks (in order): explicit ``strategy`` / ``crispr_strategy``,
    approach text from ``key_variants`` and ``clinical_trials``,
    and ``inheritance`` pattern to determine the correct CRISPR modality.
    """
    strategy = target_info.get("strategy", "") or target_info.get("crispr_strategy", "")
    if strategy:
        return strategy

    # Collect approach descriptions from key_variants and clinical_trials
    approaches = []
    for var in target_info.get("key_variants", []):
        a = var.get("approach", "")
        if a:
            approaches.append(a if isinstance(a, str) else " ".join(a))
    for trial in target_info.get("clinical_trials", []):
        a = trial.get("approach", "")
        if a:
            approaches.append(a if isinstance(a, str) else " ".join(a))

    for approach in approaches:
        approach_lower = approach.lower()
        if "base edit" in approach_lower or "base-edit" in approach_lower or "prime edit" in approach_lower:
            return "base-editing correction"
        if "knockout" in approach_lower or "disrupt" in approach_lower or "ablat" in approach_lower:
            return "knockout"
        if "gene addition" in approach_lower or "gene replacement" in approach_lower:
            return "HDR gene correction"
        if "crispri" in approach_lower or "silenc" in approach_lower:
            return "CRISPRi"
        if "crispra" in approach_lower or "activat" in approach_lower:
            return "CRISPRa"

    # Infer from inheritance pattern
    inheritance = target_info.get("inheritance", "").lower()
    if "gain-of-function" in inheritance or "gain of function" in inheritance:
        return "CRISPRi silencing of gain-of-function allele"
    if "dominant-negative" in inheritance or "dominant negative" in inheritance:
        return "knockout of dominant-negative allele"
    if "haploinsufficiency" in inheritance:
        return "CRISPRa to enhance expression from remaining allele"
    # AR, XLR = loss-of-function → correction
    if inheritance.startswith("ar") or "xlr" in inheritance or "xld" in inheritance:
        return "HDR correction of pathogenic variant"

    return ""


def _determine_construct_type(gene: str, target_info: dict) -> dict:
    """
    Decide the construct type based on the gene's therapeutic strategy.

    First checks the explicit ``strategy`` field, then infers from
    ``key_variants[].approach`` and ``inheritance`` if strategy is absent.
    """
    strategy = _infer_strategy(target_info)
    role = target_info.get("role", "")
    neuronal = "neuro" in role.lower()

    strategy_lower = strategy.lower()

    # Check explicit modality keywords before fuzzy text matches to avoid
    # e.g. "CRISPRa to reactivate silenced allele" matching CRISPRi on "silenc".
    if "CRISPRi" in strategy or ("silenc" in strategy_lower and "CRISPRa" not in strategy):
        return {
            "type": "CRISPRi",
            "effector": "dCas9-KRAB",
            "promoter": "CamKII" if neuronal else "CMV",
            "purpose": f"Transcriptional silencing of {gene}",
        }
    elif "CRISPRa" in strategy or "reactivat" in strategy_lower or "enhance" in strategy_lower:
        return {
            "type": "CRISPRa",
            "effector": "dCas9-VPR",
            "promoter": "CamKII" if neuronal else "CMV",
            "purpose": f"Transcriptional activation of {gene}",
        }
    # Check base editing BEFORE HDR — both may contain "correction"
    elif "base-edit" in strategy_lower or "base edit" in strategy_lower or "prime edit" in strategy_lower:
        return {
            "type": "base_editing",
            "effector": "ABE8e",
            "promoter": "CMV",
            "purpose": f"Base editing of pathogenic variants in {gene}",
        }
    elif "HDR" in strategy or "correction" in strategy_lower or "correct" in strategy_lower:
        return {
            "type": "HDR_correction",
            "effector": "SpCas9",
            "promoter": "CMV",
            "purpose": f"Precise mutation correction in {gene} via HDR",
        }
    elif "knockout" in strategy_lower or "disrupt" in strategy_lower:
        return {
            "type": "knockout",
            "effector": "SpCas9",
            "promoter": "CMV",
            "purpose": f"Gene disruption of {gene}",
        }
    else:
        # Default: HDR correction is safer than CRISPRi for unknown genes —
        # most rare-disease targets are loss-of-function.
        return {
            "type": "HDR_correction",
            "effector": "SpCas9",
            "promoter": "CMV",
            "purpose": f"Targeted correction of {gene}",
        }


def generate_construct(gene: str, target_info: dict,
                       guide_data: dict) -> dict:
    """
    Generate a complete construct map for one gene target.
    """
    ct = _determine_construct_type(gene, target_info)

    top_guides = guide_data.get("safe_guides") or guide_data.get("top_guides", [])
    selected_guides = top_guides[:2]  # up to 2 guides for redundancy

    elements = []

    # 5' ITR
    elements.append({"name": "5' ITR", "size_bp": 145, "type": "regulatory"})

    # Promoter
    prom = ct["promoter"]
    elements.append({
        "name": prom,
        "size_bp": PROMOTER_SIZES.get(prom, 800),
        "type": "promoter",
    })

    # Effector
    eff = ct["effector"]
    elements.append({
        "name": eff,
        "size_bp": EFFECTOR_SIZES.get(eff, 4200),
        "type": "effector",
    })

    # sgRNA cassettes
    for i, g in enumerate(selected_guides):
        cassette_seq = g.get("sequence", "") + SGRNA_SCAFFOLD
        elements.append({
            "name": f"U6::sgRNA_{gene}_guide{i + 1}",
            "size_bp": 250 + len(g.get("sequence", "")),
            "type": "sgRNA_cassette",
            "guide_sequence": g.get("sequence", ""),
            "guide_score": g.get("score", 0),
            "patient_specific": True,
        })

    # Poly-A
    elements.append({"name": "bGH polyA", "size_bp": 225, "type": "regulatory"})

    # 3' ITR
    elements.append({"name": "3' ITR", "size_bp": 145, "type": "regulatory"})

    total_bp = sum(e["size_bp"] for e in elements)

    construct = {
        "name": f"construct_{ct['type']}_{gene}",
        "gene": gene,
        "construct_type": ct["type"],
        "purpose": ct["purpose"],
        "backbone": "AAV9",
        "promoter": prom,
        "effector": eff,
        "guides_included": len(selected_guides),
        "elements": elements,
        "total_size_bp": total_bp,
        "total_size_kb": round(total_bp / 1000, 2),
    }

    if total_bp > AAV_PACKAGING_LIMIT_BP:
        construct["warning"] = (
            f"Construct is {total_bp} bp, exceeding AAV packaging limit "
            f"of ~{AAV_PACKAGING_LIMIT_BP} bp. Consider split-intein "
            f"dual-AAV approach or lentiviral delivery."
        )
        construct["recommendation"] = "DUAL_AAV_SPLIT_INTEIN"

    return construct


def generate_dreadd_construct(name: str, dreadd_info: dict,
                              safe_harbor: dict | None = None) -> dict:
    """
    Generate a DREADD chemogenetic construct map.

    DREADD constructs use Cas9 + HDR to insert the receptor transgene at a
    safe harbor locus, paired with a tissue-specific promoter.
    """
    receptor = dreadd_info["receptor"]
    reporter = dreadd_info.get("reporter", "IRES-GFP")
    promoter = dreadd_info.get("promoter", "CamKII")

    elements = []

    # 5' ITR
    elements.append({"name": "5' ITR", "size_bp": 145, "type": "regulatory"})

    # Promoter
    elements.append({
        "name": promoter,
        "size_bp": PROMOTER_SIZES.get(promoter, 800),
        "type": "promoter",
    })

    # DREADD receptor
    elements.append({
        "name": receptor,
        "size_bp": DREADD_SIZES.get(receptor, 1400),
        "type": "dreadd_receptor",
    })

    # Reporter
    elements.append({
        "name": reporter,
        "size_bp": REPORTER_SIZES.get(reporter, 1600),
        "type": "reporter",
    })

    # Poly-A
    elements.append({"name": "bGH polyA", "size_bp": 225, "type": "regulatory"})

    # 3' ITR
    elements.append({"name": "3' ITR", "size_bp": 145, "type": "regulatory"})

    total_bp = sum(e["size_bp"] for e in elements)

    construct = {
        "name": f"construct_{name}",
        "gene": name,
        "construct_type": "DREADD_insertion",
        "purpose": dreadd_info.get("purpose", f"Chemogenetic control via {receptor}"),
        "backbone": dreadd_info.get("backbone", "AAV9"),
        "promoter": promoter,
        "effector": receptor,
        "receptor": receptor,
        "reporter": reporter,
        "activation_ligand": dreadd_info.get("activation_ligand",
                                              "deschloroclozapine (DCZ)"),
        "mechanism": dreadd_info.get("mechanism", ""),
        "target_regions": dreadd_info.get("target_regions", []),
        "guides_included": 0,
        "elements": elements,
        "total_size_bp": total_bp,
        "total_size_kb": round(total_bp / 1000, 2),
    }

    if safe_harbor:
        construct["insertion_site"] = {
            "locus": safe_harbor.get("role", ""),
            "chrom": safe_harbor.get("chrom", ""),
            "start": safe_harbor.get("start", 0),
            "end": safe_harbor.get("end", 0),
        }
        construct["crispr_role"] = dreadd_info.get(
            "crispr_role",
            f"Cas9 HDR inserts transgene at safe harbor"
        )

    if total_bp > AAV_PACKAGING_LIMIT_BP:
        construct["warning"] = (
            f"Construct is {total_bp} bp, exceeding AAV packaging limit "
            f"of ~{AAV_PACKAGING_LIMIT_BP} bp. Consider split-intein "
            f"dual-AAV approach or lentiviral delivery."
        )
        construct["recommendation"] = "DUAL_AAV_SPLIT_INTEIN"

    return construct


def run(targets: dict, guide_designs: dict, output_dir: Path,
        patient_id: str,
        dreadd_constructs: dict | None = None,
        safe_harbors: dict | None = None) -> list[dict]:
    """Design constructs for all targets and write output."""
    log.info("[Step 5] Designing therapeutic constructs")

    constructs = []
    construct_dir = output_dir / "construct_maps"
    construct_dir.mkdir(exist_ok=True)

    # Standard CRISPR constructs for gene targets
    for gene, coords in targets.items():
        gd = guide_designs.get(gene, {})
        if gd.get("error"):
            log.warning("  Skipping %s: %s", gene, gd["error"])
            continue

        construct = generate_construct(gene, coords, gd)
        constructs.append(construct)

        cpath = construct_dir / f"{construct['name']}.json"
        with open(cpath, "w") as f:
            json.dump(construct, f, indent=2)

        warn = f" [!] {construct['warning']}" if "warning" in construct else ""
        log.info("  %s: %s (%d bp)%s",
                 gene, construct["construct_type"],
                 construct["total_size_bp"], warn)

    # DREADD chemogenetic constructs
    if dreadd_constructs:
        harbor_map = {
            "DREADD_inhibitory": (safe_harbors or {}).get("AAVS1"),
            "DREADD_excitatory": (safe_harbors or {}).get("ROSA26_human"),
        }
        for name, dreadd_info in dreadd_constructs.items():
            harbor = harbor_map.get(name)
            construct = generate_dreadd_construct(name, dreadd_info, harbor)
            constructs.append(construct)

            cpath = construct_dir / f"{construct['name']}.json"
            with open(cpath, "w") as f:
                json.dump(construct, f, indent=2)

            log.info("  %s: DREADD %s (%d bp) -> %s",
                     name, construct["receptor"],
                     construct["total_size_bp"],
                     ", ".join(construct.get("target_regions", [])))

    return constructs
