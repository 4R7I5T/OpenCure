#!/usr/bin/env python3
"""
Export OpenCure's curated disease–gene–CRISPR knowledge as LLM training data.

Generates instruction-tuning examples in JSONL (one JSON object per line)
from the 780+ gene targets across 60 diagnosis categories.  Each gene
produces multiple training examples covering:

  1. Gene function & disease association
  2. Variant interpretation
  3. CRISPR strategy selection & rationale
  4. Guide design constraints
  5. Construct/vector design
  6. Clinical trial status
  7. Differential diagnosis (where applicable)

Output formats
--------------
- ``alpaca``  — ``{"instruction", "input", "output"}`` triples
                 (compatible with Stanford Alpaca, Axolotl, LLaMA-Factory)
- ``sharegpt`` — ``{"conversations": [{"from": "human", ...}, ...]}``
                 (compatible with ShareGPT, OpenAI fine-tuning)

Usage::

    python -m pipeline.training.export_training_data \\
        --format alpaca --output training_data.jsonl

    python -m pipeline.training.export_training_data \\
        --format sharegpt --output training_data_sharegpt.jsonl
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Iterator

# ---------------------------------------------------------------------------
# Lazy import of all target databases (same as run_pipeline.py)
# ---------------------------------------------------------------------------

def _load_all_targets() -> dict[str, dict]:
    """Return {diagnosis: targets_dict} for every diagnosis."""
    from ..run_pipeline import _TARGETS_MAP, select_targets
    return {dx: select_targets(dx) for dx in _TARGETS_MAP}


# ---------------------------------------------------------------------------
# Example generators
# ---------------------------------------------------------------------------

def _str_field(info: dict, *keys: str) -> str:
    """Return the first non-empty string value for *keys*, or ''."""
    for k in keys:
        v = info.get(k, "")
        if isinstance(v, tuple):
            v = " ".join(v)
        if v:
            return str(v)
    return ""


def _variant_summary(info: dict) -> list[str]:
    """Summarise key_variants into short strings."""
    summaries = []
    for v in info.get("key_variants", []):
        if isinstance(v, str):
            summaries.append(v)
            continue
        name = v.get("name", "")
        rsid = v.get("rsid", "")
        conseq = v.get("consequence", "")
        sig = v.get("clinical_significance", "")
        parts = [p for p in [name, rsid, conseq, sig] if p]
        if parts:
            summaries.append(" — ".join(parts))
    return summaries


def _trial_summary(info: dict) -> list[str]:
    """Summarise clinical_trials into short strings."""
    summaries = []
    for t in info.get("clinical_trials", []):
        nct = t.get("nct", "")
        phase = t.get("phase", "")
        approach = t.get("approach", "")
        if isinstance(approach, tuple):
            approach = " ".join(approach)
        parts = [p for p in [nct, f"Phase {phase}" if phase else "", approach] if p]
        if parts:
            summaries.append("; ".join(parts))
    return summaries


# ---- Individual example generators ----

def gene_function_examples(gene: str, info: dict, diagnosis: str) -> Iterator[dict]:
    """Q&A about gene function and disease role."""
    role = _str_field(info, "role")
    disease = _str_field(info, "disease")
    conditions = info.get("conditions", [])
    inheritance = _str_field(info, "inheritance")

    if role:
        yield {
            "instruction": f"What is the function of the {gene} gene and what diseases is it associated with?",
            "input": "",
            "output": (
                f"{gene} encodes {role} "
                f"{'It is associated with ' + disease + '. ' if disease else ''}"
                f"{'Inheritance pattern: ' + inheritance + '. ' if inheritance else ''}"
                f"{'Relevant conditions: ' + ', '.join(conditions) + '.' if conditions else ''}"
            ).strip(),
        }

    if disease and inheritance:
        yield {
            "instruction": f"What is the inheritance pattern of {disease}?",
            "input": f"Gene: {gene}",
            "output": f"{disease} caused by variants in {gene} follows {inheritance} inheritance.",
        }


def variant_examples(gene: str, info: dict, diagnosis: str) -> Iterator[dict]:
    """Q&A about known pathogenic variants."""
    variants = _variant_summary(info)
    disease = _str_field(info, "disease")

    if variants:
        yield {
            "instruction": f"What are the key pathogenic variants in {gene}?",
            "input": f"Disease context: {disease}" if disease else "",
            "output": (
                f"Key pathogenic variants in {gene}:\n"
                + "\n".join(f"- {v}" for v in variants[:10])
            ),
        }

    # Per-variant Q&A for detailed entries
    for v in info.get("key_variants", []):
        if not isinstance(v, dict):
            continue
        notes = v.get("notes", "")
        if isinstance(notes, tuple):
            notes = " ".join(notes)
        name = v.get("name", "")
        if name and notes and len(notes) > 50:
            yield {
                "instruction": f"Describe the clinical significance of the {gene} variant {name}.",
                "input": "",
                "output": notes.strip(),
            }


def strategy_examples(gene: str, info: dict, diagnosis: str) -> Iterator[dict]:
    """Q&A about CRISPR therapeutic strategy."""
    strategy = _str_field(info, "strategy", "crispr_strategy")
    role = _str_field(info, "role")
    disease = _str_field(info, "disease")
    inheritance = _str_field(info, "inheritance")

    if strategy:
        yield {
            "instruction": (
                f"Design a CRISPR gene therapy strategy for {gene}-associated disease."
            ),
            "input": (
                f"Gene: {gene}\n"
                f"Disease: {disease}\n"
                f"Inheritance: {inheritance}\n"
                f"Gene function: {role[:200]}"
            ).strip(),
            "output": strategy.strip(),
        }

    # Why this modality?
    from ..steps.step5_construct import _infer_strategy, _determine_construct_type
    inferred = _infer_strategy(info)
    ct = _determine_construct_type(gene, info)
    if inferred:
        yield {
            "instruction": (
                f"What CRISPR modality should be used for {gene} and why?"
            ),
            "input": f"Inheritance: {inheritance}" if inheritance else "",
            "output": (
                f"For {gene}, the recommended modality is {ct['type']} "
                f"({ct['effector']}). {ct['purpose']}. "
                f"Rationale: {inferred}"
            ),
        }


def construct_examples(gene: str, info: dict, diagnosis: str) -> Iterator[dict]:
    """Q&A about AAV construct design."""
    from ..steps.step5_construct import _determine_construct_type, PROMOTER_SIZES, EFFECTOR_SIZES
    from ..config import AAV_PACKAGING_LIMIT_BP

    ct = _determine_construct_type(gene, info)
    role = _str_field(info, "role")
    neuronal = "neuro" in role.lower()

    prom_size = PROMOTER_SIZES.get(ct["promoter"], 800)
    eff_size = EFFECTOR_SIZES.get(ct["effector"], 4200)
    # 5'ITR + promoter + effector + 2 sgRNA cassettes + polyA + 3'ITR
    est_size = 145 + prom_size + eff_size + 2 * 270 + 225 + 145
    exceeds = est_size > AAV_PACKAGING_LIMIT_BP

    yield {
        "instruction": f"Design an AAV delivery construct for CRISPR therapy targeting {gene}.",
        "input": (
            f"Construct type: {ct['type']}\n"
            f"Effector: {ct['effector']}\n"
            f"Tissue: {'neuronal' if neuronal else 'systemic'}"
        ),
        "output": (
            f"Recommended construct for {gene} ({ct['type']}):\n"
            f"- Backbone: AAV9\n"
            f"- Promoter: {ct['promoter']} ({prom_size} bp)\n"
            f"- Effector: {ct['effector']} ({eff_size} bp)\n"
            f"- Includes 2 sgRNA cassettes (~540 bp total)\n"
            f"- Regulatory: 5'/3' ITR (290 bp) + bGH polyA (225 bp)\n"
            f"- Estimated total: ~{est_size} bp\n"
            + (f"- WARNING: Exceeds AAV packaging limit ({AAV_PACKAGING_LIMIT_BP} bp). "
               f"Consider dual-AAV split-intein approach.\n" if exceeds else "")
            + f"- Purpose: {ct['purpose']}"
        ),
    }


def guide_design_examples(gene: str, info: dict, diagnosis: str) -> Iterator[dict]:
    """Q&A about guide RNA design constraints."""
    from ..config import (GUIDE_LENGTH, PAM_SEQUENCE, SYNTHEGO_MIN_GC,
                          SYNTHEGO_MAX_GC, SYNTHEGO_MAX_HOMOPOLYMER)

    chrom = info.get("chrom", "")
    start = info.get("start", 0)

    yield {
        "instruction": (
            f"What are the design constraints for CRISPR guide RNAs targeting {gene}?"
        ),
        "input": f"Genomic locus: {chrom}:{start}",
        "output": (
            f"Guide RNA design constraints for {gene} at {chrom}:{start}:\n"
            f"- Guide length: {GUIDE_LENGTH} nt\n"
            f"- PAM: {PAM_SEQUENCE} (SpCas9)\n"
            f"- GC content: {SYNTHEGO_MIN_GC*100:.0f}–{SYNTHEGO_MAX_GC*100:.0f}%\n"
            f"- No homopolymer runs ≥{SYNTHEGO_MAX_HOMOPOLYMER} nt\n"
            f"- No TTTT (Pol III terminator signal)\n"
            f"- No BsmBI (CGTCTC) or BsaI (GGTCTC) restriction sites\n"
            f"- Target region: 1500 bp upstream to 500 bp downstream of gene start\n"
            f"- Guides are designed against the patient's variant-bearing allele\n"
            f"- Top 10 guides ranked by Doench-inspired positional scoring (0–100)\n"
            f"- Off-target analysis via Cas-OFFinder or heuristic seed assessment"
        ),
    }


def clinical_trial_examples(gene: str, info: dict, diagnosis: str) -> Iterator[dict]:
    """Q&A about clinical trials."""
    trials = _trial_summary(info)
    disease = _str_field(info, "disease")
    programs = _str_field(info, "clinical_programs")

    if trials:
        yield {
            "instruction": f"What clinical trials exist for CRISPR therapy targeting {gene}?",
            "input": f"Disease: {disease}" if disease else "",
            "output": (
                f"Clinical trials for {gene}-targeted therapy:\n"
                + "\n".join(f"- {t}" for t in trials)
                + (f"\n\nAdditional context: {programs}" if programs else "")
            ),
        }


def diagnosis_routing_examples(diagnosis: str, targets: dict) -> Iterator[dict]:
    """Q&A about which genes to target for a diagnosis."""
    gene_list = sorted(targets.keys())
    if not gene_list:
        return

    yield {
        "instruction": f"Which genes should be targeted for CRISPR therapy in {diagnosis.replace('_', ' ')} patients?",
        "input": "",
        "output": (
            f"For {diagnosis.replace('_', ' ')}, the OpenCure pipeline targets "
            f"{len(gene_list)} genes: {', '.join(gene_list)}. "
            f"Each gene is analyzed for patient-specific variants, and personalized "
            f"guide RNAs are designed against the mutant allele."
        ),
    }


# ---------------------------------------------------------------------------
# Format converters
# ---------------------------------------------------------------------------

def to_alpaca(example: dict) -> dict:
    """Return Alpaca-format dict."""
    return {
        "instruction": example["instruction"],
        "input": example.get("input", ""),
        "output": example["output"],
    }


def to_sharegpt(example: dict) -> dict:
    """Return ShareGPT-format dict."""
    human_msg = example["instruction"]
    if example.get("input"):
        human_msg += "\n\n" + example["input"]
    return {
        "conversations": [
            {"from": "human", "value": human_msg},
            {"from": "gpt", "value": example["output"]},
        ]
    }


# ---------------------------------------------------------------------------
# Main export
# ---------------------------------------------------------------------------

GENERATORS = [
    gene_function_examples,
    variant_examples,
    strategy_examples,
    construct_examples,
    guide_design_examples,
    clinical_trial_examples,
]


def generate_all() -> Iterator[dict]:
    """Yield all training examples from all databases."""
    all_targets = _load_all_targets()

    for diagnosis, targets in sorted(all_targets.items()):
        # Diagnosis-level routing example
        yield from diagnosis_routing_examples(diagnosis, targets)

        # Per-gene examples
        for gene, info in sorted(targets.items()):
            for gen_fn in GENERATORS:
                yield from gen_fn(gene, info, diagnosis)


def export(output_path: str, fmt: str = "alpaca") -> int:
    """Export training data to JSONL file. Returns example count."""
    converter = to_alpaca if fmt == "alpaca" else to_sharegpt
    count = 0

    with open(output_path, "w") as f:
        for example in generate_all():
            f.write(json.dumps(converter(example), ensure_ascii=False) + "\n")
            count += 1

    return count


def main():
    parser = argparse.ArgumentParser(
        description="Export OpenCure knowledge as LLM training data",
    )
    parser.add_argument(
        "--format", choices=["alpaca", "sharegpt"], default="alpaca",
        help="Output format (default: alpaca)",
    )
    parser.add_argument(
        "--output", "-o", default="opencure_training_data.jsonl",
        help="Output JSONL file path",
    )
    args = parser.parse_args()

    count = export(args.output, args.format)
    print(f"Exported {count} training examples to {args.output} ({args.format} format)")


if __name__ == "__main__":
    main()
