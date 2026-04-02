#!/usr/bin/env python3
"""
OpenCure Pipeline v2 — Personalized CRISPR therapy design.

Usage:
    python -m pipeline.run_pipeline \\
        --vcf Genome/AR6GE3BF6RQ_vcf.vcf.gz \\
        --bam Genome/AR6GE3BF6RQ_bam.bam \\
        --bai Genome/AR6GE3BF6RQ_bai.bam.bai \\
        --patient PATIENT_001 \\
        --diagnosis cancer

    python -m pipeline.run_pipeline \\
        --vcf patient.vcf.gz \\
        --bam patient.bam \\
        --patient PATIENT_002 \\
        --diagnosis mental_health \\
        --conditions depression anxiety

    python -m pipeline.run_pipeline \\
        --vcf patient.vcf.gz \\
        --bam patient.bam \\
        --patient PATIENT_003 \\
        --diagnosis neuropsychiatric \\
        --conditions schizophrenia cognitive_impairment

    python -m pipeline.run_pipeline \\
        --vcf patient.vcf.gz \\
        --bam patient.bam \\
        --patient PATIENT_004 \\
        --diagnosis autoimmune \\
        --conditions yao_syndrome
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

from .config import (DIAGNOSIS_CANCER, DIAGNOSIS_MENTAL_HEALTH,
                     DIAGNOSIS_NEUROPSYCH, DIAGNOSIS_AUTOIMMUNE,
                     DIAGNOSIS_HEMATOLOGICAL, DIAGNOSIS_IMMUNODEFICIENCY,
                     DIAGNOSIS_HEREDITARY, DIAGNOSIS_CARDIAC,
                     DIAGNOSIS_RARE_DISEASE, DIAGNOSIS_INFECTIOUS,
                     DIAGNOSIS_IMMUNE_DYSREGULATION, DIAGNOSIS_TELOMERE_AGING,
                     DIAGNOSIS_PULMONARY, DIAGNOSIS_HEPATIC,
                     DIAGNOSIS_SKELETAL, DIAGNOSIS_DERMATOLOGIC,
                     DIAGNOSIS_MITOCHONDRIAL, DIAGNOSIS_COMPLEMENT,
                     DIAGNOSIS_COAGULATION, DIAGNOSIS_PAIN,
                     DIAGNOSIS_EPILEPSY, DIAGNOSIS_CONGENITAL_HEART,
                     DIAGNOSIS_MONOGENIC_DIABETES, DIAGNOSIS_OBESITY,
                     DIAGNOSIS_OPHTHALMOLOGY, DIAGNOSIS_NCL,
                     DIAGNOSIS_BONE_MARROW_FAILURE, DIAGNOSIS_IMMUNE_CHECKPOINT,
                     DIAGNOSIS_CRANIOFACIAL, VALID_DIAGNOSES)
from .databases.cancer_genes import CANCER_TARGETS
from .databases.therapeutic_targets import MENTAL_HEALTH_TARGETS
from .databases.neuropsych_targets import (NEUROPSYCH_TARGETS,
                                           DREADD_CONSTRUCTS,
                                           SAFE_HARBOR_LOCI)
from .databases.autoimmune_targets import AUTOIMMUNE_TARGETS
from .databases.hematological_targets import ALL_HEMATOLOGICAL_TARGETS
from .databases.immunodeficiency_metabolic_targets import ALL_GENE_THERAPY_TARGETS
from .databases.hereditary_disease_targets import ALL_HEREDITARY_TARGETS
from .databases.cardiac_targets import ALL_CARDIAC_TARGETS
from .databases.lysosomal_kidney_endocrine_hearing_connective_targets import (
    ALL_TARGETS as ALL_RARE_DISEASE_TARGETS,
)
from .databases.additional_rare_disease_targets import (
    ALL_ADDITIONAL_RARE_TARGETS,
)
from .databases.infectious_immune_aging_targets import (
    INFECTIOUS_DISEASE_TARGETS,
    IMMUNE_DYSREGULATION_TARGETS,
    TELOMERE_AGING_TARGETS,
    PULMONARY_TARGETS,
    HEPATIC_TARGETS,
    SKELETAL_TARGETS,
    DERMATOLOGIC_TARGETS,
)
from .databases.remaining_disease_targets import (
    MITOCHONDRIAL_TARGETS,
    COMPLEMENT_TARGETS,
    COAGULATION_TARGETS,
    PAIN_TARGETS,
    EPILEPSY_TARGETS,
    CONGENITAL_HEART_TARGETS,
)
from .databases.expanded_disease_targets import (
    MONOGENIC_DIABETES_TARGETS,
    OBESITY_TARGETS,
    OPHTHALMOLOGY_TARGETS,
    NCL_TARGETS,
    BONE_MARROW_FAILURE_TARGETS,
    IMMUNE_CHECKPOINT_TARGETS,
    CRANIOFACIAL_TARGETS,
)

from .steps import (
    step1_validate,
    step2_extract_variants,
    step2b_vep_annotation,
    step2c_autoinflammatory_annotation,
    step3_design_guides,
    step4_offtarget,
    step5_construct,
    step6_report,
)

log = logging.getLogger("opencure")


def select_targets(diagnosis: str, conditions: list[str] | None = None) -> dict:
    """
    Return the gene-target dict appropriate for the diagnosis.
    Optionally filter mental-health / neuropsych targets by condition.
    """
    if diagnosis == DIAGNOSIS_CANCER:
        return dict(CANCER_TARGETS)

    if diagnosis == DIAGNOSIS_MENTAL_HEALTH:
        targets = dict(MENTAL_HEALTH_TARGETS)
        if conditions:
            filtered = {}
            for gene, info in targets.items():
                gene_conditions = info.get("conditions", [])
                if any(c in gene_conditions for c in conditions):
                    filtered[gene] = info
            return filtered if filtered else targets
        return targets

    if diagnosis == DIAGNOSIS_NEUROPSYCH:
        targets = dict(NEUROPSYCH_TARGETS)
        if conditions:
            filtered = {}
            for gene, info in targets.items():
                gene_conditions = info.get("conditions", [])
                if any(c in gene_conditions for c in conditions):
                    filtered[gene] = info
            return filtered if filtered else targets
        return targets

    if diagnosis == DIAGNOSIS_AUTOIMMUNE:
        targets = dict(AUTOIMMUNE_TARGETS)
        if conditions:
            filtered = {}
            for gene, info in targets.items():
                gene_conditions = info.get("conditions", [])
                if any(c in gene_conditions for c in conditions):
                    filtered[gene] = info
            return filtered if filtered else targets
        return targets

    if diagnosis == DIAGNOSIS_HEMATOLOGICAL:
        targets = dict(ALL_HEMATOLOGICAL_TARGETS)
        if conditions:
            filtered = {}
            for gene, info in targets.items():
                gene_conditions = info.get("conditions", [])
                if any(c in gene_conditions for c in conditions):
                    filtered[gene] = info
            return filtered if filtered else targets
        return targets

    if diagnosis == DIAGNOSIS_IMMUNODEFICIENCY:
        targets = dict(ALL_GENE_THERAPY_TARGETS)
        if conditions:
            filtered = {}
            for gene, info in targets.items():
                gene_conditions = info.get("conditions", [])
                if any(c in gene_conditions for c in conditions):
                    filtered[gene] = info
            return filtered if filtered else targets
        return targets

    if diagnosis == DIAGNOSIS_HEREDITARY:
        targets = dict(ALL_HEREDITARY_TARGETS)
        if conditions:
            filtered = {}
            for gene, info in targets.items():
                gene_conditions = info.get("conditions", [])
                if any(c in gene_conditions for c in conditions):
                    filtered[gene] = info
            return filtered if filtered else targets
        return targets

    if diagnosis == DIAGNOSIS_CARDIAC:
        targets = dict(ALL_CARDIAC_TARGETS)
        if conditions:
            filtered = {}
            for gene, info in targets.items():
                gene_conditions = info.get("conditions", [])
                if any(c in gene_conditions for c in conditions):
                    filtered[gene] = info
            return filtered if filtered else targets
        return targets

    if diagnosis == DIAGNOSIS_RARE_DISEASE:
        targets = dict(ALL_RARE_DISEASE_TARGETS)
        targets.update(ALL_ADDITIONAL_RARE_TARGETS)
        if conditions:
            filtered = {}
            for gene, info in targets.items():
                gene_conditions = info.get("conditions", [])
                if any(c in gene_conditions for c in conditions):
                    filtered[gene] = info
            return filtered if filtered else targets
        return targets

    # Map remaining diagnosis categories to their target dicts
    _CATEGORY_MAP = {
        DIAGNOSIS_INFECTIOUS: INFECTIOUS_DISEASE_TARGETS,
        DIAGNOSIS_IMMUNE_DYSREGULATION: IMMUNE_DYSREGULATION_TARGETS,
        DIAGNOSIS_TELOMERE_AGING: TELOMERE_AGING_TARGETS,
        DIAGNOSIS_PULMONARY: PULMONARY_TARGETS,
        DIAGNOSIS_HEPATIC: HEPATIC_TARGETS,
        DIAGNOSIS_SKELETAL: SKELETAL_TARGETS,
        DIAGNOSIS_DERMATOLOGIC: DERMATOLOGIC_TARGETS,
        DIAGNOSIS_MITOCHONDRIAL: MITOCHONDRIAL_TARGETS,
        DIAGNOSIS_COMPLEMENT: COMPLEMENT_TARGETS,
        DIAGNOSIS_COAGULATION: COAGULATION_TARGETS,
        DIAGNOSIS_PAIN: PAIN_TARGETS,
        DIAGNOSIS_EPILEPSY: EPILEPSY_TARGETS,
        DIAGNOSIS_CONGENITAL_HEART: CONGENITAL_HEART_TARGETS,
        DIAGNOSIS_MONOGENIC_DIABETES: MONOGENIC_DIABETES_TARGETS,
        DIAGNOSIS_OBESITY: OBESITY_TARGETS,
        DIAGNOSIS_OPHTHALMOLOGY: OPHTHALMOLOGY_TARGETS,
        DIAGNOSIS_NCL: NCL_TARGETS,
        DIAGNOSIS_BONE_MARROW_FAILURE: BONE_MARROW_FAILURE_TARGETS,
        DIAGNOSIS_IMMUNE_CHECKPOINT: IMMUNE_CHECKPOINT_TARGETS,
        DIAGNOSIS_CRANIOFACIAL: CRANIOFACIAL_TARGETS,
    }
    if diagnosis in _CATEGORY_MAP:
        targets = dict(_CATEGORY_MAP[diagnosis])
        if conditions:
            filtered = {}
            for gene, info in targets.items():
                gene_conditions = info.get("conditions", [])
                if any(c in gene_conditions for c in conditions):
                    filtered[gene] = info
            return filtered if filtered else targets
        return targets

    raise ValueError(f"Unknown diagnosis: {diagnosis}")


def run(vcf_path: str, bam_path: str, bai_path: str | None,
        patient_id: str, diagnosis: str,
        conditions: list[str] | None = None,
        output_root: str = ".") -> Path:
    """Execute the full pipeline."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(output_root) / f"{patient_id}_{timestamp}_opencure"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Banner
    print(f"""
    ╔══════════════════════════════════════════════════════╗
    ║        OpenCure — Personalized CRISPR Pipeline       ║
    ║  Patient:   {patient_id:<41s} ║
    ║  Diagnosis: {diagnosis:<41s} ║
    ║  Date:      {timestamp:<41s} ║
    ╚══════════════════════════════════════════════════════╝
    """)

    # Select gene targets
    targets = select_targets(diagnosis, conditions)
    log.info("Selected %d target genes for %s", len(targets), diagnosis)

    # Step 1: Validate inputs
    validation = step1_validate.run(vcf_path, bam_path, bai_path, targets)
    if not validation["passed"]:
        log.warning("Input validation found issues: %s", validation["errors"])
        print("[WARN] Validation issues (continuing with available data):")
        for err in validation["errors"]:
            print(f"  - {err}")

    # Step 2: Extract patient variants
    profile, strategy = step2_extract_variants.run(
        vcf_path, targets, output_dir, patient_id,
    )

    # Step 2b: VEP variant impact assessment (all non-cancer modes)
    vep_annotations = []
    if diagnosis != DIAGNOSIS_CANCER:
        vep_annotations = step2b_vep_annotation.run(
            vcf_path, targets, profile, output_dir, patient_id,
        )

    # Step 2c: Autoinflammatory variant annotation (autoimmune mode)
    autoinflammatory_assessment = {}
    if diagnosis == DIAGNOSIS_AUTOIMMUNE:
        autoinflammatory_assessment = step2c_autoinflammatory_annotation.run(
            vcf_path, targets, profile, output_dir, patient_id,
        )

    # Step 3: Design guide RNAs
    guide_designs = step3_design_guides.run(
        targets, profile, output_dir, patient_id,
        bam_path=bam_path, bai_path=bai_path,
    )

    # Step 4: Off-target analysis
    guide_designs = step4_offtarget.run(
        guide_designs, output_dir, patient_id,
    )

    # Step 5: Construct design (with DREADD for neuropsych)
    dreadd = DREADD_CONSTRUCTS if diagnosis == DIAGNOSIS_NEUROPSYCH else None
    harbors = SAFE_HARBOR_LOCI if diagnosis == DIAGNOSIS_NEUROPSYCH else None
    constructs = step5_construct.run(
        targets, guide_designs, output_dir, patient_id,
        dreadd_constructs=dreadd,
        safe_harbors=harbors,
    )

    # Step 6: Report
    report_path = step6_report.run(
        patient_id, diagnosis, output_dir,
        validation, strategy, guide_designs, constructs,
        vep_annotations=vep_annotations,
        autoinflammatory=autoinflammatory_assessment,
    )

    print(f"\n[COMPLETE] Pipeline finished. Outputs in: {output_dir}/")
    print(f"           Report: {report_path}")

    return output_dir


def main():
    parser = argparse.ArgumentParser(
        description="OpenCure Pipeline — Personalized CRISPR therapy design",
    )
    parser.add_argument("--vcf", required=True,
                        help="Path to patient VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--bam", required=True,
                        help="Path to patient BAM file")
    parser.add_argument("--bai", default=None,
                        help="Path to BAM index (.bai); auto-detected if omitted")
    parser.add_argument("--patient", required=True,
                        help="Patient identifier")
    parser.add_argument("--diagnosis", required=True,
                        choices=sorted(VALID_DIAGNOSES),
                        help="Diagnosis category")
    parser.add_argument("--conditions", nargs="*", default=None,
                        help="(mental_health only) Filter targets by condition "
                             "e.g. depression anxiety PTSD")
    parser.add_argument("--output", default=".",
                        help="Output root directory (default: cwd)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose logging")

    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    try:
        run(
            vcf_path=args.vcf,
            bam_path=args.bam,
            bai_path=args.bai,
            patient_id=args.patient,
            diagnosis=args.diagnosis,
            conditions=args.conditions,
            output_root=args.output,
        )
    except Exception as exc:
        log.error("Pipeline failed: %s", exc, exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
