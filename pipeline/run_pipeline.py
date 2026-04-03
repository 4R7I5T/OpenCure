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

from .config import (DIAGNOSIS_CANCER, DIAGNOSIS_NEUROPSYCH,
                     DIAGNOSIS_AUTOIMMUNE, DIAGNOSIS_RARE_DISEASE,
                     VALID_DIAGNOSES)
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
from .databases.inborn_errors_and_rare_targets import (
    ALL_INBORN_ERRORS_TARGETS,
)
from .databases.neurodegenerative_targets import ALL_NEURODEGENERATIVE_TARGETS
from .databases.renal_targets import ALL_RENAL_TARGETS
from .databases.metabolic_liver_targets import ALL_METABOLIC_LIVER_TARGETS
from .databases.hearing_loss_targets import HEARING_LOSS_TARGETS
from .databases.lipid_targets import LIPID_TARGETS
from .databases.neuromuscular_targets import NEUROMUSCULAR_TARGETS
from .databases.connective_tissue_targets import CONNECTIVE_TISSUE_TARGETS
from .databases.endocrine_targets import ENDOCRINE_TARGETS
from .databases.vascular_targets import VASCULAR_TARGETS
from .databases.lymphatic_targets import LYMPHATIC_TARGETS
from .databases.fertility_targets import FERTILITY_TARGETS
from .databases.leukodystrophy_targets import LEUKODYSTROPHY_TARGETS
from .databases.neurodevelopmental_targets import NEURODEVELOPMENTAL_TARGETS
from .databases.peripheral_neuropathy_targets import PERIPHERAL_NEUROPATHY_TARGETS
from .databases.movement_disorder_targets import MOVEMENT_DISORDER_TARGETS
from .databases.porphyria_targets import PORPHYRIA_TARGETS
from .databases.ciliopathy_targets import CILIOPATHY_TARGETS
from .databases.skeletal_channelopathy_targets import SKELETAL_CHANNELOPATHY_TARGETS
from .databases.platelet_targets import PLATELET_TARGETS
from .databases.prion_sleep_targets import PRION_SLEEP_TARGETS
from .databases.dental_targets import DENTAL_TARGETS
from .databases.spinocerebellar_ataxia_targets import SPINOCEREBELLAR_ATAXIA_TARGETS
from .databases.spastic_paraplegia_targets import SPASTIC_PARAPLEGIA_TARGETS
from .databases.glycosylation_targets import GLYCOSYLATION_TARGETS
from .databases.peroxisomal_targets import PEROXISOMAL_TARGETS
from .databases.organic_acidemia_targets import ORGANIC_ACIDEMIA_TARGETS
from .databases.congenital_myopathy_targets import CONGENITAL_MYOPATHY_TARGETS
from .databases.mucopolysaccharidosis_targets import MUCOPOLYSACCHARIDOSIS_TARGETS
from .databases.dna_repair_targets import DNA_REPAIR_TARGETS
from .databases.bone_dysplasia_targets import BONE_DYSPLASIA_TARGETS
from .databases.enteric_neural_crest_targets import ENTERIC_NEURAL_CREST_TARGETS

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


def _filter_by_conditions(targets: dict, conditions: list[str] | None) -> dict:
    """Filter targets by condition; returns all targets if no match."""
    if not conditions:
        return targets
    filtered = {
        gene: info for gene, info in targets.items()
        if any(c in info.get("conditions", []) for c in conditions)
    }
    return filtered if filtered else targets


# Dispatch table: diagnosis → target dict(s)
_TARGETS_MAP: dict[str, dict] = {
    "cancer":                CANCER_TARGETS,
    "mental_health":         MENTAL_HEALTH_TARGETS,
    "neuropsychiatric":      NEUROPSYCH_TARGETS,
    "autoimmune":            AUTOIMMUNE_TARGETS,
    "hematological":         ALL_HEMATOLOGICAL_TARGETS,
    "immunodeficiency":      ALL_GENE_THERAPY_TARGETS,
    "hereditary":            ALL_HEREDITARY_TARGETS,
    "cardiac":               ALL_CARDIAC_TARGETS,
    "rare_disease":          ALL_RARE_DISEASE_TARGETS,
    "infectious":            INFECTIOUS_DISEASE_TARGETS,
    "immune_dysregulation":  IMMUNE_DYSREGULATION_TARGETS,
    "telomere_aging":        TELOMERE_AGING_TARGETS,
    "pulmonary":             PULMONARY_TARGETS,
    "hepatic":               HEPATIC_TARGETS,
    "skeletal":              SKELETAL_TARGETS,
    "dermatologic":          DERMATOLOGIC_TARGETS,
    "mitochondrial":         MITOCHONDRIAL_TARGETS,
    "complement":            COMPLEMENT_TARGETS,
    "coagulation":           COAGULATION_TARGETS,
    "pain":                  PAIN_TARGETS,
    "epilepsy":              EPILEPSY_TARGETS,
    "congenital_heart":      CONGENITAL_HEART_TARGETS,
    "monogenic_diabetes":    MONOGENIC_DIABETES_TARGETS,
    "obesity":               OBESITY_TARGETS,
    "ophthalmology":         OPHTHALMOLOGY_TARGETS,
    "ncl":                   NCL_TARGETS,
    "bone_marrow_failure":   BONE_MARROW_FAILURE_TARGETS,
    "immune_checkpoint":     IMMUNE_CHECKPOINT_TARGETS,
    "craniofacial":          CRANIOFACIAL_TARGETS,
    "neurodegenerative":     ALL_NEURODEGENERATIVE_TARGETS,
    "renal":                 ALL_RENAL_TARGETS,
    "metabolic_liver":       ALL_METABOLIC_LIVER_TARGETS,
    "hearing_loss":          HEARING_LOSS_TARGETS,
    "lipid":                 LIPID_TARGETS,
    "neuromuscular":         NEUROMUSCULAR_TARGETS,
    "connective_tissue":     CONNECTIVE_TISSUE_TARGETS,
    "endocrine":             ENDOCRINE_TARGETS,
    "vascular":              VASCULAR_TARGETS,
    "lymphatic":             LYMPHATIC_TARGETS,
    "fertility":             FERTILITY_TARGETS,
    "leukodystrophy":        LEUKODYSTROPHY_TARGETS,
    "neurodevelopmental":    NEURODEVELOPMENTAL_TARGETS,
    "peripheral_neuropathy": PERIPHERAL_NEUROPATHY_TARGETS,
    "movement_disorder":     MOVEMENT_DISORDER_TARGETS,
    "porphyria":             PORPHYRIA_TARGETS,
    "ciliopathy":            CILIOPATHY_TARGETS,
    "skeletal_channelopathy": SKELETAL_CHANNELOPATHY_TARGETS,
    "platelet":              PLATELET_TARGETS,
    "prion_sleep":           PRION_SLEEP_TARGETS,
    "dental":                DENTAL_TARGETS,
    "spinocerebellar_ataxia": SPINOCEREBELLAR_ATAXIA_TARGETS,
    "spastic_paraplegia":    SPASTIC_PARAPLEGIA_TARGETS,
    "glycosylation":         GLYCOSYLATION_TARGETS,
    "peroxisomal":           PEROXISOMAL_TARGETS,
    "organic_acidemia":      ORGANIC_ACIDEMIA_TARGETS,
    "congenital_myopathy":   CONGENITAL_MYOPATHY_TARGETS,
    "mucopolysaccharidosis": MUCOPOLYSACCHARIDOSIS_TARGETS,
    "dna_repair":            DNA_REPAIR_TARGETS,
    "bone_dysplasia":        BONE_DYSPLASIA_TARGETS,
    "enteric_neural_crest":  ENTERIC_NEURAL_CREST_TARGETS,
}


def select_targets(diagnosis: str, conditions: list[str] | None = None) -> dict:
    """
    Return the gene-target dict appropriate for the diagnosis.
    Optionally filter targets by condition.
    """
    if diagnosis not in _TARGETS_MAP:
        raise ValueError(f"Unknown diagnosis: {diagnosis}")

    targets = dict(_TARGETS_MAP[diagnosis])

    # Rare disease aggregates multiple target databases
    if diagnosis == DIAGNOSIS_RARE_DISEASE:
        targets.update(ALL_ADDITIONAL_RARE_TARGETS)
        targets.update(ALL_INBORN_ERRORS_TARGETS)

    # Exclude entries without valid human genomic coordinates (e.g. viral
    # targets like HBV_cccDNA / HPV_E6E7 that need a separate workflow).
    targets = {
        gene: info for gene, info in targets.items()
        if isinstance(info.get("start"), int) and isinstance(info.get("end"), int)
    }

    return _filter_by_conditions(targets, conditions)


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
