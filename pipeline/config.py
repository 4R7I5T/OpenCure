"""
Pipeline configuration: paths, tool settings, and shared constants.
"""

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
PIPELINE_DIR = Path(__file__).resolve().parent
GENOME_DIR = PROJECT_ROOT / "Genome"

REFERENCE_GENOME = os.getenv(
    "OPENCURE_REFERENCE",
    "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
)

# ---------------------------------------------------------------------------
# External tool paths (override via env vars)
# ---------------------------------------------------------------------------
SAMTOOLS = os.getenv("SAMTOOLS", "samtools")
BCFTOOLS = os.getenv("BCFTOOLS", "bcftools")
CAS_OFFINDER = os.getenv("CAS_OFFINDER", "cas-offinder")
FLASHFRY_JAR = os.getenv("FLASHFRY_JAR", "./FlashFry-assembly-1.15.jar")
VEP_EXECUTABLE = os.getenv("VEP_EXECUTABLE", "vep")
VEP_CACHE_DIR = os.getenv("VEP_CACHE_DIR", "")

# ---------------------------------------------------------------------------
# API keys
# ---------------------------------------------------------------------------
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "")
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "")
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")

# ---------------------------------------------------------------------------
# Guide design defaults
# ---------------------------------------------------------------------------
GUIDE_LENGTH = 20
PAM_SEQUENCE = "NGG"  # SpCas9
MAX_OFFTARGET_MISMATCHES = 4
FLANKING_BASES = 125  # bases upstream/downstream for sequence context

# ---------------------------------------------------------------------------
# AAV packaging
# ---------------------------------------------------------------------------
AAV_PACKAGING_LIMIT_BP = 4700

# ---------------------------------------------------------------------------
# Diagnosis categories
# ---------------------------------------------------------------------------
DIAGNOSIS_CANCER = "cancer"
DIAGNOSIS_MENTAL_HEALTH = "mental_health"
DIAGNOSIS_NEUROPSYCH = "neuropsychiatric"
DIAGNOSIS_AUTOIMMUNE = "autoimmune"
DIAGNOSIS_HEMATOLOGICAL = "hematological"
DIAGNOSIS_IMMUNODEFICIENCY = "immunodeficiency"
DIAGNOSIS_HEREDITARY = "hereditary"
DIAGNOSIS_CARDIAC = "cardiac"
DIAGNOSIS_RARE_DISEASE = "rare_disease"
DIAGNOSIS_INFECTIOUS = "infectious"
DIAGNOSIS_IMMUNE_DYSREGULATION = "immune_dysregulation"
DIAGNOSIS_TELOMERE_AGING = "telomere_aging"
DIAGNOSIS_PULMONARY = "pulmonary"
DIAGNOSIS_HEPATIC = "hepatic"
DIAGNOSIS_SKELETAL = "skeletal"
DIAGNOSIS_DERMATOLOGIC = "dermatologic"
DIAGNOSIS_MITOCHONDRIAL = "mitochondrial"
DIAGNOSIS_COMPLEMENT = "complement"
DIAGNOSIS_COAGULATION = "coagulation"
DIAGNOSIS_PAIN = "pain"
DIAGNOSIS_EPILEPSY = "epilepsy"
DIAGNOSIS_CONGENITAL_HEART = "congenital_heart"
DIAGNOSIS_MONOGENIC_DIABETES = "monogenic_diabetes"
DIAGNOSIS_OBESITY = "obesity"
DIAGNOSIS_OPHTHALMOLOGY = "ophthalmology"
DIAGNOSIS_NCL = "ncl"
DIAGNOSIS_BONE_MARROW_FAILURE = "bone_marrow_failure"
DIAGNOSIS_IMMUNE_CHECKPOINT = "immune_checkpoint"
DIAGNOSIS_CRANIOFACIAL = "craniofacial"

VALID_DIAGNOSES = {DIAGNOSIS_CANCER, DIAGNOSIS_MENTAL_HEALTH,
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
                   DIAGNOSIS_CRANIOFACIAL}
