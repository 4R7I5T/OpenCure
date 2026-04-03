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
# Synthego sgRNA manufacturing constraints
# ---------------------------------------------------------------------------
# Synthego synthesizes chemically modified sgRNAs (2'-O-methyl + 3x PS bonds
# at each terminus).  These filters ensure guides can be manufactured and
# will perform well in vitro/in vivo.
SYNTHEGO_MIN_GC = 0.40          # hard minimum GC content
SYNTHEGO_MAX_GC = 0.70          # hard maximum GC content
SYNTHEGO_MAX_HOMOPOLYMER = 4    # reject guides with ≥N identical consecutive bases
SYNTHEGO_REJECT_POLY_T = True   # TTTT is a Pol III terminator signal -- always reject
SYNTHEGO_SEED_LENGTH = 12       # PAM-proximal seed region (critical for specificity)

# ---------------------------------------------------------------------------
# 30x WGS variant calling thresholds
# ---------------------------------------------------------------------------
# At 30x mean coverage, heterozygous variants require ≥3 alt reads and ≥10x
# total depth for reliable calling.  These thresholds filter VCF variants
# before guide design to avoid designing guides against false-positive variants.
MIN_VARIANT_DEPTH = 10          # minimum total read depth (DP) at variant site
MIN_VARIANT_ALT_DEPTH = 3       # minimum alternate allele depth (AD)
MIN_VARIANT_QUAL = 30           # minimum QUAL score
MIN_GENOTYPE_QUAL = 20          # minimum genotype quality (GQ)
LOW_COVERAGE_WARN_THRESHOLD = 15  # warn if mean coverage at target locus < this

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

# --- v2.2.0 pipelines ---
DIAGNOSIS_NEURODEGENERATIVE = "neurodegenerative"
DIAGNOSIS_RENAL = "renal"
DIAGNOSIS_METABOLIC_LIVER = "metabolic_liver"
DIAGNOSIS_HEARING_LOSS = "hearing_loss"
DIAGNOSIS_LIPID = "lipid"
DIAGNOSIS_NEUROMUSCULAR = "neuromuscular"
DIAGNOSIS_CONNECTIVE_TISSUE = "connective_tissue"
DIAGNOSIS_ENDOCRINE = "endocrine"
DIAGNOSIS_VASCULAR = "vascular"
DIAGNOSIS_LYMPHATIC = "lymphatic"
DIAGNOSIS_FERTILITY = "fertility"

# --- v2.3.0 pipelines ---
DIAGNOSIS_LEUKODYSTROPHY = "leukodystrophy"
DIAGNOSIS_NEURODEVELOPMENTAL = "neurodevelopmental"
DIAGNOSIS_PERIPHERAL_NEUROPATHY = "peripheral_neuropathy"
DIAGNOSIS_MOVEMENT_DISORDER = "movement_disorder"
DIAGNOSIS_PORPHYRIA = "porphyria"
DIAGNOSIS_CILIOPATHY = "ciliopathy"
DIAGNOSIS_SKELETAL_CHANNELOPATHY = "skeletal_channelopathy"
DIAGNOSIS_PLATELET = "platelet"
DIAGNOSIS_PRION_SLEEP = "prion_sleep"
DIAGNOSIS_DENTAL = "dental"

# --- v2.4.0 pipelines ---
DIAGNOSIS_SPINOCEREBELLAR_ATAXIA = "spinocerebellar_ataxia"
DIAGNOSIS_SPASTIC_PARAPLEGIA = "spastic_paraplegia"
DIAGNOSIS_GLYCOSYLATION = "glycosylation"
DIAGNOSIS_PEROXISOMAL = "peroxisomal"
DIAGNOSIS_ORGANIC_ACIDEMIA = "organic_acidemia"
DIAGNOSIS_CONGENITAL_MYOPATHY = "congenital_myopathy"
DIAGNOSIS_MUCOPOLYSACCHARIDOSIS = "mucopolysaccharidosis"
DIAGNOSIS_DNA_REPAIR = "dna_repair"
DIAGNOSIS_BONE_DYSPLASIA = "bone_dysplasia"
DIAGNOSIS_ENTERIC_NEURAL_CREST = "enteric_neural_crest"

VALID_DIAGNOSES = {
    v for k, v in globals().items()
    if k.startswith("DIAGNOSIS_") and isinstance(v, str)
}
