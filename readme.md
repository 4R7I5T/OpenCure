# OpenCure

**Personalized CRISPR Therapy Design Pipeline** -- From patient genome to therapeutic construct in one command.

OpenCure analyzes patient whole-genome sequencing data (VCF + BAM) against curated disease gene databases, designs patient-specific CRISPR guide RNAs, performs off-target safety analysis, and generates complete AAV delivery construct maps.

**29 diagnosis pipelines | 325 unique gene targets | 105+ rare diseases**

---

## Quick Start

```bash
# Install dependencies
conda create -n opencure python=3.10
conda activate opencure
pip install -r pipeline/requirements.txt
conda install -c bioconda samtools bcftools ensembl-vep cas-offinder

# Run the pipeline
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz \
    --bam patient.bam \
    --patient PATIENT_001 \
    --diagnosis cancer
```

## Supported Diagnosis Modes

### Core Pipelines

| Mode | Flag | Genes | Description |
|------|------|-------|-------------|
| **Cancer** | `--diagnosis cancer` | 17 | TP53, BRCA1/2, KRAS, EGFR, BRAF, PTEN, MYC, and more |
| **Mental Health** | `--diagnosis mental_health` | 13 | SLC6A4, BDNF, COMT, OPRM1, FKBP5 -- pharmacogenomics and gene therapy |
| **Neuropsychiatric** | `--diagnosis neuropsychiatric` | 9 | Circuit-level neuromodulation with DREADD chemogenetic constructs |
| **Autoimmune / Yao Syndrome** | `--diagnosis autoimmune` | 19 | NOD2, MEFV, NLRP3 -- includes Yao syndrome diagnostic scoring engine |
| **Hematological** | `--diagnosis hematological` | 11 | Sickle cell (Casgevy-target BCL11A), thalassemia, hemophilia A/B |
| **Immunodeficiency** | `--diagnosis immunodeficiency` | 19 | SCID, CGD, cystic fibrosis, PKU, familial hypercholesterolemia |
| **Hereditary** | `--diagnosis hereditary` | 18 | Huntington, SMA, DMD, retinal dystrophies, epidermolysis bullosa, TTR amyloidosis |
| **Cardiac** | `--diagnosis cardiac` | 16 | HCM, DCM, ACM, Long QT, Brugada, CPVT, Marfan, vEDS |

### Specialized Pipelines

| Mode | Flag | Genes | Description |
|------|------|-------|-------------|
| **Rare Disease** | `--diagnosis rare_disease` | 105 | Lysosomal storage, kidney, hearing loss, metabolic, neurological, and 50+ inborn errors |
| **Infectious** | `--diagnosis infectious` | 4 | HIV (CCR5 knockout), Hepatitis B, HPV |
| **Immune Dysregulation** | `--diagnosis immune_dysregulation` | 11 | IPEX, familial HLH, ALPS, CTLA-4, STAT3/1 GOF |
| **Telomere / Aging** | `--diagnosis telomere_aging` | 5 | Dyskeratosis congenita, Werner syndrome, progeria |
| **Pulmonary** | `--diagnosis pulmonary` | 3 | Pulmonary arterial hypertension, surfactant deficiency |
| **Hepatic** | `--diagnosis hepatic` | 4 | Crigler-Najjar, PFIC types 1-3 |
| **Skeletal** | `--diagnosis skeletal` | 2 | Achondroplasia, hypophosphatasia |
| **Dermatologic** | `--diagnosis dermatologic` | 2 | Ichthyosis |
| **Mitochondrial** | `--diagnosis mitochondrial` | 9 | POLG, Leigh syndrome, mtDNA depletion syndromes |
| **Complement** | `--diagnosis complement` | 6 | aHUS, complement deficiencies, C3 glomerulopathy |
| **Coagulation** | `--diagnosis coagulation` | 6 | von Willebrand, thrombophilia, Factor V Leiden |
| **Pain** | `--diagnosis pain` | 3 | Erythromelalgia, congenital insensitivity to pain |
| **Epilepsy** | `--diagnosis epilepsy` | 6 | SCN2A, KCNQ2, CDKL5, STXBP1, SCN8A, KCNT1 |
| **Congenital Heart** | `--diagnosis congenital_heart` | 4 | ASD/VSD, Holt-Oram, structural heart defects |
| **Monogenic Diabetes** | `--diagnosis monogenic_diabetes` | 6 | MODY (GCK, HNF1A, HNF4A), neonatal diabetes |
| **Obesity** | `--diagnosis obesity` | 5 | MC4R, leptin/leptin receptor, POMC deficiency |
| **Ophthalmology** | `--diagnosis ophthalmology` | 7 | Stargardt, achromatopsia, aniridia, retinoschisis |
| **NCL / Batten** | `--diagnosis ncl` | 7 | Neuronal ceroid lipofuscinoses CLN1-8 |
| **Bone Marrow Failure** | `--diagnosis bone_marrow_failure` | 4 | RUNX1, GATA2, Shwachman-Diamond |
| **Immune Checkpoint** | `--diagnosis immune_checkpoint` | 5 | PD-1, PD-L1, LAG-3, TIM-3, TIGIT |
| **Craniofacial** | `--diagnosis craniofacial` | 4 | Treacher Collins, Crouzon/Apert, Saethre-Chotzen |

## Usage Examples

### Basic Usage

```bash
# Cancer genomics -- all 17 oncogene/tumor suppressor targets
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz \
    --bam patient.bam \
    --patient PATIENT_001 \
    --diagnosis cancer

# Sickle cell disease -- BCL11A, HBB, HBG1/2, and HbF regulators
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz \
    --bam patient.bam \
    --patient PATIENT_002 \
    --diagnosis hematological \
    --conditions sickle_cell_disease

# Yao syndrome -- NOD2 with diagnostic scoring and treatment guidance
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz \
    --bam patient.bam \
    --patient PATIENT_003 \
    --diagnosis autoimmune \
    --conditions yao_syndrome
```

### Condition-Based Filtering

Use `--conditions` to focus on specific diseases within a pipeline:

```bash
# Cardiac: only Long QT syndrome genes
python -m pipeline.run_pipeline \
    --vcf p.vcf.gz --bam p.bam --patient P001 \
    --diagnosis cardiac --conditions LQTS

# Rare disease: only lysosomal storage diseases
python -m pipeline.run_pipeline \
    --vcf p.vcf.gz --bam p.bam --patient P001 \
    --diagnosis rare_disease --conditions lysosomal_storage

# Rare disease: only hearing loss genes
python -m pipeline.run_pipeline \
    --vcf p.vcf.gz --bam p.bam --patient P001 \
    --diagnosis rare_disease --conditions hearing_loss

# Immunodeficiency: only SCID genes
python -m pipeline.run_pipeline \
    --vcf p.vcf.gz --bam p.bam --patient P001 \
    --diagnosis immunodeficiency --conditions SCID

# Mitochondrial: only mtDNA depletion syndromes
python -m pipeline.run_pipeline \
    --vcf p.vcf.gz --bam p.bam --patient P001 \
    --diagnosis mitochondrial --conditions mtDNA_depletion
```

### Verbose Logging

```bash
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz --bam patient.bam \
    --patient PATIENT_001 --diagnosis autoimmune \
    --conditions yao_syndrome -v
```

## Pipeline Steps

The pipeline executes up to 8 steps automatically:

```
Step 1:   Validate inputs (VCF integrity, BAM index, coverage at target loci)
Step 2:   Extract patient variants at target gene loci + promoter regions
Step 2b:  VEP variant impact assessment (SIFT, PolyPhen, consequence)
Step 2c:  Autoinflammatory annotation (autoimmune mode only)
            -- Known variant matching, compound het detection,
               digenic screening, Yao syndrome scoring
Step 3:   Design patient-specific guide RNAs (PAM sites, scoring)
Step 4:   Off-target analysis (Cas-OFFinder or heuristic fallback)
Step 5:   Construct design (AAV9 vectors, DREADD for neuropsych)
Step 6:   Generate JSON summary + human-readable report
```

## Output Structure

```
PATIENT_001_20260402_143022_opencure/
|-- pipeline_summary.json
|-- pipeline_report.txt
|-- PATIENT_001_gene_profile.json
|-- PATIENT_001_targeting_strategy.json
|-- PATIENT_001_guide_designs.json
|-- PATIENT_001_vep_annotated.vcf
|-- PATIENT_001_vep_high_impact.json
|-- PATIENT_001_autoinflammatory_assessment.json   # autoimmune mode
|-- offtarget_analysis/
|   |-- offtarget_results.json
|   +-- filtered_safe_guides.json
+-- construct_maps/
    |-- construct_CRISPRi_GENE.json
    |-- construct_CRISPRa_GENE.json
    +-- construct_DREADD_inhibitory.json            # neuropsych mode
```

## Specialized Features

### Yao Syndrome Diagnostic Engine

When running `--diagnosis autoimmune --conditions yao_syndrome`:

- Matches patient NOD2 variants against 8 cataloged Yao-associated variants
- Detects compound heterozygosity (49% of patients carry 2+ NOD2 variants)
- Screens digenic combinations (NOD2 + MEFV/NLRP3/NLRP12/TNFRSF1A)
- Scores diagnostic likelihood (0-100)
- Generates differential diagnosis (Blau, Crohn's, FMF, CAPS)
- Provides genotype-guided treatment recommendations

### DREADD Chemogenetic Constructs

The neuropsychiatric pipeline designs DREADD constructs for on-demand neural circuit control:

- hM4Di inhibitory DREADD at AAVS1 safe harbor
- hM3Dq excitatory DREADD at ROSA26 safe harbor
- Activated by deschloroclozapine (DCZ)

### VEP Variant Impact Assessment

All non-cancer modes run Ensembl VEP for functional annotation. Falls back to heuristic scoring if VEP is not installed.

## Installation

### Prerequisites

- Python 3.10+
- conda (recommended)

### Full Installation

```bash
conda create -n opencure python=3.10
conda activate opencure

pip install -r pipeline/requirements.txt

# Bioinformatics tools (optional but recommended)
conda install -c bioconda samtools bcftools
conda install -c bioconda ensembl-vep
conda install -c bioconda cas-offinder

# Reference genome
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Environment Variables (optional)

```bash
export SAMTOOLS=samtools
export BCFTOOLS=bcftools
export VEP_EXECUTABLE=vep
export VEP_CACHE_DIR=/path/to/vep/cache
export CAS_OFFINDER=cas-offinder
export NCBI_EMAIL=your@email.com
export NCBI_API_KEY=your_ncbi_key
```

## Project Structure

```
OpenCure/
|-- pipeline/
|   |-- __init__.py              # Version 2.1.0
|   |-- __main__.py
|   |-- config.py                # 29 diagnosis constants
|   |-- run_pipeline.py          # Orchestration
|   |-- databases/               # 15 gene target databases (325 genes)
|   |-- steps/                   # 8 pipeline stages
|   +-- utils/                   # VCF, guide RNA, sequence utilities
+-- OpenCure-main/               # Legacy v1 (reference only)
```

## References

- [Casgevy](https://www.casgevyhcp.com/) -- FDA-approved CRISPR for SCD/thalassemia
- [Nomani et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11449693/) -- Yao syndrome genotype/phenotype
- NCBI Gene, ClinVar, dbSNP, gnomAD, COSMIC, OncoKB
- Ensembl VEP (Variant Effect Predictor)

## License

See [LICENSE](LICENSE) for details.

## Disclaimer

This software is for **research purposes only**. Not a medical device. All interventions require informed consent and IRB/ethics approval.
