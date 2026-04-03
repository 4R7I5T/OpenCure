# OpenCure

**Personalized CRISPR Therapy Design Pipeline** -- From patient genome to therapeutic construct in one command.

OpenCure analyzes patient whole-genome sequencing data (VCF + BAM) against curated disease gene databases, designs patient-specific CRISPR guide RNAs, performs off-target safety analysis, and generates complete AAV delivery construct maps.

**60 diagnosis pipelines | 534 gene targets | 400+ diseases**

---

## Quick Start

```bash
# Install
conda create -n opencure python=3.10
conda activate opencure
pip install -r pipeline/requirements.txt
conda install -c bioconda samtools bcftools ensembl-vep cas-offinder

# Run
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz \
    --bam patient.bam \
    --patient PATIENT_001 \
    --diagnosis cancer
```

## Supported Diagnoses

### Core Pipelines

| Mode | Flag | Genes | Description |
|------|------|-------|-------------|
| **Cancer** | `--diagnosis cancer` | 21 | TP53, BRCA1/2, KRAS, EGFR, BRAF, PTEN, MYC, MLH1/MSH2/MSH6/PMS2 (Lynch), and more |
| **Mental Health** | `--diagnosis mental_health` | 13 | SLC6A4, BDNF, COMT, OPRM1, FKBP5 -- pharmacogenomics and gene therapy |
| **Neuropsychiatric** | `--diagnosis neuropsychiatric` | 9 | Circuit-level neuromodulation with DREADD chemogenetic constructs |
| **Autoimmune / Yao Syndrome** | `--diagnosis autoimmune` | 19 | NOD2, MEFV, NLRP3 -- includes Yao syndrome diagnostic scoring engine |
| **Hematological** | `--diagnosis hematological` | 11 | Sickle cell (Casgevy-target BCL11A), thalassemia, hemophilia A/B |
| **Immunodeficiency** | `--diagnosis immunodeficiency` | 19 | SCID, CGD, cystic fibrosis, PKU, familial hypercholesterolemia |
| **Hereditary** | `--diagnosis hereditary` | 18 | Huntington, SMA, DMD, retinal dystrophies, epidermolysis bullosa, TTR amyloidosis |
| **Cardiac** | `--diagnosis cardiac` | 16 | HCM, DCM, ACM, Long QT, Brugada, CPVT, Marfan, vEDS |
| **Rare Disease** | `--diagnosis rare_disease` | 105 | Lysosomal storage, kidney, hearing loss, metabolic, neurological, and 50+ inborn errors |

### Specialized Pipelines

| Mode | Flag | Genes | Description |
|------|------|-------|-------------|
| **Infectious** | `--diagnosis infectious` | 2 | HIV (CCR5 knockout), CXCR4 |
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
| **Neurodegenerative** | `--diagnosis neurodegenerative` | 17 | Alzheimer (APP, PSEN1/2, APOE, TREM2), Parkinson (SNCA, LRRK2, GBA1, PRKN, PINK1), FTD, ALS |
| **Renal** | `--diagnosis renal` | 10 | Alport (COL4A3/4/5), ADPKD (PKD1/2), APOL1 nephropathy, cystinosis |
| **Metabolic Liver** | `--diagnosis metabolic_liver` | 16 | Urea cycle, tyrosinemia, galactosemia, GSD, MSUD, Wilson, hemochromatosis |
| **Hearing Loss** | `--diagnosis hearing_loss` | 6 | GJB2 (connexin 26), OTOF (otoferlin), TMC1, SLC26A4 |
| **Lipid Disorders** | `--diagnosis lipid` | 7 | FH (LDLR, PCSK9, APOB), ANGPTL3 (Verve base editing), Lp(a) |
| **Neuromuscular** | `--diagnosis neuromuscular` | 10 | LGMD types (CAPN3, DYSF, SGCA/B/G, FKRP), myotubular myopathy, congenital myasthenia |
| **Connective Tissue** | `--diagnosis connective_tissue` | 7 | Osteogenesis imperfecta (COL1A1/2), classical EDS, kyphoscoliotic EDS, cutis laxa |
| **Endocrine** | `--diagnosis endocrine` | 6 | Congenital adrenal hyperplasia, MEN2 (RET), MEN1, GH deficiency |
| **Vascular** | `--diagnosis vascular` | 7 | CADASIL (NOTCH3), HHT (ENG, ACVRL1), cerebral cavernous malformations, moyamoya |
| **Lymphatic** | `--diagnosis lymphatic` | 5 | Milroy disease, lymphedema-distichiasis, Hennekam syndrome |
| **Fertility** | `--diagnosis fertility` | 5 | Kallmann syndrome, primary ovarian insufficiency, hypogonadotropic hypogonadism |
| **Leukodystrophy** | `--diagnosis leukodystrophy` | 6 | MLD (ARSA, Libmeldy approved), Krabbe, X-ALD (Skysona approved), Canavan, Alexander, PMD |
| **Neurodevelopmental** | `--diagnosis neurodevelopmental` | 7 | Rett (MECP2), Fragile X (FMR1), Angelman (UBE3A), Phelan-McDermid, SYNGAP1, Kabuki |
| **Peripheral Neuropathy** | `--diagnosis peripheral_neuropathy` | 6 | CMT1A (PMP22), CMT2A (MFN2), CMTX1 (GJB1), CMT1B, CMT4C, HSAN1 |
| **Movement Disorder** | `--diagnosis movement_disorder` | 7 | DYT1 dystonia, dopa-responsive dystonia, AHC/RDP, PKD, myoclonus-dystonia |
| **Porphyria** | `--diagnosis porphyria` | 5 | AIP (givosiran approved), EPP, CEP, variegate porphyria, XLSA |
| **Ciliopathy** | `--diagnosis ciliopathy` | 7 | Primary ciliary dyskinesia, Bardet-Biedl, Joubert, nephronophthisis |
| **Skeletal Channelopathy** | `--diagnosis skeletal_channelopathy` | 5 | Malignant hyperthermia, periodic paralysis, myotonia congenita, Andersen-Tawil |
| **Platelet Disorder** | `--diagnosis platelet` | 5 | Glanzmann thrombasthenia, Bernard-Soulier, gray platelet, Hermansky-Pudlak |
| **Prion / Sleep** | `--diagnosis prion_sleep` | 2 | Fatal familial insomnia / genetic CJD (PRNP), narcolepsy (HCRT) |
| **Dental** | `--diagnosis dental` | 5 | Amelogenesis imperfecta, dentinogenesis imperfecta |
| **Spinocerebellar Ataxia** | `--diagnosis spinocerebellar_ataxia` | 6 | SCA1, SCA3/MJD, SCA6, SCA7, SCA13, SCA15 |
| **Spastic Paraplegia** | `--diagnosis spastic_paraplegia` | 6 | SPG4, SPG3A, SPG7, SPG11, SPG10, SPG31 |
| **Glycosylation Disorder** | `--diagnosis glycosylation` | 5 | PMM2-CDG, ALG1-CDG, ALG6-CDG, DPAGT1-CDG, SRD5A3-CDG |
| **Peroxisomal** | `--diagnosis peroxisomal` | 5 | Zellweger spectrum (PEX1/6/12), RCDP, ACOX1 deficiency |
| **Organic Acidemia** | `--diagnosis organic_acidemia` | 5 | MMA (active gene therapy trials), glutaric aciduria, alkaptonuria, VLCAD |
| **Congenital Myopathy** | `--diagnosis congenital_myopathy` | 5 | Nemaline myopathy (NEB, ACTA1, TPM3), centronuclear myopathy |
| **Mucopolysaccharidosis** | `--diagnosis mucopolysaccharidosis` | 5 | Sanfilippo A/B, Morquio A, GM1, MPS VI |
| **DNA Repair** | `--diagnosis dna_repair` | 5 | Xeroderma pigmentosum, Cockayne syndrome, trichothiodystrophy |
| **Bone Dysplasia** | `--diagnosis bone_dysplasia` | 5 | FOP (palovarotene approved), XLH (burosumab approved), hereditary exostoses |
| **Enteric / Neural Crest** | `--diagnosis enteric_neural_crest` | 4 | Hirschsprung, Ondine's curse/CCHS, Waardenburg-Shah |

## Pipeline Steps

```
Step 1    Validate inputs (VCF integrity, BAM index, per-gene coverage)
Step 2    Extract patient variants at target gene loci + promoter regions
            30x WGS quality filters (DP >= 10, AD >= 3, QUAL >= 30)
Step 2b   VEP variant impact assessment (SIFT, PolyPhen, consequence)
Step 2c   Autoinflammatory annotation (autoimmune mode only)
            Known-variant matching, compound het detection, Yao scoring
Step 3    Design patient-specific guide RNAs
            3-tier sequence fetch: local FASTA -> BAM consensus -> NCBI Entrez
            Doench-inspired scoring, Synthego manufacturing filters
Step 4    Off-target analysis (Cas-OFFinder or heuristic fallback)
Step 5    Construct design
            Infers CRISPR modality from strategy/inheritance/clinical trials
            CRISPRi, CRISPRa, HDR, base editing, knockout, DREADD
Step 6    Generate JSON summary + human-readable report
```

## Output Structure

```
PATIENT_001_20260403_opencure/
|-- pipeline_summary.json
|-- pipeline_report.txt
|-- PATIENT_001_gene_profile.json
|-- PATIENT_001_targeting_strategy.json
|-- PATIENT_001_guide_designs.json
|-- PATIENT_001_vep_annotated.vcf          # non-cancer modes
|-- PATIENT_001_vep_high_impact.json       # non-cancer modes
|-- PATIENT_001_autoinflammatory.json      # autoimmune mode
|-- offtarget_analysis/
|   |-- offtarget_results.json
|   +-- filtered_safe_guides.json
+-- construct_maps/
    |-- construct_CRISPRi_GENE.json
    |-- construct_CRISPRa_GENE.json
    |-- construct_HDR_correction_GENE.json
    |-- construct_base_editing_GENE.json
    +-- construct_DREADD_inhibitory.json   # neuropsych mode
```

## Construct Type Selection

Step 5 automatically selects the correct CRISPR modality for each gene:

| Modality | Effector | When Used |
|----------|----------|-----------|
| **CRISPRi** | dCas9-KRAB | Gain-of-function mutations, oncogene silencing |
| **CRISPRa** | dCas9-VPR | Haploinsufficiency, tumor suppressor reactivation |
| **HDR Correction** | SpCas9 | Autosomal recessive loss-of-function, gene replacement |
| **Base Editing** | ABE8e | Point mutation correction (clinical trial approaches) |
| **Knockout** | SpCas9 | Dominant-negative alleles, viral targets (CCR5) |
| **DREADD** | hM4Di/hM3Dq | Neuropsychiatric circuit modulation |

The modality is inferred from (in priority order):
1. Explicit `strategy` field in the gene database
2. `crispr_strategy` field
3. `clinical_trials[].approach` keywords
4. `key_variants[].approach` keywords
5. `inheritance` pattern (AR -> HDR, AD GOF -> CRISPRi, haploinsufficiency -> CRISPRa)

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

All non-cancer modes run Ensembl VEP for functional annotation (SIFT, PolyPhen, consequence classification). Falls back to heuristic scoring if VEP is not installed.

### Guide RNA Scoring

Guides are scored 0-100 using a Doench 2014/2016-inspired algorithm:

- Position-specific nucleotide weights (PAM-proximal seed weighted most heavily)
- GC content optimization (40-60% ideal)
- Homopolymer and self-complementarity penalties
- Synthego manufacturing compatibility filters (GC range, Pol III terminators, restriction sites)
- Top 10 guides per gene, ranked by score

## LLM Training Data Export

OpenCure can export its curated knowledge as training data for fine-tuning biomedical LLMs:

```bash
# Generate 3,900+ instruction-tuning examples
python -m pipeline.training.export_training_data \
    --format alpaca -o training_data.jsonl

# ShareGPT format (for OpenAI fine-tuning)
python -m pipeline.training.export_training_data \
    --format sharegpt -o training_data.jsonl
```

Covers gene function, variant interpretation, CRISPR strategy selection, construct design, guide constraints, and clinical trials across all 60 disease categories.

### Fine-Tuning on RunPod

```bash
# On a RunPod pod with A100 80GB (~$3-10 total)
pip install "unsloth[colab-new] @ git+https://github.com/unslothai/unsloth.git"
pip install --no-deps "trl<0.9.0" peft accelerate bitsandbytes

python pipeline/training/train_runpod.py \
    --data training_data.jsonl \
    --output opencure-bio-llm \
    --merge
```

Uses QLoRA (4-bit) on LLaMA 3.1 8B Instruct by default. Trains in 45-90 minutes.

## Usage Examples

```bash
# Cancer genomics
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz --bam patient.bam \
    --patient P001 --diagnosis cancer

# Sickle cell disease
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz --bam patient.bam \
    --patient P002 --diagnosis hematological \
    --conditions sickle_cell_disease

# Yao syndrome with diagnostic scoring
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz --bam patient.bam \
    --patient P003 --diagnosis autoimmune \
    --conditions yao_syndrome

# Huntington's / Parkinson's / ALS
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz --bam patient.bam \
    --patient P004 --diagnosis neurodegenerative \
    --conditions parkinsons_disease

# PCSK9 base editing (Verve approach)
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz --bam patient.bam \
    --patient P005 --diagnosis lipid \
    --conditions familial_hypercholesterolemia

# OTOF hearing loss gene therapy
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz --bam patient.bam \
    --patient P006 --diagnosis hearing_loss \
    --conditions DFNB9

# Rett syndrome (MECP2 gene therapy trials active)
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz --bam patient.bam \
    --patient P007 --diagnosis neurodevelopmental \
    --conditions rett_syndrome

# Verbose logging
python -m pipeline.run_pipeline \
    --vcf patient.vcf.gz --bam patient.bam \
    --patient P001 --diagnosis cancer -v
```

### Condition-Based Filtering

Use `--conditions` to focus on specific diseases within any pipeline:

```bash
--diagnosis cardiac --conditions LQTS
--diagnosis rare_disease --conditions lysosomal_storage
--diagnosis immunodeficiency --conditions SCID
--diagnosis mitochondrial --conditions mtDNA_depletion
```

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

### Environment Variables

```bash
export SAMTOOLS=samtools
export BCFTOOLS=bcftools
export VEP_EXECUTABLE=vep
export VEP_CACHE_DIR=/path/to/vep/cache
export CAS_OFFINDER=cas-offinder
export NCBI_EMAIL=your@email.com      # enables remote sequence fetch
export NCBI_API_KEY=your_ncbi_key     # optional, increases rate limit
```

## Project Structure

```
OpenCure/
|-- pipeline/
|   |-- __init__.py              # v3.0.0
|   |-- __main__.py              # CLI entry point
|   |-- config.py                # Constants, thresholds, diagnosis categories
|   |-- run_pipeline.py          # Orchestrator, target selection, condition filtering
|   |-- databases/               # 45 curated gene target databases (534 genes)
|   |   |-- cancer_genes.py
|   |   |-- therapeutic_targets.py
|   |   |-- neuropsych_targets.py
|   |   |-- autoimmune_targets.py
|   |   |-- hematological_targets.py
|   |   |-- cardiac_targets.py
|   |   |-- hereditary_disease_targets.py
|   |   |-- immunodeficiency_metabolic_targets.py
|   |   |-- infectious_immune_aging_targets.py
|   |   |-- neurodegenerative_targets.py
|   |   |-- ... (35 more disease-specific databases)
|   |   +-- yao_variants.py      # NOD2 variant catalog for Yao syndrome
|   |-- steps/
|   |   |-- step1_validate.py
|   |   |-- step2_extract_variants.py
|   |   |-- step2b_vep_annotation.py
|   |   |-- step2c_autoinflammatory_annotation.py
|   |   |-- step3_design_guides.py
|   |   |-- step4_offtarget.py
|   |   |-- step5_construct.py
|   |   +-- step6_report.py
|   |-- utils/
|   |   |-- guide_utils.py       # PAM finding, Doench scoring, Synthego filters
|   |   |-- sequence_utils.py    # Reverse complement, variant application
|   |   +-- vcf_utils.py         # VCF parsing, pysam integration
|   +-- training/
|       |-- export_training_data.py  # LLM training data generator
|       +-- train_runpod.py          # QLoRA fine-tuning script
+-- pipeline/requirements.txt
```

## Gene Database Format

Each gene entry contains:

```python
"BRCA1": {
    "chrom": "chr17",
    "start": 43044295,
    "end": 43170245,
    "role": "DNA double-strand break repair via homologous recombination",
    "strategy": "HDR-mediated correction of pathogenic variants",
    "conditions": ["breast_cancer", "ovarian_cancer"],
    "inheritance": "AD",               # some databases
    "key_variants": [...],             # rsIDs, consequences, notes
    "clinical_trials": [...],          # NCT IDs, phase, approach
}
```

All coordinates are GRCh38/hg38 verified against NCBI Gene.

## References

- [Casgevy](https://www.casgevyhcp.com/) -- FDA-approved CRISPR for SCD/thalassemia (Dec 2023)
- [Nomani et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11449693/) -- Yao syndrome genotype/phenotype
- NCBI Gene, ClinVar, dbSNP, gnomAD, COSMIC, OncoKB
- Ensembl VEP (Variant Effect Predictor)
- Doench et al. 2014, 2016 -- sgRNA activity scoring
- Synthego sgRNA design guidelines

## License

See [LICENSE](LICENSE) for details.

## Disclaimer

This software is for **research purposes only**. Not a medical device. Not FDA-approved. All therapeutic interventions require informed patient consent and IRB/ethics approval. Verify all genomic coordinates against current NCBI/Ensembl releases before any clinical application.
