"""
Additional rare disease gene-therapy / CRISPR targets -- GRCh38 (hg38) coordinates.

Covers neurological, metabolic, hematologic, and endocrine rare diseases
NOT already present in the existing OpenCure databases (which cover
HTT, SMN1, SOD1, PKU/PAH, Wilson/ATP7B, AATD/SERPINA1, SCD, hemophilia,
and lysosomal/kidney/endocrine/hearing/connective tissue targets).

Categories:
  1. Neurological rare diseases
     - Rett syndrome (MECP2)
     - Fragile X syndrome (FMR1)
     - Angelman syndrome (UBE3A)
     - Dravet syndrome (SCN1A)
     - Tuberous sclerosis (TSC1, TSC2)
     - Neurofibromatosis (NF1, NF2)
     - Charcot-Marie-Tooth disease (PMP22, GJB1, MFN2)
     - Canavan disease (ASPA)
  2. Metabolic rare diseases
     - Maple syrup urine disease (BCKDHA, BCKDHB, DBT)
     - Tyrosinemia type 1 (FAH)
     - Glycogen storage disease type Ia (G6PC1)
     - Classical galactosemia (GALT)
     - Homocystinuria (CBS)
     - OTC deficiency (OTC)
     - Citrullinemia type I (ASS1)
  3. Hematologic rare diseases
     - HBS1L-MYB intergenic (thalassemia modifier)
     - Congenital neutropenia (ELANE)
     - Hereditary spherocytosis (ANK1, SLC4A1)
  4. Endocrine rare diseases
     - Congenital hypothyroidism (TSHR, PAX8)
     - Multiple endocrine neoplasia (MEN1, RET)

Coordinates are from NCBI Gene, GRCh38.p14 (RefSeq annotation RS_2025_08).
Variant positions from NCBI dbSNP / ClinVar.

IMPORTANT -- verify every coordinate against current NCBI / Ensembl releases
before production use.  Numbering can shift between patch levels.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14)
  - OMIM
  - ClinVar
  - ClinicalTrials.gov
  - CRISPR Medicine News
  - Published literature through early 2026
"""


# ============================================================================
# 1. NEUROLOGICAL RARE DISEASES
# ============================================================================

NEUROLOGICAL_RARE_TARGETS = {

    # -------------------------------------------------------------------
    # 1a.  Rett Syndrome -- MECP2
    # -------------------------------------------------------------------
    "MECP2": {
        "gene_id": 4204,
        "chrom": "chrX",
        "start": 154_021_573,
        "end": 154_097_717,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq28",
        "exon_count": 4,
        "role": (
            "Methyl-CpG binding protein 2 -- transcriptional regulator that "
            "binds methylated DNA to modulate chromatin architecture and gene "
            "expression.  Highly expressed in mature neurons.  Acts as both "
            "transcriptional activator and repressor depending on context.  "
            "Loss-of-function causes Rett syndrome (progressive "
            "neurodevelopmental disorder, almost exclusively in females)."
        ),
        "disease": "Rett Syndrome",
        "omim_disease": 312750,
        "omim_gene": 300005,
        "inheritance": "XLD (de novo in >99%)",
        "mutation_type": "point_mutations_and_small_indels",
        "key_variants": [
            {
                "name": "p.Thr158Met (T158M)",
                "hgvs_coding": "NM_004992.4:c.473C>T",
                "rsid": "rs28934907",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "frequency": "~10% of all Rett cases",
                "notes": (
                    "MBD domain; most common individual MECP2 mutation.  "
                    "Disrupts methyl-CpG binding.  Generally associated "
                    "with classic Rett phenotype."
                ),
            },
            {
                "name": "p.Arg168Ter (R168X)",
                "hgvs_coding": "NM_004992.4:c.502C>T",
                "rsid": "rs28934906",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "TRD domain; truncating variant.  One of eight recurrent "
                    "mutations that account for ~70% of Rett cases.  More "
                    "severe phenotype than missense mutations."
                ),
            },
            {
                "name": "p.Arg255Ter (R255X)",
                "hgvs_coding": "NM_004992.4:c.763C>T",
                "rsid": "rs28934908",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Truncating variant; severe phenotype.",
            },
            {
                "name": "p.Arg306Cys (R306C)",
                "hgvs_coding": "NM_004992.4:c.916C>T",
                "rsid": "rs28934909",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "NCoR/SMRT interaction domain; milder phenotype.  "
                    "Preserved speech variant more common."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 900 pathogenic/likely-pathogenic variants in ClinVar.  Eight "
            "recurrent point mutations (R106W, R133C, T158M, R168X, R255X, "
            "R270X, R294X, R306C) account for ~65-70% of cases.  Large "
            "deletions account for ~10%.  Males with MECP2 loss-of-function "
            "typically have severe neonatal encephalopathy (lethal or near-"
            "lethal).  Females are heterozygous with skewed X-inactivation "
            "influencing severity."
        ),
        "crispr_strategy": (
            "Multiple approaches under investigation: "
            "(1) AAV-mediated MECP2 gene replacement -- challenging because "
            "MECP2 overexpression is toxic (MECP2 duplication syndrome).  "
            "Requires tight dosage control via endogenous promoter or "
            "miRNA-regulated transgene.  Neurogene/Taysha NGN-401 uses "
            "self-complementary AAV9 with engineered MECP2 promoter fragment "
            "and miR-responsive element for dose control. "
            "(2) X-chromosome reactivation -- CRISPRa or epigenetic editing "
            "to reactivate the silenced wild-type MECP2 on the inactive X "
            "chromosome (Xi).  Targets XIST RNA or Xi-specific epigenetic "
            "marks.  Broad Institute/MIT work on CRISPR-mediated Xi "
            "reactivation published 2020. "
            "(3) Base editing -- adenine or cytosine base editors to correct "
            "the 8 recurrent missense/nonsense mutations, which are "
            "overwhelmingly C>T transitions at CpG sites (ideal ABE targets). "
            "(4) RNA editing -- ADAR-based RNA editing to correct G>A "
            "mutations at the transcript level."
        ),
        "clinical_trials": [
            {
                "product": "NGN-401 (Neurogene)",
                "sponsor": "Neurogene Inc. (formerly Taysha Gene Therapies)",
                "phase": "I/II",
                "nct": "NCT05898620",
                "approach": (
                    "Intrathecal scAAV9 delivering MECP2 cDNA with engineered "
                    "miniMECP2 promoter + 3'UTR miRNA-responsive element "
                    "(miR-Responsive Auto-Regulatory Element, miRARE) to "
                    "prevent overexpression toxicity."
                ),
                "status_2026": (
                    "Phase I/II initiated mid-2023 for girls age 4-10 with "
                    "classic Rett and specific MECP2 mutations.  Dose "
                    "escalation ongoing.  FDA granted Rare Pediatric Disease "
                    "Designation.  Interim data expected H2 2025.  Key "
                    "challenge is achieving sufficient transduction of "
                    "cortical neurons while preventing overexpression in "
                    "cells where Xi carries the normal allele."
                ),
            },
        ],
        "conditions": ["rett_syndrome", "neurological", "neurodevelopmental"],
    },

    # -------------------------------------------------------------------
    # 1b.  Fragile X Syndrome -- FMR1
    # -------------------------------------------------------------------
    "FMR1": {
        "gene_id": 2332,
        "chrom": "chrX",
        "start": 147_911_919,
        "end": 147_951_125,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq27.3",
        "exon_count": 17,
        "role": (
            "Fragile X messenger ribonucleoprotein 1 (FMRP) -- RNA-binding "
            "protein that regulates translation of ~800 mRNAs at synapses.  "
            "Essential for synaptic plasticity, dendritic spine maturation, "
            "and mGluR-dependent long-term depression.  CGG trinucleotide "
            "repeat expansion in the 5'UTR causes hypermethylation and "
            "transcriptional silencing (full mutation >=200 repeats)."
        ),
        "disease": "Fragile X Syndrome",
        "omim_disease": 300624,
        "omim_gene": 309550,
        "inheritance": "XLD (trinucleotide repeat expansion)",
        "mutation_type": "trinucleotide_repeat_expansion",
        "repeat_details": {
            "repeat_unit": "CGG",
            "location": "5'UTR of FMR1 exon 1",
            "normal_range": "5-44 repeats",
            "intermediate": "45-54 repeats (gray zone)",
            "premutation": "55-200 repeats (FXTAS, FXPOI risk; unstable expansion)",
            "full_mutation": ">=200 repeats (hypermethylation, gene silencing)",
            "protein_effect": (
                "Full mutation causes CpG methylation of the CGG repeat and "
                "upstream FMR1 promoter, leading to transcriptional silencing "
                "and absence of FMRP.  Premutation alleles produce elevated "
                "FMR1 mRNA with toxic RNA gain-of-function (FXTAS)."
            ),
        },
        "key_variants": [
            {
                "name": "CGG repeat expansion (>=200 repeats, full mutation)",
                "consequence": "transcriptional_silencing",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Not a SNP -- trinucleotide repeat expansion with "
                    "epigenetic silencing.  Expansion from premutation to "
                    "full mutation occurs during maternal transmission.  "
                    "Mosaicism (methylation or repeat length) present in "
                    "~20-40% of affected males."
                ),
            },
        ],
        "mutational_landscape": (
            "Unlike most genetic diseases, Fragile X is almost exclusively "
            "caused by a single mutational mechanism: CGG repeat expansion "
            "in the FMR1 5'UTR.  Rare cases involve point mutations or "
            "deletions within FMR1 coding sequence.  Premutation carriers "
            "(55-200 repeats) number ~1 in 200 females and ~1 in 400 males."
        ),
        "crispr_strategy": (
            "Two main CRISPR approaches: "
            "(1) Epigenetic reactivation -- dCas9 fused to TET1 demethylase "
            "(or p300/VP64) targeted to the methylated FMR1 promoter/CGG "
            "repeat to remove DNA methylation and restore FMR1 transcription.  "
            "Liu et al. (Cell 2018) demonstrated dCas9-TET1 reactivation of "
            "FMR1 in iPSC-derived neurons with sustained expression.  Key "
            "challenge: delivery to sufficient neurons in vivo, and ensuring "
            "sustained demethylation without re-silencing. "
            "(2) CGG repeat excision -- dual Cas9 guides flanking the expanded "
            "CGG repeat to excise the expansion, leaving a normal-length "
            "repeat.  Park et al. (Cell Reports 2015) demonstrated excision "
            "in iPSCs.  Risk: large deletion may disrupt promoter elements. "
            "(3) Gene replacement -- AAV-delivered FMR1 cDNA (~2.1 kb for "
            "major isoform).  Challenge: FMRP is expressed in most neurons "
            "and requires precise dosage regulation."
        ),
        "clinical_trials": [
            {
                "product": "None as of early 2026 for gene/CRISPR therapy",
                "status_2026": (
                    "No clinical trials for CRISPR or gene therapy in Fragile X "
                    "as of early 2026.  Preclinical work on dCas9-TET1 "
                    "epigenetic reactivation is the most advanced CRISPR "
                    "approach (Bhatt lab, UC Davis).  Small molecule trials "
                    "targeting mGluR5 pathway (mavoglurant, basimglurant) "
                    "failed in Phase II/III.  Current clinical focus: "
                    "gaboxadol (GABA-A agonist) Phase III, zatolmilast "
                    "(PDE4D inhibitor) Phase II."
                ),
            },
        ],
        "conditions": ["fragile_x", "neurological", "neurodevelopmental",
                       "intellectual_disability"],
    },

    # -------------------------------------------------------------------
    # 1c.  Angelman Syndrome -- UBE3A
    # -------------------------------------------------------------------
    "UBE3A": {
        "gene_id": 7337,
        "chrom": "chr15",
        "start": 25_333_728,
        "end": 25_439_008,
        "strand": "-",
        "refseq": "NC_000015.10",
        "cytoband": "15q11.2",
        "exon_count": 16,
        "role": (
            "Ubiquitin protein ligase E3A -- HECT domain E3 ubiquitin ligase "
            "with dual function: (1) protein ubiquitination/degradation and "
            "(2) transcriptional coactivation of steroid hormone receptors.  "
            "Imprinted gene: expressed exclusively from the maternal allele "
            "in neurons (paternal allele silenced by UBE3A-ATS long non-coding "
            "RNA antisense transcript).  Loss of maternal UBE3A causes "
            "Angelman syndrome."
        ),
        "disease": "Angelman Syndrome",
        "omim_disease": 105830,
        "omim_gene": 601623,
        "inheritance": "Imprinted (maternal expression in neurons)",
        "mutation_type": "deletion_or_point_mutation_or_imprinting_defect",
        "key_variants": [
            {
                "name": "15q11.2-q13 maternal deletion (~5-7 Mb)",
                "consequence": "gene_deletion",
                "clinical_significance": "pathogenic",
                "frequency": "~70% of Angelman cases",
                "notes": (
                    "Large deletion of maternal 15q11.2-q13 encompassing "
                    "UBE3A and surrounding imprinted genes.  Most common "
                    "mechanism.  Deletion class I (~6-7 Mb, BP1-BP3) or "
                    "class II (~5-6 Mb, BP2-BP3).  More severe phenotype "
                    "than other mechanisms due to haploinsufficiency of "
                    "additional genes (GABRB3, GABRA5, GABRG3)."
                ),
            },
            {
                "name": "p.Arg506Cys (R506C)",
                "rsid": "rs587784338",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Catalytic HECT domain; abolishes E3 ligase activity.",
            },
        ],
        "mutational_landscape": (
            "Angelman syndrome has four molecular classes: "
            "(1) Maternal 15q11.2-q13 deletion (~70%); "
            "(2) UBE3A point mutations (~11%); "
            "(3) Paternal uniparental disomy (UPD) (~3-7%); "
            "(4) Imprinting center defect (~3%).  "
            "~10% of clinically diagnosed cases have no identifiable "
            "molecular mechanism."
        ),
        "crispr_strategy": (
            "Primary CRISPR approach: Unsilencing the intact paternal UBE3A "
            "allele by disrupting the UBE3A antisense transcript (UBE3A-ATS).  "
            "The paternal UBE3A gene is structurally intact but silenced by "
            "the overlapping antisense transcript.  CRISPR disruption of "
            "UBE3A-ATS or its regulatory elements can derepress paternal "
            "UBE3A expression in neurons.  "
            "Approaches: (1) CRISPRi targeting UBE3A-ATS promoter (Snord116/ "
            "IC region); (2) Cas9 cleavage of UBE3A-ATS to terminate "
            "transcription; (3) ASO-mediated UBE3A-ATS knockdown (GeneTx/ "
            "Ultragenyx GTX-102, clinical).  "
            "For deletion cases: gene replacement with AAV-UBE3A (~2.6 kb "
            "cDNA fits AAV) under neuron-specific promoter."
        ),
        "clinical_trials": [
            {
                "product": "GTX-102 (Ultragenyx/GeneTx)",
                "sponsor": "Ultragenyx Pharmaceutical",
                "phase": "I/II",
                "nct": "NCT04259281",
                "approach": (
                    "Intrathecal ASO targeting UBE3A-ATS to unsilence "
                    "paternal UBE3A.  Not CRISPR-based but validates the "
                    "antisense knockdown mechanism."
                ),
                "status_2026": (
                    "Phase I/II reinitiated after initial safety pause (lower "
                    "extremity weakness, attributed to drug distribution).  "
                    "Lumbar puncture dosing with modified protocol.  Dose "
                    "escalation ongoing.  Preliminary data suggest "
                    "UBE3A protein restoration in CSF biomarkers."
                ),
            },
        ],
        "conditions": ["angelman_syndrome", "neurological", "neurodevelopmental",
                       "imprinting_disorder"],
    },

    # -------------------------------------------------------------------
    # 1d.  Dravet Syndrome -- SCN1A
    # -------------------------------------------------------------------
    "SCN1A": {
        "gene_id": 6323,
        "chrom": "chr2",
        "start": 165_984_641,
        "end": 166_149_160,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q24.3",
        "exon_count": 26,
        "role": (
            "Sodium voltage-gated channel alpha subunit 1 (Nav1.1) -- "
            "primary sodium channel in GABAergic inhibitory interneurons.  "
            "Loss-of-function in interneurons reduces inhibitory tone, "
            "causing network hyperexcitability and seizures.  Paradoxically, "
            "sodium channel blockers (carbamazepine, phenytoin) worsen Dravet "
            "because they further suppress inhibitory interneuron firing."
        ),
        "disease": "Dravet Syndrome (Severe Myoclonic Epilepsy of Infancy)",
        "omim_disease": 607208,
        "omim_gene": 182389,
        "inheritance": "AD (de novo in >90%)",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "p.Arg1648His (R1648H)",
                "hgvs_coding": "NM_001165963.4:c.4943G>A",
                "rsid": "rs121917911",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Domain IV S4 voltage sensor; reduces channel activation.",
            },
            {
                "name": "p.Arg1657Cys (R1657C)",
                "rsid": "rs121917910",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Domain IV S4; recurrent de novo mutation.",
            },
        ],
        "mutational_landscape": (
            ">1,800 pathogenic/likely-pathogenic SCN1A variants in ClinVar.  "
            "~80% are loss-of-function (truncating, splice, frameshift).  "
            "Extreme allelic heterogeneity -- most families carry private "
            "mutations.  De novo in >90% of cases.  Genotype-phenotype: "
            "truncating variants generally cause more severe Dravet; missense "
            "variants can cause Dravet or milder GEFS+."
        ),
        "crispr_strategy": (
            "Primary approach: CRISPRa upregulation of the remaining wild-type "
            "SCN1A allele to compensate for haploinsufficiency.  Dravet is "
            "caused by ~50% reduction in Nav1.1 (haploinsufficiency), so "
            "boosting expression from the intact allele can theoretically "
            "restore normal inhibitory neuron function.  "
            "Stoke Therapeutics / Encoded Therapeutics approach: Targeted "
            "Augmentation of Nuclear Gene Output (TANGO) -- ASO targeting "
            "a naturally occurring non-productive exon (poison exon 20N) in "
            "SCN1A pre-mRNA to boost productive mRNA.  ETX101 (Encoded/"
            "Ultragenyx) uses AAV-delivered CRISPRa (dCas9-VPR) or "
            "engineered transcription factor (eTF) to upregulate SCN1A "
            "specifically in interneurons via Dlx5/6 enhancer.  "
            "Gene replacement difficult: SCN1A cDNA is ~6.0 kb, at the "
            "limit for AAV.  Dual AAV or oversized AAV approaches possible."
        ),
        "clinical_trials": [
            {
                "product": "ETX101 (Encoded Therapeutics / Ultragenyx)",
                "sponsor": "Ultragenyx",
                "phase": "I/II",
                "nct": "NCT05901350",
                "approach": (
                    "AAV9-delivered engineered transcription factor (eTF) "
                    "under GABAergic interneuron-specific Dlx5/6 enhancer "
                    "to upregulate endogenous SCN1A.  Intracerebroventricular "
                    "(ICV) delivery."
                ),
                "status_2026": (
                    "Phase I/II initiated late 2023 for children 2-18 years "
                    "with genetically confirmed Dravet (SCN1A loss-of-function).  "
                    "Dose escalation ongoing.  FDA Fast Track and Orphan Drug "
                    "designations.  First engineered transcription factor gene "
                    "therapy to enter clinical trials for any disease."
                ),
            },
            {
                "product": "STK-001 (Stoke Therapeutics)",
                "sponsor": "Stoke Therapeutics",
                "phase": "II",
                "nct": "NCT04740476",
                "approach": (
                    "Intrathecal ASO targeting SCN1A poison exon (TANGO "
                    "mechanism) to increase productive SCN1A mRNA."
                ),
                "status_2026": (
                    "Phase II MONARCH study ongoing for children/adolescents "
                    "2-18 years.  Phase I/II SWALLOWTAIL data showed "
                    "reduction in convulsive seizure frequency.  Not "
                    "CRISPR-based but validates upregulation approach."
                ),
            },
        ],
        "conditions": ["dravet_syndrome", "epilepsy", "neurological",
                       "channelopathy"],
    },

    # -------------------------------------------------------------------
    # 1e.  Tuberous Sclerosis Complex -- TSC1
    # -------------------------------------------------------------------
    "TSC1": {
        "gene_id": 7248,
        "chrom": "chr9",
        "start": 132_891_348,
        "end": 132_945_370,
        "strand": "-",
        "refseq": "NC_000009.12",
        "cytoband": "9q34.13",
        "exon_count": 23,
        "role": (
            "TSC complex subunit 1 (hamartin) -- forms a complex with "
            "tuberin (TSC2) that acts as a GAP (GTPase-activating protein) "
            "for the small GTPase Rheb, thereby inhibiting mTORC1 signalling.  "
            "Loss-of-function leads to constitutive mTORC1 activation, "
            "causing hamartomatous growths in brain, skin, kidneys, heart, "
            "and lungs (lymphangioleiomyomatosis)."
        ),
        "disease": "Tuberous Sclerosis Complex (TSC)",
        "omim_disease": 191100,
        "omim_gene": 605284,
        "inheritance": "AD",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "Various truncating mutations throughout TSC1",
                "consequence": "truncating",
                "clinical_significance": "pathogenic",
                "notes": (
                    "TSC1 mutations account for ~30% of TSC cases.  "
                    "Predominantly nonsense and frameshift mutations "
                    "(>90% truncating).  Missense variants rare in TSC1.  "
                    "Generally milder phenotype than TSC2 mutations."
                ),
            },
        ],
        "crispr_strategy": (
            "mTOR pathway modulation rather than direct gene correction is "
            "the primary therapeutic approach.  mTOR inhibitors (everolimus, "
            "sirolimus) are FDA-approved for TSC-associated subependymal "
            "giant cell astrocytomas (SEGA), renal angiomyolipomas, and "
            "epilepsy.  CRISPR approaches are preclinical: "
            "(1) Base editing or prime editing to correct specific TSC1 "
            "truncating mutations in somatic mosaicism cases; "
            "(2) CRISPRa upregulation of remaining wild-type allele; "
            "(3) Gene replacement (TSC1 cDNA ~3.5 kb fits AAV) with tissue-"
            "specific promoters.  Main challenge: TSC is a multisystem "
            "disease requiring systemic delivery."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": (
                    "Everolimus (Afinitor) and sirolimus are standard of care.  "
                    "CRISPR approaches remain preclinical.  Academic programs "
                    "at Boston Children's and Cincinnati Children's exploring "
                    "AAV-TSC1/TSC2 gene replacement."
                ),
            },
        ],
        "conditions": ["tuberous_sclerosis", "TSC", "neurological", "epilepsy",
                       "tumor_predisposition"],
    },

    # -------------------------------------------------------------------
    # 1f.  Tuberous Sclerosis Complex -- TSC2
    # -------------------------------------------------------------------
    "TSC2": {
        "gene_id": 7249,
        "chrom": "chr16",
        "start": 2_047_936,
        "end": 2_089_491,
        "strand": "+",
        "refseq": "NC_000016.10",
        "cytoband": "16p13.3",
        "exon_count": 42,
        "role": (
            "TSC complex subunit 2 (tuberin) -- GTPase-activating protein "
            "for Rheb that inhibits mTORC1 signalling.  Contains the GAP "
            "domain (C-terminal).  Forms obligate complex with hamartin "
            "(TSC1) and TBC1D7.  TSC2 mutations account for ~70% of TSC "
            "cases and are associated with more severe phenotype than TSC1."
        ),
        "disease": "Tuberous Sclerosis Complex (TSC)",
        "omim_disease": 613254,
        "omim_gene": 191092,
        "inheritance": "AD",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "Various (truncating, missense, splice, large deletions)",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "TSC2 mutations account for ~70% of TSC cases.  More "
                    "diverse mutation spectrum than TSC1: truncating (~60%), "
                    "missense (~20%), splice-site, in-frame deletions, and "
                    "large genomic rearrangements.  Contiguous gene deletion "
                    "of TSC2 + PKD1 (adjacent on 16p13.3) causes TSC with "
                    "early-onset polycystic kidney disease."
                ),
            },
        ],
        "crispr_strategy": (
            "Similar to TSC1.  TSC2 cDNA (~5.4 kb) is at the upper limit "
            "for single AAV packaging.  Dual AAV split-intein approach may "
            "be necessary.  Base/prime editing for specific missense "
            "mutations in the GAP domain.  Somatic gene correction "
            "approaches under preclinical investigation."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": (
                    "Preclinical only.  mTOR inhibitors remain standard of "
                    "care.  PKD1-TSC2 contiguous deletion cases require "
                    "special consideration."
                ),
            },
        ],
        "conditions": ["tuberous_sclerosis", "TSC", "neurological", "epilepsy",
                       "tumor_predisposition"],
    },

    # -------------------------------------------------------------------
    # 1g.  Neurofibromatosis Type 1 -- NF1
    # -------------------------------------------------------------------
    "NF1": {
        "gene_id": 4763,
        "chrom": "chr17",
        "start": 31_094_927,
        "end": 31_377_677,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q11.2",
        "exon_count": 58,
        "role": (
            "Neurofibromin 1 -- RasGAP protein that negatively regulates "
            "RAS/MAPK signalling by catalysing GTP hydrolysis of RAS.  "
            "Tumor suppressor: loss-of-function leads to constitutive RAS "
            "activation.  Very large gene (~283 kb, one of the largest "
            "human genes).  NF1 is the most common genetic tumor "
            "predisposition syndrome (~1 in 3,000)."
        ),
        "disease": "Neurofibromatosis Type 1 (NF1 / von Recklinghausen disease)",
        "omim_disease": 162200,
        "omim_gene": 613113,
        "inheritance": "AD",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "Extreme allelic heterogeneity (>3,000 variants reported)",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "No single hotspot mutation.  ~5% have whole-gene "
                    "microdeletion (17q11.2 microdeletion syndrome) with "
                    "more severe phenotype.  Truncating mutations predominate "
                    "(nonsense, frameshift, splice).  High de novo mutation "
                    "rate (~50%)."
                ),
            },
        ],
        "crispr_strategy": (
            "NF1 cDNA is ~8.5 kb, exceeding single AAV capacity.  "
            "Approaches: (1) Dual AAV split-intein for gene replacement; "
            "(2) MEK inhibitors (selumetinib, FDA-approved for NF1 "
            "plexiform neurofibromas) are current standard of care; "
            "(3) Base/prime editing for specific truncating variants; "
            "(4) Exon skipping for select in-frame deletion targets.  "
            "CRISPR approaches remain early preclinical."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": (
                    "Selumetinib (Koselugo) FDA-approved 2020 for NF1 "
                    "plexiform neurofibromas in children >=2 years.  "
                    "Gene therapy approaches in early preclinical stage.  "
                    "University of Florida program exploring dual AAV for "
                    "NF1 Schwann cell-targeted gene replacement."
                ),
            },
        ],
        "conditions": ["neurofibromatosis", "NF1", "neurological",
                       "tumor_predisposition", "RASopathy"],
    },

    # -------------------------------------------------------------------
    # 1h.  Neurofibromatosis Type 2 -- NF2 (MERLIN/Schwannomin)
    # -------------------------------------------------------------------
    "NF2": {
        "gene_id": 4771,
        "chrom": "chr22",
        "start": 29_603_556,
        "end": 29_698_600,
        "strand": "+",
        "refseq": "NC_000022.11",
        "cytoband": "22q12.2",
        "exon_count": 17,
        "role": (
            "NF2 (merlin/schwannomin) -- cytoskeletal linker protein of "
            "the ERM (ezrin-radixin-moesin) family.  Functions as tumor "
            "suppressor by regulating Hippo/YAP, mTOR, and Rac/PAK "
            "signalling pathways.  Loss-of-function causes bilateral "
            "vestibular schwannomas (hallmark), meningiomas, and "
            "ependymomas."
        ),
        "disease": "Neurofibromatosis Type 2 (NF2-related schwannomatosis)",
        "omim_disease": 101000,
        "omim_gene": 607379,
        "inheritance": "AD",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "Various truncating mutations",
                "consequence": "truncating",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Truncating mutations (nonsense, frameshift) associated "
                    "with severe phenotype (early onset, more tumors).  "
                    "Missense/in-frame mutations and splice-site variants "
                    "associated with milder disease.  ~50% de novo."
                ),
            },
        ],
        "crispr_strategy": (
            "NF2/merlin cDNA is ~1.8 kb (fits AAV).  Approaches: "
            "(1) AAV-mediated NF2 gene replacement with Schwann cell "
            "promoter (P0/MPZ) for targeted expression; "
            "(2) Base editing for specific truncating variants; "
            "(3) CRISPRa upregulation of wild-type allele.  "
            "Early preclinical only.  Current treatment is surgery, "
            "radiation, and bevacizumab."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": (
                    "Preclinical AAV-NF2 gene replacement studies at "
                    "Massachusetts Eye and Ear / Harvard.  "
                    "Bevacizumab (anti-VEGF) used off-label for "
                    "vestibular schwannoma volume reduction."
                ),
            },
        ],
        "conditions": ["neurofibromatosis", "NF2", "neurological",
                       "tumor_predisposition", "schwannomatosis"],
    },

    # -------------------------------------------------------------------
    # 1i.  Charcot-Marie-Tooth Type 1A -- PMP22
    # -------------------------------------------------------------------
    "PMP22": {
        "gene_id": 5376,
        "chrom": "chr17",
        "start": 15_230_987,
        "end": 15_265_328,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17p12",
        "exon_count": 5,
        "role": (
            "Peripheral myelin protein 22 -- integral membrane glycoprotein "
            "comprising ~2-5% of peripheral nerve myelin.  Essential for "
            "Schwann cell proliferation, differentiation, and myelin "
            "compaction.  CMT1A caused by 1.5 Mb duplication at 17p12 "
            "containing PMP22 (gene dosage: 3 copies instead of 2).  "
            "PMP22 deletion causes HNPP (hereditary neuropathy with "
            "liability to pressure palsies)."
        ),
        "disease": "Charcot-Marie-Tooth Type 1A (CMT1A)",
        "omim_disease": 118220,
        "omim_gene": 601097,
        "inheritance": "AD (gene duplication)",
        "mutation_type": "gene_duplication",
        "key_variants": [
            {
                "name": "17p12 duplication (1.5 Mb, PMP22 gene dosage)",
                "consequence": "gene_duplication",
                "clinical_significance": "pathogenic",
                "frequency": "~70-80% of CMT1A",
                "notes": (
                    "Most common inherited neuropathy (~1 in 2,500).  "
                    "Duplication arises from unequal crossing-over between "
                    "flanking CMT1A-REP repeats.  Rarely, PMP22 point "
                    "mutations cause CMT1A (severe) or Dejerine-Sottas (DSN)."
                ),
            },
        ],
        "crispr_strategy": (
            "Gene dosage reduction approaches: "
            "(1) CRISPRi silencing of one PMP22 copy to normalize gene "
            "dosage from 3 copies to ~2.  dCas9-KRAB targeted to PMP22 "
            "promoter.  Challenge: delivering to Schwann cells in "
            "peripheral nerves. "
            "(2) Allele-specific silencing of one duplicated copy using "
            "SNP-based discrimination (if heterozygous). "
            "(3) CRISPR excision of the duplicated segment (deletion of "
            "one 1.5 Mb copy).  Technically challenging. "
            "Preclinical: AAV-mediated delivery of CRISPRi or shRNA "
            "targeting PMP22 in Schwann cells using P0/MPZ promoter. "
            "PMP22 dosage must be precisely balanced: too little causes "
            "HNPP, too much causes CMT1A."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": (
                    "Preclinical.  Small molecule approaches: PXT3003 "
                    "(combination baclofen/naltrexone/sorbitol) Phase III "
                    "(Pharnext/Pherecydes).  AAV-based PMP22 regulation "
                    "at early preclinical stage (University of Florida, "
                    "Nationwide Children's Hospital)."
                ),
            },
        ],
        "conditions": ["charcot_marie_tooth", "CMT1A", "neurological",
                       "neuropathy", "peripheral_neuropathy"],
    },

    # -------------------------------------------------------------------
    # 1j.  Charcot-Marie-Tooth Type 1X -- GJB1 (Connexin 32)
    # -------------------------------------------------------------------
    "GJB1": {
        "gene_id": 2705,
        "chrom": "chrX",
        "start": 71_223_251,
        "end": 71_233_378,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq13.1",
        "exon_count": 2,
        "role": (
            "Gap junction protein beta 1 (connexin 32, Cx32) -- forms gap "
            "junctions in Schwann cells (and hepatocytes).  In peripheral "
            "nerve, Cx32 gap junctions provide a radial pathway across "
            "myelin layers.  Loss-of-function causes CMTX1 (second most "
            "common CMT subtype)."
        ),
        "disease": "Charcot-Marie-Tooth Type X1 (CMTX1)",
        "omim_disease": 302800,
        "omim_gene": 304040,
        "inheritance": "XLD",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": ">460 different GJB1 mutations reported",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Allelic heterogeneity.  Missense, nonsense, frameshift, "
                    "and non-coding regulatory mutations.  Males more severely "
                    "affected.  Females variably affected due to X-inactivation."
                ),
            },
        ],
        "crispr_strategy": (
            "AAV-mediated GJB1 gene replacement (cDNA ~0.85 kb, ideal AAV "
            "target).  Preclinical AAV9 gene replacement in Gjb1-null mice "
            "demonstrated rescue of neuropathy phenotype (Kleopa lab, Cyprus "
            "Institute of Neurology).  Schwann cell-specific expression via "
            "MPZ promoter.  CRISPR correction of individual mutations also "
            "feasible given small gene size.  Gene replacement preferred "
            "due to extreme allelic heterogeneity."
        ),
        "clinical_trials": [
            {
                "product": "No clinical trials as of early 2026",
                "status_2026": (
                    "Preclinical AAV gene replacement is the most advanced "
                    "approach.  IND-enabling studies ongoing.  Anticipated "
                    "first-in-human within 2-3 years."
                ),
            },
        ],
        "conditions": ["charcot_marie_tooth", "CMTX1", "neurological",
                       "neuropathy", "peripheral_neuropathy"],
    },

    # -------------------------------------------------------------------
    # 1k.  Charcot-Marie-Tooth Type 2A -- MFN2
    # -------------------------------------------------------------------
    "MFN2": {
        "gene_id": 9927,
        "chrom": "chr1",
        "start": 11_980_181,
        "end": 12_013_491,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1p36.22",
        "exon_count": 19,
        "role": (
            "Mitofusin 2 -- mitochondrial outer membrane GTPase required for "
            "mitochondrial fusion.  Also mediates ER-mitochondria tethering "
            "and mitochondrial transport along axons.  Dominant mutations "
            "cause CMT2A, the most common axonal CMT subtype.  Some "
            "recessive mutations cause a severe, early-onset neuropathy."
        ),
        "disease": "Charcot-Marie-Tooth Type 2A (CMT2A)",
        "omim_disease": 609260,
        "omim_gene": 608507,
        "inheritance": "AD (occasionally AR for severe forms)",
        "mutation_type": "dominant_negative_or_loss_of_function",
        "key_variants": [
            {
                "name": "p.Arg94Gln (R94Q)",
                "rsid": "rs119103262",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common MFN2 mutation (~10% of CMT2A).  GTPase "
                    "domain.  Early onset, severe axonal neuropathy."
                ),
            },
        ],
        "crispr_strategy": (
            "Allele-specific silencing of dominant-negative MFN2 mutant "
            "allele while preserving wild-type.  MFN2 cDNA (~2.3 kb) fits "
            "AAV for suppression-and-replacement approach.  Preclinical "
            "only.  Challenge: peripheral nerve/axonal delivery.  "
            "Alternative: gene replacement with rescue of mitochondrial "
            "transport."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": "Preclinical only.",
            },
        ],
        "conditions": ["charcot_marie_tooth", "CMT2A", "neurological",
                       "neuropathy", "peripheral_neuropathy",
                       "mitochondrial"],
    },

    # -------------------------------------------------------------------
    # 1l.  Canavan Disease -- ASPA
    # -------------------------------------------------------------------
    "ASPA": {
        "gene_id": 443,
        "chrom": "chr17",
        "start": 3_474_937,
        "end": 3_504_377,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17p13.2",
        "exon_count": 6,
        "role": (
            "Aspartoacylase -- enzyme that hydrolyses N-acetylaspartate (NAA) "
            "to aspartate and acetate in oligodendrocytes.  The acetate is "
            "used for myelin lipid synthesis.  Deficiency causes NAA "
            "accumulation and spongiform leukodystrophy (Canavan disease) "
            "with progressive white matter vacuolation."
        ),
        "disease": "Canavan Disease (spongiform leukodystrophy)",
        "omim_disease": 271900,
        "omim_gene": 608034,
        "inheritance": "AR",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "p.Glu285Ala (E285A)",
                "hgvs_coding": "NM_000049.4:c.854A>C",
                "rsid": "rs28940279",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "frequency": "Most common worldwide; >98% in Ashkenazi Jewish",
                "notes": (
                    "Founder mutation in Ashkenazi Jewish population.  "
                    "Carrier frequency ~1/40 in AJ.  Virtually no residual "
                    "enzyme activity."
                ),
            },
            {
                "name": "p.Tyr231Ter (Y231X)",
                "hgvs_coding": "NM_000049.4:c.693C>A",
                "rsid": "rs28940280",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Second most common; pan-ethnic.",
            },
        ],
        "crispr_strategy": (
            "Gene replacement is the primary approach.  ASPA cDNA is ~0.9 kb "
            "(ideal AAV target).  AAV-ASPA gene therapy is one of the most "
            "advanced gene therapies for leukodystrophies: "
            "(1) Myrtelle rAAV-Olig001-ASPA (MYT-01) Phase I/II -- "
            "AAV-Olig001 (oligodendrocyte-tropic capsid) delivering ASPA "
            "under CBh promoter.  Intracerebroventricular delivery. "
            "(2) Base editing to correct E285A (A>C) could use CBE. "
            "Preclinical studies show significant NAA reduction and "
            "myelination improvement after AAV-ASPA delivery."
        ),
        "clinical_trials": [
            {
                "product": "rAAV-Olig001-ASPA (MYT-01, Myrtelle)",
                "sponsor": "Myrtelle Inc.",
                "phase": "I/II",
                "nct": "NCT04833907",
                "approach": (
                    "AAV-Olig001 serotype (engineered for oligodendrocyte "
                    "tropism) delivering ASPA cDNA.  ICV administration."
                ),
                "status_2026": (
                    "Phase I/II ongoing for children 6 months to 6 years.  "
                    "Dose escalation complete.  Interim data (2025): NAA "
                    "levels reduced on MR spectroscopy; myelination "
                    "improvement on brain MRI.  FDA Fast Track and Rare "
                    "Pediatric Disease designations."
                ),
            },
        ],
        "conditions": ["canavan_disease", "leukodystrophy", "neurological",
                       "metabolic", "lysosomal"],
    },
}


# ============================================================================
# 2. METABOLIC RARE DISEASES
# ============================================================================

METABOLIC_RARE_TARGETS = {

    # -------------------------------------------------------------------
    # 2a.  Maple Syrup Urine Disease -- BCKDHA
    # -------------------------------------------------------------------
    "BCKDHA": {
        "gene_id": 593,
        "chrom": "chr19",
        "start": 41_395_063,
        "end": 41_416_270,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19q13.2",
        "exon_count": 9,
        "role": (
            "Branched chain keto acid dehydrogenase E1 subunit alpha -- "
            "alpha subunit of the E1 component of the mitochondrial "
            "branched-chain alpha-keto acid dehydrogenase complex (BCKDH).  "
            "This complex catalyses the rate-limiting oxidative "
            "decarboxylation of branched-chain alpha-keto acids derived "
            "from leucine, isoleucine, and valine.  Deficiency causes "
            "MSUD type Ia."
        ),
        "disease": "Maple Syrup Urine Disease Type Ia (MSUD)",
        "omim_disease": 248600,
        "omim_gene": 608348,
        "inheritance": "AR",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "p.Tyr393Asn (Y393N)",
                "rsid": "rs121964987",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Thiamine-responsive MSUD variant.  ~25% residual "
                    "enzyme activity."
                ),
            },
        ],
        "crispr_strategy": (
            "Liver-directed gene replacement or editing.  BCKDH complex is "
            "primarily active in liver.  Approaches: "
            "(1) AAV-mediated liver gene replacement (cDNA ~1.3 kb, ideal "
            "AAV target); "
            "(2) LNP-delivered base editing for recurrent missense variants; "
            "(3) Liver-directed mRNA therapy (Moderna has explored).  "
            "Liver transplantation is curative for MSUD (the liver alone "
            "provides ~10% of BCKDH activity, sufficient to prevent crises), "
            "supporting liver-directed gene therapy rationale."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": (
                    "Preclinical AAV gene therapy (Horizon Therapeutics, "
                    "academic programs).  mRNA therapy in preclinical "
                    "development.  Liver transplant remains definitive "
                    "treatment for severe classic MSUD."
                ),
            },
        ],
        "conditions": ["MSUD", "maple_syrup_urine_disease", "metabolic",
                       "amino_acid_disorder"],
    },

    # -------------------------------------------------------------------
    # 2b.  Maple Syrup Urine Disease -- BCKDHB
    # -------------------------------------------------------------------
    "BCKDHB": {
        "gene_id": 594,
        "chrom": "chr6",
        "start": 80_818_586,
        "end": 80_888_622,
        "strand": "-",
        "refseq": "NC_000006.12",
        "cytoband": "6q14.1",
        "exon_count": 11,
        "role": (
            "Branched chain keto acid dehydrogenase E1 subunit beta -- "
            "beta subunit of the E1 component of the BCKDH complex.  "
            "E1 alpha and beta subunits form an alpha2-beta2 heterotetramer.  "
            "Mutations cause MSUD type Ib."
        ),
        "disease": "Maple Syrup Urine Disease Type Ib (MSUD)",
        "omim_disease": 248600,
        "omim_gene": 248611,
        "inheritance": "AR",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "Various loss-of-function mutations",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Less common than BCKDHA mutations.  Same clinical "
                    "presentation as type Ia."
                ),
            },
        ],
        "crispr_strategy": (
            "Same liver-directed approaches as BCKDHA.  BCKDHB cDNA ~1.2 kb "
            "(ideal AAV target).  LNP-delivered base editing feasible for "
            "recurrent variants."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": "Preclinical only.",
            },
        ],
        "conditions": ["MSUD", "maple_syrup_urine_disease", "metabolic",
                       "amino_acid_disorder"],
    },

    # -------------------------------------------------------------------
    # 2c.  Maple Syrup Urine Disease -- DBT (E2 subunit)
    # -------------------------------------------------------------------
    "DBT": {
        "gene_id": 1629,
        "chrom": "chr1",
        "start": 100_189_037,
        "end": 100_237_635,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p21.2",
        "exon_count": 11,
        "role": (
            "Dihydrolipoamide branched chain transacylase E2 -- core "
            "structural and catalytic subunit of the BCKDH complex.  "
            "24 copies form the cubic inner core.  Mutations cause MSUD "
            "type II."
        ),
        "disease": "Maple Syrup Urine Disease Type II (MSUD)",
        "omim_disease": 248600,
        "omim_gene": 248610,
        "inheritance": "AR",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "Various loss-of-function mutations",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Least common MSUD subtype.  Includes thiamine-responsive "
                    "forms."
                ),
            },
        ],
        "crispr_strategy": (
            "Liver-directed gene replacement.  DBT cDNA ~1.4 kb (fits AAV).  "
            "Same delivery rationale as BCKDHA/B: liver transplant is curative, "
            "supporting liver-directed gene therapy."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": "Preclinical only.",
            },
        ],
        "conditions": ["MSUD", "maple_syrup_urine_disease", "metabolic",
                       "amino_acid_disorder"],
    },

    # -------------------------------------------------------------------
    # 2d.  Tyrosinemia Type 1 -- FAH
    # -------------------------------------------------------------------
    "FAH": {
        "gene_id": 2184,
        "chrom": "chr15",
        "start": 80_154_851,
        "end": 80_187_567,
        "strand": "+",
        "refseq": "NC_000015.10",
        "cytoband": "15q25.1",
        "exon_count": 14,
        "role": (
            "Fumarylacetoacetate hydrolase -- terminal enzyme in tyrosine "
            "catabolism pathway.  Catalyses hydrolysis of fumarylacetoacetate "
            "(FAA) to fumarate and acetoacetate.  Deficiency causes "
            "accumulation of FAA and its metabolite succinylacetone (SA), "
            "which are hepatotoxic and nephrotoxic.  SA also inhibits "
            "porphobilinogen synthase causing porphyria-like crises."
        ),
        "disease": "Tyrosinemia Type 1 (Hepatorenal Tyrosinemia, HT1)",
        "omim_disease": 276700,
        "omim_gene": 613871,
        "inheritance": "AR",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "IVS12+5G>A (c.1062+5G>A, Quebec founder mutation)",
                "hgvs_coding": "NM_000137.4:c.1062+5G>A",
                "rsid": "rs121964992",
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "frequency": "Most common in French-Canadian (Quebec) population",
                "notes": (
                    "Founder effect in Saguenay-Lac-Saint-Jean region of "
                    "Quebec.  Carrier frequency ~1/20 in this region.  "
                    "Splice site variant causing exon 12 skipping."
                ),
            },
            {
                "name": "p.Pro261Leu (P261L)",
                "rsid": "rs121964990",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Northern European / Scandinavian populations.",
            },
        ],
        "mutational_landscape": (
            ">100 pathogenic FAH variants reported.  Splice-site and missense "
            "predominate.  Strong founder effects: IVS12+5G>A in Quebec, "
            "W262X in Finland, P261L in Scandinavia.  Spontaneous reversion "
            "to normal hepatocytes occurs in ~50% of patients (natural gene "
            "therapy via somatic mosaicism in liver)."
        ),
        "crispr_strategy": (
            "FAH is a prime CRISPR gene therapy target because: "
            "(1) Corrected hepatocytes have a selective growth advantage "
            "(self-corrected nodules observed naturally); "
            "(2) Liver-directed delivery is well-established (LNP or AAV); "
            "(3) FAH cDNA ~1.3 kb (ideal AAV target); "
            "(4) Nitisinone (NTBC) blocks the upstream enzyme (4-HPPD) and "
            "prevents toxic metabolite accumulation, serving as a bridge "
            "to gene therapy.  "
            "Approaches: (1) In vivo liver-directed base editing via LNP; "
            "(2) AAV-FAH gene replacement; (3) CRISPR-Cas9 HDR correction "
            "in hepatocytes (demonstrated in Fah-/- mouse model by Yin et "
            "al., Nature Biotechnology 2014 -- landmark paper for in vivo "
            "CRISPR gene correction).  "
            "NOTE: Tyrosinemia type 1 was used in the original Yin et al. "
            "2014 proof-of-concept for in vivo CRISPR gene correction."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR clinical trials as of early 2026",
                "status_2026": (
                    "Despite being a landmark preclinical CRISPR target "
                    "(Yin et al. 2014), no clinical trials have been "
                    "initiated.  Nitisinone + dietary restriction is highly "
                    "effective, reducing urgency.  Preclinical programs "
                    "ongoing at several academic centers."
                ),
            },
        ],
        "conditions": ["tyrosinemia", "HT1", "metabolic", "liver",
                       "amino_acid_disorder"],
    },

    # -------------------------------------------------------------------
    # 2e.  Glycogen Storage Disease Type Ia -- G6PC1
    # -------------------------------------------------------------------
    "G6PC1": {
        "gene_id": 2538,
        "chrom": "chr17",
        "start": 42_900_796,
        "end": 42_914_437,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.31",
        "exon_count": 5,
        "role": (
            "Glucose-6-phosphatase catalytic subunit 1 -- ER-membrane-bound "
            "enzyme that catalyses the terminal step of gluconeogenesis and "
            "glycogenolysis (hydrolysis of glucose-6-phosphate to glucose + "
            "phosphate).  Expressed primarily in liver, kidney, and intestine.  "
            "Deficiency causes GSD Ia (von Gierke disease): severe fasting "
            "hypoglycemia, hepatomegaly, lactic acidosis, hyperuricemia, "
            "hyperlipidemia, and hepatocellular adenoma/carcinoma risk."
        ),
        "disease": "Glycogen Storage Disease Type Ia (GSD Ia, von Gierke disease)",
        "omim_disease": 232200,
        "omim_gene": 613742,
        "inheritance": "AR",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "p.Arg83Cys (R83C)",
                "rsid": "rs80356484",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "frequency": "Most common in Caucasian populations (~30% of alleles)",
                "notes": (
                    "Active site region.  <3% residual enzyme activity.  "
                    "Severe classic GSD Ia phenotype."
                ),
            },
            {
                "name": "p.Gln347Ter (Q347X)",
                "rsid": "rs80356492",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Common in Ashkenazi Jewish population (~100% of alleles).",
            },
            {
                "name": "c.648G>T (p.Leu216Leu / aberrant splicing)",
                "rsid": "rs80356486",
                "consequence": "synonymous_splice",
                "clinical_significance": "pathogenic",
                "notes": "Most common in East Asian populations (Chinese, Japanese, Korean).",
            },
        ],
        "crispr_strategy": (
            "AAV gene replacement is the leading approach: "
            "(1) Ultragenyx DTX401 (AAV8-G6PC) Phase III trial -- "
            "liver-directed AAV8 delivering human G6PC1 cDNA (~1.1 kb, "
            "ideal AAV target) under hybrid liver promoter. "
            "(2) Base editing for recurrent R83C (C>T, ideal for ABE). "
            "(3) mRNA therapy (Moderna mRNA-3745) in clinical development. "
            "G6PC1 gene therapy has strong clinical rationale: disease is "
            "liver-centric, gene is small, and dietary management is "
            "burdensome (continuous glucose feeds to prevent hypoglycemia)."
        ),
        "clinical_trials": [
            {
                "product": "DTX401 (Ultragenyx)",
                "sponsor": "Ultragenyx Pharmaceutical",
                "phase": "III",
                "nct": "NCT05139316",
                "approach": (
                    "AAV8 vector delivering human G6PC1 cDNA under hybrid "
                    "liver-specific promoter.  Single IV infusion."
                ),
                "status_2026": (
                    "Phase III pivotal trial ongoing for adults with GSD Ia.  "
                    "Phase I/II data showed dose-dependent improvements in "
                    "fasting tolerance, lactate levels, and reduced dietary "
                    "starch requirements.  One of the most advanced gene "
                    "therapy programs for a metabolic liver disease.  "
                    "BLA submission anticipated 2026-2027."
                ),
            },
            {
                "product": "mRNA-3745 (Moderna)",
                "sponsor": "Moderna",
                "phase": "I/II",
                "nct": "NCT05095727",
                "approach": (
                    "IV LNP-encapsulated mRNA encoding G6PC1.  Repeat "
                    "dosing required (mRNA is transient)."
                ),
                "status_2026": (
                    "Phase I/II ongoing.  Dose escalation in adults.  "
                    "Interim data showed dose-dependent glucose response.  "
                    "Different paradigm from AAV (repeat dosing vs one-time)."
                ),
            },
        ],
        "conditions": ["GSD_Ia", "glycogen_storage_disease", "metabolic",
                       "liver", "hypoglycemia"],
    },

    # -------------------------------------------------------------------
    # 2f.  Classical Galactosemia -- GALT
    # -------------------------------------------------------------------
    "GALT": {
        "gene_id": 2592,
        "chrom": "chr9",
        "start": 34_636_559,
        "end": 34_640_807,
        "strand": "-",
        "refseq": "NC_000009.12",
        "cytoband": "9p13.3",
        "exon_count": 11,
        "role": (
            "Galactose-1-phosphate uridylyltransferase -- enzyme in the "
            "Leloir pathway of galactose metabolism.  Catalyses conversion "
            "of galactose-1-phosphate + UDP-glucose to glucose-1-phosphate "
            "+ UDP-galactose.  Deficiency causes accumulation of "
            "galactose-1-phosphate and galactitol, leading to liver "
            "disease, cataracts, and long-term cognitive/neurological "
            "complications despite dietary galactose restriction."
        ),
        "disease": "Classical Galactosemia (Galactosemia Type I)",
        "omim_disease": 230400,
        "omim_gene": 606999,
        "inheritance": "AR",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "p.Gln188Arg (Q188R)",
                "hgvs_coding": "NM_000155.4:c.563A>G",
                "rsid": "rs75391579",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "frequency": "~70% of alleles in European populations",
                "notes": (
                    "Most common galactosemia mutation.  No residual enzyme "
                    "activity.  Not responsive to pharmacological chaperones."
                ),
            },
            {
                "name": "p.Lys285Asn (K285N / Duarte variant)",
                "rsid": "rs2070074",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Common in many populations.  Reduced but not absent "
                    "activity (~25% of normal when homozygous).  Duarte "
                    "galactosemia (D/G compound heterozygote) usually "
                    "benign."
                ),
            },
        ],
        "crispr_strategy": (
            "Liver-directed gene replacement or base editing.  GALT cDNA "
            "~1.1 kb (ideal AAV target).  Despite effective newborn screening "
            "and dietary restriction, long-term complications persist "
            "(speech/language delays, cognitive difficulties, premature "
            "ovarian insufficiency), indicating need for definitive therapy.  "
            "AAV-GALT gene replacement in GALT-null mice restored enzyme "
            "activity and normalized galactose-1-phosphate.  Base editing "
            "for Q188R (A>G, ideal for ABE) is feasible."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy clinical trials as of early 2026",
                "status_2026": (
                    "Preclinical.  Academic AAV-GALT programs at UT "
                    "Southwestern and University of Utah.  GALT-null mouse "
                    "model gene therapy data published.  No IND filed."
                ),
            },
        ],
        "conditions": ["galactosemia", "metabolic", "liver",
                       "carbohydrate_disorder"],
    },

    # -------------------------------------------------------------------
    # 2g.  Homocystinuria -- CBS
    # -------------------------------------------------------------------
    "CBS": {
        "gene_id": 875,
        "chrom": "chr21",
        "start": 43_053_189,
        "end": 43_076_943,
        "strand": "-",
        "refseq": "NC_000021.9",
        "cytoband": "21q22.3",
        "exon_count": 17,
        "role": (
            "Cystathionine beta-synthase -- pyridoxal 5'-phosphate (vitamin "
            "B6)-dependent enzyme that catalyses the first step of the "
            "transsulfuration pathway: condensation of homocysteine + serine "
            "to form cystathionine.  Deficiency causes classic homocystinuria: "
            "severe hyperhomocysteinemia with lens subluxation, skeletal "
            "abnormalities (marfanoid habitus), thromboembolism, and "
            "intellectual disability."
        ),
        "disease": "Classic Homocystinuria (CBS deficiency)",
        "omim_disease": 236200,
        "omim_gene": 613381,
        "inheritance": "AR",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "p.Ile278Thr (I278T)",
                "rsid": "rs5742905",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "frequency": "Most common in pan-ethnic; B6-responsive",
                "notes": (
                    "Pyridoxine (B6)-responsive variant.  Most common "
                    "mutation worldwide.  ~50% of patients are B6-responsive."
                ),
            },
            {
                "name": "p.Gly307Ser (G307S)",
                "rsid": "rs121964962",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "B6-nonresponsive.  Common in Celtic/Irish populations.  "
                    "Carrier frequency ~1/50 in Ireland."
                ),
            },
        ],
        "crispr_strategy": (
            "Liver-directed gene therapy.  CBS is expressed primarily in "
            "liver.  CBS cDNA ~1.7 kb (fits AAV).  Approaches: "
            "(1) AAV-CBS liver gene replacement; "
            "(2) mRNA therapy (Moderna preclinical program); "
            "(3) Base editing for I278T (T>C, ideal for ABE).  "
            "Preclinical AAV-CBS in Cbs-/- mice restored homocysteine "
            "levels to near-normal.  B6-responsive patients managed "
            "medically, but B6-nonresponsive patients need alternative."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy clinical trials as of early 2026",
                "status_2026": (
                    "Preclinical.  OT-601 (Orphan Technologies) enzyme "
                    "replacement therapy (PEGylated CBS) entered Phase I/II.  "
                    "Pegtibatinase Phase III.  Gene therapy in preclinical stage."
                ),
            },
        ],
        "conditions": ["homocystinuria", "metabolic", "amino_acid_disorder",
                       "liver"],
    },

    # -------------------------------------------------------------------
    # 2h.  OTC Deficiency -- OTC (Ornithine Transcarbamylase)
    # -------------------------------------------------------------------
    "OTC": {
        "gene_id": 5009,
        "chrom": "chrX",
        "start": 38_352_984,
        "end": 38_421_450,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xp11.4",
        "exon_count": 10,
        "role": (
            "Ornithine transcarbamylase -- mitochondrial matrix enzyme that "
            "catalyses the second step of the urea cycle: condensation of "
            "ornithine + carbamoyl phosphate to form citrulline.  Most "
            "common urea cycle defect.  X-linked: hemizygous males have "
            "severe neonatal hyperammonemia; heterozygous females have "
            "variable phenotype (mild to severe depending on X-inactivation)."
        ),
        "disease": "OTC Deficiency (Ornithine Transcarbamylase Deficiency)",
        "omim_disease": 311250,
        "omim_gene": 300461,
        "inheritance": "XL",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "p.Arg40His (R40H, late-onset)",
                "rsid": "rs121434385",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Partial enzyme deficiency; late-onset form.",
            },
            {
                "name": "p.Arg277Trp (R277W)",
                "rsid": "rs121434386",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Severe neonatal-onset form.",
            },
        ],
        "mutational_landscape": (
            ">400 pathogenic OTC variants reported.  Extreme allelic "
            "heterogeneity.  ~80% are point mutations (missense predominate).  "
            "~20% are small deletions/insertions.  De novo rate ~10% in "
            "severe male cases.  No clear genotype-phenotype correlation "
            "for many variants."
        ),
        "crispr_strategy": (
            "Liver-directed gene therapy -- one of the most clinically "
            "advanced metabolic gene therapy targets: "
            "(1) Ultragenyx/GeneTx DTX301 (AAV8-OTC) Phase I/II completed "
            "with mixed results -- dose-dependent OTC expression but "
            "variable durability. "
            "(2) LOGICA program: LNP-CRISPR for liver-directed OTC "
            "correction (preclinical). "
            "(3) OTC cDNA ~1.1 kb (ideal AAV target). "
            "(4) Base editing for recurrent missense mutations. "
            "HISTORICAL NOTE: OTC deficiency was the site of the Jesse "
            "Gelsinger tragedy (1999), the first gene therapy death, "
            "which used an adenoviral vector.  Modern AAV and LNP "
            "approaches have much better safety profiles."
        ),
        "clinical_trials": [
            {
                "product": "DTX301 (Ultragenyx)",
                "sponsor": "Ultragenyx Pharmaceutical",
                "phase": "I/II",
                "nct": "NCT02991144",
                "approach": (
                    "AAV8 vector delivering human OTC cDNA under liver-"
                    "specific promoter.  Single IV infusion."
                ),
                "status_2026": (
                    "Phase I/II completed.  Dose-dependent reduction in "
                    "ammonia levels and increase in ureagenesis.  However, "
                    "durability was variable.  Ultragenyx evaluating next-"
                    "generation vectors.  Phase III planning."
                ),
            },
        ],
        "conditions": ["OTC_deficiency", "urea_cycle_defect", "metabolic",
                       "liver", "hyperammonemia"],
    },

    # -------------------------------------------------------------------
    # 2i.  Citrullinemia Type I -- ASS1
    # -------------------------------------------------------------------
    "ASS1": {
        "gene_id": 445,
        "chrom": "chr9",
        "start": 130_444_880,
        "end": 130_501_274,
        "strand": "+",
        "refseq": "NC_000009.12",
        "cytoband": "9q34.11",
        "exon_count": 16,
        "role": (
            "Argininosuccinate synthase 1 -- cytosolic enzyme catalysing the "
            "third step of the urea cycle: condensation of citrulline + "
            "aspartate to form argininosuccinate (ATP-dependent).  Deficiency "
            "causes citrullinemia type I: hyperammonemia with elevated "
            "citrulline."
        ),
        "disease": "Citrullinemia Type I (CTLN1, Argininosuccinate Synthetase Deficiency)",
        "omim_disease": 215700,
        "omim_gene": 603470,
        "inheritance": "AR",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "p.Gly390Arg (G390R)",
                "rsid": "rs121908641",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Recurrent in East Asian populations.",
            },
            {
                "name": "IVS6-2A>G (c.535-2A>G)",
                "rsid": "rs121908642",
                "consequence": "splice_acceptor",
                "clinical_significance": "pathogenic",
                "notes": "Common in Japanese population (neonatal form).",
            },
        ],
        "crispr_strategy": (
            "Liver-directed gene replacement.  ASS1 cDNA ~1.2 kb (ideal "
            "AAV target).  Liver is the primary site of urea cycle activity.  "
            "Same rationale as OTC deficiency: liver transplant is curative, "
            "supporting liver-directed gene therapy.  Preclinical AAV-ASS1 "
            "in Ass1-deficient mice showed correction.  Base editing for "
            "recurrent variants is feasible."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy clinical trials as of early 2026",
                "status_2026": (
                    "Preclinical.  Academic programs exploring AAV-ASS1 "
                    "liver gene replacement.  Liver transplant remains "
                    "definitive treatment for severe cases."
                ),
            },
        ],
        "conditions": ["citrullinemia", "urea_cycle_defect", "metabolic",
                       "liver", "hyperammonemia"],
    },
}


# ============================================================================
# 3. HEMATOLOGIC RARE DISEASES
# ============================================================================

HEMATOLOGIC_RARE_TARGETS = {

    # -------------------------------------------------------------------
    # 3a.  Thalassemia Intermedia Modifier -- HBS1L-MYB intergenic region
    # -------------------------------------------------------------------
    "HBS1L_MYB_INTERGENIC": {
        "gene_id": None,  # intergenic regulatory region
        "chrom": "chr6",
        "start": 135_324_854,
        "end": 135_426_882,
        "strand": "N/A (intergenic)",
        "refseq": "NC_000006.12",
        "cytoband": "6q23.3",
        "exon_count": None,
        "role": (
            "HBS1L-MYB intergenic region -- contains multiple erythroid "
            "enhancers that regulate MYB expression in erythroid cells.  "
            "MYB is a key transcription factor controlling erythropoiesis "
            "and fetal hemoglobin (HbF) silencing.  SNPs in this region "
            "are the second strongest genetic modifier of HbF levels "
            "(after BCL11A).  Reduced MYB expression increases HbF "
            "production, ameliorating beta-thalassemia and SCD severity."
        ),
        "disease": "Beta-Thalassemia / SCD modifier (HbF regulation)",
        "omim_gene": None,
        "inheritance": "Complex trait modifier",
        "mutation_type": "regulatory_variants",
        "key_variants": [
            {
                "name": "rs9399137 (HMIP-2 core element)",
                "rsid": "rs9399137",
                "consequence": "regulatory",
                "clinical_significance": "modifier",
                "notes": (
                    "GWAS hit for HbF levels.  T allele associated with "
                    "higher HbF and milder beta-thalassemia phenotype.  "
                    "Located in HMIP-2 (HBS1L-MYB intergenic polymorphism "
                    "block 2), an erythroid enhancer ~84 kb upstream of MYB."
                ),
            },
            {
                "name": "rs66650371 (3-bp deletion in HMIP-2)",
                "rsid": "rs66650371",
                "consequence": "regulatory",
                "clinical_significance": "modifier",
                "notes": (
                    "3-bp deletion disrupting GATA1 binding motif in HMIP-2 "
                    "enhancer.  Reduces MYB expression, increases HbF.  "
                    "Strong HbF QTL."
                ),
            },
        ],
        "crispr_strategy": (
            "CRISPR disruption of the HMIP-2 enhancer to reduce MYB "
            "expression and elevate HbF.  Similar conceptual approach to "
            "BCL11A enhancer editing (Casgevy/exa-cel).  Ex vivo editing "
            "in CD34+ HSPCs.  Sankaran lab (Boston Children's/Broad) has "
            "demonstrated that CRISPR disruption of specific GATA1 binding "
            "sites within HMIP-2 phenocopies the naturally occurring "
            "HbF-boosting variants.  Could be combined with BCL11A "
            "enhancer editing for additive HbF induction."
        ),
        "clinical_trials": [
            {
                "product": "No clinical trials specifically targeting HBS1L-MYB region",
                "status_2026": (
                    "Preclinical.  Most clinical programs for HbF induction "
                    "target BCL11A (Casgevy) or HBG1/HBG2 promoters "
                    "(BEAM-101 / risto-cel).  HBS1L-MYB is an attractive "
                    "second target for combination editing strategies."
                ),
            },
        ],
        "conditions": ["beta_thalassemia", "sickle_cell_disease",
                       "hemoglobin_disorder", "HbF_modifier"],
    },

    # -------------------------------------------------------------------
    # 3b.  Congenital Neutropenia -- ELANE
    # -------------------------------------------------------------------
    "ELANE": {
        "gene_id": 1991,
        "chrom": "chr19",
        "start": 851_094,
        "end": 855_326,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19p13.3",
        "exon_count": 5,
        "role": (
            "Neutrophil elastase (elastase, neutrophil expressed) -- serine "
            "protease stored in azurophil granules of neutrophils.  Functions "
            "in innate immunity (microbial killing) and inflammation.  "
            "Mutations cause misfolded protein accumulation in the ER of "
            "neutrophil precursors, triggering unfolded protein response "
            "(UPR) and apoptosis.  Most common genetic cause of severe "
            "congenital neutropenia (SCN1, Kostmann syndrome)."
        ),
        "disease": "Severe Congenital Neutropenia Type 1 (SCN1 / Kostmann syndrome)",
        "omim_disease": 202700,
        "omim_gene": 130130,
        "inheritance": "AD",
        "mutation_type": "gain_of_toxic_function",
        "key_variants": [
            {
                "name": "Various (>100 mutations, extreme allelic heterogeneity)",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Mutations cluster in all 5 exons.  Most are missense.  "
                    "Some splice-site variants cause cyclic neutropenia "
                    "(milder phenotype) rather than SCN.  Genotype-phenotype "
                    "correlation imperfect."
                ),
            },
        ],
        "crispr_strategy": (
            "Ex vivo gene editing in autologous CD34+ HSPCs.  Approaches: "
            "(1) CRISPR-Cas9 knockout of mutant ELANE allele -- since "
            "disease is due to toxic gain-of-function (misfolded protein), "
            "eliminating the mutant allele should be therapeutic even "
            "without replacement (haploinsufficiency of ELANE is well "
            "tolerated). "
            "(2) Base editing to correct specific mutations. "
            "(3) Gene disruption combined with granulocyte colony-"
            "stimulating factor (G-CSF) if needed. "
            "IMPORTANT: ELANE-SCN patients are at high risk (~20% lifetime) "
            "for MDS/AML due to acquired CSF3R mutations.  Definitive "
            "cure via gene therapy may reduce this risk."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR clinical trials as of early 2026",
                "status_2026": (
                    "Preclinical.  Academic programs at NIH and Dana-Farber "
                    "exploring CRISPR knockout of mutant ELANE in HSPCs.  "
                    "G-CSF (filgrastim) remains standard of care.  HSCT "
                    "curative but carries transplant-related mortality."
                ),
            },
        ],
        "conditions": ["congenital_neutropenia", "SCN", "hematological",
                       "immunodeficiency"],
    },

    # -------------------------------------------------------------------
    # 3c.  Hereditary Spherocytosis -- ANK1
    # -------------------------------------------------------------------
    "ANK1": {
        "gene_id": 286,
        "chrom": "chr8",
        "start": 41_510_082,
        "end": 41_705_903,
        "strand": "-",
        "refseq": "NC_000008.11",
        "cytoband": "8p11.21",
        "exon_count": 43,
        "role": (
            "Ankyrin 1 (erythrocytic) -- adaptor protein that links the "
            "spectrin-actin cytoskeleton to band 3 (SLC4A1) and other "
            "integral membrane proteins in the red cell membrane.  "
            "Mutations cause hereditary spherocytosis (HS), the most "
            "common inherited hemolytic anemia in Northern European "
            "populations.  ANK1 mutations account for ~40-65% of HS."
        ),
        "disease": "Hereditary Spherocytosis (HS)",
        "omim_disease": 182900,
        "omim_gene": 612641,
        "inheritance": "AD (most common) or AR (rare, severe)",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "Extreme allelic heterogeneity (>60 mutations reported)",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Frameshift, nonsense, splice, and missense mutations "
                    "throughout the gene.  Most AD cases are de novo or "
                    "inherited haploinsufficiency.  AR cases with homozygous "
                    "or compound heterozygous mutations cause near-lethal "
                    "hydrops fetalis or transfusion-dependent anemia."
                ),
            },
        ],
        "crispr_strategy": (
            "Gene therapy for HS is not actively pursued because "
            "splenectomy effectively manages most cases.  However, for "
            "severe AR forms: "
            "(1) Ex vivo gene replacement in erythroid progenitors; "
            "(2) ANK1 cDNA is ~5.4 kb (at upper AAV limit, dual AAV may "
            "be needed). "
            "Not a priority for CRISPR given effective surgical management."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": (
                    "Splenectomy and supportive care remain standard.  "
                    "No gene therapy programs in development."
                ),
            },
        ],
        "conditions": ["hereditary_spherocytosis", "hemolytic_anemia",
                       "hematological"],
    },

    # -------------------------------------------------------------------
    # 3d.  Hereditary Spherocytosis -- SLC4A1 (Band 3)
    # -------------------------------------------------------------------
    "SLC4A1": {
        "gene_id": 6521,
        "chrom": "chr17",
        "start": 44_244_852,
        "end": 44_264_618,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.31",
        "exon_count": 20,
        "role": (
            "Solute carrier family 4 member 1 (anion exchanger 1, band 3) -- "
            "major integral membrane protein of the red blood cell membrane.  "
            "Functions as a chloride/bicarbonate exchanger (CO2 transport) "
            "and as an anchor for the spectrin-based cytoskeleton via "
            "ankyrin binding.  Mutations cause hereditary spherocytosis "
            "(~20-25% of HS), Southeast Asian ovalocytosis (SAO), and "
            "distal renal tubular acidosis (dRTA)."
        ),
        "disease": "Hereditary Spherocytosis (HS) / Distal RTA",
        "omim_disease": 612653,
        "omim_gene": 109270,
        "inheritance": "AD",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "Band 3 variants (various)",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Membrane-spanning domain mutations: HS.  Cytoplasmic "
                    "domain mutations: dRTA.  SAO: 27-bp deletion in exon "
                    "11 (common in Southeast Asia, protective against "
                    "P. falciparum malaria)."
                ),
            },
        ],
        "crispr_strategy": (
            "Not a priority target for CRISPR gene therapy.  Same rationale "
            "as ANK1: splenectomy manages HS effectively.  SLC4A1 cDNA "
            "~2.7 kb (fits AAV if needed for severe AR cases or dRTA)."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials as of early 2026",
                "status_2026": "No gene therapy programs in development.",
            },
        ],
        "conditions": ["hereditary_spherocytosis", "hemolytic_anemia",
                       "hematological", "renal_tubular_acidosis"],
    },
}


# ============================================================================
# 4. ENDOCRINE RARE DISEASES
# ============================================================================

ENDOCRINE_RARE_TARGETS = {

    # -------------------------------------------------------------------
    # 4a.  Congenital Hypothyroidism -- TSHR
    # -------------------------------------------------------------------
    "TSHR": {
        "gene_id": 7253,
        "chrom": "chr14",
        "start": 81_421_152,
        "end": 81_612_908,
        "strand": "+",
        "refseq": "NC_000014.9",
        "cytoband": "14q31.1",
        "exon_count": 10,
        "role": (
            "Thyroid stimulating hormone receptor -- G protein-coupled "
            "receptor (GPCR) that binds TSH and activates Gs/cAMP and "
            "Gq/PLC signalling in thyroid follicular cells.  Drives "
            "thyroid hormone synthesis and secretion.  Loss-of-function "
            "mutations cause congenital hypothyroidism with thyroid "
            "hypoplasia.  Gain-of-function causes familial non-autoimmune "
            "hyperthyroidism."
        ),
        "disease": "Congenital Hypothyroidism (TSH resistance)",
        "omim_disease": 275200,
        "omim_gene": 603372,
        "inheritance": "AR (LOF) or AD (GOF)",
        "mutation_type": "loss_of_function_or_gain_of_function",
        "key_variants": [
            {
                "name": "Various loss-of-function mutations",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "LOF: inactivating mutations in extracellular or "
                    "transmembrane domains.  Homozygous/compound het causes "
                    "TSH resistance.  Heterozygous carriers mildly elevated "
                    "TSH.  GOF: constitutively activating mutations cause "
                    "non-autoimmune hyperthyroidism."
                ),
            },
        ],
        "crispr_strategy": (
            "Gene therapy is NOT a priority for TSHR-related congenital "
            "hypothyroidism because levothyroxine (T4) replacement is "
            "safe, effective, and inexpensive.  However, for thyroid "
            "dysgenesis cases requiring lifelong replacement: "
            "(1) Thyroid organoid/cell therapy with gene-corrected cells; "
            "(2) TSHR cDNA ~2.3 kb (fits AAV) for potential thyroid-"
            "directed gene replacement in organoid context."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials",
                "status_2026": (
                    "Levothyroxine replacement is standard of care.  No "
                    "gene therapy programs.  Not a viable CRISPR target "
                    "given effective pharmacological management."
                ),
            },
        ],
        "conditions": ["congenital_hypothyroidism", "endocrine",
                       "thyroid"],
    },

    # -------------------------------------------------------------------
    # 4b.  Congenital Hypothyroidism -- PAX8
    # -------------------------------------------------------------------
    "PAX8": {
        "gene_id": 7849,
        "chrom": "chr2",
        "start": 113_215_996,
        "end": 113_278_954,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2q14.1",
        "exon_count": 12,
        "role": (
            "Paired box 8 -- transcription factor essential for thyroid "
            "gland development and differentiation.  Regulates expression "
            "of thyroglobulin (TG), thyroid peroxidase (TPO), and sodium-"
            "iodide symporter (NIS/SLC5A5).  Haploinsufficiency causes "
            "thyroid dysgenesis (hypoplasia/ectopia) with congenital "
            "hypothyroidism."
        ),
        "disease": "Congenital Hypothyroidism (thyroid dysgenesis)",
        "omim_disease": 218700,
        "omim_gene": 167415,
        "inheritance": "AD (haploinsufficiency)",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": "Various (paired domain mutations most common)",
                "consequence": "mixed",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~2-5% of congenital hypothyroidism cases.  Mutations "
                    "cluster in the paired DNA-binding domain.  Variable "
                    "expressivity even within families."
                ),
            },
        ],
        "crispr_strategy": (
            "Not a practical CRISPR target.  PAX8 is a developmental "
            "transcription factor -- gene therapy would need to be "
            "delivered during thyroid organogenesis (prenatal).  "
            "Levothyroxine replacement is the definitive management.  "
            "Research interest: iPSC-derived thyroid organoids with "
            "gene-corrected PAX8 for potential regenerative medicine."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials",
                "status_2026": (
                    "Not applicable.  Levothyroxine replacement is safe "
                    "and effective.  Research focus on thyroid organoid "
                    "generation from PAX8-corrected iPSCs."
                ),
            },
        ],
        "conditions": ["congenital_hypothyroidism", "endocrine",
                       "thyroid", "developmental"],
    },

    # -------------------------------------------------------------------
    # 4c.  Multiple Endocrine Neoplasia Type 1 -- MEN1
    # -------------------------------------------------------------------
    "MEN1": {
        "gene_id": 4221,
        "chrom": "chr11",
        "start": 64_572_557,
        "end": 64_578_765,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11q13.1",
        "exon_count": 10,
        "role": (
            "Menin 1 -- nuclear scaffold protein and tumor suppressor.  "
            "Interacts with >40 proteins including JunD, NF-kB, Smads, "
            "FANCD2, and the MLL/SET1 histone methyltransferase complex.  "
            "Functions in transcription regulation, DNA repair, chromatin "
            "remodeling, and cell cycle control.  Loss-of-function causes "
            "MEN1 syndrome: tumors of parathyroid (~95%), enteropancreatic "
            "neuroendocrine (~40-70%), and anterior pituitary (~30-40%)."
        ),
        "disease": "Multiple Endocrine Neoplasia Type 1 (MEN1)",
        "omim_disease": 131100,
        "omim_gene": 613733,
        "inheritance": "AD (tumor suppressor, two-hit model)",
        "mutation_type": "loss_of_function",
        "key_variants": [
            {
                "name": ">1,500 unique MEN1 mutations reported",
                "consequence": "mixed (>70% truncating)",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Extreme allelic heterogeneity.  No genotype-phenotype "
                    "correlation.  >70% are truncating (frameshift, "
                    "nonsense, splice).  Mutations distributed throughout "
                    "coding region with no hotspots.  ~10% de novo."
                ),
            },
        ],
        "crispr_strategy": (
            "Tumor suppressor gene -- CRISPR approaches limited: "
            "(1) Correcting inherited mutations would not address "
            "second-hit somatic inactivation in tumors; "
            "(2) Potential for prophylactic gene replacement in at-risk "
            "endocrine tissues -- extremely challenging given multisystem "
            "involvement; "
            "(3) MEN1 cDNA ~1.8 kb (fits AAV) but tissue targeting to "
            "parathyroid, islets, and pituitary simultaneously is not "
            "feasible with current technology.  "
            "Standard of care: surveillance and surgery for tumors."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials",
                "status_2026": (
                    "No gene therapy programs.  Surveillance-and-surgery "
                    "paradigm.  Not a viable current CRISPR target due to "
                    "tumor suppressor biology and multisystem involvement."
                ),
            },
        ],
        "conditions": ["MEN1", "multiple_endocrine_neoplasia", "endocrine",
                       "tumor_predisposition"],
    },

    # -------------------------------------------------------------------
    # 4d.  Multiple Endocrine Neoplasia Type 2 -- RET
    # -------------------------------------------------------------------
    "RET": {
        "gene_id": 5979,
        "chrom": "chr10",
        "start": 43_077_027,
        "end": 43_130_351,
        "strand": "+",
        "refseq": "NC_000010.11",
        "cytoband": "10q11.21",
        "exon_count": 20,
        "role": (
            "RET proto-oncogene -- receptor tyrosine kinase that binds GDNF "
            "family ligands (GDNF, neurturin, artemin, persephin) via "
            "GFRalpha co-receptors.  Activates RAS/MAPK, PI3K/AKT, and "
            "PLCgamma pathways.  Gain-of-function mutations cause MEN2: "
            "MEN2A (medullary thyroid carcinoma + pheochromocytoma + "
            "hyperparathyroidism) and MEN2B (MTC + pheo + mucosal neuromas).  "
            "Loss-of-function causes Hirschsprung disease."
        ),
        "disease": "Multiple Endocrine Neoplasia Type 2 (MEN2A, MEN2B)",
        "omim_disease": 171400,
        "omim_gene": 164761,
        "inheritance": "AD (gain-of-function)",
        "mutation_type": "gain_of_function",
        "key_variants": [
            {
                "name": "p.Cys634Arg (C634R)",
                "hgvs_coding": "NM_020975.6:c.1900T>C",
                "rsid": "rs77724903",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common MEN2A mutation (~80% of MEN2A families).  "
                    "Extracellular cysteine-rich domain.  Causes ligand-"
                    "independent dimerization and constitutive kinase "
                    "activation."
                ),
            },
            {
                "name": "p.Met918Thr (M918T)",
                "hgvs_coding": "NM_020975.6:c.2753T>C",
                "rsid": "rs74799832",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "MEN2B (~95% of cases).  Intracellular tyrosine kinase "
                    "domain.  Most aggressive phenotype -- prophylactic "
                    "thyroidectomy recommended in first 6 months of life.  "
                    "De novo in ~50%."
                ),
            },
        ],
        "crispr_strategy": (
            "Allele-specific disruption of gain-of-function RET mutant "
            "allele.  Approaches: "
            "(1) CRISPR-Cas9 allele-specific knockout using SNP-guided "
            "targeting to distinguish mutant from wild-type allele; "
            "(2) Base editing to revert specific GOF mutations (e.g., "
            "M918T: T>C, could use ABE); "
            "(3) For MEN2B: given severity, prophylactic thyroidectomy "
            "in infancy is current standard -- gene editing would need "
            "to be delivered in utero or neonatally.  "
            "Selpercatinib and pralsetinib (RET kinase inhibitors) are "
            "FDA-approved for RET-altered cancers, providing pharmacological "
            "alternative.  CRISPR approaches are preclinical only."
        ),
        "clinical_trials": [
            {
                "product": "No CRISPR/gene therapy trials",
                "status_2026": (
                    "RET kinase inhibitors (selpercatinib, pralsetinib) "
                    "FDA-approved for advanced MTC.  Prophylactic "
                    "thyroidectomy for MEN2.  No gene therapy programs "
                    "for germline RET correction.  CRISPR approaches "
                    "remain theoretical."
                ),
            },
        ],
        "conditions": ["MEN2", "multiple_endocrine_neoplasia", "endocrine",
                       "tumor_predisposition", "medullary_thyroid_carcinoma"],
    },
}


# ============================================================================
# UNIFIED ACCESS: all additional rare disease targets in a single dictionary
# ============================================================================

ALL_ADDITIONAL_RARE_TARGETS = {
    **NEUROLOGICAL_RARE_TARGETS,
    **METABOLIC_RARE_TARGETS,
    **HEMATOLOGIC_RARE_TARGETS,
    **ENDOCRINE_RARE_TARGETS,
}


# ============================================================================
# Helper functions
# ============================================================================

def get_target(gene_symbol: str) -> dict | None:
    """Return the target dict for *gene_symbol*, or None."""
    return ALL_ADDITIONAL_RARE_TARGETS.get(gene_symbol.upper())


def get_targets_by_condition(condition: str) -> dict[str, dict]:
    """Return all targets matching a given condition string."""
    return {
        gene: info
        for gene, info in ALL_ADDITIONAL_RARE_TARGETS.items()
        if condition in info.get("conditions", [])
    }


def get_targets_by_category(category: str) -> dict:
    """Return the appropriate sub-dictionary by disease category."""
    categories = {
        "neurological": NEUROLOGICAL_RARE_TARGETS,
        "metabolic": METABOLIC_RARE_TARGETS,
        "hematologic": HEMATOLOGIC_RARE_TARGETS,
        "endocrine": ENDOCRINE_RARE_TARGETS,
    }
    return categories.get(category, {})


def get_bed_regions() -> list[tuple[str, int, int, str]]:
    """Return a BED-format list of (chrom, start, end, name) for all targets."""
    regions = []
    for name, info in ALL_ADDITIONAL_RARE_TARGETS.items():
        if info.get("start") and info.get("end"):
            regions.append((info["chrom"], info["start"], info["end"], name))
    return sorted(regions, key=lambda r: (r[0], r[1]))


def get_active_clinical_programs() -> dict[str, dict]:
    """Return only genes with active (non-preclinical) clinical programs."""
    active = {}
    for gene, info in ALL_ADDITIONAL_RARE_TARGETS.items():
        trials = info.get("clinical_trials", [])
        for trial in trials:
            phase = trial.get("phase", "")
            if phase and "No" not in trial.get("product", ""):
                active[gene] = info
                break
    return active


def get_crispr_amenable_targets() -> dict[str, dict]:
    """Return targets with active CRISPR or gene therapy clinical programs."""
    return {
        gene: info
        for gene, info in ALL_ADDITIONAL_RARE_TARGETS.items()
        if any(
            "Phase" in t.get("phase", "") and "No" not in t.get("product", "")
            for t in info.get("clinical_trials", [])
        )
    }


# ============================================================================
# Summary table for quick reference
# ============================================================================

GENE_COORDINATE_SUMMARY = {
    gene: {
        "chrom": info["chrom"],
        "start": info["start"],
        "end": info["end"],
        "strand": info.get("strand", "N/A"),
        "cytoband": info.get("cytoband", "N/A"),
        "disease": info.get("disease", "N/A"),
    }
    for gene, info in ALL_ADDITIONAL_RARE_TARGETS.items()
    if info.get("start") is not None
}
