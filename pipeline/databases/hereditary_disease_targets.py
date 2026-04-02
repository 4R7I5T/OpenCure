"""
Hereditary disease gene-therapy and CRISPR targets -- GRCh38 coordinates.

Comprehensive database covering neurodegenerative diseases, muscular
dystrophies, hereditary blindness/retinal dystrophies, skin diseases,
and systemic hereditary conditions.  Each entry includes verified NCBI
GRCh38.p14 coordinates, key pathogenic variants, CRISPR/gene-therapy
strategies, and clinical-trial status as of early 2026.

Categories:
  1. Neurodegenerative diseases (HD, SMA, ALS-SOD1, Friedreich ataxia)
  2. Muscular dystrophies (DMD/BMD, myotonic dystrophy type 1)
  3. Hereditary blindness (LCA-RPE65, LCA10-CEP290, RP-RPGR/RHO,
     Usher syndrome MYO7A/USH2A)
  4. Skin diseases (epidermolysis bullosa COL7A1/KRT5/KRT14/LAMB3)
  5. Systemic hereditary (transthyretin amyloidosis - TTR)

All coordinates are GRCh38/hg38 (NCBI Annotation Release RS_2025_08).

Sources:
  - NCBI Gene (GRCh38.p14)
  - OMIM
  - ClinVar
  - ClinicalTrials.gov
  - Published literature through early 2026
"""


# ==========================================================================
# 1. NEURODEGENERATIVE DISEASES
# ==========================================================================

NEURODEGENERATIVE_TARGETS = {

    # ------------------------------------------------------------------
    # 1a. Huntington's Disease
    # ------------------------------------------------------------------
    "HTT": {
        "chrom": "chr4",
        "start": 3074681,
        "end": 3243960,
        "strand": "+",
        "cytoband": "4p16.3",
        "ncbi_gene_id": 3064,
        "refseq": "NC_000004.12",
        "role": "Huntingtin protein -- cytoplasmic/nuclear scaffolding, vesicle "
                "transport, transcription regulation, BDNF transport",
        "disease": "Huntington's disease (HD)",
        "omim_disease": 143100,
        "omim_gene": 613004,
        "inheritance": "AD",
        "mutation_type": "trinucleotide_repeat_expansion",
        "repeat_details": {
            "repeat_unit": "CAG",
            "location": "exon 1 (5' coding region)",
            "normal_range": "6-26 repeats",
            "intermediate": "27-35 repeats (no disease, may expand in offspring)",
            "reduced_penetrance": "36-39 repeats",
            "full_penetrance": ">=40 repeats (disease certain)",
            "juvenile_onset": ">=60 repeats (onset <20 years)",
            "protein_effect": "Expanded polyglutamine tract in N-terminal huntingtin "
                             "causes toxic gain-of-function: aggregation, "
                             "transcriptional dysregulation, mitochondrial "
                             "dysfunction, excitotoxicity",
        },
        "key_variants": [
            "CAG_expansion_exon1 (>=36 repeats pathogenic)",
        ],
        "strategy": (
            "Multiple CRISPR approaches: (1) Allele-specific silencing -- use "
            "SNP-guided Cas9 or CRISPRi to selectively silence the mutant HTT "
            "allele while preserving wild-type; (2) CAG repeat contraction -- "
            "paired Cas9 nickases or dual-guide excision flanking the expanded "
            "CAG tract to contract repeats below pathogenic threshold; "
            "(3) RNA-targeting -- Cas13d/CasRx for mRNA knockdown of mutant "
            "HTT with allele-specificity via SNP targeting; (4) Base editing -- "
            "adenine/cytosine base editors to reduce somatic repeat expansion "
            "(Nature Genetics 2025); (5) Total HTT knockdown -- CRISPRi or "
            "miRNA if partial HTT loss is tolerated"
        ),
        "clinical_programs": {
            "AMT-130_uniQure": {
                "type": "AAV5-miRNA (HTT silencing)",
                "route": "intrastriatal",
                "status": "Phase I/II -- 36-month data showed 75% slowing of "
                          "disease progression (cUHDRS, p=0.003). Pre-BLA "
                          "meeting with FDA held Oct 2025; Type A meeting "
                          "Jan 2026. Potential first disease-modifying therapy.",
            },
            "LETI-101_LifeEdit": {
                "type": "CRISPR allele-specific editing (AAV-delivered)",
                "route": "intrastriatal",
                "status": "Preclinical -- up to 80% mutant HTT reduction in "
                          "striatum, with cortex/thalamus effects exceeding "
                          "40% threshold. Human trials planned 2026.",
            },
            "WVE-003_Wave": {
                "type": "Allele-selective ASO (targets SNP rs362307)",
                "route": "intrathecal",
                "status": "Phase I/II -- selectively silences mutant HTT "
                          "mRNA while sparing wild-type allele.",
            },
        },
        "conditions": ["huntington_disease", "neurodegenerative"],
    },

    # ------------------------------------------------------------------
    # 1b. Spinal Muscular Atrophy -- SMN1
    # ------------------------------------------------------------------
    "SMN1": {
        "chrom": "chr5",
        "start": 70924941,
        "end": 70966375,
        "strand": "+",
        "cytoband": "5q13.2",
        "ncbi_gene_id": 6606,
        "refseq": "NC_000005.10",
        "role": "Survival of motor neuron 1 -- snRNP assembly, mRNA splicing, "
                "motor neuron survival",
        "disease": "Spinal muscular atrophy (SMA types I-IV)",
        "omim_disease": 253300,
        "omim_gene": 600354,
        "inheritance": "AR",
        "mutation_type": "deletion/conversion",
        "mutation_details": (
            "Homozygous deletion or gene conversion of SMN1 exon 7 causes ~95% "
            "of SMA. Exon 7 deletion eliminates full-length SMN protein "
            "production from SMN1. SMN2 partially compensates but produces "
            "mostly truncated protein (Delta7-SMN) due to C-to-T change at "
            "position c.840 in exon 7 disrupting an exonic splicing enhancer."
        ),
        "key_variants": [
            "exon_7_deletion (homozygous -- 95% of SMA)",
            "exon_7_8_deletion (homozygous)",
            "intragenic_point_mutations (compound het with deletion, ~5%)",
        ],
        "strategy": (
            "Gene replacement: Zolgensma (onasemnogene abeparvovec) delivers "
            "functional SMN1 via AAV9 (see clinical_programs). "
            "CRISPR approaches: (1) SMN2-to-SMN1 conversion via adenine base "
            "editing of c.840T>C in exon 7, restoring exonic splicing enhancer "
            "and full-length SMN production (up to 99% reversion in "
            "fibroblasts, Science 2023); (2) Prime editing of c.840 T-to-C "
            "plus c.859 G-to-C to fully convert SMN2 splicing; "
            "(3) CRISPRa to upregulate SMN2 expression from existing copies"
        ),
        "clinical_programs": {
            "Zolgensma_Novartis": {
                "type": "AAV9-SMN1 gene replacement",
                "route": "IV (single dose)",
                "status": "FDA-approved 2019 for SMA patients <2 years. "
                          "First gene therapy for a neuromuscular disease. "
                          "One-time IV infusion delivers functional SMN1 "
                          "transgene. ~$2.1M per dose.",
            },
            "nusinersen_Spinraza": {
                "type": "ASO splice-switching (promotes SMN2 exon 7 inclusion)",
                "route": "intrathecal",
                "status": "FDA-approved 2016. Chronic dosing required.",
            },
            "risdiplam_Evrysdi": {
                "type": "Small molecule splice modifier (oral)",
                "route": "oral",
                "status": "FDA-approved 2020. Daily oral dosing.",
            },
        },
        "conditions": ["spinal_muscular_atrophy", "neurodegenerative",
                        "motor_neuron_disease"],
    },

    # ------------------------------------------------------------------
    # 1b2. Spinal Muscular Atrophy -- SMN2 (therapeutic target)
    # ------------------------------------------------------------------
    "SMN2": {
        "chrom": "chr5",
        "start": 70049523,
        "end": 70090528,
        "strand": "+",
        "cytoband": "5q13.2",
        "ncbi_gene_id": 6607,
        "refseq": "NC_000005.10",
        "role": "Survival of motor neuron 2 -- paralog of SMN1, produces ~10-20% "
                "full-length SMN protein due to alternative splicing of exon 7",
        "disease": "SMA modifier (copy number determines severity)",
        "omim_gene": 601627,
        "inheritance": "modifier",
        "mutation_type": "critical_nucleotide_difference",
        "mutation_details": (
            "Single C-to-T transition at c.840 (exon 7, position 6) converts "
            "an exonic splicing enhancer (ESE) to an exonic splicing silencer "
            "(ESS). This causes ~80-90% of SMN2 transcripts to skip exon 7, "
            "producing truncated, unstable Delta7-SMN protein. SMN2 copy number "
            "inversely correlates with severity: 1 copy = type I (severe), "
            "2 copies = type II, 3+ copies = type III/IV (milder)."
        ),
        "key_variants": [
            "c.840C>T (exon 7 ESE disruption -- THE critical difference)",
            "c.859G>C (exon 7 -- additional splicing modifier)",
        ],
        "crispr_target_region": {
            "description": "Exon 7 region for base editing SMN2->SMN1 conversion",
            "exon7_start_approx": 70076874,
            "exon7_end_approx": 70076928,
            "critical_nucleotide_c840": "Position within exon 7 for A-to-G "
                                        "base edit (on antisense strand) to "
                                        "revert T>C and restore ESE",
        },
        "strategy": (
            "Primary CRISPR target for SMA: adenine base editing (ABE) to "
            "convert c.840T back to C in SMN2 exon 7, restoring full-length "
            "SMN protein production from all SMN2 copies. Nature Biomedical "
            "Engineering 2023 optimized ABE8e delivery achieved near-complete "
            "exon 7 inclusion. Also: CRISPR/Cpf1 + ssODN for gene conversion "
            "(4/36 efficiency). Prime editing for dual c.840+c.859 correction."
        ),
        "conditions": ["spinal_muscular_atrophy", "therapeutic_target"],
    },

    # ------------------------------------------------------------------
    # 1c. ALS (SOD1-linked)
    # ------------------------------------------------------------------
    "SOD1": {
        "chrom": "chr21",
        "start": 31659693,
        "end": 31668931,
        "strand": "+",
        "cytoband": "21q22.11",
        "ncbi_gene_id": 6647,
        "refseq": "NC_000021.9",
        "role": "Superoxide dismutase 1 -- cytoplasmic Cu/Zn SOD, free radical "
                "detoxification. Toxic gain-of-function mutations cause motor "
                "neuron death via protein aggregation, oxidative stress, "
                "mitochondrial dysfunction",
        "disease": "Amyotrophic lateral sclerosis (ALS, SOD1-linked)",
        "omim_disease": 105400,
        "omim_gene": 147450,
        "inheritance": "AD (most variants); AR (D90A homozygous in Scandinavian)",
        "mutation_type": "missense (>180 known pathogenic variants)",
        "key_variants": [
            "A4V (A5V per mature protein; rs121912442) -- most common in "
            "North America, ~50% of familial SOD1-ALS, aggressive course "
            "(median survival ~1.2 years from onset)",
            "D90A (D91A per mature protein; rs80265967) -- most common "
            "worldwide; AR in Scandinavian families (mild, slow progression), "
            "AD elsewhere (more severe)",
            "G93A (rs121912438) -- used in SOD1-G93A transgenic mouse model",
            "H46R (rs121912440) -- common in Japan",
            "E100G, I113T, L144F -- other recurrent pathogenic variants",
        ],
        "strategy": (
            "Goal: reduce toxic mutant SOD1 protein. (1) CRISPR-Cas9 knockout "
            "of mutant SOD1 allele -- AAV-delivered Cas9+sgRNA targeting SOD1 "
            "exons reduced protein >2.5-fold in spinal cord of G93A mice, "
            "extending survival 25% (Science Advances 2018); (2) Allele-specific "
            "CRISPR using PAM-proximal SNPs to distinguish mutant from WT; "
            "(3) CRISPRi/dCas9-KRAB for transcriptional silencing; "
            "(4) Cas13-based mRNA knockdown for reversible suppression; "
            "(5) Total SOD1 knockdown may be tolerated as SOD1-null mice "
            "are viable (though with some axonopathy)"
        ),
        "clinical_programs": {
            "tofersen_Qalsody_Biogen": {
                "type": "ASO (SOD1 mRNA degradation via RNase H)",
                "route": "intrathecal",
                "status": "FDA accelerated approval April 2023. First treatment "
                          "targeting genetic cause of ALS. Long-term data (2025) "
                          "show slowed progression; ~25% of patients show "
                          "stabilization or improvement. Phase 3 ATLAS trial "
                          "in presymptomatic SOD1-ALS expected results 2027.",
            },
        },
        "conditions": ["ALS", "neurodegenerative", "motor_neuron_disease"],
    },

    # ------------------------------------------------------------------
    # 1d. Friedreich's Ataxia
    # ------------------------------------------------------------------
    "FXN": {
        "chrom": "chr9",
        "start": 69035752,
        "end": 69079076,
        "strand": "-",
        "cytoband": "9q21.11",
        "ncbi_gene_id": 2395,
        "refseq": "NC_000009.12",
        "role": "Frataxin -- mitochondrial iron-sulfur cluster biogenesis, "
                "iron homeostasis, oxidative stress protection",
        "disease": "Friedreich's ataxia (FRDA)",
        "omim_disease": 229300,
        "omim_gene": 606829,
        "inheritance": "AR",
        "mutation_type": "trinucleotide_repeat_expansion",
        "repeat_details": {
            "repeat_unit": "GAA",
            "location": "intron 1 of FXN",
            "normal_range": "5-33 repeats (85% of alleles have <12)",
            "premutation": "34-65 repeats",
            "pathogenic": "66-1700 repeats (most patients 600-900)",
            "mechanism": "Expanded GAA repeat forms unusual DNA structures "
                        "(triplex, sticky DNA, R-loops) that silence FXN "
                        "transcription via heterochromatin spreading. "
                        "Frataxin protein reduced to 5-30% of normal levels. "
                        "Repeat length inversely correlates with age of onset.",
        },
        "key_variants": [
            "GAA expansion in intron 1 (homozygous in ~96% of patients)",
            "GAA expansion + point mutation (compound het, ~4%)",
            "Point mutations: I154F, W155R, G130V (rare compound hets)",
        ],
        "strategy": (
            "Multiple approaches: (1) GAA repeat excision -- paired CRISPR-Cas9 "
            "guides flanking the expanded GAA tract to delete the expansion "
            "and restore FXN transcription (demonstrated in YG8R mice and "
            "patient iPSC-derived cells); (2) Base editing to reduce repeat "
            "expansions -- adenine/cytosine base editors reduce somatic "
            "expansion in CNS of FA mouse models (Nature Genetics 2025); "
            "(3) CRISPRa/epigenetic reactivation -- dCas9-VPR or dCas9-p300 "
            "to counteract heterochromatin silencing at FXN locus; "
            "(4) AAV-frataxin gene replacement (see clinical programs); "
            "(5) Anti-GAA oligonucleotides to disrupt triplex structures "
            "and increase FXN expression (Mol Ther Nucleic Acids 2025)"
        ),
        "clinical_programs": {
            "SGT-212_SolidBio": {
                "type": "AAV-frataxin gene replacement (dual-route delivery)",
                "route": "intradentate nucleus (MRI-guided) + IV",
                "status": "Phase Ib FALCON trial -- first patient dosed "
                          "Jan 2026. FDA Fast Track designation Jan 2025. "
                          "Designed to deliver frataxin to both cerebellar "
                          "dentate nuclei and cardiomyocytes. First-in-class "
                          "dual-route gene therapy for FA.",
            },
            "Voyager_Neurocrine": {
                "type": "AAV-FXN with BBB-penetrant novel capsid",
                "route": "IV",
                "status": "Development candidate selected 2024; FIH trials "
                          "expected 2025-2026.",
            },
        },
        "conditions": ["friedreich_ataxia", "neurodegenerative",
                        "trinucleotide_repeat"],
    },
}


# ==========================================================================
# 2. MUSCULAR DYSTROPHIES
# ==========================================================================

MUSCULAR_DYSTROPHY_TARGETS = {

    # ------------------------------------------------------------------
    # 2a. Duchenne/Becker Muscular Dystrophy
    # ------------------------------------------------------------------
    "DMD": {
        "chrom": "chrX",
        "start": 31119222,
        "end": 33339388,
        "strand": "-",
        "cytoband": "Xp21.2-p21.1",
        "ncbi_gene_id": 1756,
        "refseq": "NC_000023.11",
        "gene_size_bp": 2220166,
        "role": "Dystrophin -- cytoskeletal linker protein connecting actin "
                "cytoskeleton to extracellular matrix via dystroglycan complex. "
                "Largest known human gene (~2.2 Mb, 79 exons). Protein: 427 kDa "
                "(3685 aa), with N-terminal actin-binding domain, central rod "
                "domain (24 spectrin-like repeats), and C-terminal dystroglycan "
                "binding domain",
        "disease": "Duchenne muscular dystrophy (DMD) / Becker muscular "
                   "dystrophy (BMD)",
        "omim_disease_dmd": 310200,
        "omim_disease_bmd": 300376,
        "omim_gene": 300377,
        "inheritance": "XLR",
        "mutation_type": "large deletions (60-70%), duplications (5-15%), "
                         "point mutations (25-30%)",
        "mutation_details": (
            "DMD: out-of-frame mutations abolish dystrophin production. "
            "BMD: in-frame mutations produce truncated but partially "
            "functional dystrophin. Mutation hotspots: exons 45-55 (rod "
            "domain) and exons 2-10 (actin-binding domain). The central "
            "rod domain is partially redundant, so in-frame deletions "
            "produce the milder Becker phenotype."
        ),
        "key_variants": [
            "exon_45-55_deletions (most common hotspot region)",
            "exon_2-10_deletions (second hotspot)",
            "exon_51_amenable_deletions (~13% of DMD patients)",
            "exon_53_amenable_deletions (~8% of DMD patients)",
            "exon_45_amenable_deletions (~8% of DMD patients)",
        ],
        "exon_skipping_targets": {
            "exon_51": {
                "patient_coverage": "~13% of DMD patients",
                "approved_ASO": "eteplirsen (Exondys 51, FDA 2016)",
                "crispr_approach": "CRISPR disruption of exon 51 splice sites "
                                   "to induce permanent exon skipping and "
                                   "restore reading frame. HG302 (HuidaGene) "
                                   "uses hfCas12Max to edit exon 51 splice-donor.",
            },
            "exon_53": {
                "patient_coverage": "~8% of DMD patients",
                "approved_ASO": "golodirsen (Vyondys 53, FDA 2019)",
                "crispr_approach": "CRISPR disruption of exon 53 splice "
                                   "acceptor/donor for permanent skipping",
            },
            "exon_45": {
                "patient_coverage": "~8% of DMD patients",
                "approved_ASO": "casimersen (Amondys 45, FDA 2021)",
                "crispr_approach": "Base editing of exon 45 splice sites "
                                   "(demonstrated preclinically, Mol Ther "
                                   "Nucleic Acids 2023)",
            },
        },
        "strategy": (
            "CRISPR approaches applicable to ~85% of DMD patients: "
            "(1) Exon skipping -- single-guide disruption of splice donor/acceptor "
            "sites flanking deletable exons to restore reading frame (converts "
            "DMD to milder Becker phenotype); (2) Exon excision -- dual-guide "
            "deletion of one or more exons; (3) Exon reframing -- small indels "
            "at exon junctions to restore frame; (4) Base editing -- precise "
            "splice-site modification without DSBs; (5) Prime editing -- precise "
            "correction of point mutations; (6) Micro-dystrophin gene replacement "
            "via AAV (see Elevidys). Multi-exon skipping of exons 45-55 could "
            "treat ~63% of all DMD patients."
        ),
        "clinical_programs": {
            "HG302_HuidaGene": {
                "type": "CRISPR-Cas12 exon 51 splice-site editing (AAV)",
                "route": "IV",
                "status": "Phase I MUSCLE trial (Shanghai) -- first patient "
                          "dosed Dec 2024. First-in-human CRISPR gene editing "
                          "for DMD. Early data (ASGCT May 2025) show dystrophin "
                          "restoration with no adverse effects. Uses AI-optimized "
                          "hfCas12Max nuclease.",
            },
            "Elevidys_Sarepta": {
                "type": "AAV-rh74 micro-dystrophin gene replacement",
                "route": "IV",
                "status": "FDA accelerated approval June 2023 for ambulatory "
                          "DMD patients 4-5 years. NOTE: two non-ambulatory "
                          "patients died of acute liver failure ~2 months "
                          "post-dosing (2025); enrollment of non-ambulatory "
                          "patients halted.",
            },
        },
        "conditions": ["duchenne_muscular_dystrophy", "becker_muscular_dystrophy",
                        "muscular_dystrophy"],
    },

    # ------------------------------------------------------------------
    # 2b. Myotonic Dystrophy Type 1
    # ------------------------------------------------------------------
    "DMPK": {
        "chrom": "chr19",
        "start": 45769717,
        "end": 45782490,
        "strand": "-",
        "cytoband": "19q13.32",
        "ncbi_gene_id": 1760,
        "refseq": "NC_000019.10",
        "role": "DM1 protein kinase (myotonin) -- serine/threonine kinase. "
                "CTG expansion in 3'UTR causes toxic RNA gain-of-function: "
                "expanded CUG-repeat RNA sequesters MBNL1 splicing factor "
                "into nuclear foci, causing widespread mis-splicing",
        "disease": "Myotonic dystrophy type 1 (DM1, Steinert disease)",
        "omim_disease": 160900,
        "omim_gene": 605377,
        "inheritance": "AD",
        "mutation_type": "trinucleotide_repeat_expansion",
        "repeat_details": {
            "repeat_unit": "CTG",
            "location": "3' untranslated region (3'UTR) of DMPK",
            "normal_range": "5-34 repeats",
            "premutation": "35-49 repeats (may expand in offspring)",
            "mild_DM1": "50-150 repeats (late onset, cataracts, mild myotonia)",
            "classic_DM1": "100-1000 repeats (adult onset, progressive weakness)",
            "congenital_DM1": ">1000 repeats (severe neonatal presentation)",
            "mechanism": "Expanded CUG-repeat RNA forms nuclear foci that "
                        "sequester MBNL1/MBNL2 splicing regulators, causing "
                        "mis-splicing of hundreds of transcripts (CLCN1, "
                        "insulin receptor, cardiac troponin T, etc.). Also "
                        "upregulates CELF1. Additionally, antisense "
                        "transcription through the locus is affected.",
        },
        "key_variants": [
            "CTG expansion in 3'UTR (50 to >3000 repeats in patients)",
        ],
        "strategy": (
            "(1) CRISPRi silencing -- dCas9-KRAB targeted to DMPK promoter "
            "reduces CUG-expanded transcripts up to 80%, eliminating nuclear "
            "foci and correcting mis-splicing (Mol Ther Nucleic Acids 2023); "
            "(2) Direct CTG repeat deletion -- dual SaCas9 guides flanking "
            "the CTG expansion excise the repeat tract; demonstrated in DM1 "
            "patient myoblasts carrying 2600 CTG repeats, eliminating foci "
            "and correcting splicing; (3) Cas9 nickase repeat contraction -- "
            "repeat-targeted nickase induces contraction of expanded CTGs; "
            "(4) Allele-specific approaches -- target expanded allele via "
            "flanking SNPs or repeat-length-dependent accessibility; "
            "(5) RNA-targeting -- Cas13/CasRx to degrade CUG-expanded "
            "transcripts without genomic editing"
        ),
        "clinical_programs": {
            "note": "No CRISPR clinical trials yet for DM1 as of early 2026. "
                    "DYNE-101 (Dyne Therapeutics) is an antibody-conjugated "
                    "ASO targeting DMPK in Phase I/II (2024-2025). Multiple "
                    "preclinical CRISPR programs active.",
        },
        "conditions": ["myotonic_dystrophy", "muscular_dystrophy",
                        "trinucleotide_repeat"],
    },
}


# ==========================================================================
# 3. HEREDITARY BLINDNESS / RETINAL DYSTROPHIES
# ==========================================================================

HEREDITARY_BLINDNESS_TARGETS = {

    # ------------------------------------------------------------------
    # 3a. Leber Congenital Amaurosis (RPE65)
    # ------------------------------------------------------------------
    "RPE65": {
        "chrom": "chr1",
        "start": 68428822,
        "end": 68449954,
        "strand": "-",
        "cytoband": "1p31.2",
        "ncbi_gene_id": 6121,
        "refseq": "NC_000001.11",
        "role": "Retinoid isomerohydrolase RPE65 -- essential enzyme in retinal "
                "pigment epithelium (RPE) that converts all-trans-retinyl esters "
                "to 11-cis-retinol, a critical step in the visual (retinoid) "
                "cycle for photoreceptor function",
        "disease": "Leber congenital amaurosis type 2 (LCA2) / "
                   "RPE65-associated retinal dystrophy",
        "omim_disease": 204100,
        "omim_gene": 180069,
        "inheritance": "AR",
        "mutation_type": "missense, nonsense, splice_site, frameshift",
        "key_variants": [
            "Various loss-of-function mutations (>100 reported)",
            "Most frequent: missense mutations disrupting iron-binding site "
            "or protein folding",
        ],
        "strategy": (
            "Gene replacement (Luxturna) is the established standard. CRISPR "
            "approaches: (1) Base editing to correct specific missense "
            "mutations in RPE65; (2) HDR-mediated correction of pathogenic "
            "variants in RPE cells. For most patients, AAV-mediated gene "
            "replacement remains optimal due to the diversity of mutations "
            "and RPE65's small size (1602 bp CDS) fitting easily in AAV."
        ),
        "clinical_programs": {
            "Luxturna_Spark_Novartis": {
                "type": "AAV2-RPE65 gene replacement",
                "route": "subretinal injection",
                "status": "FDA-approved December 2017. FIRST FDA-approved gene "
                          "therapy for a genetic disease. Delivers functional "
                          "RPE65 to RPE cells. Demonstrated durable improvement "
                          "in navigational vision and light sensitivity. "
                          "~$425K per eye. Requires biallelic RPE65 mutations "
                          "with viable retinal cells.",
            },
        },
        "conditions": ["leber_congenital_amaurosis", "retinal_dystrophy",
                        "hereditary_blindness"],
    },

    # ------------------------------------------------------------------
    # 3b. LCA10 (CEP290)
    # ------------------------------------------------------------------
    "CEP290": {
        "chrom": "chr12",
        "start": 88049016,
        "end": 88142088,
        "strand": "-",
        "cytoband": "12q21.32",
        "ncbi_gene_id": 80184,
        "refseq": "NC_000012.12",
        "role": "Centrosomal protein 290 kDa -- essential for cilium assembly "
                "and maintenance in photoreceptors. Required for transition "
                "zone formation in connecting cilium of rods and cones",
        "disease": "Leber congenital amaurosis type 10 (LCA10)",
        "omim_disease": 611755,
        "omim_gene": 610142,
        "inheritance": "AR",
        "mutation_type": "deep_intronic (most common), plus coding mutations",
        "key_variants": [
            "c.2991+1655A>G (IVS26 mutation) -- deep intronic mutation in "
            "intron 26 creating a cryptic splice donor site. Inserts a 128-bp "
            "cryptic exon containing a premature stop codon (p.Cys998X). "
            "Accounts for ~15% of all LCA cases in Western populations. "
            "Most common single LCA mutation worldwide.",
        ],
        "ivs26_mutation_details": {
            "hgvs": "c.2991+1655A>G",
            "location": "intron 26, 1655 bp downstream of exon 26",
            "effect": "Creates cryptic splice donor site, leading to inclusion "
                     "of 128-bp pseudoexon with premature stop codon",
            "protein_effect": "p.Cys998X (truncation)",
            "allele_frequency": "~1% in European populations",
            "approximate_grch38_position": 88071000,
        },
        "strategy": (
            "(1) CRISPR excision of IVS26 cryptic exon -- EDIT-101 used paired "
            "SaCas9 guides to excise the region containing the aberrant splice "
            "donor site created by c.2991+1655A>G, restoring normal CEP290 "
            "splicing; (2) Single-cut disruption of cryptic splice site; "
            "(3) Base editing to revert A>G mutation; "
            "(4) Gene replacement is challenging because CEP290 CDS is ~7.4 kb, "
            "exceeding AAV packaging capacity (~4.7 kb), requiring dual-AAV or "
            "alternative vectors. This makes CRISPR correction of the IVS26 "
            "mutation particularly attractive as it bypasses the size constraint."
        ),
        "clinical_programs": {
            "EDIT-101_Editas": {
                "type": "AAV5-SaCas9 + dual sgRNAs (in vivo retinal CRISPR)",
                "route": "subretinal injection",
                "status": "BRILLIANCE Phase I/II -- 14 patients treated. "
                          "3/14 showed clinically meaningful BCVA improvement. "
                          "Well-tolerated, no ocular SAEs. HOWEVER: program "
                          "paused by Editas (2023-2024) due to small eligible "
                          "patient population (~300 homozygous IVS26 in US). "
                          "Editas seeking collaboration partner. FIRST in vivo "
                          "CRISPR medicine trial in humans (dosed 2020).",
            },
        },
        "conditions": ["leber_congenital_amaurosis", "LCA10",
                        "retinal_dystrophy", "hereditary_blindness",
                        "ciliopathy"],
    },

    # ------------------------------------------------------------------
    # 3c. Retinitis Pigmentosa -- RPGR (X-linked)
    # ------------------------------------------------------------------
    "RPGR": {
        "chrom": "chrX",
        "start": 38269163,
        "end": 38327509,
        "strand": "-",
        "cytoband": "Xp11.4",
        "ncbi_gene_id": 6103,
        "refseq": "NC_000023.11",
        "role": "Retinitis pigmentosa GTPase regulator -- ciliary protein "
                "essential for photoreceptor outer segment disc morphogenesis "
                "and protein trafficking. Multiple isoforms: constitutive "
                "RPGR(1-19), retina-specific RPGR-ORF15 (predominant isoform "
                "in photoreceptors), and human-specific RPGR(s14/15)",
        "disease": "X-linked retinitis pigmentosa (XLRP, RP3)",
        "omim_disease": 300029,
        "omim_gene": 312610,
        "inheritance": "XLR",
        "mutation_type": "frameshift, missense, nonsense (majority in ORF15 "
                         "exon -- a repetitive, purine-rich region prone to "
                         "mutations)",
        "key_variants": [
            "ORF15 mutations (~60-80% of RPGR-XLRP)",
            "Exon 1-14 mutations (remaining cases)",
            "RPGR mutations cause >70% of all X-linked RP",
        ],
        "strategy": (
            "(1) AAV-RPGR gene replacement -- multiple clinical trials active; "
            "delivers codon-optimized RPGR-ORF15 via AAV2/5 or AAV8; "
            "(2) CRISPR frame restoration in ORF15 -- single-guide indels "
            "to restore reading frame in the repetitive ORF15 region; "
            "(3) Base editing for point mutations in exons 1-14"
        ),
        "clinical_programs": {
            "botaretigene_sparoparvovec_Janssen_MeiraGTx": {
                "type": "AAV-RPGR-ORF15 gene replacement",
                "route": "subretinal injection",
                "status": "Phase III trial active (NCT04671433). Phase I/II "
                          "data showed significant vision improvement at "
                          "6-month timepoint. Efficacy/safety data presented "
                          "at FFB Retinal Therapeutics Summit 2025. Leading "
                          "RPGR gene therapy program.",
            },
            "AGTC-501_Beacon": {
                "type": "AAV-RPGR gene replacement",
                "route": "subretinal injection",
                "status": "Phase I/II HORIZON trial -- 24-month data show "
                          "50% of patients at highest dose maintained >=7 dB "
                          "improvement in >=5 loci.",
            },
        },
        "conditions": ["retinitis_pigmentosa", "XLRP", "retinal_dystrophy",
                        "hereditary_blindness"],
    },

    # ------------------------------------------------------------------
    # 3d. Retinitis Pigmentosa -- RHO (autosomal dominant)
    # ------------------------------------------------------------------
    "RHO": {
        "chrom": "chr3",
        "start": 129528639,
        "end": 129535344,
        "strand": "+",
        "cytoband": "3q22.1",
        "ncbi_gene_id": 6010,
        "refseq": "NC_000003.12",
        "role": "Rhodopsin -- the visual pigment of rod photoreceptors. "
                "G-protein-coupled receptor that initiates phototransduction "
                "cascade upon photon absorption. Most commonly mutated gene "
                "in autosomal dominant RP (~25-30% of adRP cases)",
        "disease": "Autosomal dominant retinitis pigmentosa (adRP, RP4)",
        "omim_disease": 613731,
        "omim_gene": 180380,
        "inheritance": "AD (most mutations); AR (rare)",
        "mutation_type": "missense (majority), nonsense, small indels",
        "key_variants": [
            "P23H (rs28933963) -- most common RHO mutation in North America "
            "(~10% of all adRP), causes rhodopsin misfolding and ER stress",
            "T58R -- severe, early onset",
            "P347L, P347S, P347R -- C-terminal mutations affecting trafficking",
            "G106R -- structural mutation in transmembrane domain",
            ">150 pathogenic variants reported",
        ],
        "strategy": (
            "(1) Allele-specific CRISPR knockout -- selectively disrupt "
            "mutant RHO allele using SNP-guided PAM sites (one functional "
            "RHO allele is sufficient for rod function); "
            "(2) Ablate-and-replace -- CRISPR knockout of both endogenous "
            "RHO alleles + AAV delivery of codon-optimized, guide-resistant "
            "RHO transgene; (3) Base editing for P23H and other specific "
            "missense mutations; (4) CRISPRi of mutant allele; "
            "(5) Dominant-negative mutations make simple gene augmentation "
            "insufficient -- must silence the mutant allele"
        ),
        "clinical_programs": {
            "note": "No CRISPR clinical trials for RHO-RP as of early 2026. "
                    "ProQR (sepofarsen/QR-421a) ASO programs for various RP "
                    "genes. Multiple preclinical CRISPR ablate-and-replace "
                    "studies showing efficacy in RHO-P23H mouse models.",
        },
        "conditions": ["retinitis_pigmentosa", "adRP", "retinal_dystrophy",
                        "hereditary_blindness"],
    },

    # ------------------------------------------------------------------
    # 3e. Usher Syndrome Type 1B -- MYO7A
    # ------------------------------------------------------------------
    "MYO7A": {
        "chrom": "chr11",
        "start": 77128246,
        "end": 77215241,
        "strand": "+",
        "cytoband": "11q13.5",
        "ncbi_gene_id": 4647,
        "refseq": "NC_000011.10",
        "role": "Myosin VIIA -- unconventional myosin motor protein essential "
                "for stereocilia organization in hair cells and RPE "
                "melanosome transport in retina. CDS ~6.6 kb (exceeds AAV "
                "packaging limit)",
        "disease": "Usher syndrome type 1B (USH1B) -- congenital profound "
                   "deafness + vestibular dysfunction + progressive RP",
        "omim_disease": 276900,
        "omim_gene": 276903,
        "inheritance": "AR",
        "mutation_type": "missense, nonsense, frameshift, splice_site "
                         "(>400 pathogenic variants reported)",
        "key_variants": [
            "Widely heterogeneous -- no single common mutation",
            "R634X, R1240X (nonsense -- common in some populations)",
            "c.3503G>A (p.G1168D), c.5648G>A (p.R1883Q)",
        ],
        "strategy": (
            "(1) Dual-AAV gene replacement -- MYO7A CDS (~6.6 kb) exceeds "
            "AAV packaging limit, requiring split-vector strategy. "
            "AAVantgarde Bio uses dual-AAV intein-mediated protein trans-splicing; "
            "(2) CRISPR correction of specific mutations in patient iPSC-derived "
            "retinal organoids; (3) Base editing for prevalent missense mutations; "
            "(4) Exon skipping for in-frame deletable exons"
        ),
        "clinical_programs": {
            "AAVantgarde_Bio": {
                "type": "Dual-AAV MYO7A gene replacement (split intein)",
                "route": "subretinal injection",
                "status": "Phase I/II completed dosing (2025). Dual-AAV "
                          "delivers MYO7A gene split across two AAV vectors "
                          "that reconstitute full-length protein in retinal "
                          "cells via intein-mediated trans-splicing.",
            },
        },
        "conditions": ["usher_syndrome", "USH1B", "retinal_dystrophy",
                        "hereditary_blindness", "deafness"],
    },

    # ------------------------------------------------------------------
    # 3f. Usher Syndrome Type 2A / Non-syndromic RP -- USH2A
    # ------------------------------------------------------------------
    "USH2A": {
        "chrom": "chr1",
        "start": 215622891,
        "end": 216423448,
        "strand": "-",
        "cytoband": "1q41",
        "ncbi_gene_id": 7399,
        "refseq": "NC_000001.11",
        "gene_size_bp": 800557,
        "role": "Usherin -- large transmembrane protein (5202 aa) in the "
                "extracellular matrix of photoreceptors and hair cells. "
                "Essential for long-term maintenance of photoreceptor and "
                "cochlear hair cell structure. CDS ~15.6 kb (far exceeds "
                "AAV packaging limit)",
        "disease": "Usher syndrome type 2A (USH2A) -- moderate-to-severe "
                   "congenital hearing loss + progressive RP, OR non-syndromic "
                   "autosomal recessive RP (arRP)",
        "omim_disease": 276901,
        "omim_gene": 608400,
        "inheritance": "AR",
        "mutation_type": "missense, nonsense, frameshift, splice_site",
        "key_variants": [
            "c.2299delG (p.Glu767Serfs*21) -- most common worldwide, "
            "especially in European populations",
            "c.2276G>T (p.Cys759Phe) -- common in East Asian populations",
            "Exon 13 mutations -- therapeutic target for ASO approaches",
            ">400 pathogenic variants across 72 exons",
        ],
        "strategy": (
            "USH2A CDS is ~15.6 kb -- far too large for AAV gene replacement. "
            "Strategies: (1) Exon skipping (ASO or CRISPR) -- particularly "
            "exon 13 skip to treat exon 13 mutations; ultevursen (ASO) is "
            "in Phase 2b LUNA trial; (2) CRISPR correction of common point "
            "mutations (c.2299delG) via HDR or prime editing; "
            "(3) Truncated mini-usherin gene replacement (identifying minimal "
            "functional domains); (4) Dual/triple-AAV strategies for partial "
            "gene delivery; (5) CRISPR base editing for prevalent missense "
            "variants"
        ),
        "clinical_programs": {
            "ultevursen_Sepul_Bio": {
                "type": "ASO (exon 13 skipping)",
                "route": "intravitreal injection",
                "status": "Phase 2b LUNA trial -- first patients dosed 2025. "
                          "Targets USH2A exon 13 mutations specifically. "
                          "Originally ProQR, now Sepul Bio (Thea subsidiary).",
            },
            "FFB_Natural_History": {
                "type": "Natural history study",
                "status": "4-year data released Nov 2025. Critical for "
                          "understanding progression and defining endpoints "
                          "for future gene therapy trials.",
            },
        },
        "conditions": ["usher_syndrome", "USH2A", "retinitis_pigmentosa",
                        "retinal_dystrophy", "hereditary_blindness", "deafness"],
    },
}


# ==========================================================================
# 4. SKIN DISEASES -- EPIDERMOLYSIS BULLOSA
# ==========================================================================

SKIN_DISEASE_TARGETS = {

    # ------------------------------------------------------------------
    # 4a. Dystrophic Epidermolysis Bullosa -- COL7A1
    # ------------------------------------------------------------------
    "COL7A1": {
        "chrom": "chr3",
        "start": 48564073,
        "end": 48595329,
        "strand": "-",
        "cytoband": "3p21.31",
        "ncbi_gene_id": 1294,
        "refseq": "NC_000003.12",
        "role": "Collagen type VII alpha 1 chain -- major component of anchoring "
                "fibrils that attach epidermis to dermis at the dermal-epidermal "
                "junction basement membrane zone",
        "disease": "Dystrophic epidermolysis bullosa (DEB) -- dominant (DDEB) "
                   "and recessive (RDEB) forms",
        "omim_disease_rdeb": 226600,
        "omim_disease_ddeb": 131750,
        "omim_gene": 120120,
        "inheritance": "AD (DDEB) or AR (RDEB)",
        "mutation_type": "missense, nonsense, splice_site, frameshift",
        "key_variants": [
            "Glycine substitutions in collagenous domain (dominant)",
            "Premature stop codons throughout (recessive)",
            ">800 pathogenic variants reported across 118 exons",
            "107 of 118 exons can be individually removed without disrupting "
            "the reading frame, making exon skipping broadly applicable",
        ],
        "strategy": (
            "(1) Topical HSV-1 gene replacement (Vyjuvek) -- approved therapy; "
            "(2) CRISPR exon skipping -- dual-sgRNA/Cas9 RNP deletion of "
            "mutated exons (up to 95% deletion efficiency demonstrated 2025); "
            "107/118 exons are skippable in-frame, making this broadly "
            "applicable; (3) HDR correction of specific point mutations in "
            "patient keratinocytes ex vivo, followed by autologous skin graft; "
            "(4) Base editing for common missense/nonsense mutations; "
            "(5) In vivo CRISPR delivery via adenoviral vectors (preclinical)"
        ),
        "clinical_programs": {
            "Vyjuvek_Krystal_Bio": {
                "type": "HSV-1 vector delivering COL7A1 transgene (topical)",
                "route": "topical application to wounds",
                "status": "FDA-approved May 2023. FIRST topical gene therapy. "
                          "First redosable gene therapy. HSV-1 vector delivers "
                          "functional COL7A1 to wound keratinocytes and "
                          "fibroblasts. Label expanded Sept 2025 to include "
                          "patients from birth. Promotes C7 deposition at "
                          "basement membrane zone and wound healing.",
            },
            "CRISPR_exon_skipping_preclinical": {
                "type": "CRISPR-Cas9 RNP exon removal (ex vivo)",
                "status": "Preclinical 2025 -- demonstrated deletion of exons "
                          "73 and 105 restores C7 production in RDEB patient "
                          "cells; functional C7 deposited in regenerated skin "
                          "grafts on immunocompromised mice (Human Gene Therapy "
                          "2025).",
            },
        },
        "conditions": ["epidermolysis_bullosa", "DEB", "RDEB", "skin_disease"],
    },

    # ------------------------------------------------------------------
    # 4b. Epidermolysis Bullosa Simplex -- KRT5
    # ------------------------------------------------------------------
    "KRT5": {
        "chrom": "chr12",
        "start": 52514575,
        "end": 52520394,
        "strand": "-",
        "cytoband": "12q13.13",
        "ncbi_gene_id": 3852,
        "refseq": "NC_000012.12",
        "role": "Keratin 5 -- type II intermediate filament protein expressed "
                "in basal keratinocytes. Forms heterodimers with KRT14 to "
                "create the cytoskeletal network providing mechanical "
                "resilience to basal epidermal cells",
        "disease": "Epidermolysis bullosa simplex (EBS) -- intraepidermal "
                   "blistering due to basal keratinocyte fragility",
        "omim_disease": 131900,
        "omim_gene": 148040,
        "inheritance": "AD (most forms), AR (rare severe forms)",
        "mutation_type": "missense (hotspots in helix boundary motifs)",
        "key_variants": [
            "p.Glu475Gly (E475G) -- rod domain helix initiation motif, "
            "severe generalized EBS (Dowling-Meara)",
            "p.Leu469Pro -- severe generalized EBS",
            "Mutations cluster in helix initiation (1A) and termination "
            "(2B) peptides of the rod domain",
        ],
        "strategy": (
            "(1) Allele-specific CRISPR knockout of dominant-negative mutant "
            "allele (single functional KRT5 allele is sufficient); "
            "(2) Ablate-and-replace -- knock out both alleles + deliver "
            "guide-resistant KRT5 transgene; (3) Base editing for specific "
            "missense mutations in helix boundary motifs; "
            "(4) Ex vivo correction of autologous keratinocytes + skin graft"
        ),
        "clinical_programs": {
            "note": "No approved gene therapies for EBS-KRT5 as of early 2026. "
                    "Preclinical CRISPR studies demonstrate allele-specific "
                    "knockout feasibility.",
        },
        "conditions": ["epidermolysis_bullosa", "EBS", "skin_disease"],
    },

    # ------------------------------------------------------------------
    # 4c. Epidermolysis Bullosa Simplex -- KRT14
    # ------------------------------------------------------------------
    "KRT14": {
        "chrom": "chr17",
        "start": 41582279,
        "end": 41586895,
        "strand": "-",
        "cytoband": "17q21.2",
        "ncbi_gene_id": 3861,
        "refseq": "NC_000017.11",
        "role": "Keratin 14 -- type I intermediate filament protein that "
                "heterodimerizes with KRT5 in basal keratinocytes. Mutations "
                "disrupt cytoskeletal integrity",
        "disease": "Epidermolysis bullosa simplex (EBS)",
        "omim_disease": 131900,
        "omim_gene": 148066,
        "inheritance": "AD (most), AR (severe generalized)",
        "mutation_type": "missense (helix boundary motifs), nonsense",
        "key_variants": [
            "p.Arg125Cys (R125C) -- helix initiation motif, severe EBS",
            "p.Arg125His (R125H) -- same residue, milder phenotype",
            "p.Leu384Pro -- rod domain, Dowling-Meara EBS",
            "p.Val133_Leu136del -- in-frame deletion, AR severe",
        ],
        "strategy": (
            "(1) Allele-specific CRISPR knockout of dominant-negative KRT14 "
            "mutant allele; (2) Base editing for recurrent missense mutations "
            "(R125C, R125H); (3) Ex vivo HDR correction in patient "
            "keratinocytes for autologous grafting; (4) CRISPRi silencing "
            "of mutant allele using allele-specific guides"
        ),
        "clinical_programs": {
            "note": "No approved gene therapies for EBS-KRT14 as of early 2026. "
                    "Preclinical CRISPR allele-specific strategies in development.",
        },
        "conditions": ["epidermolysis_bullosa", "EBS", "skin_disease"],
    },

    # ------------------------------------------------------------------
    # 4d. Junctional Epidermolysis Bullosa -- LAMB3
    # ------------------------------------------------------------------
    "LAMB3": {
        "chrom": "chr1",
        "start": 209614870,
        "end": 209652425,
        "strand": "-",
        "cytoband": "1q32.2",
        "ncbi_gene_id": 3914,
        "refseq": "NC_000001.11",
        "role": "Laminin subunit beta-3 -- component of laminin-332 "
                "(previously laminin-5), a key adhesion ligand in the "
                "basement membrane that anchors epidermis to dermis via "
                "hemidesmosomes and anchoring filaments",
        "disease": "Junctional epidermolysis bullosa (JEB) -- blistering at "
                   "the lamina lucida of the basement membrane",
        "omim_disease": 226700,
        "omim_gene": 150310,
        "inheritance": "AR",
        "mutation_type": "nonsense, frameshift, splice_site (mostly "
                         "loss-of-function)",
        "key_variants": [
            "c.1903C>T (p.Arg635Ter) -- recurrent nonsense mutation",
            "c.628G>A (p.Glu210Lys) -- missense",
            "Various truncating mutations across 23 exons",
        ],
        "strategy": (
            "(1) Ex vivo gene replacement + autologous skin graft -- landmark "
            "case: patient with JEB-LAMB3 received transgenic epidermal "
            "grafts covering ~80% body surface area using retroviral-corrected "
            "autologous keratinocytes (Bochum, 2017; Nature 2017), with "
            "sustained skin regeneration >5 years; (2) CRISPR-HDR correction "
            "in patient keratinocyte stem cells ex vivo; (3) Base editing "
            "for nonsense mutations (e.g., R635X); (4) In vivo HSV-1 vector "
            "delivery (following Vyjuvek model)"
        ),
        "clinical_programs": {
            "Holostem_EB-101": {
                "type": "Retroviral LAMB3 gene replacement in autologous "
                        "epidermal stem cells",
                "route": "ex vivo corrected skin grafts",
                "status": "Holostem/Chiesi -- clinical development for JEB. "
                          "Based on the landmark Bochum whole-body skin "
                          "regeneration case (2017).",
            },
        },
        "conditions": ["epidermolysis_bullosa", "JEB", "skin_disease"],
    },
}


# ==========================================================================
# 5. SYSTEMIC HEREDITARY DISEASES
# ==========================================================================

SYSTEMIC_HEREDITARY_TARGETS = {

    # ------------------------------------------------------------------
    # 5a. Transthyretin Amyloidosis -- TTR
    # ------------------------------------------------------------------
    "TTR": {
        "chrom": "chr18",
        "start": 31591877,
        "end": 31598821,
        "strand": "+",
        "cytoband": "18q12.1",
        "ncbi_gene_id": 7276,
        "refseq": "NC_000018.10",
        "role": "Transthyretin -- liver-produced serum protein that transports "
                "thyroxine (T4) and retinol-binding protein. Tetramer that "
                "becomes unstable with pathogenic mutations, dissociates into "
                "monomers that misfold and aggregate into amyloid fibrils "
                "depositing in heart, nerves, and other tissues",
        "disease": "Hereditary transthyretin amyloidosis (hATTR) -- "
                   "polyneuropathy (ATTRv-PN) and/or cardiomyopathy (ATTR-CM)",
        "omim_disease": 105210,
        "omim_gene": 176300,
        "inheritance": "AD",
        "mutation_type": "missense (>130 amyloidogenic variants)",
        "key_variants": [
            "V30M (Val50Met per mature protein; rs28933979) -- most common "
            "worldwide, endemic in Portugal/Sweden/Japan. Causes FAP "
            "(familial amyloid polyneuropathy)",
            "V122I (Val142Ile; rs76992529) -- most common in African Americans "
            "(3.4% carrier frequency), primarily cardiac phenotype",
            "T60A (Thr80Ala) -- common in Ireland/UK, mixed phenotype",
            "E89Q, S77Y, I84S -- other common amyloidogenic variants",
        ],
        "strategy": (
            "(1) In vivo CRISPR gene knockout (NTLA-2001/nex-z) -- LNP-delivered "
            "Cas9+sgRNA targeting TTR in hepatocytes to permanently inactivate "
            "TTR production. Single IV infusion. Landmark: first in vivo CRISPR "
            "therapy in humans (2021); (2) Base editing to correct specific "
            "pathogenic variants while preserving TTR function; "
            "(3) siRNA (patisiran/Onpattro) and ASO (inotersen/Tegsedi) for "
            "chronic TTR knockdown (approved therapies); "
            "(4) TTR stabilizers (tafamidis/Vyndaqel, diflunisal) to prevent "
            "tetramer dissociation"
        ),
        "clinical_programs": {
            "nexiguran_ziclumeran_nex_z_Intellia": {
                "type": "LNP-delivered Cas9+sgRNA (in vivo hepatic CRISPR)",
                "route": "single IV infusion",
                "status": "Phase III MAGNITUDE (ATTR-CM) and MAGNITUDE-2 "
                          "(ATTRv-PN) trials. Demonstrated rapid, deep, durable "
                          "TTR reduction (~90%) from single dose. HOWEVER: "
                          "clinical hold placed Oct 2025 after grade 4 hepatic "
                          "AE in one ATTR-CM patient who subsequently died "
                          "(Nov 2025). Trials paused pending safety review. "
                          "FIRST in vivo CRISPR medicine administered to "
                          "humans (June 2021). Phase I showed TTR knockdown "
                          "of 87% at 0.3 mg/kg dose.",
            },
            "patisiran_Onpattro_Alnylam": {
                "type": "siRNA (lipid nanoparticle)",
                "route": "IV every 3 weeks",
                "status": "FDA-approved 2018 for ATTRv-PN. First RNAi therapy.",
            },
            "vutrisiran_Amvuttra_Alnylam": {
                "type": "siRNA (GalNAc-conjugated, subcutaneous)",
                "route": "SC every 3 months",
                "status": "FDA-approved 2022 for ATTRv-PN.",
            },
        },
        "conditions": ["transthyretin_amyloidosis", "hATTR", "ATTR_CM",
                        "ATTR_PN", "amyloidosis", "cardiomyopathy",
                        "polyneuropathy"],
    },
}


# ==========================================================================
# UNIFIED ACCESS: all targets in a single dictionary
# ==========================================================================

ALL_HEREDITARY_TARGETS = {
    **NEURODEGENERATIVE_TARGETS,
    **MUSCULAR_DYSTROPHY_TARGETS,
    **HEREDITARY_BLINDNESS_TARGETS,
    **SKIN_DISEASE_TARGETS,
    **SYSTEMIC_HEREDITARY_TARGETS,
}


# ==========================================================================
# Disease-to-gene mapping for pipeline routing
# ==========================================================================

HEREDITARY_DISEASES = {
    "huntington_disease": {
        "primary_gene": "HTT",
        "omim": 143100,
        "inheritance": "AD",
        "key_pathway": "protein_aggregation, transcription_dysregulation",
        "mutation_class": "trinucleotide_repeat",
    },
    "spinal_muscular_atrophy": {
        "primary_gene": "SMN1",
        "modifier_genes": ["SMN2"],
        "omim": 253300,
        "inheritance": "AR",
        "key_pathway": "snRNP_assembly, motor_neuron_survival",
        "mutation_class": "deletion/conversion",
    },
    "ALS_SOD1": {
        "primary_gene": "SOD1",
        "omim": 105400,
        "inheritance": "AD/AR",
        "key_pathway": "oxidative_stress, protein_aggregation",
        "mutation_class": "missense_gain_of_function",
    },
    "friedreich_ataxia": {
        "primary_gene": "FXN",
        "omim": 229300,
        "inheritance": "AR",
        "key_pathway": "mitochondrial_iron_sulfur_cluster, oxidative_stress",
        "mutation_class": "trinucleotide_repeat",
    },
    "duchenne_muscular_dystrophy": {
        "primary_gene": "DMD",
        "omim": 310200,
        "inheritance": "XLR",
        "key_pathway": "dystrophin_glycoprotein_complex, membrane_stability",
        "mutation_class": "large_deletion/duplication/point",
    },
    "becker_muscular_dystrophy": {
        "primary_gene": "DMD",
        "omim": 300376,
        "inheritance": "XLR",
        "key_pathway": "dystrophin_glycoprotein_complex",
        "mutation_class": "in_frame_deletion",
    },
    "myotonic_dystrophy_1": {
        "primary_gene": "DMPK",
        "omim": 160900,
        "inheritance": "AD",
        "key_pathway": "RNA_toxicity, MBNL1_sequestration, mis_splicing",
        "mutation_class": "trinucleotide_repeat",
    },
    "LCA2_RPE65": {
        "primary_gene": "RPE65",
        "omim": 204100,
        "inheritance": "AR",
        "key_pathway": "retinoid_cycle, visual_transduction",
        "mutation_class": "loss_of_function",
        "approved_therapy": "Luxturna",
    },
    "LCA10_CEP290": {
        "primary_gene": "CEP290",
        "omim": 611755,
        "inheritance": "AR",
        "key_pathway": "ciliogenesis, photoreceptor_outer_segment",
        "mutation_class": "deep_intronic_splice",
    },
    "XLRP_RPGR": {
        "primary_gene": "RPGR",
        "omim": 300029,
        "inheritance": "XLR",
        "key_pathway": "ciliary_transport, photoreceptor_maintenance",
        "mutation_class": "frameshift_ORF15",
    },
    "adRP_RHO": {
        "primary_gene": "RHO",
        "omim": 613731,
        "inheritance": "AD",
        "key_pathway": "phototransduction, rhodopsin_folding",
        "mutation_class": "dominant_negative_missense",
    },
    "usher_1B": {
        "primary_gene": "MYO7A",
        "omim": 276900,
        "inheritance": "AR",
        "key_pathway": "stereocilia_organization, RPE_transport",
        "mutation_class": "loss_of_function",
    },
    "usher_2A": {
        "primary_gene": "USH2A",
        "omim": 276901,
        "inheritance": "AR",
        "key_pathway": "photoreceptor_cilium, hair_cell_maintenance",
        "mutation_class": "loss_of_function",
    },
    "dystrophic_EB": {
        "primary_gene": "COL7A1",
        "omim": 226600,
        "inheritance": "AR/AD",
        "key_pathway": "anchoring_fibrils, dermal_epidermal_junction",
        "mutation_class": "loss_of_function",
        "approved_therapy": "Vyjuvek",
    },
    "EB_simplex_KRT5": {
        "primary_gene": "KRT5",
        "omim": 131900,
        "inheritance": "AD",
        "key_pathway": "keratin_cytoskeleton, basal_keratinocyte_integrity",
        "mutation_class": "dominant_negative_missense",
    },
    "EB_simplex_KRT14": {
        "primary_gene": "KRT14",
        "omim": 131900,
        "inheritance": "AD",
        "key_pathway": "keratin_cytoskeleton, basal_keratinocyte_integrity",
        "mutation_class": "dominant_negative_missense",
    },
    "junctional_EB_LAMB3": {
        "primary_gene": "LAMB3",
        "omim": 226700,
        "inheritance": "AR",
        "key_pathway": "laminin_332, hemidesmosome_adhesion",
        "mutation_class": "loss_of_function",
    },
    "hATTR_amyloidosis": {
        "primary_gene": "TTR",
        "omim": 105210,
        "inheritance": "AD",
        "key_pathway": "protein_misfolding, amyloid_deposition",
        "mutation_class": "missense_destabilizing",
    },
}


# ==========================================================================
# Trinucleotide repeat reference table
# ==========================================================================

TRINUCLEOTIDE_REPEAT_DISEASES = {
    "HTT": {
        "repeat_unit": "CAG",
        "location": "exon 1 (coding)",
        "normal": "6-26",
        "pathogenic": ">=40",
        "protein_effect": "polyglutamine expansion",
    },
    "FXN": {
        "repeat_unit": "GAA",
        "location": "intron 1 (non-coding)",
        "normal": "5-33",
        "pathogenic": "66-1700",
        "protein_effect": "transcriptional silencing (loss of frataxin)",
    },
    "DMPK": {
        "repeat_unit": "CTG",
        "location": "3'UTR (non-coding)",
        "normal": "5-34",
        "pathogenic": "50-3000+",
        "protein_effect": "RNA gain-of-function (MBNL1 sequestration)",
    },
}


# ==========================================================================
# FDA-approved gene therapies referenced in this database
# ==========================================================================

APPROVED_GENE_THERAPIES = {
    "Luxturna": {
        "generic_name": "voretigene neparvovec-rzyl",
        "target_gene": "RPE65",
        "disease": "RPE65-associated retinal dystrophy (LCA2/RP)",
        "vector": "AAV2",
        "route": "subretinal",
        "approval_year": 2017,
        "manufacturer": "Spark/Novartis",
        "mechanism": "Gene replacement -- delivers functional RPE65",
        "note": "First gene therapy for a genetic disease (US)",
    },
    "Zolgensma": {
        "generic_name": "onasemnogene abeparvovec-xioi",
        "target_gene": "SMN1",
        "disease": "Spinal muscular atrophy (SMA)",
        "vector": "AAV9",
        "route": "IV",
        "approval_year": 2019,
        "manufacturer": "Novartis",
        "mechanism": "Gene replacement -- delivers functional SMN1",
        "note": "One-time IV infusion, most expensive drug at launch (~$2.1M)",
    },
    "Elevidys": {
        "generic_name": "delandistrogene moxeparvovec-rokl",
        "target_gene": "DMD (micro-dystrophin)",
        "disease": "Duchenne muscular dystrophy",
        "vector": "AAVrh74",
        "route": "IV",
        "approval_year": 2023,
        "manufacturer": "Sarepta",
        "mechanism": "Micro-dystrophin gene replacement (truncated but "
                     "functional dystrophin)",
        "note": "Accelerated approval for ambulatory DMD patients 4-5 yrs. "
                "Safety signal in non-ambulatory patients (2025).",
    },
    "Vyjuvek": {
        "generic_name": "beremagene geperpavec-svdt",
        "target_gene": "COL7A1",
        "disease": "Dystrophic epidermolysis bullosa (DEB)",
        "vector": "HSV-1 (replication-deficient)",
        "route": "topical",
        "approval_year": 2023,
        "manufacturer": "Krystal Biotech",
        "mechanism": "Gene delivery -- HSV-1 vector delivers COL7A1 to "
                     "wound keratinocytes/fibroblasts",
        "note": "First topical gene therapy. First redosable gene therapy. "
                "Label expanded to from-birth Sept 2025.",
    },
    "Qalsody": {
        "generic_name": "tofersen",
        "target_gene": "SOD1",
        "disease": "SOD1-ALS",
        "vector": "N/A (antisense oligonucleotide)",
        "route": "intrathecal",
        "approval_year": 2023,
        "manufacturer": "Biogen",
        "mechanism": "ASO -- degrades SOD1 mRNA via RNase H",
        "note": "Accelerated approval. First treatment targeting a genetic "
                "cause of ALS. Not gene therapy per se but RNA therapeutic.",
    },
}


# ==========================================================================
# Helper functions
# ==========================================================================

def get_target(gene_name: str) -> dict | None:
    """Look up a gene target by name from the unified dictionary."""
    return ALL_HEREDITARY_TARGETS.get(gene_name)


def get_targets_by_condition(condition: str) -> dict[str, dict]:
    """Return all targets matching a given condition string."""
    return {
        gene: info
        for gene, info in ALL_HEREDITARY_TARGETS.items()
        if condition in info.get("conditions", [])
    }


def get_targets_by_category(category: str) -> dict:
    """Return the appropriate sub-dictionary by disease category."""
    categories = {
        "neurodegenerative": NEURODEGENERATIVE_TARGETS,
        "muscular_dystrophy": MUSCULAR_DYSTROPHY_TARGETS,
        "hereditary_blindness": HEREDITARY_BLINDNESS_TARGETS,
        "skin_disease": SKIN_DISEASE_TARGETS,
        "systemic_hereditary": SYSTEMIC_HEREDITARY_TARGETS,
    }
    return categories.get(category, {})


def get_repeat_expansion_targets() -> dict[str, dict]:
    """Return only trinucleotide repeat expansion disease targets."""
    return {
        gene: info
        for gene, info in ALL_HEREDITARY_TARGETS.items()
        if info.get("mutation_type", "").startswith("trinucleotide")
    }


def get_approved_therapy_genes() -> list[str]:
    """Return gene names that have FDA-approved gene therapies."""
    return [t["target_gene"] for t in APPROVED_GENE_THERAPIES.values()]
