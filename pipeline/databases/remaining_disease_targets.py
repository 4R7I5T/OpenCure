"""
Remaining disease gene-therapy / CRISPR targets -- GRCh38 (hg38) coordinates.

Covers mitochondrial diseases, complement deficiencies, coagulation/thrombotic
disorders, pain disorders, epilepsy (beyond SCN1A), and congenital heart
defects NOT already present in the existing OpenCure databases.

Categories:
  1. Mitochondrial Diseases
     - POLG, SURF1, NDUFV1, NDUFS1, SDHA, TWNK, DGUOK, TK2, SCO2
  2. Complement Deficiencies
     - C1QA, C2, C3, CFH, CFI, CD46
  3. Coagulation / Thrombotic Disorders
     - VWF, FGA, SERPINC1, PROC, PROS1, F5
  4. Pain Disorders
     - SCN9A, SCN11A, NTRK1
  5. Epilepsy Panel (beyond SCN1A)
     - SCN2A, KCNQ2, CDKL5, STXBP1, SCN8A, KCNT1
  6. Congenital Heart Defects
     - NKX2-5, TBX5, GATA4, MYH6

Coordinates are from NCBI Gene, GRCh38.p14 (RefSeq annotation RS_2025_08).
Variant positions from NCBI dbSNP / ClinVar.

IMPORTANT -- verify every coordinate against current NCBI / Ensembl releases
before production use.  Numbering can shift between patch levels.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14)
  - NCBI Datasets API (accessed April 2026)
  - OMIM
  - ClinVar / dbSNP
  - ClinicalTrials.gov
  - CRISPR Medicine News
  - Innovative Genomics Institute clinical-trial tracker
  - Published literature through early 2026
"""


# ============================================================================
# 1. MITOCHONDRIAL DISEASES
# ============================================================================

MITOCHONDRIAL_TARGETS = {

    # -------------------------------------------------------------------
    # 1a.  POLG -- mitochondrial DNA polymerase gamma, catalytic subunit
    # -------------------------------------------------------------------
    "POLG": {
        "gene_id": 5428,
        "chrom": "chr15",
        "start": 89_316_320,
        "end": 89_334_824,
        "strand": "-",
        "refseq": "NC_000015.10",
        "cytoband": "15q26.1",
        "exon_count": 23,
        "role": (
            "DNA polymerase gamma, catalytic subunit -- sole polymerase "
            "responsible for mitochondrial DNA replication and repair.  "
            "Forms a heterotrimer with POLG2 accessory subunit homodimer.  "
            "Most common nuclear gene cause of mitochondrial disease.  "
            "Loss-of-function causes progressive external ophthalmoplegia "
            "(PEO), Alpers-Huttenlocher syndrome, SANDO, myoclonic epilepsy "
            "myopathy sensory ataxia (MEMSA), and mitochondrial DNA "
            "depletion syndrome 4A/4B."
        ),
        "disease": "POLG-related mitochondrial disorders",
        "omim_gene": 174763,
        "inheritance": "AR and AD (PEO)",
        "key_variants": [
            {
                "name": "p.Ala467Thr (A467T)",
                "hgvs_coding": "NM_002693.3:c.1399G>A",
                "rsid": "rs113994095",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "frequency": "Most common POLG mutation (~31% of mutant alleles)",
                "notes": (
                    "Linker domain; reduces polymerase activity to 4-20% of "
                    "wild-type and impairs interaction with POLG2 accessory "
                    "subunit.  Homozygosity or compound heterozygosity causes "
                    "Alpers syndrome, SANDO, or myoclonic epilepsy."
                ),
            },
            {
                "name": "p.Trp748Ser (W748S)",
                "hgvs_coding": "NM_002693.3:c.2243G>C",
                "rsid": "rs113994096",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Polymerase domain; second most common European POLG "
                    "mutation.  Typically found in cis with E1143G.  "
                    "Associated with ataxia-neuropathy spectrum."
                ),
            },
            {
                "name": "p.Gly848Ser (G848S)",
                "hgvs_coding": "NM_002693.3:c.2542G>A",
                "rsid": "rs113994097",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Polymerase domain; compound het with A467T in Alpers.",
            },
        ],
        "strategy": (
            "Gene replacement therapy via AAV9 delivering POLG cDNA (~3.7 kb, "
            "fits AAV capacity) to affected tissues (CNS, muscle, liver).  "
            "Challenges: achieving mitochondrial import of functional protein "
            "and correcting mtDNA depletion/deletions already present.  "
            "Alternative approaches: (1) base editing for the common A467T "
            "(G>A, ideal ABE target); (2) allele-specific silencing for "
            "dominant PEO mutations; (3) nucleoside bypass therapy "
            "(deoxynucleoside supplementation) as a non-genetic approach.  "
            "No CRISPR clinical trials for POLG as of 2026; all preclinical."
        ),
        "clinical_programs": (
            "No gene therapy clinical trials.  Supportive management only.  "
            "Preclinical AAV-POLG studies in mouse models show partial rescue "
            "of mtDNA depletion."
        ),
        "conditions": [
            "mitochondrial_disease", "POLG", "alpers_syndrome", "PEO",
            "SANDO", "mtDNA_depletion",
        ],
    },

    # -------------------------------------------------------------------
    # 1b.  SURF1 -- cytochrome c oxidase assembly factor
    # -------------------------------------------------------------------
    "SURF1": {
        "gene_id": 6834,
        "chrom": "chr9",
        "start": 133_351_758,
        "end": 133_356_487,
        "strand": "-",
        "refseq": "NC_000009.12",
        "cytoband": "9q34.2",
        "exon_count": 9,
        "role": (
            "Cytochrome c oxidase (Complex IV) assembly factor -- inner "
            "mitochondrial membrane protein required for proper assembly of "
            "the COX holoenzyme.  Loss-of-function causes isolated COX "
            "deficiency leading to Leigh syndrome, the most common "
            "pediatric mitochondrial disease presentation."
        ),
        "disease": "Leigh syndrome (mitochondrial Complex IV deficiency)",
        "omim_gene": 185620,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.312_321del10insAT",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common SURF1 mutation in European populations.  "
                    "Leads to truncated, non-functional protein."
                ),
            },
            {
                "name": "p.Arg192Trp (R192W)",
                "hgvs_coding": "NM_003172.4:c.574C>T",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Associated with attenuated Leigh syndrome phenotype.",
            },
            {
                "name": "c.588+1G>A",
                "rsid": None,
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "notes": "Splice site mutation leading to exon skipping.",
            },
        ],
        "strategy": (
            "AAV9-mediated gene replacement (cDNA ~1.0 kb, easily fits AAV).  "
            "Preclinical studies: scAAV9-hSURF1 delivered intrathecally in "
            "Surf1-knockout mice shows partial rescue of COX assembly and "
            "improved motor function.  Challenges: need to target CNS broadly "
            "(Leigh syndrome is a CNS disease) and achieve mitochondrial "
            "localization.  No clinical trials as of 2026."
        ),
        "clinical_programs": (
            "Preclinical only.  AAV9 gene replacement in Surf1-KO mice "
            "published 2020-2022.  IND-enabling studies not yet reported."
        ),
        "conditions": [
            "mitochondrial_disease", "leigh_syndrome", "COX_deficiency",
            "complex_IV_deficiency",
        ],
    },

    # -------------------------------------------------------------------
    # 1c.  NDUFV1 -- Complex I subunit (NADH:ubiquinone oxidoreductase)
    # -------------------------------------------------------------------
    "NDUFV1": {
        "gene_id": 4723,
        "chrom": "chr11",
        "start": 67_606_936,
        "end": 67_612_554,
        "strand": "+",
        "refseq": "NC_000011.10",
        "cytoband": "11q13.2",
        "exon_count": 10,
        "role": (
            "NADH:ubiquinone oxidoreductase core subunit V1 -- 51-kDa "
            "subunit of mitochondrial Complex I that contains the NADH "
            "binding site and FMN prosthetic group.  Catalyses the first "
            "step of the mitochondrial electron transport chain.  Biallelic "
            "mutations cause mitochondrial Complex I deficiency with Leigh "
            "syndrome or leukodystrophy."
        ),
        "disease": "Mitochondrial Complex I deficiency",
        "omim_gene": 161015,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Ala341Val (A341V)",
                "hgvs_coding": "NM_007103.4:c.1022C>T",
                "rsid": "rs137852789",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts NADH-binding pocket; Leigh syndrome.",
            },
            {
                "name": "p.Thr423Met (T423M)",
                "hgvs_coding": "NM_007103.4:c.1268C>T",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "FMN-binding domain.",
            },
        ],
        "strategy": (
            "AAV-mediated gene replacement (cDNA ~1.5 kb, fits AAV).  "
            "Target CNS and skeletal muscle.  Alternatively, base editing "
            "for recurrent missense variants.  Preclinical only."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "mitochondrial_disease", "complex_I_deficiency",
            "leigh_syndrome", "leukodystrophy",
        ],
    },

    # -------------------------------------------------------------------
    # 1d.  NDUFS1 -- Complex I subunit (75 kDa)
    # -------------------------------------------------------------------
    "NDUFS1": {
        "gene_id": 4719,
        "chrom": "chr2",
        "start": 206_114_817,
        "end": 206_159_444,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q33.3",
        "exon_count": 20,
        "role": (
            "NADH:ubiquinone oxidoreductase 75-kDa Fe-S subunit -- largest "
            "subunit of Complex I containing multiple iron-sulfur clusters "
            "essential for electron transfer.  Biallelic mutations cause "
            "severe mitochondrial Complex I deficiency with leukoencephalopathy "
            "or fatal infantile lactic acidosis."
        ),
        "disease": "Mitochondrial Complex I deficiency",
        "omim_gene": 157655,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg557Ter (R557X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Truncating variant; severe phenotype.",
            },
            {
                "name": "p.Asp252Gly (D252G)",
                "rsid": "rs267607203",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Iron-sulfur cluster domain.",
            },
        ],
        "strategy": (
            "AAV gene replacement (cDNA ~2.2 kb, fits AAV).  CNS-targeted "
            "delivery required.  Challenges include achieving sufficient "
            "transduction of affected neurons and glia.  Preclinical only."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "mitochondrial_disease", "complex_I_deficiency",
            "leukoencephalopathy", "lactic_acidosis",
        ],
    },

    # -------------------------------------------------------------------
    # 1e.  SDHA -- Complex II (succinate dehydrogenase)
    # -------------------------------------------------------------------
    "SDHA": {
        "gene_id": 6389,
        "chrom": "chr5",
        "start": 218_320,
        "end": 268_746,
        "strand": "+",
        "refseq": "NC_000005.10",
        "cytoband": "5p15.33",
        "exon_count": 15,
        "role": (
            "Succinate dehydrogenase complex flavoprotein subunit A -- "
            "catalytic subunit of mitochondrial Complex II (succinate: "
            "ubiquinone oxidoreductase), the only enzyme complex that "
            "participates in both the TCA cycle and the electron transport "
            "chain.  Biallelic loss-of-function causes Leigh syndrome or "
            "mitochondrial Complex II deficiency.  Heterozygous variants "
            "also predispose to paraganglioma-pheochromocytoma syndrome."
        ),
        "disease": "Mitochondrial Complex II deficiency / Leigh syndrome",
        "omim_gene": 600857,
        "inheritance": "AR (mito disease); AD (paraganglioma/pheo)",
        "key_variants": [
            {
                "name": "p.Arg554Trp (R554W)",
                "rsid": "rs121908998",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Late-onset neurodegenerative phenotype.",
            },
            {
                "name": "p.Ala524Val (A524V)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Leigh syndrome in compound heterozygotes.",
            },
        ],
        "strategy": (
            "AAV gene replacement (cDNA ~2.0 kb, fits AAV).  Must target "
            "CNS for Leigh syndrome.  Alternative: base editing for "
            "recurrent missense mutations.  Preclinical only."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "mitochondrial_disease", "complex_II_deficiency",
            "leigh_syndrome", "paraganglioma",
        ],
    },

    # -------------------------------------------------------------------
    # 1f.  TWNK -- mitochondrial helicase (formerly PEO1/Twinkle)
    # -------------------------------------------------------------------
    "TWNK": {
        "gene_id": 56652,
        "chrom": "chr10",
        "start": 100_987_543,
        "end": 100_994_403,
        "strand": "+",
        "refseq": "NC_000010.11",
        "cytoband": "10q24.31",
        "exon_count": 5,
        "role": (
            "Twinkle mtDNA helicase -- hexameric helicase essential for "
            "mitochondrial DNA replication; unwinds mtDNA at the "
            "replication fork.  Dominant mutations cause progressive "
            "external ophthalmoplegia with mtDNA deletions (adPEO).  "
            "Recessive mutations cause infantile-onset spinocerebellar "
            "ataxia (IOSCA) or hepatocerebral mtDNA depletion syndrome."
        ),
        "disease": "PEO / IOSCA / mtDNA depletion syndrome 7",
        "omim_gene": 606075,
        "inheritance": "AD (PEO); AR (IOSCA, depletion)",
        "key_variants": [
            {
                "name": "p.Arg374Gln (R374Q)",
                "rsid": "rs137852562",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common PEO-causing TWNK variant; dominant.",
            },
            {
                "name": "p.Tyr508Cys (Y508C)",
                "rsid": "rs137852563",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Finnish IOSCA founder mutation; recessive.",
            },
        ],
        "strategy": (
            "AAV gene replacement (cDNA ~2.0 kb, fits AAV).  For dominant "
            "PEO: allele-specific silencing of the dominant-negative allele "
            "combined with replacement.  For recessive IOSCA: straightforward "
            "gene supplementation.  Preclinical only."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "mitochondrial_disease", "PEO", "IOSCA",
            "mtDNA_depletion", "spinocerebellar_ataxia",
        ],
    },

    # -------------------------------------------------------------------
    # 1g.  DGUOK -- deoxyguanosine kinase
    # -------------------------------------------------------------------
    "DGUOK": {
        "gene_id": 1716,
        "chrom": "chr2",
        "start": 73_926_880,
        "end": 73_958_946,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2p13.1",
        "exon_count": 7,
        "role": (
            "Deoxyguanosine kinase -- mitochondrial enzyme that "
            "phosphorylates purine deoxyribonucleosides (dGuo and dAdo) "
            "for mitochondrial dNTP pool maintenance.  Biallelic loss-of-"
            "function causes hepatocerebral mtDNA depletion syndrome "
            "(MTDPS3), typically presenting with neonatal liver failure "
            "and progressive neurological decline."
        ),
        "disease": "Hepatocerebral mtDNA depletion syndrome (MTDPS3)",
        "omim_gene": 601465,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.444-62C>A (splice)",
                "rsid": None,
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": "Deep intronic variant causing pseudo-exon inclusion.",
            },
            {
                "name": "p.Arg142Lys (R142K)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Active site region; severe hepatocerebral phenotype.",
            },
        ],
        "strategy": (
            "AAV gene replacement (cDNA ~0.9 kb, small, fits AAV easily).  "
            "Liver-targeted AAV8 or AAVrh10 for hepatic phenotype; CNS "
            "targeting for neurological component.  Nucleoside bypass "
            "therapy (deoxyguanosine/deoxyadenosine supplementation) as "
            "alternative non-genetic approach.  Preclinical only."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "mitochondrial_disease", "mtDNA_depletion",
            "hepatocerebral_syndrome", "liver_failure",
        ],
    },

    # -------------------------------------------------------------------
    # 1h.  TK2 -- thymidine kinase 2
    # -------------------------------------------------------------------
    "TK2": {
        "gene_id": 7084,
        "chrom": "chr16",
        "start": 66_508_003,
        "end": 66_550_291,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16q21",
        "exon_count": 10,
        "role": (
            "Thymidine kinase 2 -- mitochondrial enzyme that phosphorylates "
            "pyrimidine deoxyribonucleosides (thymidine and deoxycytidine) "
            "to maintain the mitochondrial dNTP pool.  Biallelic mutations "
            "cause myopathic mtDNA depletion syndrome (MTDPS2) with "
            "progressive myopathy, often fatal in infancy or childhood."
        ),
        "disease": "Myopathic mtDNA depletion syndrome (MTDPS2)",
        "omim_gene": 188250,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg183Gly (R183G)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Active site mutation; severe infantile myopathy.",
            },
            {
                "name": "p.Thr108Met (T108M)",
                "rsid": "rs121434584",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Milder, later-onset SMA-like presentation.",
            },
        ],
        "strategy": (
            "AAV gene replacement (cDNA ~0.8 kb, small, fits AAV easily); "
            "muscle-targeted AAV serotypes (AAV9, AAVrh74).  Nucleoside "
            "bypass therapy with deoxythymidine/deoxycytidine is already "
            "in compassionate use and clinical trials (FDA granted IND to "
            "Modis Therapeutics / Zogenix for oral dThd+dCtd).  Gene "
            "therapy remains preclinical."
        ),
        "clinical_programs": (
            "Nucleoside bypass therapy: FDA expanded access.  "
            "Gene therapy: preclinical only."
        ),
        "conditions": [
            "mitochondrial_disease", "mtDNA_depletion",
            "myopathy", "TK2_deficiency",
        ],
    },

    # -------------------------------------------------------------------
    # 1i.  SCO2 -- cytochrome c oxidase assembly factor
    # -------------------------------------------------------------------
    "SCO2": {
        "gene_id": 9997,
        "chrom": "chr22",
        "start": 50_523_568,
        "end": 50_526_442,
        "strand": "-",
        "refseq": "NC_000022.11",
        "cytoband": "22q13.33",
        "exon_count": 2,
        "role": (
            "SCO2, cytochrome c oxidase assembly factor -- copper chaperone "
            "that delivers copper to the CuA site of COX subunit II.  "
            "Essential for Complex IV assembly.  Biallelic mutations cause "
            "fatal infantile cardioencephalomyopathy with severe COX "
            "deficiency, typically presenting with hypertrophic "
            "cardiomyopathy and respiratory failure in the first weeks "
            "of life."
        ),
        "disease": "Fatal infantile cardioencephalomyopathy (COX deficiency)",
        "omim_gene": 604272,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Glu140Lys (E140K)",
                "hgvs_coding": "NM_005138.3:c.418G>A",
                "rsid": "rs74315294",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "frequency": "Present on ~80% of pathogenic SCO2 alleles",
                "notes": (
                    "Most common SCO2 mutation; disrupts copper binding.  "
                    "Found on virtually all pathogenic alleles (in cis or "
                    "trans with a more severe allele)."
                ),
            },
        ],
        "strategy": (
            "AAV gene replacement (cDNA ~0.8 kb, small, fits AAV).  "
            "Cardiac-targeted AAV9 plus CNS delivery.  Challenge: extremely "
            "early onset (neonatal) limits window for intervention -- "
            "prenatal or immediate postnatal delivery may be required.  "
            "Copper supplementation as adjunctive therapy.  Preclinical only."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "mitochondrial_disease", "COX_deficiency",
            "complex_IV_deficiency", "cardiomyopathy",
            "fatal_infantile_cardioencephalomyopathy",
        ],
    },
}


# ============================================================================
# 2. COMPLEMENT DEFICIENCIES
# ============================================================================

COMPLEMENT_TARGETS = {

    # -------------------------------------------------------------------
    # 2a.  C1QA -- complement C1q A chain
    # -------------------------------------------------------------------
    "C1QA": {
        "gene_id": 712,
        "chrom": "chr1",
        "start": 22_636_463,
        "end": 22_639_678,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p36.12",
        "exon_count": 2,
        "role": (
            "Complement C1q A chain -- subcomponent of the C1 complex that "
            "initiates the classical complement pathway.  C1q binds IgG/IgM "
            "immune complexes and apoptotic cells.  Homozygous C1q deficiency "
            "is the strongest single-gene risk factor for systemic lupus "
            "erythematosus (SLE), with >90% of C1q-deficient individuals "
            "developing SLE due to impaired clearance of apoptotic debris."
        ),
        "disease": "C1q deficiency / SLE",
        "omim_gene": 120550,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gly6Arg (G6R)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts collagenous domain triple helix.",
            },
            {
                "name": "p.Gly34Asp (G34D)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Gly-X-Y collagen repeat disruption.",
            },
        ],
        "strategy": (
            "Liver-directed AAV gene replacement for C1q chains (C1QA "
            "cDNA ~0.8 kb).  Requires co-expression of C1QB and C1QC for "
            "functional C1q heterotrimer assembly -- triple AAV or single "
            "polycistronic vector approach.  Alternatively, recombinant "
            "C1q protein replacement.  Preclinical only."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "complement_deficiency", "SLE", "lupus",
            "C1q_deficiency", "classical_pathway",
        ],
    },

    # -------------------------------------------------------------------
    # 2b.  C2 -- complement component 2
    # -------------------------------------------------------------------
    "C2": {
        "gene_id": 717,
        "chrom": "chr6",
        "start": 31_897_783,
        "end": 31_945_672,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6p21.33",
        "exon_count": 18,
        "role": (
            "Complement component 2 -- serine protease in the classical "
            "and lectin complement pathways.  C2 is cleaved by C1s to "
            "generate C2a (the catalytic subunit of the C3 convertase "
            "C4b2a).  C2 deficiency is the most common complement "
            "deficiency (~1:20,000 in European populations), associated "
            "with increased risk of SLE and recurrent infections with "
            "encapsulated bacteria."
        ),
        "disease": "C2 deficiency (most common complement deficiency)",
        "omim_gene": 613927,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Type I: 28-bp deletion in exon 6",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "frequency": "~90% of C2-deficient alleles in Europeans",
                "notes": (
                    "Linked to the MHC haplotype HLA-A25, B18, DR2.  "
                    "28-bp deletion creates premature stop codon."
                ),
            },
            {
                "name": "Type II: c.954G>C (selective deficiency)",
                "rsid": None,
                "consequence": "missense_affecting_secretion",
                "clinical_significance": "pathogenic",
                "notes": "C2 protein is synthesized but not secreted.",
            },
        ],
        "strategy": (
            "Liver-directed AAV gene replacement (cDNA ~2.4 kb, fits AAV).  "
            "C2 is primarily hepatocyte-synthesized, making liver-targeted "
            "AAV8/AAV5 an appropriate vector.  Alternatively, recombinant "
            "C2 protein replacement during acute infections.  Preclinical."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Managed with vaccination against "
            "encapsulated organisms and antibiotic prophylaxis."
        ),
        "conditions": [
            "complement_deficiency", "SLE", "recurrent_infections",
            "C2_deficiency", "classical_pathway",
        ],
    },

    # -------------------------------------------------------------------
    # 2c.  C3 -- central complement component
    # -------------------------------------------------------------------
    "C3": {
        "gene_id": 718,
        "chrom": "chr19",
        "start": 6_677_704,
        "end": 6_720_650,
        "strand": "-",
        "refseq": "NC_000019.10",
        "cytoband": "19p13.3",
        "exon_count": 41,
        "role": (
            "Complement component 3 -- the central molecule of the "
            "complement system, at the convergence of all three activation "
            "pathways (classical, lectin, alternative).  Cleaved to C3a "
            "(anaphylatoxin) and C3b (opsonin / C5 convertase component).  "
            "Complete C3 deficiency causes severe recurrent pyogenic "
            "infections and membranoproliferative glomerulonephritis.  "
            "Gain-of-function mutations cause atypical HUS (aHUS) and "
            "C3 glomerulopathy."
        ),
        "disease": "C3 deficiency / C3 glomerulopathy / aHUS",
        "omim_gene": 120700,
        "inheritance": "AR (deficiency); AD (aHUS, GOF)",
        "key_variants": [
            {
                "name": "p.Arg102Gly (C3 slow/fast polymorphism)",
                "rsid": "rs2230199",
                "consequence": "missense",
                "clinical_significance": "risk_factor",
                "notes": "C3F allele; risk factor for AMD and C3 glomerulopathy.",
            },
            {
                "name": "p.Ile1157Thr (aHUS-associated GOF)",
                "rsid": "rs121909582",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Gain-of-function; resistance to CFH-mediated regulation.",
            },
        ],
        "strategy": (
            "Liver-directed AAV gene replacement for C3 deficiency (cDNA "
            "~5.0 kb, near AAV limit but feasible with compact promoters).  "
            "C3 is primarily liver-synthesized (~90%).  For GOF mutations "
            "causing aHUS: (1) allele-specific silencing; (2) base editing "
            "to correct specific GOF variants; (3) anti-C5 mAb therapy "
            "(eculizumab/ravulizumab) is approved and effective.  Gene "
            "therapy for C3 deficiency remains preclinical."
        ),
        "clinical_programs": (
            "No gene therapy trials for C3 deficiency.  aHUS managed with "
            "eculizumab (Soliris) / ravulizumab (Ultomiris) anti-C5 therapy."
        ),
        "conditions": [
            "complement_deficiency", "aHUS", "C3_glomerulopathy",
            "C3_deficiency", "recurrent_infections",
        ],
    },

    # -------------------------------------------------------------------
    # 2d.  CFH -- complement factor H
    # -------------------------------------------------------------------
    "CFH": {
        "gene_id": 3075,
        "chrom": "chr1",
        "start": 196_652_043,
        "end": 196_747_504,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q31.3",
        "exon_count": 23,
        "role": (
            "Complement factor H -- major soluble regulator of the "
            "alternative complement pathway.  Binds C3b and accelerates "
            "decay of the C3 convertase (C3bBb) while serving as a "
            "cofactor for CFI-mediated C3b cleavage.  CFH mutations cause "
            "atypical hemolytic uremic syndrome (aHUS), C3 glomerulopathy, "
            "and increased risk of age-related macular degeneration (AMD).  "
            "The CFH-CFHR gene cluster on 1q31.3 is a hotspot for genomic "
            "rearrangements."
        ),
        "disease": "Atypical HUS / C3 glomerulopathy / AMD risk",
        "omim_gene": 134370,
        "inheritance": "AD (aHUS, incomplete penetrance)",
        "key_variants": [
            {
                "name": "p.Tyr402His (Y402H)",
                "rsid": "rs1061170",
                "consequence": "missense",
                "clinical_significance": "risk_factor",
                "frequency": "~35% in Europeans",
                "notes": (
                    "SCR7 domain; major AMD risk variant (OR ~2.5 per allele).  "
                    "Reduces CFH binding to CRP and heparin on Bruch membrane."
                ),
            },
            {
                "name": "p.Arg1210Cys (R1210C)",
                "rsid": "rs121913059",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "SCR20 domain; aHUS-associated.  Creates aberrant disulfide "
                    "bond with albumin, trapping CFH in circulation."
                ),
            },
            {
                "name": "CFH-CFHR3/1 deletion (deltaH3-1)",
                "rsid": None,
                "consequence": "structural_variant",
                "clinical_significance": "risk_modifier",
                "notes": (
                    "Common 84-kb deletion removing CFHR3 and CFHR1; protective "
                    "against AMD but can predispose to anti-FH autoantibody "
                    "aHUS in some contexts."
                ),
            },
        ],
        "strategy": (
            "Liver-directed AAV gene replacement (cDNA ~3.7 kb, fits AAV).  "
            "CFH is primarily hepatocyte-synthesized.  For aHUS patients "
            "with CFH mutations, liver-directed gene therapy could restore "
            "circulating CFH levels.  Anti-C5 mAb therapy "
            "(eculizumab/ravulizumab) is the current standard of care.  "
            "CRISPR approaches: (1) base editing for recurrent variants like "
            "R1210C; (2) for AMD, retinal delivery of CFH via subretinal AAV.  "
            "All gene therapy approaches remain preclinical."
        ),
        "clinical_programs": (
            "No gene therapy trials.  aHUS managed with eculizumab/ravulizumab.  "
            "Gemini Therapeutics explored recombinant CFH for AMD "
            "(discontinued 2022)."
        ),
        "conditions": [
            "complement_deficiency", "aHUS", "C3_glomerulopathy",
            "AMD", "macular_degeneration",
        ],
    },

    # -------------------------------------------------------------------
    # 2e.  CFI -- complement factor I
    # -------------------------------------------------------------------
    "CFI": {
        "gene_id": 3426,
        "chrom": "chr4",
        "start": 109_730_982,
        "end": 109_801_999,
        "strand": "-",
        "refseq": "NC_000004.12",
        "cytoband": "4q25",
        "exon_count": 13,
        "role": (
            "Complement factor I -- serine protease that cleaves and "
            "inactivates C3b (to iC3b) and C4b in the presence of cofactors "
            "(CFH, MCP/CD46, C4BP).  Essential negative regulator of "
            "complement.  CFI haploinsufficiency causes atypical HUS, C3 "
            "glomerulopathy, and increased susceptibility to infections.  "
            "CFI variants also associated with AMD risk."
        ),
        "disease": "Atypical HUS / C3 glomerulopathy / AMD",
        "omim_gene": 217030,
        "inheritance": "AD (aHUS, haploinsufficiency)",
        "key_variants": [
            {
                "name": "p.Gly261Asp (G261D)",
                "rsid": "rs141853578",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Serine protease domain; aHUS-associated.",
            },
            {
                "name": "p.Ile340Thr (I340T)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "likely_pathogenic",
                "notes": "Reduced C3b cleavage activity.",
            },
        ],
        "strategy": (
            "Liver-directed AAV gene replacement (cDNA ~1.8 kb, fits AAV).  "
            "CFI is liver-synthesized.  Ionis/Roche GT005 (AAV2-CFI) was "
            "in Phase I/II for geographic atrophy (dry AMD) via subretinal "
            "injection (NCT03846193) but program paused.  For aHUS, "
            "systemic liver-targeted approach.  Anti-C5 therapy is current "
            "standard for aHUS."
        ),
        "clinical_programs": (
            "Gyroscope Therapeutics (Novartis) GT005: Phase I/II for dry AMD "
            "(subretinal AAV2-CFI); program paused 2023.  No aHUS gene "
            "therapy trials."
        ),
        "conditions": [
            "complement_deficiency", "aHUS", "C3_glomerulopathy",
            "AMD", "macular_degeneration",
        ],
    },

    # -------------------------------------------------------------------
    # 2f.  CD46 (MCP) -- membrane cofactor protein
    # -------------------------------------------------------------------
    "CD46": {
        "gene_id": 4179,
        "chrom": "chr1",
        "start": 207_752_038,
        "end": 207_795_516,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q32.2",
        "exon_count": 14,
        "role": (
            "Membrane cofactor protein (CD46/MCP) -- type I transmembrane "
            "glycoprotein that serves as a cofactor for CFI-mediated "
            "cleavage of C3b and C4b deposited on host cell surfaces.  "
            "Expressed on all nucleated cells.  Also serves as a receptor "
            "for measles virus, HHV-6, and several bacterial pathogens.  "
            "Heterozygous loss-of-function mutations cause ~10-15% of "
            "aHUS cases, with generally better prognosis than CFH-aHUS."
        ),
        "disease": "Atypical HUS",
        "omim_gene": 120920,
        "inheritance": "AD (incomplete penetrance ~50%)",
        "key_variants": [
            {
                "name": "c.286+2T>G (splice donor exon 2-3)",
                "rsid": None,
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "notes": "Exon skipping; reduced MCP surface expression.",
            },
            {
                "name": "p.Ser206Pro (S206P)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "CCP3-4 domain; impairs C3b cofactor activity.",
            },
        ],
        "strategy": (
            "CD46 is a membrane-bound protein expressed on all nucleated "
            "cells, making systemic gene replacement challenging.  "
            "Approaches: (1) ex vivo CRISPR correction in HSPCs for "
            "kidney-targeted rescue (renal endothelium); (2) base editing "
            "for recurrent variants; (3) in vivo targeted delivery to "
            "renal endothelium.  aHUS from CD46 mutations has better "
            "prognosis and often remits spontaneously.  Anti-C5 therapy "
            "available.  All gene therapy preclinical."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "complement_deficiency", "aHUS",
            "CD46_deficiency", "MCP_deficiency",
        ],
    },
}


# ============================================================================
# 3. COAGULATION / THROMBOTIC DISORDERS
# ============================================================================

COAGULATION_TARGETS = {

    # -------------------------------------------------------------------
    # 3a.  VWF -- von Willebrand factor
    # -------------------------------------------------------------------
    "VWF": {
        "gene_id": 7450,
        "chrom": "chr12",
        "start": 5_948_877,
        "end": 6_124_670,
        "strand": "-",
        "refseq": "NC_000012.12",
        "cytoband": "12p13.31",
        "exon_count": 52,
        "role": (
            "Von Willebrand factor -- large multimeric glycoprotein "
            "essential for primary hemostasis.  Mediates platelet adhesion "
            "to exposed subendothelium at sites of vascular injury and "
            "serves as a carrier protein for coagulation factor VIII, "
            "protecting it from premature clearance.  Synthesized by "
            "endothelial cells and megakaryocytes.  VWF mutations cause "
            "von Willebrand disease (VWD), the most common inherited "
            "bleeding disorder (prevalence ~1%)."
        ),
        "disease": "Von Willebrand disease (types 1, 2, 3)",
        "omim_gene": 613160,
        "inheritance": "AD (types 1, 2A/B/M); AR (type 3)",
        "key_variants": [
            {
                "name": "p.Arg1306Trp (type 2B, Malmoe)",
                "rsid": "rs61748466",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "A1 domain gain-of-function; spontaneous platelet binding.  "
                    "Type 2B VWD with thrombocytopenia."
                ),
            },
            {
                "name": "p.Arg854Gln (type 2N, Normandy)",
                "rsid": "rs41276738",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "D' domain; impairs FVIII binding.  Mimics mild hemophilia A."
                ),
            },
            {
                "name": "p.Cys1190Arg (type 2A)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts disulfide bond critical for multimerization.",
            },
        ],
        "strategy": (
            "VWF gene replacement is extremely challenging: cDNA ~8.4 kb "
            "exceeds AAV packaging capacity, and VWF requires endothelial "
            "cell-specific post-translational processing (multimerization, "
            "Weibel-Palade body storage).  Approaches: (1) dual AAV "
            "split-intein for type 3 VWD; (2) lentiviral delivery to "
            "endothelial cells; (3) CRISPR correction of specific variants "
            "in endothelial cells or iPSC-derived endothelium; (4) base "
            "editing for recurrent missense variants.  All approaches "
            "preclinical.  Standard care: DDAVP, VWF/FVIII concentrates, "
            "emicizumab for type 2N."
        ),
        "clinical_programs": (
            "No gene therapy clinical trials.  Preclinical: dual hybrid "
            "AAV endothelial-targeted VWF expression (Bhatt et al. 2020).  "
            "CRISPR proof-of-concept in endothelial cells (Bar et al. 2025)."
        ),
        "conditions": [
            "coagulation_disorder", "von_willebrand_disease", "VWD",
            "bleeding_disorder",
        ],
    },

    # -------------------------------------------------------------------
    # 3b.  FGA -- fibrinogen alpha chain
    # -------------------------------------------------------------------
    "FGA": {
        "gene_id": 2243,
        "chrom": "chr4",
        "start": 154_583_126,
        "end": 154_590_742,
        "strand": "-",
        "refseq": "NC_000004.12",
        "cytoband": "4q31.3",
        "exon_count": 6,
        "role": (
            "Fibrinogen alpha chain -- one of three polypeptide chains "
            "(alpha, beta, gamma) that form fibrinogen, the soluble "
            "precursor to fibrin.  Thrombin cleaves fibrinopeptide A from "
            "the alpha chain to initiate polymerization.  FGA mutations "
            "cause congenital afibrinogenemia (null), hypofibrinogenemia, "
            "or dysfibrinogenemia (qualitative defect).  Also associated "
            "with hereditary renal amyloidosis (Ostertag type)."
        ),
        "disease": "Congenital afibrinogenemia / dysfibrinogenemia / "
                   "hereditary renal amyloidosis",
        "omim_gene": 134820,
        "inheritance": "AR (afibrinogenemia); AD (dysfibrinogenemia, amyloidosis)",
        "key_variants": [
            {
                "name": "p.Arg35His (fibrinogen Detroit)",
                "rsid": "rs121909597",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Thrombin cleavage site; impairs fibrinopeptide A "
                    "release.  Classic dysfibrinogenemia."
                ),
            },
            {
                "name": "p.Glu526Val (amyloidosis)",
                "rsid": "rs121909599",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Causes hereditary renal amyloidosis (Ostertag type).",
            },
        ],
        "strategy": (
            "Liver-directed AAV gene replacement for afibrinogenemia (FGA "
            "cDNA ~2.7 kb).  Fibrinogen is exclusively liver-synthesized.  "
            "Challenges: restoring the heterohexameric (alpha2-beta2-gamma2) "
            "complex may require co-expression of all three chains, or "
            "supplementation of the missing chain.  For dysfibrinogenemia: "
            "allele-specific silencing of the dominant-negative allele.  "
            "For amyloidosis: liver transplant is curative (removes source "
            "of amyloidogenic fibrinogen).  All gene therapy preclinical."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "coagulation_disorder", "afibrinogenemia",
            "dysfibrinogenemia", "bleeding_disorder", "amyloidosis",
        ],
    },

    # -------------------------------------------------------------------
    # 3c.  SERPINC1 -- antithrombin III
    # -------------------------------------------------------------------
    "SERPINC1": {
        "gene_id": 462,
        "chrom": "chr1",
        "start": 173_903_800,
        "end": 173_917_327,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1q25.1",
        "exon_count": 7,
        "role": (
            "Serpin family C member 1 (antithrombin III) -- major "
            "physiological inhibitor of thrombin (IIa), factor Xa, and "
            "other coagulation serine proteases.  Heparin dramatically "
            "accelerates AT-mediated inhibition.  Heterozygous deficiency "
            "causes hereditary thrombophilia with high risk of venous "
            "thromboembolism (VTE).  Homozygous deficiency of type I "
            "(quantitative) is embryonic lethal."
        ),
        "disease": "Hereditary antithrombin deficiency / thrombophilia",
        "omim_gene": 107300,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "Type II HBS: p.Arg47Cys (heparin-binding site)",
                "rsid": "rs121909548",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Disrupts heparin binding; heparin-resistant thrombophilia."
                ),
            },
            {
                "name": "Type II RS: p.Arg393Cys (reactive site)",
                "rsid": "rs121909550",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Reactive center loop; impairs protease inhibition.",
            },
        ],
        "strategy": (
            "Liver-directed AAV gene replacement (cDNA ~1.4 kb, fits AAV "
            "easily).  AT is liver-synthesized.  Goal: restore circulating "
            "AT to >50% of normal to prevent thrombosis.  Base editing "
            "for specific missense variants.  Current management: lifelong "
            "anticoagulation (heparin, warfarin, DOACs).  AT concentrate "
            "(Thrombate III, ATryn) available for acute situations.  "
            "Gene therapy preclinical."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "coagulation_disorder", "thrombophilia",
            "antithrombin_deficiency", "VTE",
        ],
    },

    # -------------------------------------------------------------------
    # 3d.  PROC -- protein C
    # -------------------------------------------------------------------
    "PROC": {
        "gene_id": 5624,
        "chrom": "chr2",
        "start": 127_418_427,
        "end": 127_429_242,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2q14.3",
        "exon_count": 9,
        "role": (
            "Protein C -- vitamin K-dependent serine protease zymogen.  "
            "Activated protein C (APC) inactivates factors Va and VIIIa, "
            "providing a critical negative feedback loop in coagulation.  "
            "Heterozygous protein C deficiency causes thrombophilia; "
            "homozygous deficiency causes neonatal purpura fulminans, "
            "a life-threatening condition requiring immediate protein C "
            "replacement."
        ),
        "disease": "Protein C deficiency / neonatal purpura fulminans",
        "omim_gene": 612283,
        "inheritance": "AD (heterozygous thrombophilia); AR (purpura fulminans)",
        "key_variants": [
            {
                "name": "p.Arg169Trp",
                "rsid": "rs121918150",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Type I deficiency; reduced protein C levels.",
            },
            {
                "name": "p.Arg306Cys (type II)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Catalytic domain; normal levels but reduced activity.",
            },
        ],
        "strategy": (
            "Liver-directed AAV gene replacement (cDNA ~1.4 kb, fits AAV).  "
            "Protein C is liver-synthesized and vitamin K-dependent.  "
            "Particularly compelling for homozygous neonatal purpura "
            "fulminans, where lifelong protein C concentrate infusions "
            "are required.  Gene therapy could provide sustained "
            "endogenous production.  Preclinical only."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Managed with protein C concentrate "
            "(Ceprotin) and anticoagulation."
        ),
        "conditions": [
            "coagulation_disorder", "thrombophilia",
            "protein_C_deficiency", "purpura_fulminans",
        ],
    },

    # -------------------------------------------------------------------
    # 3e.  PROS1 -- protein S
    # -------------------------------------------------------------------
    "PROS1": {
        "gene_id": 5627,
        "chrom": "chr3",
        "start": 93_873_051,
        "end": 93_973_896,
        "strand": "-",
        "refseq": "NC_000003.12",
        "cytoband": "3q11.1",
        "exon_count": 15,
        "role": (
            "Protein S -- vitamin K-dependent glycoprotein that serves as "
            "a non-enzymatic cofactor for activated protein C (APC) in the "
            "inactivation of factors Va and VIIIa.  Also has direct "
            "APC-independent anticoagulant activity via interaction with "
            "TFPI.  Synthesized primarily by hepatocytes, endothelial "
            "cells, and megakaryocytes.  Deficiency causes thrombophilia "
            "with increased VTE risk."
        ),
        "disease": "Protein S deficiency / thrombophilia",
        "omim_gene": 176880,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "p.Lys155Glu (type I)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Reduced total and free protein S levels.",
            },
            {
                "name": "Large deletions (PROS1/PROS2 pseudogene conversion)",
                "rsid": None,
                "consequence": "structural_variant",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Gene conversion with PROS2 pseudogene causes ~10% "
                    "of protein S deficiency."
                ),
            },
        ],
        "strategy": (
            "Liver-directed AAV gene replacement (cDNA ~2.1 kb, fits AAV).  "
            "Challenges: PROS1 has a nearby pseudogene (PROS2) with high "
            "homology, complicating CRISPR targeting.  Codon-optimized "
            "transgene design to avoid pseudogene cross-reactivity.  "
            "Management: anticoagulation therapy.  Gene therapy preclinical."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "coagulation_disorder", "thrombophilia",
            "protein_S_deficiency", "VTE",
        ],
    },

    # -------------------------------------------------------------------
    # 3f.  F5 -- coagulation factor V (Factor V Leiden)
    # -------------------------------------------------------------------
    "F5": {
        "gene_id": 2153,
        "chrom": "chr1",
        "start": 169_511_951,
        "end": 169_586_481,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1q24.2",
        "exon_count": 25,
        "role": (
            "Coagulation factor V -- large glycoprotein cofactor in the "
            "prothrombinase complex (FVa-FXa) that converts prothrombin "
            "to thrombin.  Also has anticoagulant function as a cofactor "
            "for APC-mediated inactivation of FVIIIa.  The Factor V "
            "Leiden mutation (p.Arg534Gln, formerly R506Q) abolishes one "
            "of three APC cleavage sites, causing activated protein C "
            "resistance -- the most common inherited thrombophilia "
            "(~5% of Europeans are heterozygous)."
        ),
        "disease": "Factor V Leiden thrombophilia / APC resistance",
        "omim_gene": 612309,
        "inheritance": "AD (incomplete penetrance; dose-dependent risk)",
        "key_variants": [
            {
                "name": "p.Arg534Gln (Factor V Leiden, R506Q legacy)",
                "hgvs_coding": "NM_000130.4:c.1601G>A",
                "hgvs_genomic": "NC_000001.11:g.169549811C>T",
                "rsid": "rs6025",
                "consequence": "missense",
                "clinical_significance": "pathogenic / risk_factor",
                "frequency": "~5% heterozygotes in Europeans",
                "notes": (
                    "Abolishes Arg534 APC cleavage site in factor Va.  "
                    "VTE risk increased 3-8x in heterozygotes, 9-80x in "
                    "homozygotes.  Incomplete penetrance: many carriers "
                    "never develop thrombosis."
                ),
            },
            {
                "name": "p.Arg334Gln (Factor V Cambridge)",
                "rsid": "rs118203906",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Second APC cleavage site; rare, APC resistance.",
            },
        ],
        "strategy": (
            "Base editing to correct the rs6025 G>A Leiden variant (ideal "
            "ABE target: restoring Gln534 back to Arg534).  Liver-directed "
            "delivery (LNP or AAV) since FV is hepatocyte-synthesized.  "
            "However, clinical need is debatable as most carriers never "
            "develop thrombosis and anticoagulation is effective.  Gene "
            "therapy would be most justified for homozygous individuals "
            "with recurrent severe VTE despite anticoagulation.  "
            "All approaches preclinical."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Managed with anticoagulation "
            "(heparin, warfarin, DOACs) during/after thrombotic events."
        ),
        "conditions": [
            "coagulation_disorder", "thrombophilia",
            "factor_V_Leiden", "APC_resistance", "VTE",
        ],
    },
}


# ============================================================================
# 4. PAIN DISORDERS
# ============================================================================

PAIN_TARGETS = {

    # -------------------------------------------------------------------
    # 4a.  SCN9A -- Nav1.7 sodium channel
    # -------------------------------------------------------------------
    "SCN9A": {
        "gene_id": 6335,
        "chrom": "chr2",
        "start": 166_195_185,
        "end": 166_375_987,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q24.3",
        "exon_count": 27,
        "role": (
            "Sodium voltage-gated channel alpha subunit 9 (Nav1.7) -- "
            "peripheral sodium channel expressed in dorsal root ganglion "
            "(DRG) nociceptors and sympathetic neurons.  Acts as a "
            "threshold channel that amplifies subthreshold stimuli to "
            "trigger action potentials.  Gain-of-function mutations cause "
            "inherited erythromelalgia (IEM) and paroxysmal extreme pain "
            "disorder (PEPD).  Loss-of-function mutations cause congenital "
            "insensitivity to pain (CIP) without other neurological "
            "deficits -- validating Nav1.7 as a pain target."
        ),
        "disease": "Erythromelalgia / PEPD / congenital insensitivity to pain",
        "omim_gene": 603415,
        "inheritance": "AD (IEM, PEPD); AR (CIP)",
        "key_variants": [
            {
                "name": "p.Ile848Thr (I848T, IEM)",
                "rsid": "rs121908915",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Gain-of-function; hyperpolarizing shift in activation.  "
                    "Primary erythromelalgia with severe burning pain."
                ),
            },
            {
                "name": "p.Trp897Ter (W897X, CIP)",
                "rsid": "rs121908916",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Loss-of-function; complete insensitivity to pain.  "
                    "Pakistani kindred, first CIP-SCN9A family described."
                ),
            },
            {
                "name": "p.Phe1449Val (F1449V, IEM)",
                "rsid": "rs121908919",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Gain-of-function; erythromelalgia.",
            },
        ],
        "strategy": (
            "Epigenetic silencing of SCN9A in DRG neurons for chronic pain "
            "management.  Navega Therapeutics LATER approach: AAV-delivered "
            "dCas9 fused to KRAB repressor targeting SCN9A promoter in DRG "
            "neurons, providing long-lasting analgesia without DNA cutting.  "
            "NT-Z001 (zinc finger-KRAB repressor variant) received $4M CIRM "
            "grant (Feb 2025) for IND-enabling studies.  Alternative: "
            "(1) base editing to correct specific GOF mutations causing "
            "IEM/PEPD; (2) small molecule Nav1.7 blockers (multiple pharma "
            "programs, most failed clinical trials due to selectivity).  "
            "Preclinical for gene therapy; no human trials yet."
        ),
        "clinical_programs": (
            "Navega Therapeutics NT-Z001: preclinical (CIRM-funded).  "
            "Multiple Nav1.7 small molecule programs in Phase I/II.  "
            "No CRISPR/gene therapy clinical trials."
        ),
        "conditions": [
            "pain_disorder", "erythromelalgia", "CIP",
            "congenital_insensitivity_to_pain", "PEPD", "neuropathic_pain",
        ],
    },

    # -------------------------------------------------------------------
    # 4b.  SCN11A -- Nav1.9 sodium channel
    # -------------------------------------------------------------------
    "SCN11A": {
        "gene_id": 11280,
        "chrom": "chr3",
        "start": 38_845_764,
        "end": 39_051_944,
        "strand": "-",
        "refseq": "NC_000003.12",
        "cytoband": "3p22.2",
        "exon_count": 27,
        "role": (
            "Sodium voltage-gated channel alpha subunit 11 (Nav1.9) -- "
            "persistent sodium current channel expressed in DRG nociceptors "
            "and myenteric neurons.  Sets the resting membrane potential "
            "of nociceptors close to action potential threshold.  "
            "Gain-of-function mutations cause familial episodic pain "
            "syndrome type 3 (FEPS3) and painful small fiber neuropathy.  "
            "Loss-of-function mutations cause congenital insensitivity "
            "to pain with GI dysmotility."
        ),
        "disease": "Familial episodic pain syndrome 3 / painful neuropathy",
        "omim_gene": 604385,
        "inheritance": "AD (GOF: episodic pain); AD (LOF: insensitivity)",
        "key_variants": [
            {
                "name": "p.Arg225Cys (R225C, FEPS3)",
                "rsid": "rs587777037",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Gain-of-function; enhances persistent Na+ current.  "
                    "Episodic pain triggered by cold or fasting."
                ),
            },
            {
                "name": "p.Leu811Pro (L811P, insensitivity to pain)",
                "rsid": "rs587777038",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Paradoxical: GOF at channel level but causes "
                    "depolarization block leading to insensitivity to pain."
                ),
            },
        ],
        "strategy": (
            "Allele-specific silencing of the GOF allele for episodic pain.  "
            "CRISPRi/dCas9-KRAB approach similar to SCN9A strategy.  "
            "Base editing for specific GOF missense variants.  Less "
            "explored than SCN9A; earlier stage preclinical."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "pain_disorder", "familial_episodic_pain",
            "small_fiber_neuropathy", "FEPS3",
        ],
    },

    # -------------------------------------------------------------------
    # 4c.  NTRK1 -- neurotrophic receptor tyrosine kinase 1
    # -------------------------------------------------------------------
    "NTRK1": {
        "gene_id": 4914,
        "chrom": "chr1",
        "start": 156_815_750,
        "end": 156_881_850,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q23.1",
        "exon_count": 17,
        "role": (
            "Neurotrophic receptor tyrosine kinase 1 (TrkA) -- high-"
            "affinity receptor for nerve growth factor (NGF).  Essential "
            "for survival and differentiation of nociceptive sensory "
            "neurons and sympathetic neurons during development.  Biallelic "
            "loss-of-function mutations cause CIPA (congenital insensitivity "
            "to pain with anhidrosis, HSAN IV) -- complete absence of "
            "pain sensation, temperature sensation, and sweating.  "
            "Associated with intellectual disability, self-mutilation, "
            "recurrent fractures, and hyperthermia."
        ),
        "disease": "CIPA (HSAN type IV)",
        "omim_gene": 191315,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg774Pro (R774P)",
                "rsid": "rs80356677",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Tyrosine kinase domain; abolishes catalytic activity.",
            },
            {
                "name": "c.851-33T>A (splice)",
                "rsid": None,
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Founder mutation in Israeli-Bedouin population.  "
                    "Creates cryptic splice site."
                ),
            },
            {
                "name": "p.Leu213Pro (L213P)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Leucine-rich repeat domain; disrupts NGF binding.",
            },
        ],
        "strategy": (
            "AAV gene replacement targeting DRG neurons (cDNA ~2.4 kb, "
            "fits AAV).  Challenges: (1) restoring pain sensation in "
            "patients who have never had it requires establishing new "
            "neural circuits; (2) therapeutic window unclear -- treatment "
            "in neonates/infants may preserve developing nociceptive "
            "neurons before irreversible loss; (3) risk of neuropathic "
            "pain from reinnervation.  Base editing for recurrent "
            "variants.  Earlier-stage than SCN9A programs."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "pain_disorder", "CIPA", "HSAN_IV",
            "congenital_insensitivity_to_pain", "anhidrosis",
        ],
    },
}


# ============================================================================
# 5. EPILEPSY PANEL (beyond SCN1A already covered elsewhere)
# ============================================================================

EPILEPSY_TARGETS = {

    # -------------------------------------------------------------------
    # 5a.  SCN2A -- Nav1.2 sodium channel
    # -------------------------------------------------------------------
    "SCN2A": {
        "gene_id": 6326,
        "chrom": "chr2",
        "start": 165_239_414,
        "end": 165_392_304,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2q24.3",
        "exon_count": 27,
        "role": (
            "Sodium voltage-gated channel alpha subunit 2 (Nav1.2) -- "
            "predominant sodium channel in excitatory neurons, especially "
            "at the axon initial segment.  Critical for action potential "
            "initiation and propagation in cortical pyramidal neurons.  "
            "Gain-of-function mutations cause early infantile epileptic "
            "encephalopathy (Ohtahara/West syndrome).  Loss-of-function "
            "mutations cause autism spectrum disorder and late-onset "
            "epilepsy.  Second most common genetic cause of developmental "
            "and epileptic encephalopathy (DEE) after SCN1A."
        ),
        "disease": "SCN2A-related DEE / autism / BFNIE",
        "omim_gene": 182390,
        "inheritance": "AD (de novo in DEE)",
        "key_variants": [
            {
                "name": "p.Arg1882Gln (R1882Q, GOF)",
                "rsid": "rs796052984",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Gain-of-function; enhanced persistent current.  "
                    "Severe neonatal epileptic encephalopathy."
                ),
            },
            {
                "name": "p.Ala263Val (A263V, GOF)",
                "rsid": "rs387906669",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Gain-of-function; neonatal/infantile seizures.",
            },
            {
                "name": "p.Arg937Cys (LOF-associated)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Loss-of-function; associated with ASD and late epilepsy.",
            },
        ],
        "strategy": (
            "Strategy depends on variant mechanism: "
            "(1) GOF mutations: antisense oligonucleotide (ASO) knockdown "
            "of SCN2A.  Praxis Precision Medicine elsunersen (PRAX-222) is "
            "a gapmer ASO in clinical trials for GOF-DEE.  First preterm "
            "infant treated (Nature Medicine 2025) with 60% seizure reduction.  "
            "(2) LOF mutations: CRISPR activation (CRISPRa) of the remaining "
            "wild-type allele.  UCSF team demonstrated CRISPRa rescue of "
            "SCN2A haploinsufficiency in mice (Nature, Sept 2025).  "
            "(3) Gene replacement: cDNA ~6.0 kb, near AAV limit.  "
            "Regel Therapeutics RT102 (preclinical) for LOF."
        ),
        "clinical_programs": (
            "Praxis elsunersen (PRAX-222): Phase I/II ASO for GOF-DEE.  "
            "Praxis relutrigene (PRAX-562): Phase II, 53% seizure reduction.  "
            "Regel RT102: preclinical gene therapy for LOF.  "
            "UCSF CRISPRa: preclinical (Nature 2025)."
        ),
        "conditions": [
            "epilepsy", "DEE", "epileptic_encephalopathy",
            "SCN2A", "autism", "BFNIE",
        ],
    },

    # -------------------------------------------------------------------
    # 5b.  KCNQ2 -- Kv7.2 potassium channel
    # -------------------------------------------------------------------
    "KCNQ2": {
        "gene_id": 3785,
        "chrom": "chr20",
        "start": 63_400_208,
        "end": 63_472_655,
        "strand": "-",
        "refseq": "NC_000020.11",
        "cytoband": "20q13.33",
        "exon_count": 17,
        "role": (
            "Potassium voltage-gated channel subfamily Q member 2 (Kv7.2) "
            "-- forms heteromeric channels with KCNQ3 (Kv7.3) that "
            "generate the M-current, a slowly activating/deactivating "
            "potassium current that stabilizes neuronal resting potential "
            "and controls excitability.  GOF and LOF mutations both cause "
            "epilepsy: benign familial neonatal epilepsy (BFNE, milder "
            "LOF) or KCNQ2 developmental and epileptic encephalopathy "
            "(KCNQ2-DEE, dominant-negative LOF or GOF)."
        ),
        "disease": "KCNQ2-DEE / benign familial neonatal epilepsy",
        "omim_gene": 602235,
        "inheritance": "AD (mostly de novo for DEE)",
        "key_variants": [
            {
                "name": "p.Ala294Val (A294V, DEE)",
                "rsid": "rs587783292",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Pore domain; dominant-negative loss of M-current.  "
                    "Severe neonatal DEE."
                ),
            },
            {
                "name": "p.Arg213Trp (R213W, DEE)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "S4 voltage-sensor domain; severe DEE.",
            },
            {
                "name": "p.Arg207Trp (R207W, GOF)",
                "rsid": "rs387906682",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Gain-of-function variant; neonatal encephalopathy.",
            },
        ],
        "strategy": (
            "Gene replacement or augmentation for LOF-KCNQ2 (cDNA ~2.6 kb, "
            "fits AAV).  For dominant-negative variants: allele-specific "
            "silencing of the mutant allele + supplementation with "
            "wild-type.  Kv7 channel openers (retigabine/ezogabine) are "
            "available but systemic side effects limit use.  XEN1101 "
            "(Xenon) is a selective Kv7.2/7.3 opener in Phase III for "
            "focal epilepsy.  Gene therapy preclinical."
        ),
        "clinical_programs": (
            "No gene therapy clinical trials.  Xenon XEN1101 Phase III "
            "(epilepsy, but broad population).  Managed with sodium channel "
            "blockers (carbamazepine) for BFNE."
        ),
        "conditions": [
            "epilepsy", "DEE", "epileptic_encephalopathy",
            "KCNQ2", "neonatal_epilepsy", "BFNE",
        ],
    },

    # -------------------------------------------------------------------
    # 5c.  CDKL5 -- cyclin-dependent kinase-like 5
    # -------------------------------------------------------------------
    "CDKL5": {
        "gene_id": 6792,
        "chrom": "chrX",
        "start": 18_425_608,
        "end": 18_653_629,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xp22.13",
        "exon_count": 24,
        "role": (
            "Cyclin-dependent kinase-like 5 -- serine/threonine kinase "
            "highly expressed in neurons.  Phosphorylates multiple "
            "substrates involved in neuronal morphogenesis, synapse "
            "formation, and microtubule dynamics.  Loss-of-function "
            "causes CDKL5 deficiency disorder (CDD), formerly classified "
            "as an atypical Rett syndrome variant.  CDD presents with "
            "early-onset refractory seizures (often within first 3 "
            "months), severe intellectual disability, cortical visual "
            "impairment, and stereotypic hand movements.  Predominantly "
            "affects females; severe in males."
        ),
        "disease": "CDKL5 deficiency disorder (CDD)",
        "omim_gene": 300203,
        "inheritance": "XLD (mostly de novo)",
        "key_variants": [
            {
                "name": "p.Arg59Ter (R59X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Kinase domain truncation; severe phenotype.",
            },
            {
                "name": "p.Arg178Pro (R178P)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Kinase domain; abolishes catalytic activity.",
            },
            {
                "name": "p.Leu220Pro (L220P)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts ATP-binding pocket.",
            },
        ],
        "strategy": (
            "AAV9-mediated gene replacement delivered intrathecally or "
            "intracerebroventricularly (cDNA ~3.2 kb, fits AAV).  "
            "Loulou Foundation / Elaaj Bio leading CDD gene therapy "
            "program: dose-finding in mouse models complete; toxicology "
            "in NHP complete.  Partnership with Gemma Bio for vector "
            "manufacturing.  Clinical trial authorization anticipated "
            "late 2026 / early 2027.  Challenges: X-inactivation in "
            "females (mosaic expression), need for broad CNS transduction, "
            "overexpression toxicity risk.  Ganaxolone (Ztalmy) approved "
            "2022 for CDD seizures but does not address underlying cause."
        ),
        "clinical_programs": (
            "Loulou Foundation / Elaaj Bio: AAV9-CDKL5 gene therapy, "
            "IND-enabling stage.  Clinical trial anticipated 2027.  "
            "Ganaxolone (Ztalmy) approved for CDD seizures (2022)."
        ),
        "conditions": [
            "epilepsy", "CDKL5_deficiency", "CDD",
            "epileptic_encephalopathy", "DEE",
        ],
    },

    # -------------------------------------------------------------------
    # 5d.  STXBP1 -- syntaxin-binding protein 1
    # -------------------------------------------------------------------
    "STXBP1": {
        "gene_id": 6812,
        "chrom": "chr9",
        "start": 127_611_912,
        "end": 127_696_029,
        "strand": "+",
        "refseq": "NC_000009.12",
        "cytoband": "9q34.11",
        "exon_count": 20,
        "role": (
            "Syntaxin-binding protein 1 (Munc18-1) -- essential component "
            "of the synaptic vesicle fusion machinery.  Binds syntaxin-1 "
            "in a closed conformation to regulate SNARE complex assembly "
            "and neurotransmitter release.  Haploinsufficiency causes "
            "STXBP1 encephalopathy (formerly Ohtahara syndrome / EIEE4), "
            "featuring early-onset refractory seizures, severe "
            "intellectual disability, and movement disorders."
        ),
        "disease": "STXBP1 encephalopathy (Ohtahara syndrome)",
        "omim_gene": 602926,
        "inheritance": "AD (de novo, haploinsufficiency)",
        "key_variants": [
            {
                "name": "p.Arg406Cys (R406C)",
                "rsid": "rs587783786",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Domain 3a; destabilizes syntaxin-1 binding.",
            },
            {
                "name": "p.Arg292His (R292H)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Recurrent de novo variant; Ohtahara syndrome.",
            },
            {
                "name": "p.Arg406His (R406H)",
                "rsid": "rs587783787",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Same residue as R406C; severe phenotype.",
            },
        ],
        "strategy": (
            "AAV9 gene replacement for haploinsufficiency (cDNA ~1.8 kb, "
            "fits AAV easily).  Restoring STXBP1 levels to >50% of normal "
            "may rescue synaptic transmission.  Challenges: ubiquitous "
            "neuronal expression required; overexpression may be toxic.  "
            "Tight promoter regulation needed (e.g., synapsin-1 promoter "
            "with miRNA-regulated elements).  Preclinical only."
        ),
        "clinical_programs": "Preclinical; no clinical trials.",
        "conditions": [
            "epilepsy", "STXBP1_encephalopathy", "ohtahara_syndrome",
            "DEE", "epileptic_encephalopathy",
        ],
    },

    # -------------------------------------------------------------------
    # 5e.  SCN8A -- Nav1.6 sodium channel
    # -------------------------------------------------------------------
    "SCN8A": {
        "gene_id": 6334,
        "chrom": "chr12",
        "start": 51_591_233,
        "end": 51_812_864,
        "strand": "+",
        "refseq": "NC_000012.12",
        "cytoband": "12q13.13",
        "exon_count": 27,
        "role": (
            "Sodium voltage-gated channel alpha subunit 8 (Nav1.6) -- "
            "the predominant sodium channel at nodes of Ranvier and axon "
            "initial segments in mature neurons.  Responsible for the "
            "resurgent sodium current.  Gain-of-function mutations cause "
            "SCN8A-related epileptic encephalopathy (EIEE13), often with "
            "SUDEP risk.  Also causes movement disorders and intellectual "
            "disability.  Loss-of-function is associated with cognitive "
            "impairment without epilepsy."
        ),
        "disease": "SCN8A-related epileptic encephalopathy (EIEE13)",
        "omim_gene": 600702,
        "inheritance": "AD (de novo GOF for DEE)",
        "key_variants": [
            {
                "name": "p.Arg1617Gln (R1617Q)",
                "rsid": "rs796052826",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Gain-of-function; enhanced persistent/resurgent current.  "
                    "Severe infantile-onset DEE with SUDEP risk."
                ),
            },
            {
                "name": "p.Asn1768Asp (N1768D)",
                "rsid": "rs397515754",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "First SCN8A-DEE mutation identified.  GOF with "
                    "increased persistent current."
                ),
            },
        ],
        "strategy": (
            "ASO-mediated knockdown of the GOF allele (analogous to "
            "SCN2A ASO approach).  Allele-specific silencing.  Encoded "
            "Therapeutics exploring engineered regulatory elements.  "
            "Precision medicine: sodium channel blockers (phenytoin, "
            "quinidine) partially effective for GOF variants.  Gene "
            "therapy preclinical."
        ),
        "clinical_programs": (
            "No gene therapy clinical trials.  ASO approaches in "
            "preclinical development.  Managed with high-dose phenytoin "
            "or quinidine."
        ),
        "conditions": [
            "epilepsy", "SCN8A", "epileptic_encephalopathy",
            "DEE", "SUDEP",
        ],
    },

    # -------------------------------------------------------------------
    # 5f.  KCNT1 -- KNa1.1 sodium-activated potassium channel
    # -------------------------------------------------------------------
    "KCNT1": {
        "gene_id": 57582,
        "chrom": "chr9",
        "start": 135_702_185,
        "end": 135_795_502,
        "strand": "+",
        "refseq": "NC_000009.12",
        "cytoband": "9q34.3",
        "exon_count": 31,
        "role": (
            "Potassium sodium-activated channel subfamily T member 1 "
            "(KNa1.1, Slack) -- large-conductance sodium-activated "
            "potassium channel expressed in cortical neurons.  Modulates "
            "neuronal firing patterns and adaptation.  Gain-of-function "
            "mutations cause epilepsy of infancy with migrating focal "
            "seizures (EIMFS/EIEE14) and autosomal dominant nocturnal "
            "frontal lobe epilepsy (ADNFLE).  GOF increases K+ current "
            "but paradoxically increases neuronal excitability through "
            "altered firing patterns."
        ),
        "disease": "EIMFS / ADNFLE / KCNT1-related epilepsy",
        "omim_gene": 608167,
        "inheritance": "AD (de novo for EIMFS; inherited for ADNFLE)",
        "key_variants": [
            {
                "name": "p.Arg474His (R474H)",
                "rsid": "rs587784235",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "RCK1 domain; recurrent GOF mutation in EIMFS.",
            },
            {
                "name": "p.Arg928Cys (R928C)",
                "rsid": "rs587784236",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "RCK2 domain; GOF, ADNFLE and EIMFS.",
            },
            {
                "name": "p.Gly288Ser (G288S)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "S5-S6 pore region; severe EIMFS.",
            },
        ],
        "strategy": (
            "Quinidine (off-label) as a KCNT1 channel blocker has shown "
            "variable efficacy.  Gene therapy approaches: (1) ASO-mediated "
            "knockdown of the GOF allele; (2) allele-specific CRISPR "
            "silencing; (3) dominant-negative peptide expression.  "
            "Praxis Precision Medicine was developing KCNT1-targeted "
            "therapeutics.  All gene therapy preclinical."
        ),
        "clinical_programs": (
            "No gene therapy clinical trials.  Off-label quinidine with "
            "variable results.  Praxis pipeline (early stage)."
        ),
        "conditions": [
            "epilepsy", "EIMFS", "ADNFLE", "KCNT1",
            "epileptic_encephalopathy", "DEE",
        ],
    },
}


# ============================================================================
# 6. CONGENITAL HEART DEFECTS
# ============================================================================

CONGENITAL_HEART_TARGETS = {

    # -------------------------------------------------------------------
    # 6a.  NKX2-5 -- cardiac homeobox transcription factor
    # -------------------------------------------------------------------
    "NKX2-5": {
        "gene_id": 1482,
        "chrom": "chr5",
        "start": 173_232_109,
        "end": 173_235_206,
        "strand": "-",
        "refseq": "NC_000005.10",
        "cytoband": "5q35.1",
        "exon_count": 2,
        "role": (
            "NK2 homeobox 5 -- master cardiac transcription factor "
            "essential for heart morphogenesis, septation, and conduction "
            "system development.  One of the earliest markers of cardiac "
            "progenitor cells.  Heterozygous loss-of-function mutations "
            "cause atrial septal defects (ASD), ventricular septal defects "
            "(VSD), atrioventricular conduction block, and tetralogy of "
            "Fallot.  Also associated with non-goitrous congenital "
            "hypothyroidism type 5."
        ),
        "disease": "Congenital heart defects (ASD, VSD, conduction defects)",
        "omim_gene": 600584,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "p.Gln170Ter (Q170X)",
                "rsid": "rs104894398",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Homeodomain truncation; ASD + progressive AV block.  "
                    "Original kindred (Schott et al. 1998, Science)."
                ),
            },
            {
                "name": "p.Arg161Pro (R161P)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Homeodomain; disrupts DNA binding.",
            },
            {
                "name": "p.Ala119Ser (A119S)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Associated with tetralogy of Fallot.",
            },
        ],
        "strategy": (
            "Congenital heart defects are structural malformations that "
            "occur during embryonic development, making postnatal gene "
            "therapy unable to reverse the morphological defect.  However, "
            "progressive conduction disease (AV block) from NKX2-5 "
            "mutations may benefit from gene augmentation in the "
            "conduction system.  Approaches: (1) AAV9-NKX2-5 gene "
            "supplementation (cDNA ~1.0 kb, small) targeted to cardiac "
            "conduction system to slow progressive AV block; (2) CRISPR "
            "correction in iPSC-derived cardiomyocytes for research; "
            "(3) preimplantation genetic testing to prevent transmission.  "
            "Surgical repair remains standard for structural defects."
        ),
        "clinical_programs": (
            "No gene therapy trials.  CHD managed surgically.  "
            "Pacemaker for progressive AV block."
        ),
        "conditions": [
            "congenital_heart_defect", "ASD", "VSD",
            "AV_block", "tetralogy_of_fallot",
        ],
    },

    # -------------------------------------------------------------------
    # 6b.  TBX5 -- T-box transcription factor 5
    # -------------------------------------------------------------------
    "TBX5": {
        "gene_id": 6910,
        "chrom": "chr12",
        "start": 114_353_911,
        "end": 114_408_442,
        "strand": "-",
        "refseq": "NC_000012.12",
        "cytoband": "12q24.21",
        "exon_count": 9,
        "role": (
            "T-box transcription factor 5 -- essential for upper limb and "
            "cardiac development.  Interacts with NKX2-5, GATA4, and TBX20 "
            "to regulate cardiac gene expression.  Heterozygous mutations "
            "cause Holt-Oram syndrome (heart-hand syndrome): congenital "
            "heart defects (ASD, VSD) plus upper limb malformations "
            "(especially radial ray defects: thumb anomalies, absent/short "
            "radius).  Also causes isolated cardiac conduction disease."
        ),
        "disease": "Holt-Oram syndrome (heart-hand syndrome)",
        "omim_gene": 601620,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "p.Gly80Arg (G80R)",
                "rsid": "rs104894737",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "T-box DNA-binding domain; severe cardiac and limb defects.",
            },
            {
                "name": "p.Arg237Trp (R237W)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "T-box domain; associated with ASD + conduction disease.",
            },
            {
                "name": "p.Arg237Gln (R237Q)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Reduces DNA-binding affinity.",
            },
        ],
        "strategy": (
            "Structural defects formed during embryonic development cannot "
            "be reversed by postnatal gene therapy.  Potential applications: "
            "(1) gene supplementation for progressive conduction disease; "
            "(2) prenatal gene therapy (experimental concept stage only); "
            "(3) CRISPR correction in iPSC cardiomyocytes for research and "
            "drug screening.  Preimplantation genetic testing available "
            "for affected families.  Surgical repair for heart and limb "
            "defects remains standard of care."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Managed with cardiac surgery and "
            "orthopedic intervention."
        ),
        "conditions": [
            "congenital_heart_defect", "holt_oram_syndrome",
            "heart_hand_syndrome", "ASD", "limb_malformation",
        ],
    },

    # -------------------------------------------------------------------
    # 6c.  GATA4 -- GATA-binding protein 4
    # -------------------------------------------------------------------
    "GATA4": {
        "gene_id": 2626,
        "chrom": "chr8",
        "start": 11_676_935,
        "end": 11_760_002,
        "strand": "+",
        "refseq": "NC_000008.11",
        "cytoband": "8p23.1",
        "exon_count": 7,
        "role": (
            "GATA-binding protein 4 -- zinc-finger transcription factor "
            "essential for cardiac development, cardiomyocyte "
            "differentiation, and cardiac gene regulation.  Binds the "
            "consensus sequence WGATAR in cardiac promoters/enhancers "
            "and interacts with NKX2-5, TBX5, and FOG2.  Heterozygous "
            "mutations cause isolated congenital heart defects, primarily "
            "atrial septal defects (ASD, secundum type), ventricular "
            "septal defects, and pulmonary stenosis.  Also implicated "
            "in tetralogy of Fallot and atrioventricular canal defects."
        ),
        "disease": "Congenital heart defects (ASD, VSD, pulmonary stenosis)",
        "omim_gene": 600576,
        "inheritance": "AD (haploinsufficiency/dominant-negative)",
        "key_variants": [
            {
                "name": "p.Gly296Ser (G296S)",
                "rsid": "rs104894090",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Zinc finger 2 domain; disrupts GATA4-TBX5 "
                    "interaction.  Familial ASD."
                ),
            },
            {
                "name": "p.Glu359del (deletion)",
                "rsid": None,
                "consequence": "in-frame_deletion",
                "clinical_significance": "pathogenic",
                "notes": "Reduces transcriptional activation.",
            },
            {
                "name": "p.Ala411Val (A411V)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "C-terminal activation domain; ASD.",
            },
        ],
        "strategy": (
            "Similar to NKX2-5 and TBX5: structural heart defects "
            "established during embryogenesis cannot be reversed "
            "postnatally.  Research applications: (1) CRISPR correction "
            "in hiPSC-derived cardiomyocytes for disease modeling (widely "
            "published); (2) understanding GATA4-NKX2-5-TBX5 regulatory "
            "network for potential regenerative medicine approaches; "
            "(3) preimplantation genetic testing for familial cases.  "
            "AAV9-mediated somatic CRISPR has been used to correct "
            "cardiac gene mutations in neonatal mice (proof-of-concept "
            "stage).  Surgical repair is standard of care."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Extensive hiPSC-CRISPR research "
            "for disease modeling.  Managed with cardiac surgery."
        ),
        "conditions": [
            "congenital_heart_defect", "ASD", "VSD",
            "pulmonary_stenosis", "tetralogy_of_fallot",
        ],
    },

    # -------------------------------------------------------------------
    # 6d.  MYH6 -- alpha-myosin heavy chain
    # -------------------------------------------------------------------
    "MYH6": {
        "gene_id": 4624,
        "chrom": "chr14",
        "start": 23_381_987,
        "end": 23_408_273,
        "strand": "-",
        "refseq": "NC_000014.9",
        "cytoband": "14q11.2",
        "exon_count": 39,
        "role": (
            "Myosin heavy chain 6 (alpha-MHC) -- predominant sarcomeric "
            "myosin isoform in human atrial cardiomyocytes.  Functions "
            "as the fast-twitch isoform with higher ATPase activity "
            "compared to MYH7 (beta-MHC).  Mutations cause atrial septal "
            "defect type 3 (ASD3), sick sinus syndrome (SSS1), and "
            "dilated/hypertrophic cardiomyopathy.  MYH6 and MYH7 are "
            "arranged in a head-to-head tandem on 14q11.2."
        ),
        "disease": "ASD3 / sick sinus syndrome / cardiomyopathy",
        "omim_gene": 160710,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Ile820Asn (I820N)",
                "rsid": "rs267607048",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Converter domain; ASD type 3 with conduction defects."
                ),
            },
            {
                "name": "p.Arg795Gln (R795Q)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Associated with sick sinus syndrome.",
            },
            {
                "name": "p.Ala230Pro (A230P)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Head domain; dilated cardiomyopathy.",
            },
        ],
        "strategy": (
            "AAV gene replacement is challenging: MYH6 cDNA ~5.8 kb, "
            "at the limit of AAV capacity.  For ASD: surgical closure "
            "remains definitive.  For sick sinus syndrome: pacemaker "
            "implantation.  For cardiomyopathy: standard HCM/DCM "
            "management.  Research approaches: (1) dual AAV split-intein "
            "for gene replacement; (2) allele-specific silencing of "
            "dominant-negative alleles; (3) base editing for recurrent "
            "missense variants.  Preclinical only."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Managed with cardiac surgery "
            "(ASD), pacemaker (SSS), and standard HF therapy (DCM)."
        ),
        "conditions": [
            "congenital_heart_defect", "ASD", "sick_sinus_syndrome",
            "cardiomyopathy", "DCM", "HCM",
        ],
    },
}


# ============================================================================
# UNIFIED ACCESS: all remaining disease targets in a single dictionary
# ============================================================================

ALL_REMAINING_TARGETS = {
    **MITOCHONDRIAL_TARGETS,
    **COMPLEMENT_TARGETS,
    **COAGULATION_TARGETS,
    **PAIN_TARGETS,
    **EPILEPSY_TARGETS,
    **CONGENITAL_HEART_TARGETS,
}


# ============================================================================
# Helper functions
# ============================================================================

def get_target(gene_symbol: str) -> dict | None:
    """Return the target dict for *gene_symbol*, or None."""
    return ALL_REMAINING_TARGETS.get(gene_symbol.upper())


def get_targets_by_condition(condition: str) -> dict[str, dict]:
    """Return all targets matching a given condition string."""
    return {
        gene: info
        for gene, info in ALL_REMAINING_TARGETS.items()
        if condition in info.get("conditions", [])
    }


def get_targets_by_category(category: str) -> dict:
    """Return the appropriate sub-dictionary by disease category."""
    categories = {
        "mitochondrial": MITOCHONDRIAL_TARGETS,
        "complement": COMPLEMENT_TARGETS,
        "coagulation": COAGULATION_TARGETS,
        "pain": PAIN_TARGETS,
        "epilepsy": EPILEPSY_TARGETS,
        "congenital_heart": CONGENITAL_HEART_TARGETS,
    }
    return categories.get(category.lower(), {})


def get_bed_regions() -> list[tuple[str, int, int, str]]:
    """Return a BED-format list of (chrom, start, end, name) for all targets."""
    regions = []
    for name, info in ALL_REMAINING_TARGETS.items():
        if info.get("start") and info.get("end"):
            regions.append((info["chrom"], info["start"], info["end"], name))
    return sorted(regions, key=lambda r: (r[0], r[1]))


def get_active_clinical_programs() -> dict[str, dict]:
    """Return only genes with active (non-preclinical) clinical programs."""
    active = {}
    for gene, info in ALL_REMAINING_TARGETS.items():
        programs = info.get("clinical_programs", "")
        if isinstance(programs, str):
            if "Phase" in programs or "approved" in programs.lower():
                active[gene] = info
        elif isinstance(programs, list):
            for trial in programs:
                phase = trial.get("phase", "")
                if phase and "preclinical" not in phase.lower():
                    active[gene] = info
                    break
    return active


# ============================================================================
# COORDINATE SUMMARY TABLE (for quick pipeline reference)
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
    for gene, info in ALL_REMAINING_TARGETS.items()
    if info.get("start") is not None
}
