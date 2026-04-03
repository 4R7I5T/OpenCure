"""
Renal / kidney disease gene-therapy / CRISPR targets -- GRCh38 (hg38) coordinates.

Covers hereditary nephropathies, polycystic kidney disease, nephrotic syndrome,
and renal tubular disorders amenable to gene therapy or CRISPR correction.

Categories:
  1. Alport Syndrome (COL4A3, COL4A4, COL4A5)
  2. Polycystic Kidney Disease (PKD1, PKD2)
  3. APOL1 Nephropathy (APOL1)
  4. Cystinosis (CTNS)
  5. Congenital Nephrotic Syndrome (NPHS1, NPHS2)
  6. Fabry Disease -- renal manifestation (GLA)

Coordinates are from NCBI Gene, GRCh38.p14 (RefSeq annotation RS_2025_08).

IMPORTANT -- verify every coordinate against current NCBI / Ensembl releases
before production use.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, dbSNP, OMIM
  - ClinicalTrials.gov
  - Published literature through early 2026
"""


# ============================================================================
# 1. ALPORT SYNDROME
# ============================================================================

ALPORT_TARGETS = {

    "COL4A5": {
        "gene_id": 1287,
        "chrom": "chrX",
        "start": 107_683_074,
        "end": 107_940_775,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xq22.3",
        "exon_count": 51,
        "role": (
            "Type IV collagen alpha-5 chain -- essential component of the "
            "glomerular basement membrane (GBM), cochlea, and lens capsule.  "
            "Mutations cause X-linked Alport syndrome (~80% of all Alport), "
            "characterised by progressive glomerulonephritis, sensorineural "
            "hearing loss, and anterior lenticonus."
        ),
        "disease": "X-linked Alport syndrome",
        "omim_disease": 301050,
        "omim_gene": 303630,
        "inheritance": "XL",
        "key_variants": [
            {
                "name": "p.Gly624Asp (G624D)",
                "rsid": "rs104886043",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Gly-X-Y repeat disruption in collagenous domain.",
            },
            {
                "name": "Large deletions (exon-level)",
                "rsid": None,
                "consequence": "large_deletion",
                "clinical_significance": "pathogenic",
                "notes": "~15% of X-linked Alport; no residual protein.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated gene replacement -- challenging due to large cDNA "
            "(~5.1 kb) but feasible with compact promoter; "
            "(2) Base editing for Gly->X missense mutations in collagenous "
            "domain; "
            "(3) Ex vivo CRISPR correction of patient podocyte progenitors.  "
            "Delivery: kidney-tropic AAV (AAV2/9 or engineered capsids)."
        ),
        "clinical_programs": (
            "No gene therapy trials for COL4A5 as of 2026.  "
            "ACEi/ARB standard of care to slow progression.  "
            "Bardoxolone methyl (Reata) Phase 3 for Alport CKD."
        ),
        "conditions": ["alport_syndrome", "hereditary_nephritis",
                        "renal", "glomerular_disease", "hearing_loss"],
    },

    "COL4A3": {
        "gene_id": 1285,
        "chrom": "chr2",
        "start": 227_011_935,
        "end": 227_164_684,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q36.3",
        "exon_count": 52,
        "role": (
            "Type IV collagen alpha-3 chain -- GBM structural protein.  "
            "Biallelic mutations cause autosomal recessive Alport syndrome.  "
            "Heterozygous variants cause thin basement membrane nephropathy "
            "(benign familial hematuria) or autosomal dominant Alport."
        ),
        "disease": "Autosomal recessive/dominant Alport syndrome",
        "omim_disease": 203780,
        "omim_gene": 120070,
        "inheritance": "AR / AD",
        "key_variants": [
            {
                "name": "p.Gly695Arg",
                "rsid": "rs200562792",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Gly-X-Y disruption; collagenous domain.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated gene replacement (cDNA ~5.1 kb, tight AAV fit); "
            "(2) Base editing for recurrent Gly substitutions; "
            "(3) CRISPRa to upregulate wild-type allele in heterozygous patients."
        ),
        "clinical_programs": "No gene therapy or CRISPR trials as of 2026.",
        "conditions": ["alport_syndrome", "hereditary_nephritis", "renal",
                        "thin_basement_membrane"],
    },

    "COL4A4": {
        "gene_id": 1286,
        "chrom": "chr2",
        "start": 227_164_799,
        "end": 227_318_282,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2q36.3",
        "exon_count": 48,
        "role": (
            "Type IV collagen alpha-4 chain -- GBM structural protein.  "
            "Mutations cause autosomal Alport syndrome.  COL4A3 and COL4A4 "
            "are adjacent head-to-head on chromosome 2."
        ),
        "disease": "Autosomal Alport syndrome",
        "omim_disease": 203780,
        "omim_gene": 120131,
        "inheritance": "AR / AD",
        "key_variants": [
            {
                "name": "p.Gly545Glu",
                "rsid": "rs200090482",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Collagenous domain; Gly substitution disrupts triple helix.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~5.0 kb); "
            "(2) Base editing for specific Gly substitutions."
        ),
        "clinical_programs": "No gene therapy or CRISPR trials as of 2026.",
        "conditions": ["alport_syndrome", "hereditary_nephritis", "renal"],
    },
}


# ============================================================================
# 2. POLYCYSTIC KIDNEY DISEASE
# ============================================================================

PKD_TARGETS = {

    "PKD1": {
        "gene_id": 5310,
        "chrom": "chr16",
        "start": 2_088_708,
        "end": 2_135_898,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16p13.3",
        "exon_count": 46,
        "role": (
            "Polycystin-1 -- large transmembrane glycoprotein (~4303 aa), "
            "mechanosensor on primary cilia.  Mutations cause autosomal "
            "dominant polycystic kidney disease (ADPKD, ~85% of cases).  "
            "Loss of polycystin signaling -> uncontrolled tubular cell "
            "proliferation and cyst formation."
        ),
        "disease": "Autosomal dominant polycystic kidney disease (ADPKD)",
        "omim_disease": 173900,
        "omim_gene": 601313,
        "inheritance": "AD (two-hit model)",
        "key_variants": [
            {
                "name": "Truncating variants (>1000 described)",
                "rsid": None,
                "consequence": "nonsense/frameshift/splice",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Extreme allelic heterogeneity.  Truncating variants "
                    "have more severe course than non-truncating."
                ),
            },
        ],
        "mutational_landscape": (
            ">2600 pathogenic variants in the PKD Mutation Database.  "
            "No single common mutation -- most families have private "
            "mutations.  6 pseudogenes on 16p13.1 complicate genetic testing."
        ),
        "strategy": (
            "(1) CRISPRa to upregulate the wild-type PKD1 allele (two-hit "
            "model: increasing WT dosage may prevent cystogenesis); "
            "(2) Gene too large for AAV (~12.9 kb cDNA); dual-AAV or "
            "lentiviral delivery possible; "
            "(3) Base editing for specific recurrent missense mutations; "
            "(4) In vivo CRISPR activation of PKD1 in renal tubular cells.  "
            "Delivery: kidney-tropic AAV or lipid nanoparticle (LNP) renal."
        ),
        "clinical_programs": (
            "Tolvaptan (Jynarque) approved for ADPKD (vasopressin antagonist).  "
            "No gene therapy or CRISPR trials for PKD1 as of 2026.  "
            "Regulus RG-012 (anti-miR-17) Phase 2 for ADPKD."
        ),
        "conditions": ["ADPKD", "polycystic_kidney_disease", "renal",
                        "kidney_cysts"],
    },

    "PKD2": {
        "gene_id": 5311,
        "chrom": "chr4",
        "start": 88_007_647,
        "end": 88_077_759,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4q22.1",
        "exon_count": 15,
        "role": (
            "Polycystin-2 (TRPP2) -- calcium-permeable cation channel on "
            "primary cilia and ER.  Forms complex with polycystin-1.  "
            "Mutations cause ADPKD (~15% of cases), generally milder course "
            "than PKD1 (ESRD ~10 years later)."
        ),
        "disease": "ADPKD (type 2)",
        "omim_disease": 613095,
        "omim_gene": 173910,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Truncating variants throughout gene",
                "rsid": None,
                "consequence": "nonsense/frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Milder phenotype than PKD1 truncating variants.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.9 kb, fits AAV comfortably); "
            "(2) CRISPRa to upregulate wild-type allele; "
            "(3) Base editing for specific missense variants.  "
            "Delivery: kidney-tropic AAV."
        ),
        "clinical_programs": "No gene therapy trials for PKD2 as of 2026.",
        "conditions": ["ADPKD", "polycystic_kidney_disease", "renal",
                        "kidney_cysts"],
    },
}


# ============================================================================
# 3. APOL1 NEPHROPATHY
# ============================================================================

APOL1_TARGETS = {

    "APOL1": {
        "gene_id": 8542,
        "chrom": "chr22",
        "start": 36_253_070,
        "end": 36_267_530,
        "strand": "-",
        "refseq": "NC_000022.11",
        "cytoband": "22q12.3",
        "exon_count": 7,
        "role": (
            "Apolipoprotein L1 -- innate immune factor that kills "
            "Trypanosoma brucei (African sleeping sickness).  G1 and G2 "
            "risk alleles evolved under selection pressure from T. brucei "
            "rhodesiense.  Two copies of risk alleles (G1/G1, G2/G2, or "
            "G1/G2) dramatically increase risk of FSGS, HIVAN, hypertensive "
            "nephrosclerosis, and sickle cell nephropathy.  "
            "~13% of African Americans carry two risk alleles."
        ),
        "disease": "APOL1-associated nephropathy (FSGS, HIVAN)",
        "omim_disease": 612551,
        "omim_gene": 603743,
        "inheritance": "recessive risk (two risk alleles needed)",
        "key_variants": [
            {
                "name": "G1 haplotype (S342G + I384M)",
                "rsid": "rs73885319 + rs60910145",
                "consequence": "missense",
                "clinical_significance": "risk_factor",
                "notes": (
                    "Two linked missense variants.  ~22% allele frequency "
                    "in African Americans."
                ),
            },
            {
                "name": "G2 (del388N389Y)",
                "rsid": "rs71785313",
                "consequence": "in_frame_deletion",
                "clinical_significance": "risk_factor",
                "notes": (
                    "6-bp deletion removing N388/Y389.  ~13% allele frequency "
                    "in African Americans.  Alters SRA binding domain."
                ),
            },
        ],
        "strategy": (
            "(1) Base editing to revert G1 risk variants (S342G, I384M) to "
            "reference sequence in podocytes; "
            "(2) CRISPRi to reduce APOL1 expression in kidney (APOL1 "
            "reduction is tolerable -- absent in most non-African populations); "
            "(3) ASO-mediated APOL1 reduction (Ionis approach).  "
            "Delivery: kidney-tropic LNP or AAV with podocyte promoter."
        ),
        "clinical_programs": (
            "Vertex VX-147 (inaxaplin) -- small molecule APOL1 inhibitor, "
            "Phase 2/3 for APOL1-mediated FSGS (NCT05312879).  "
            "Ionis APOL1 ASO preclinical.  "
            "No CRISPR trials for APOL1 as of 2026."
        ),
        "conditions": ["APOL1_nephropathy", "FSGS", "HIVAN", "renal",
                        "kidney_disease", "nephrotic_syndrome"],
    },
}


# ============================================================================
# 4. CYSTINOSIS
# ============================================================================

CYSTINOSIS_TARGETS = {

    "CTNS": {
        "gene_id": 1497,
        "chrom": "chr17",
        "start": 3_636_116,
        "end": 3_662_038,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17p13.2",
        "exon_count": 12,
        "role": (
            "Cystinosin -- lysosomal cystine transporter.  Loss-of-function "
            "causes cystinosis (cystine accumulation in lysosomes), leading to "
            "renal Fanconi syndrome, progressive CKD, and multisystem disease.  "
            "Infantile nephropathic form is most severe."
        ),
        "disease": "Nephropathic cystinosis",
        "omim_disease": 219800,
        "omim_gene": 606272,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "57-kb deletion (European founder)",
                "rsid": None,
                "consequence": "large_deletion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mutation in Northern Europeans (~76% of "
                    "alleles).  Removes CTNS exons 1-10 and upstream genes."
                ),
            },
            {
                "name": "p.Trp138Ter (W138X)",
                "rsid": "rs121918228",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Second most common in European populations.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated CTNS gene replacement (cDNA ~1.1 kb, ideal for "
            "AAV) -- AVROBIO approach using HSC gene therapy; "
            "(2) Ex vivo CRISPR correction in autologous HSCs (gene therapy "
            "allows cross-correction via enzyme secretion); "
            "(3) In vivo AAV delivery to kidney with CTNS cDNA.  "
            "Delivery: lentiviral HSC gene therapy or AAV renal."
        ),
        "clinical_programs": (
            "AVROBIO AVR-RD-04 -- lentiviral HSC gene therapy for cystinosis "
            "Phase 1/2 (NCT03897361).  "
            "UCSD CTNS gene therapy preclinical.  "
            "Cysteamine (Cystagon/Procysbi) standard of care."
        ),
        "conditions": ["cystinosis", "renal_fanconi_syndrome", "renal",
                        "lysosomal_storage"],
    },
}


# ============================================================================
# 5. CONGENITAL NEPHROTIC SYNDROME
# ============================================================================

NEPHROTIC_TARGETS = {

    "NPHS1": {
        "gene_id": 4868,
        "chrom": "chr19",
        "start": 35_819_518,
        "end": 35_845_370,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19q13.12",
        "exon_count": 29,
        "role": (
            "Nephrin -- podocyte slit diaphragm protein, essential for "
            "glomerular filtration barrier.  Biallelic loss-of-function causes "
            "congenital nephrotic syndrome of the Finnish type (CNF), the "
            "most severe form with massive proteinuria from birth."
        ),
        "disease": "Congenital nephrotic syndrome, Finnish type (CNF)",
        "omim_disease": 256300,
        "omim_gene": 602716,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Fin-major (c.121delCT, p.Leu41fs)",
                "rsid": "rs104893768",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Finnish founder mutation; ~78% of Finnish alleles.",
            },
            {
                "name": "Fin-minor (c.3325C>T, p.Arg1109Ter)",
                "rsid": "rs104893769",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Second Finnish founder; ~16% of Finnish alleles.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated gene replacement (cDNA ~3.6 kb, fits AAV) with "
            "podocyte-specific promoter (NPHS2 promoter); "
            "(2) Ex vivo CRISPR correction in iPSC-derived podocytes.  "
            "Delivery: kidney-tropic AAV with podocyte promoter."
        ),
        "clinical_programs": (
            "No gene therapy trials as of 2026.  "
            "Management: albumin infusion, bilateral nephrectomy, "
            "and kidney transplant (definitive cure)."
        ),
        "conditions": ["congenital_nephrotic_syndrome", "CNF", "renal",
                        "nephrotic_syndrome", "podocytopathy"],
    },

    "NPHS2": {
        "gene_id": 7827,
        "chrom": "chr1",
        "start": 179_520_506,
        "end": 179_545_225,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q25.2",
        "exon_count": 8,
        "role": (
            "Podocin -- podocyte lipid raft-associated protein at the slit "
            "diaphragm.  Recruits nephrin to lipid rafts for signaling.  "
            "Mutations cause steroid-resistant nephrotic syndrome (SRNS), "
            "the most common genetic form."
        ),
        "disease": "Steroid-resistant nephrotic syndrome (SRNS)",
        "omim_disease": 600995,
        "omim_gene": 604766,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg138Gln (R138Q)",
                "rsid": "rs74315342",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common NPHS2 mutation worldwide; exon 3.",
            },
            {
                "name": "p.Val260Glu (V260E)",
                "rsid": "rs28939693",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "C-terminal domain; impairs nephrin interaction.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.1 kb, ideal for AAV); "
            "(2) Base editing for R138Q (single A>G change); "
            "(3) CRISPRa to boost residual expression.  "
            "Delivery: kidney-tropic AAV with NPHS2 promoter."
        ),
        "clinical_programs": "No gene therapy or CRISPR trials as of 2026.",
        "conditions": ["steroid_resistant_nephrotic_syndrome", "SRNS", "renal",
                        "nephrotic_syndrome", "podocytopathy", "FSGS"],
    },
}


# ============================================================================
# 6. FABRY DISEASE (RENAL MANIFESTATION)
# ============================================================================

FABRY_RENAL_TARGETS = {

    "GLA": {
        "gene_id": 2717,
        "chrom": "chrX",
        "start": 101_397_795,
        "end": 101_407_925,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq22.1",
        "exon_count": 7,
        "role": (
            "Alpha-galactosidase A -- lysosomal enzyme cleaving "
            "globotriaosylceramide (Gb3/GL-3).  Deficiency causes Fabry "
            "disease: Gb3 accumulation in vascular endothelium, podocytes, "
            "cardiomyocytes, and dorsal root ganglia.  Progressive CKD is "
            "a major cause of morbidity."
        ),
        "disease": "Fabry disease",
        "omim_disease": 301500,
        "omim_gene": 300644,
        "inheritance": "XL",
        "key_variants": [
            {
                "name": "p.Arg227Ter (R227X)",
                "rsid": "rs104894834",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Classic severe Fabry; no residual enzyme activity.",
            },
            {
                "name": "p.Asn215Ser (N215S)",
                "rsid": "rs28935486",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Cardiac variant Fabry; later onset, predominantly cardiac.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.3 kb, ideal for AAV) -- "
            "liver-directed AAV to produce secreted GLA for cross-correction; "
            "(2) Ex vivo CRISPR insertion of GLA at safe harbor (AAVS1) in "
            "autologous HSCs; "
            "(3) mRNA/LNP for transient enzyme replacement.  "
            "Delivery: AAV8/AAV6 liver-directed (Freeline approach) or "
            "lentiviral HSC."
        ),
        "clinical_programs": (
            "Freeline FLT190 (AAV-GLA liver) Phase 1/2.  "
            "4D Molecular Therapeutics 4D-310 (AAV-GLA cardiac) Phase 1/2.  "
            "Sangamo ST-920 (AAV-GLA liver) Phase 1/2.  "
            "Avrobio AVR-RD-01 (lentiviral HSC) Phase 1/2.  "
            "ERT approved: agalsidase beta (Fabrazyme), agalsidase alfa "
            "(Replagal).  Migalastat (chaperone) for amenable mutations."
        ),
        "conditions": ["fabry_disease", "renal", "lysosomal_storage",
                        "GLA_deficiency", "nephropathy"],
    },
}


# ============================================================================
# Combined export
# ============================================================================

ALL_RENAL_TARGETS = {}
ALL_RENAL_TARGETS.update(ALPORT_TARGETS)
ALL_RENAL_TARGETS.update(PKD_TARGETS)
ALL_RENAL_TARGETS.update(APOL1_TARGETS)
ALL_RENAL_TARGETS.update(CYSTINOSIS_TARGETS)
ALL_RENAL_TARGETS.update(NEPHROTIC_TARGETS)
ALL_RENAL_TARGETS.update(FABRY_RENAL_TARGETS)
