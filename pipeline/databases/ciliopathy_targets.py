"""
Ciliopathy gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers primary ciliary dyskinesia (PCD), Bardet-Biedl syndrome (BBS),
and Joubert syndrome.  Ciliopathies are a class of diseases caused by
dysfunction of primary or motile cilia.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


CILIOPATHY_TARGETS = {

    # ===================================================================
    # PCD -- DNAI1
    # ===================================================================
    "DNAI1": {
        "gene_id": 27019,
        "chrom": "chr9",
        "start": 34_460_471,
        "end": 34_520_880,
        "strand": "+",
        "refseq": "NC_000009.12",
        "cytoband": "9p13.3",
        "exon_count": 20,
        "role": (
            "Dynein axonemal intermediate chain 1 -- component of the outer "
            "dynein arm (ODA), the motor powering ciliary beating.  Biallelic "
            "mutations cause primary ciliary dyskinesia (PCD): chronic "
            "sinopulmonary disease, bronchiectasis, situs inversus (~50%), "
            "male infertility.  Incidence ~1:15,000-20,000."
        ),
        "disease": "Primary ciliary dyskinesia (PCD/Kartagener syndrome)",
        "omim_disease": 244400,
        "omim_gene": 604366,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "IVS1+2_3insT (splice, founder)",
                "rsid": "rs387907025",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": "Most common DNAI1 mutation; ~55% of DNAI1-PCD alleles.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement to airway epithelium (cDNA ~2.1 kb, "
            "fits AAV); "
            "(2) mRNA/LNP nebulized delivery to airway epithelial cells; "
            "(3) CRISPR correction in airway basal stem cells (self-renewing).  "
            "Delivery: inhaled AAV or LNP (lung-directed)."
        ),
        "clinical_programs": "Preclinical only. Airway clearance therapy SOC.",
        "conditions": ["primary_ciliary_dyskinesia", "PCD", "kartagener",
                        "ciliopathy", "bronchiectasis", "situs_inversus"],
    },

    "DNAH5": {
        "gene_id": 1767,
        "chrom": "chr5",
        "start": 13_690_445,
        "end": 13_944_451,
        "strand": "+",
        "refseq": "NC_000005.10",
        "cytoband": "5p15.2",
        "exon_count": 79,
        "role": (
            "Dynein axonemal heavy chain 5 -- ODA heavy chain, the motor "
            "ATPase.  Most common PCD gene (~28% of PCD).  Absent ODA on "
            "electron microscopy."
        ),
        "disease": "PCD (most common gene)",
        "omim_disease": 608644,
        "omim_gene": 603335,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various truncating (>100 mutations)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Extreme heterogeneity; no hotspot.",
            },
        ],
        "strategy": (
            "(1) Gene too large for AAV (~14 kb cDNA); dual-AAV or "
            "mRNA/LNP approach; "
            "(2) mRNA inhaled delivery (transient but repeatable); "
            "(3) CRISPR correction of specific variants in airway basal cells."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["primary_ciliary_dyskinesia", "PCD", "ciliopathy",
                        "ODA_defect"],
    },

    "CCDC39": {
        "gene_id": 339829,
        "chrom": "chr3",
        "start": 180_330_974,
        "end": 180_427_264,
        "strand": "+",
        "refseq": "NC_000003.12",
        "cytoband": "3q26.33",
        "exon_count": 20,
        "role": (
            "Coiled-coil domain containing 39 -- axonemal ruler protein "
            "determining 96-nm repeat periodicity.  Mutations cause PCD with "
            "inner dynein arm + microtubular disorganization defect."
        ),
        "disease": "PCD with IDA/MTD defect",
        "omim_disease": 613807,
        "omim_gene": 613798,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.2190delA (p.Glu730fs)",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Consanguineous families; severe PCD.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.8 kb, fits AAV); "
            "(2) mRNA/LNP inhaled delivery."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["primary_ciliary_dyskinesia", "PCD", "ciliopathy"],
    },

    # ===================================================================
    # Bardet-Biedl Syndrome
    # ===================================================================
    "BBS1": {
        "gene_id": 582,
        "chrom": "chr11",
        "start": 66_278_094,
        "end": 66_301_131,
        "strand": "+",
        "refseq": "NC_000011.10",
        "cytoband": "11q13.2",
        "exon_count": 17,
        "role": (
            "Bardet-Biedl syndrome 1 -- BBSome complex component for "
            "ciliary protein trafficking.  Most common BBS gene (~23% of "
            "cases).  BBS: rod-cone dystrophy, obesity, polydactyly, renal "
            "anomalies, learning difficulties, hypogonadism."
        ),
        "disease": "Bardet-Biedl syndrome (most common gene)",
        "omim_disease": 209900,
        "omim_gene": 209901,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Met390Arg (M390R)",
                "rsid": "rs113993960",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common BBS1 mutation (~80% of BBS1 alleles).",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.8 kb, fits AAV); "
            "(2) Subretinal AAV for retinal degeneration component; "
            "(3) Base editing for M390R.  "
            "Delivery: AAV for retinal component; systemic for other features."
        ),
        "clinical_programs": (
            "Setmelanotide (MC4R agonist) FDA-approved 2020 for BBS obesity.  "
            "No gene therapy trials for BBS1 as of 2026."
        ),
        "conditions": ["bardet_biedl_syndrome", "BBS", "ciliopathy",
                        "retinal_dystrophy", "syndromic_obesity"],
    },

    "BBS10": {
        "gene_id": 79738,
        "chrom": "chr12",
        "start": 76_013_075,
        "end": 76_019_172,
        "strand": "-",
        "refseq": "NC_000012.12",
        "cytoband": "12q21.2",
        "exon_count": 2,
        "role": (
            "Bardet-Biedl syndrome 10 -- chaperonin-like protein required "
            "for BBSome assembly.  Second most common BBS gene (~20%)."
        ),
        "disease": "Bardet-Biedl syndrome",
        "omim_disease": 615987,
        "omim_gene": 610148,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.271dupT (p.Cys91Leufs)",
                "rsid": "rs587776662",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Founder mutation in Faroe Islands (carrier 1/6).",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.1 kb, fits AAV); "
            "(2) Base editing for specific mutations."
        ),
        "clinical_programs": "Setmelanotide for obesity component. No gene therapy.",
        "conditions": ["bardet_biedl_syndrome", "BBS", "ciliopathy"],
    },

    # ===================================================================
    # Joubert Syndrome
    # ===================================================================
    "TMEM67": {
        "gene_id": 91147,
        "chrom": "chr8",
        "start": 93_768_566,
        "end": 93_833_638,
        "strand": "-",
        "refseq": "NC_000008.11",
        "cytoband": "8q22.1",
        "exon_count": 27,
        "role": (
            "Transmembrane protein 67 (meckelin/MKS3) -- transition zone "
            "protein of primary cilia.  Mutations cause Joubert syndrome "
            "(molar tooth sign on MRI, ataxia, abnormal breathing, ID) "
            "and Meckel-Gruber syndrome (lethal, encephalocele + polycystic "
            "kidneys + polydactyly)."
        ),
        "disease": "Joubert syndrome / Meckel-Gruber syndrome",
        "omim_disease": 610688,
        "omim_gene": 609884,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various missense and truncating",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Genotype-phenotype correlations weak.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~3.0 kb, fits AAV); "
            "(2) CRISPRa to upregulate residual expression in hypomorphic alleles."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["joubert_syndrome", "meckel_gruber", "ciliopathy",
                        "molar_tooth_sign"],
    },

    "NPHP1": {
        "gene_id": 4867,
        "chrom": "chr2",
        "start": 110_862_588,
        "end": 110_982_579,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q13",
        "exon_count": 20,
        "role": (
            "Nephrocystin-1 -- ciliary/centrosomal protein.  Homozygous "
            "deletions are the most common cause of nephronophthisis (NPHP): "
            "autosomal recessive cystic kidney disease, the most common "
            "genetic cause of ESRD in children and young adults."
        ),
        "disease": "Nephronophthisis type 1 / Joubert with NPHP",
        "omim_disease": 256100,
        "omim_gene": 607100,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Homozygous deletion (~290 kb)",
                "rsid": None,
                "consequence": "large_deletion",
                "clinical_significance": "pathogenic",
                "notes": "~21% of all NPHP; NAHR between flanking repeats.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.2 kb, fits AAV) to kidney; "
            "(2) Kidney transplant is curative for renal disease.  "
            "Delivery: kidney-tropic AAV."
        ),
        "clinical_programs": "No gene therapy trials. Kidney transplant SOC.",
        "conditions": ["nephronophthisis", "NPHP", "ciliopathy",
                        "cystic_kidney_disease", "renal"],
    },
}
