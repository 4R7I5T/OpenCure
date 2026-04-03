"""
Endocrine disease gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers congenital adrenal hyperplasia (CAH), multiple endocrine neoplasia
(MEN), growth hormone deficiency, and related monogenic endocrinopathies.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - ClinicalTrials.gov, Neurocrine/Spruce CYP21A2 programs
"""


ENDOCRINE_TARGETS = {

    # ===================================================================
    # CAH -- CYP21A2 (21-Hydroxylase)
    # ===================================================================
    "CYP21A2": {
        "gene_id": 1589,
        "chrom": "chr6",
        "start": 32_038_265,
        "end": 32_041_670,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6p21.33",
        "exon_count": 10,
        "role": (
            "Steroid 21-hydroxylase -- adrenal cortex enzyme converting "
            "17-hydroxyprogesterone to 11-deoxycortisol (cortisol pathway) "
            "and progesterone to deoxycorticosterone (aldosterone pathway).  "
            "Deficiency causes >95% of congenital adrenal hyperplasia (CAH): "
            "cortisol deficiency, aldosterone deficiency (salt-wasting), and "
            "androgen excess (virilization).  Incidence ~1:15,000."
        ),
        "disease": "Congenital adrenal hyperplasia (21-hydroxylase deficiency)",
        "omim_disease": 201910,
        "omim_gene": 613815,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "IVS2 splice (c.293-13A/C>G)",
                "rsid": "rs6467",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common severe CAH mutation (~25% of alleles).  "
                    "Creates cryptic splice site in intron 2.  Gene conversion "
                    "from CYP21A1P pseudogene."
                ),
            },
            {
                "name": "p.Ile172Asn (I172N)",
                "rsid": "rs6475",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Simple virilizing form; ~2% residual activity.",
            },
            {
                "name": "p.Val281Leu (V281L)",
                "rsid": "rs6471",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Non-classic CAH; ~20-50% residual activity.  Most common NC-CAH allele.",
            },
            {
                "name": "Large deletion / gene conversion",
                "rsid": None,
                "consequence": "large_deletion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~20% of alleles.  Unequal crossover between CYP21A2 and "
                    "CYP21A1P pseudogene (tandem in HLA region)."
                ),
            },
        ],
        "strategy": (
            "(1) AAV adrenal-directed gene replacement (cDNA ~1.5 kb, ideal "
            "for AAV) -- Neurocrine/Spruce approach; "
            "(2) Base editing to correct IVS2 splice or I172N mutations; "
            "(3) CRISPR-mediated gene conversion repair (complex due to "
            "pseudogene/gene architecture in MHC region).  "
            "Delivery: AAV5 or AAV8 with adrenal-specific promoter (MC2R "
            "promoter or ACTH-responsive element)."
        ),
        "clinical_programs": (
            "Neurocrine Biosciences/Spruce Biosciences BBP-631 (AAV5-CYP21A2) "
            "Phase 1/2 for CAH (NCT04783610).  "
            "Crinecerfont (CRF1 antagonist) Phase 3 -- small molecule, "
            "not gene therapy.  "
            "Standard of care: glucocorticoid + fludrocortisone replacement."
        ),
        "conditions": ["congenital_adrenal_hyperplasia", "CAH", "endocrine",
                        "21_hydroxylase_deficiency", "adrenal"],
    },

    # ===================================================================
    # CAH -- CYP11B1 (11-Beta-Hydroxylase)
    # ===================================================================
    "CYP11B1": {
        "gene_id": 1584,
        "chrom": "chr8",
        "start": 143_953_773,
        "end": 143_961_234,
        "strand": "+",
        "refseq": "NC_000008.11",
        "cytoband": "8q24.3",
        "exon_count": 9,
        "role": (
            "11-beta-hydroxylase -- converts 11-deoxycortisol to cortisol.  "
            "Deficiency causes CAH with hypertension (~5-8% of CAH).  "
            "Excess deoxycorticosterone -> mineralocorticoid hypertension + "
            "virilization."
        ),
        "disease": "CAH (11-beta-hydroxylase deficiency)",
        "omim_disease": 202010,
        "omim_gene": 610613,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg448His (R448H)",
                "rsid": "rs104894065",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Moroccan Jewish founder; ~40% of alleles in this population.",
            },
        ],
        "strategy": (
            "(1) AAV adrenal gene replacement (cDNA ~1.5 kb, fits AAV); "
            "(2) Base editing for R448H."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["congenital_adrenal_hyperplasia", "CAH", "endocrine",
                        "hypertensive_CAH"],
    },

    # ===================================================================
    # MEN2 -- RET Proto-Oncogene
    # ===================================================================
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
            "RET -- receptor tyrosine kinase.  Gain-of-function mutations "
            "cause MEN2A (medullary thyroid carcinoma + pheochromocytoma) "
            "and MEN2B (aggressive MTC + mucosal neuromas).  Loss-of-function "
            "causes Hirschsprung disease."
        ),
        "disease": "Multiple endocrine neoplasia type 2 (MEN2)",
        "omim_disease": 171400,
        "omim_gene": 164761,
        "inheritance": "AD (GOF)",
        "key_variants": [
            {
                "name": "p.Cys634Arg (C634R)",
                "rsid": "rs77724903",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common MEN2A mutation; cysteine-rich domain; ~85% penetrance for MTC.",
            },
            {
                "name": "p.Met918Thr (M918T)",
                "rsid": "rs74799832",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "MEN2B; most aggressive (MTC by age 1); kinase domain activation.",
            },
        ],
        "strategy": (
            "(1) Allele-specific CRISPR disruption of mutant RET allele to "
            "prevent MTC development (prophylactic gene therapy as alternative "
            "to prophylactic thyroidectomy); "
            "(2) CRISPRi to reduce mutant RET expression in thyroid C-cells; "
            "(3) Base editing to correct C634R or M918T.  "
            "Delivery: AAV with C-cell-specific promoter (calcitonin promoter)."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "Selpercatinib (RET inhibitor) FDA-approved for RET-mutant cancers.  "
            "Prophylactic thyroidectomy standard of care for MEN2."
        ),
        "conditions": ["MEN2", "MEN2A", "MEN2B", "medullary_thyroid_carcinoma",
                        "endocrine", "pheochromocytoma", "hirschsprung"],
    },

    # ===================================================================
    # MEN1 -- Menin
    # ===================================================================
    "MEN1": {
        "gene_id": 4221,
        "chrom": "chr11",
        "start": 64_803_512,
        "end": 64_810_960,
        "strand": "+",
        "refseq": "NC_000011.10",
        "cytoband": "11q13.1",
        "exon_count": 10,
        "role": (
            "Menin -- tumor suppressor, scaffolding protein regulating gene "
            "expression, DNA repair, and cell cycle.  Loss-of-function causes "
            "MEN1: parathyroid hyperplasia/adenoma, pancreatic neuroendocrine "
            "tumors, pituitary adenomas."
        ),
        "disease": "Multiple endocrine neoplasia type 1 (MEN1)",
        "omim_disease": 131100,
        "omim_gene": 613733,
        "inheritance": "AD (tumor suppressor, LOH)",
        "key_variants": [
            {
                "name": ">1500 unique mutations (no hotspot)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Truncating mutations throughout; extreme allelic heterogeneity.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type MEN1 allele in target tissues; "
            "(2) AAV gene replacement (cDNA ~1.8 kb, fits AAV); "
            "(3) CRISPR-based gene reactivation of epigenetically silenced allele.  "
            "Delivery: tissue-specific AAV to parathyroid/pancreas."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "Surveillance + surgical management standard of care."
        ),
        "conditions": ["MEN1", "endocrine", "parathyroid_adenoma",
                        "pancreatic_NET", "pituitary_adenoma"],
    },

    # ===================================================================
    # Growth Hormone Deficiency -- GH1
    # ===================================================================
    "GH1": {
        "gene_id": 2688,
        "chrom": "chr17",
        "start": 63_917_193,
        "end": 63_918_847,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q23.3",
        "exon_count": 5,
        "role": (
            "Growth hormone 1 -- pituitary somatotroph hormone.  Mutations "
            "cause isolated growth hormone deficiency (IGHD): types IA (AR, "
            "complete absence), IB (AR, partial), II (AD, splice mutations), "
            "III (X-linked)."
        ),
        "disease": "Isolated growth hormone deficiency (IGHD)",
        "omim_disease": 262400,
        "omim_gene": 139250,
        "inheritance": "AR (IA/IB) / AD (II)",
        "key_variants": [
            {
                "name": "GH1 gene deletion (6.7 kb)",
                "rsid": None,
                "consequence": "large_deletion",
                "clinical_significance": "pathogenic",
                "notes": "IGHD type IA; homozygous deletion.  Anti-GH antibodies with ERT.",
            },
            {
                "name": "IVS3+1 or +2 splice mutations",
                "rsid": None,
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": "IGHD type II (AD); exon 3 skipping -> dominant-negative.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~0.6 kb, ideal for AAV) with "
            "pituitary-specific promoter; "
            "(2) Allele-specific silencing for dominant type II; "
            "(3) AAV liver-directed GH expression for durable replacement.  "
            "Delivery: AAV to pituitary (challenging) or AAV8 liver-directed."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "Recombinant GH (somatropin) standard of care.  "
            "Long-acting GH (somapacitan, lonapegsomatropin) approved."
        ),
        "conditions": ["growth_hormone_deficiency", "IGHD", "endocrine",
                        "short_stature", "pituitary"],
    },

    # ===================================================================
    # Pseudohypoparathyroidism -- GNAS
    # ===================================================================
    "GNAS": {
        "gene_id": 2778,
        "chrom": "chr20",
        "start": 58_839_718,
        "end": 58_911_192,
        "strand": "+",
        "refseq": "NC_000020.11",
        "cytoband": "20q13.32",
        "exon_count": 16,
        "role": (
            "Gs-alpha subunit -- stimulatory G-protein alpha subunit, "
            "essential for cAMP signaling downstream of PTH, TSH, GHRH, and "
            "other GPCR hormones.  GNAS is imprinted: maternal allele active "
            "in renal proximal tubule.  Maternal mutations cause "
            "pseudohypoparathyroidism type Ia (PHP-Ia / Albright hereditary "
            "osteodystrophy): PTH resistance, obesity, brachydactyly."
        ),
        "disease": "Pseudohypoparathyroidism type Ia (PHP-Ia)",
        "omim_disease": 103580,
        "omim_gene": 139320,
        "inheritance": "AD (imprinted -- maternal allele)",
        "key_variants": [
            {
                "name": "c.565_568delGACT (4-bp deletion, exon 7)",
                "rsid": "rs387906653",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Hotspot deletion in poly-GACT tract; most common mutation.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to reactivate silenced paternal GNAS allele in "
            "kidney (imprinting reversal -- novel approach); "
            "(2) AAV gene replacement (cDNA ~1.2 kb, fits AAV) with renal "
            "proximal tubule promoter; "
            "(3) Base editing for specific mutations.  "
            "Delivery: kidney-tropic AAV."
        ),
        "clinical_programs": "No gene therapy trials. Calcium + calcitriol SOC.",
        "conditions": ["pseudohypoparathyroidism", "PHP", "endocrine",
                        "Albright_hereditary_osteodystrophy", "PTH_resistance"],
    },
}
