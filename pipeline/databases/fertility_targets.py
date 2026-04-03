"""
Reproductive / fertility disorder gene-therapy targets -- GRCh38 (hg38).

Covers Kallmann syndrome, primary ovarian insufficiency (POI), congenital
bilateral absence of vas deferens (CBAVD), and related monogenic infertility.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


FERTILITY_TARGETS = {

    # ===================================================================
    # Kallmann Syndrome -- ANOS1 (KAL1)
    # ===================================================================
    "ANOS1": {
        "gene_id": 3730,
        "chrom": "chrX",
        "start": 8_508_768,
        "end": 8_729_048,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xp22.31",
        "exon_count": 14,
        "role": (
            "Anosmin-1 -- extracellular matrix protein guiding GnRH neuron "
            "migration from olfactory placode to hypothalamus during "
            "embryogenesis.  Mutations cause X-linked Kallmann syndrome: "
            "hypogonadotropic hypogonadism + anosmia (absent smell)."
        ),
        "disease": "Kallmann syndrome (X-linked)",
        "omim_disease": 308700,
        "omim_gene": 300836,
        "inheritance": "XLR",
        "key_variants": [
            {
                "name": "Deletions and truncating mutations",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Various; gene deletions common.",
            },
        ],
        "strategy": (
            "(1) GnRH replacement (pulsatile GnRH pump) is effective and "
            "non-genetic; "
            "(2) AAV gene replacement (cDNA ~2.0 kb fits AAV) could restore "
            "endogenous GnRH neuron function if delivered early; "
            "(3) CRISPR correction in iPSC-derived GnRH neurons for research.  "
            "Note: exogenous hormone therapy is highly effective, reducing "
            "gene therapy urgency."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "Pulsatile GnRH or gonadotropin therapy standard of care."
        ),
        "conditions": ["kallmann_syndrome", "hypogonadotropic_hypogonadism",
                        "fertility", "anosmia", "reproductive"],
    },

    # ===================================================================
    # Kallmann -- FGFR1
    # ===================================================================
    "FGFR1": {
        "gene_id": 2260,
        "chrom": "chr8",
        "start": 38_411_138,
        "end": 38_468_834,
        "strand": "+",
        "refseq": "NC_000008.11",
        "cytoband": "8p11.23",
        "exon_count": 24,
        "role": (
            "Fibroblast growth factor receptor 1 -- receptor tyrosine kinase "
            "regulating GnRH neuron migration and olfactory bulb development.  "
            "Loss-of-function causes autosomal dominant Kallmann syndrome "
            "(KAL2, ~10% of cases).  Also: Pfeiffer syndrome (GOF, separate)."
        ),
        "disease": "Kallmann syndrome (KAL2, autosomal dominant)",
        "omim_disease": 147950,
        "omim_gene": 136350,
        "inheritance": "AD (incomplete penetrance)",
        "key_variants": [
            {
                "name": "Various LOF missense and truncating",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Oligogenic modification common (FGFR1 + other loci).",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type FGFR1 allele; "
            "(2) Hormone replacement therapy is effective standard of care."
        ),
        "clinical_programs": "No gene therapy trials. Hormone therapy SOC.",
        "conditions": ["kallmann_syndrome", "hypogonadotropic_hypogonadism",
                        "fertility", "reproductive"],
    },

    # ===================================================================
    # Primary Ovarian Insufficiency -- BMP15
    # ===================================================================
    "BMP15": {
        "gene_id": 9210,
        "chrom": "chrX",
        "start": 50_910_076,
        "end": 50_916_212,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xp11.22",
        "exon_count": 2,
        "role": (
            "Bone morphogenetic protein 15 -- oocyte-secreted growth factor "
            "regulating folliculogenesis and ovulation.  Mutations cause "
            "primary ovarian insufficiency (POI) / premature ovarian failure.  "
            "Important for oocyte-granulosa cell communication."
        ),
        "disease": "Primary ovarian insufficiency (POI)",
        "omim_disease": 300510,
        "omim_gene": 300247,
        "inheritance": "XL (dosage-sensitive)",
        "key_variants": [
            {
                "name": "p.Tyr235Cys (Y235C)",
                "rsid": "rs104894767",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts BMP signaling; associated with POI.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.2 kb, ideal for AAV) with "
            "oocyte-specific promoter (ZP3); "
            "(2) CRISPRa to boost BMP15 expression; "
            "(3) Recombinant BMP15 protein as non-genetic approach.  "
            "Delivery: ovarian AAV injection (challenging)."
        ),
        "clinical_programs": "No gene therapy trials. IVF with donor oocytes SOC.",
        "conditions": ["primary_ovarian_insufficiency", "POI", "fertility",
                        "premature_ovarian_failure", "reproductive"],
    },

    # ===================================================================
    # POI -- GDF9
    # ===================================================================
    "GDF9": {
        "gene_id": 2661,
        "chrom": "chr5",
        "start": 132_861_029,
        "end": 132_866_531,
        "strand": "+",
        "refseq": "NC_000005.10",
        "cytoband": "5q31.1",
        "exon_count": 2,
        "role": (
            "Growth differentiation factor 9 -- oocyte-secreted TGF-beta "
            "family member essential for primary follicle development beyond "
            "the primordial stage.  Mutations associated with POI and "
            "dizygotic twinning."
        ),
        "disease": "POI / follicular arrest",
        "omim_disease": 618014,
        "omim_gene": 601918,
        "inheritance": "AD / AR",
        "key_variants": [
            {
                "name": "p.Ser186Tyr (S186Y)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Prodomain; impairs GDF9 processing.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.4 kb, ideal for AAV) with "
            "oocyte-specific promoter; "
            "(2) Recombinant GDF9 for in vitro oocyte maturation."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["primary_ovarian_insufficiency", "POI", "fertility",
                        "reproductive", "follicular_arrest"],
    },

    # ===================================================================
    # Hypogonadotropic Hypogonadism -- GNRHR
    # ===================================================================
    "GNRHR": {
        "gene_id": 2798,
        "chrom": "chr4",
        "start": 68_281_295,
        "end": 68_298_781,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4q13.2",
        "exon_count": 3,
        "role": (
            "GnRH receptor -- G-protein coupled receptor on pituitary "
            "gonadotrophs, mediates GnRH stimulation of LH/FSH secretion.  "
            "Biallelic mutations cause normosmic hypogonadotropic "
            "hypogonadism (nHH, no anosmia)."
        ),
        "disease": "Normosmic hypogonadotropic hypogonadism",
        "omim_disease": 146110,
        "omim_gene": 138850,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gln106Arg (Q106R)",
                "rsid": "rs104893836",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common mutation; reduces GnRH binding affinity.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.0 kb, ideal for AAV) to "
            "pituitary gonadotrophs; "
            "(2) Pulsatile GnRH therapy highly effective without gene therapy."
        ),
        "clinical_programs": "No gene therapy trials. Pulsatile GnRH SOC.",
        "conditions": ["hypogonadotropic_hypogonadism", "nHH", "fertility",
                        "reproductive"],
    },
}
