"""
Congenital myopathy gene-therapy targets -- GRCh38 (hg38).

Covers nemaline myopathy, centronuclear myopathy, and related structural
congenital myopathies NOT already covered by the neuromuscular pipeline
(which covers LGMD, MTM1) or skeletal channelopathy (RYR1).

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


CONGENITAL_MYOPATHY_TARGETS = {

    "NEB": {
        "gene_id": 4703,
        "chrom": "chr2",
        "start": 151_486_799,
        "end": 151_745_202,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q23.3",
        "exon_count": 183,
        "role": (
            "Nebulin -- giant sarcomeric protein (~800 kDa) that spans the "
            "thin filament length in skeletal muscle, regulating actin filament "
            "length and cross-bridge cycling.  Biallelic mutations cause "
            "typical nemaline myopathy (NEM2), the most common form of "
            "nemaline myopathy (~50% of all NM).  Neonatal hypotonia, "
            "proximal weakness, respiratory compromise, nemaline rods on "
            "biopsy."
        ),
        "disease": "Nemaline myopathy type 2 (NEM2, most common NM)",
        "omim_disease": 256030,
        "omim_gene": 161650,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Exon 55 deletion (Ashkenazi Jewish founder)",
                "rsid": None,
                "consequence": "in_frame_deletion",
                "clinical_significance": "pathogenic",
                "notes": "Carrier rate ~1/108 in Ashkenazi Jewish; typical NM.",
            },
            {
                "name": ">200 pathogenic variants (extreme heterogeneity)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Largest exon count of any human gene (183 exons, "
                    "~26 kb cDNA).  Compound heterozygosity with missense "
                    "allele = milder; two null alleles = severe."
                ),
            },
        ],
        "strategy": (
            "(1) Gene far too large for AAV (~26 kb cDNA, 183 exons); "
            "(2) Exon skipping for specific exon deletions (analogous to "
            "DMD exon skipping); "
            "(3) CRISPRa to upregulate residual nebulin from hypomorphic "
            "alleles; "
            "(4) Troponin activators (CK-2127107, reldesemtiv) to enhance "
            "residual thin filament function.  "
            "Delivery: ASO (intramuscular/systemic) for exon skipping."
        ),
        "clinical_programs": (
            "No gene therapy trials as of 2026.  "
            "Reldesemtiv (Cytokinetics) troponin activator Phase 2/3 for "
            "NM (mechanism bypass, not gene correction).  "
            "NEB exon skipping preclinical."
        ),
        "conditions": ["nemaline_myopathy", "NEM2", "congenital_myopathy",
                        "neonatal_hypotonia", "nemaline_rods"],
    },

    "ACTA1": {
        "gene_id": 58,
        "chrom": "chr1",
        "start": 229_565_921,
        "end": 229_569_156,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q42.13",
        "exon_count": 7,
        "role": (
            "Skeletal alpha-actin -- major actin isoform in postnatal skeletal "
            "muscle thin filaments.  Mutations cause nemaline myopathy type 3 "
            "(NEM3), actin myopathy, intranuclear rod myopathy, and cap "
            "myopathy.  Most mutations are de novo dominant (severe neonatal) "
            "or recessive (milder)."
        ),
        "disease": "Nemaline myopathy type 3 / Actin myopathy",
        "omim_disease": 161800,
        "omim_gene": 102610,
        "inheritance": "AD (usually de novo) / AR",
        "key_variants": [
            {
                "name": "Various missense (dominant-negative)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    ">200 mutations described.  Dominant mutations poison "
                    "the thin filament.  Genotype-phenotype correlations "
                    "are poor."
                ),
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant mutant ACTA1 allele "
            "(cardiac alpha-actin ACTC1 can compensate in skeletal muscle); "
            "(2) CRISPRa to upregulate ACTC1 (fetal cardiac actin) as "
            "surrogate replacement; "
            "(3) AAV-mediated ACTC1 overexpression rescued ACTA1-null "
            "mice (Ravenscroft et al.).  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "ACTC1 replacement strategy validated in ACTA1-null mice "
            "(Ravenscroft et al., Hum Mol Genet 2011)."
        ),
        "conditions": ["nemaline_myopathy", "NEM3", "actin_myopathy",
                        "congenital_myopathy", "neonatal_hypotonia"],
    },

    "DNM2": {
        "gene_id": 1785,
        "chrom": "chr19",
        "start": 10_677_857,
        "end": 10_729_459,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19p13.2",
        "exon_count": 22,
        "role": (
            "Dynamin 2 -- large GTPase involved in membrane trafficking, "
            "endocytosis, and cytoskeleton remodeling.  Dominant mutations "
            "cause autosomal dominant centronuclear myopathy (CNM): "
            "progressive proximal weakness, ptosis, ophthalmoplegia, "
            "centrally located nuclei on biopsy.  Also causes CMT2M and "
            "dominant intermediate CMT (CMTDIB)."
        ),
        "disease": "Autosomal dominant centronuclear myopathy / CMT2M",
        "omim_disease": 160150,
        "omim_gene": 602556,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Arg465Trp (R465W)",
                "rsid": "rs121918308",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "PH domain; most common CNM mutation; neonatal/childhood onset.",
            },
            {
                "name": "p.Glu368Lys (E368K)",
                "rsid": "rs121918307",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Middle domain; milder, later onset.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant mutant allele; "
            "(2) DNM2 reduction (even wild-type allele) rescues phenotype "
            "in animal models -- Dynacure DYN101 ASO approach; "
            "(3) Base editing for R465W.  "
            "Delivery: subcutaneous ASO or AAV9."
        ),
        "clinical_programs": (
            "Dynacure DYN101 (anti-DNM2 ASO) Phase 1/2 for XLMTM and "
            "centronuclear myopathies (NCT04033159).  Reducing DNM2 "
            "compensates for MTM1 loss (cross-rescue)."
        ),
        "conditions": ["centronuclear_myopathy", "CNM", "congenital_myopathy",
                        "CMT2M"],
    },

    "BIN1": {
        "gene_id": 274,
        "chrom": "chr2",
        "start": 127_048_052,
        "end": 127_107_286,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q14.3",
        "exon_count": 20,
        "role": (
            "Amphiphysin 2 (BIN1/SH3P9) -- BAR domain protein involved in "
            "T-tubule biogenesis, endocytosis, and membrane remodeling.  "
            "Recessive mutations cause autosomal recessive centronuclear "
            "myopathy: rapidly progressive, severe, often fatal in infancy."
        ),
        "disease": "Autosomal recessive centronuclear myopathy",
        "omim_disease": 255200,
        "omim_gene": 601248,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various truncating and missense",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Rare; ~30 families reported.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.7 kb, fits AAV); "
            "(2) AAV-BIN1 restored T-tubule structure in BIN1-null mice.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": "Preclinical AAV-BIN1 gene therapy.",
        "conditions": ["centronuclear_myopathy", "CNM", "congenital_myopathy",
                        "neonatal_hypotonia"],
    },

    "TPM3": {
        "gene_id": 7170,
        "chrom": "chr1",
        "start": 154_127_462,
        "end": 154_165_104,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1q21.3",
        "exon_count": 13,
        "role": (
            "Tropomyosin 3 -- slow skeletal muscle tropomyosin, regulates "
            "actin-myosin interaction.  Mutations cause nemaline myopathy "
            "type 1 (NEM1, AD) and cap myopathy.  The first NM gene "
            "identified."
        ),
        "disease": "Nemaline myopathy type 1 (NEM1) / Cap myopathy",
        "omim_disease": 609284,
        "omim_gene": 191030,
        "inheritance": "AD / AR",
        "key_variants": [
            {
                "name": "p.Met9Arg (M9R)",
                "rsid": "rs121918254",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "N-terminal; dominant; severe neonatal form.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant mutant allele; "
            "(2) AAV gene replacement (cDNA ~0.8 kb, ideal for AAV) for "
            "recessive forms."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["nemaline_myopathy", "NEM1", "cap_myopathy",
                        "congenital_myopathy"],
    },
}
