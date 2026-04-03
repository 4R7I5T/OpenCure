"""
Lymphatic disease gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers hereditary primary lymphedema genes.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


LYMPHATIC_TARGETS = {

    "FLT4": {
        "gene_id": 2324,
        "chrom": "chr5",
        "start": 180_028_876,
        "end": 180_076_624,
        "strand": "+",
        "refseq": "NC_000005.10",
        "cytoband": "5q35.3",
        "exon_count": 31,
        "role": (
            "VEGFR-3 (vascular endothelial growth factor receptor 3) -- "
            "receptor tyrosine kinase essential for lymphangiogenesis.  "
            "Missense mutations in the kinase domain cause Milroy disease "
            "(hereditary lymphedema type IA): congenital lower-extremity "
            "lymphedema.  Most common genetic cause of primary lymphedema."
        ),
        "disease": "Milroy disease (hereditary lymphedema type IA)",
        "omim_disease": 153100,
        "omim_gene": 136352,
        "inheritance": "AD (incomplete penetrance)",
        "key_variants": [
            {
                "name": "Various kinase domain missense mutations",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Cluster in tyrosine kinase domain; impair signaling.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement with functional VEGFR-3 (cDNA ~4.1 kb, "
            "fits AAV with compact promoter); "
            "(2) CRISPRa to upregulate wild-type allele; "
            "(3) VEGF-C protein therapy as adjunct.  "
            "Delivery: lymphatic-tropic AAV or direct lymph node injection."
        ),
        "clinical_programs": (
            "Herantis VEGF-C gene therapy for lymphedema Phase 1 (Finland).  "
            "No CRISPR trials."
        ),
        "conditions": ["milroy_disease", "primary_lymphedema", "lymphatic",
                        "lymphedema", "congenital_lymphedema"],
    },

    "FOXC2": {
        "gene_id": 2303,
        "chrom": "chr16",
        "start": 86_567_983,
        "end": 86_569_817,
        "strand": "+",
        "refseq": "NC_000016.10",
        "cytoband": "16q24.1",
        "exon_count": 1,
        "role": (
            "Forkhead box C2 -- transcription factor essential for lymphatic "
            "valve development.  Haploinsufficiency causes lymphedema-"
            "distichiasis syndrome: late-onset lymphedema + double row of "
            "eyelashes + varicose veins."
        ),
        "disease": "Lymphedema-distichiasis syndrome",
        "omim_disease": 153400,
        "omim_gene": 602402,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Truncating mutations throughout single exon",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Single-exon gene; haploinsufficiency mechanism.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.5 kb, ideal for AAV); "
            "(2) CRISPRa to upregulate wild-type allele."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["lymphedema_distichiasis", "primary_lymphedema",
                        "lymphatic", "lymphedema"],
    },

    "GJC2": {
        "gene_id": 57165,
        "chrom": "chr1",
        "start": 228_322_093,
        "end": 228_331_592,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1q42.13",
        "exon_count": 2,
        "role": (
            "Connexin 47 (GJC2) -- gap junction protein in lymphatic "
            "endothelium and oligodendrocytes.  Mutations cause hereditary "
            "lymphedema type IC or Pelizaeus-Merzbacher-like disease (CNS)."
        ),
        "disease": "Hereditary lymphedema type IC",
        "omim_disease": 613480,
        "omim_gene": 608803,
        "inheritance": "AD (lymphedema) / AR (PMLD)",
        "key_variants": [
            {
                "name": "c.-167G>T (5'UTR, lymphedema-specific)",
                "rsid": "rs587776654",
                "consequence": "regulatory",
                "clinical_significance": "pathogenic",
                "notes": "Founder variant; creates upstream ORF reducing translation.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.3 kb, ideal for AAV); "
            "(2) Base editing to correct 5'UTR variant."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["primary_lymphedema", "lymphatic", "lymphedema"],
    },

    "CCBE1": {
        "gene_id": 147372,
        "chrom": "chr18",
        "start": 59_418_454,
        "end": 59_481_779,
        "strand": "+",
        "refseq": "NC_000018.10",
        "cytoband": "18q21.32",
        "exon_count": 12,
        "role": (
            "Collagen and calcium-binding EGF domain-containing protein 1 -- "
            "essential for VEGF-C processing and lymphangiogenesis.  Biallelic "
            "loss causes Hennekam lymphangiectasia-lymphedema syndrome: "
            "severe lymphedema, intestinal lymphangiectasia, intellectual "
            "disability, facial dysmorphism."
        ),
        "disease": "Hennekam syndrome (lymphangiectasia-lymphedema)",
        "omim_disease": 235510,
        "omim_gene": 612753,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Cys174Arg (C174R)",
                "rsid": "rs398124411",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts EGF domain; impairs VEGF-C activation.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.2 kb, ideal for AAV); "
            "(2) Base editing for C174R."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["hennekam_syndrome", "lymphangiectasia", "lymphatic",
                        "lymphedema"],
    },

    "PIEZO1": {
        "gene_id": 9780,
        "chrom": "chr16",
        "start": 88_713_548,
        "end": 88_792_021,
        "strand": "+",
        "refseq": "NC_000016.10",
        "cytoband": "16q24.3",
        "exon_count": 51,
        "role": (
            "PIEZO1 -- mechanosensitive cation channel.  Gain-of-function "
            "causes dehydrated hereditary stomatocytosis (DHS).  Biallelic "
            "loss-of-function causes generalized lymphatic dysplasia "
            "(nonimmune hydrops fetalis in severe cases)."
        ),
        "disease": "Generalized lymphatic dysplasia / DHS",
        "omim_disease": 616843,
        "omim_gene": 611184,
        "inheritance": "AR (lymphatic) / AD (DHS, GOF)",
        "key_variants": [
            {
                "name": "Various LOF mutations (lymphatic)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Biallelic LOF -> lymphatic dysplasia, hydrops fetalis.",
            },
        ],
        "strategy": (
            "(1) Gene too large for single AAV (~7.6 kb cDNA); dual-AAV or "
            "mRNA approach; "
            "(2) CRISPRa to upregulate residual expression in hypomorphic "
            "cases."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["lymphatic_dysplasia", "lymphatic", "hydrops_fetalis",
                        "stomatocytosis"],
    },
}
