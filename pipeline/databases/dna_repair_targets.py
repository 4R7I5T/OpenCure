"""
DNA repair disorder gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers xeroderma pigmentosum (XP), Cockayne syndrome (CS), and
trichothiodystrophy (TTD) -- nucleotide excision repair (NER) defects
causing UV sensitivity, neurodegeneration, and/or cancer predisposition.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


DNA_REPAIR_TARGETS = {

    "XPA": {
        "gene_id": 7507,
        "chrom": "chr9",
        "start": 97_656_217,
        "end": 97_682_172,
        "strand": "+",
        "refseq": "NC_000009.12",
        "cytoband": "9q22.33",
        "exon_count": 6,
        "role": (
            "Xeroderma pigmentosum group A -- core NER scaffolding protein "
            "that verifies DNA damage and positions the repair complex.  "
            "Deficiency causes XP-A, the most severe XP group: extreme UV "
            "sensitivity, >1000x skin cancer risk, progressive neurodegeneration "
            "(cerebellar ataxia, sensorineural deafness, cognitive decline).  "
            "Most common XP type in Japan (founder IVS3-1G>C)."
        ),
        "disease": "Xeroderma pigmentosum group A (XP-A)",
        "omim_disease": 278700,
        "omim_gene": 611153,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "IVS3-1G>C (Japanese founder)",
                "rsid": "rs104894143",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common XP mutation worldwide due to Japanese "
                    "founder effect.  Causes exon 3 skipping -> no XPA protein."
                ),
            },
            {
                "name": "p.Arg228Ter (R228X)",
                "rsid": "rs104894144",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "North African founder; severe XP-A with neurodegeneration.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~0.8 kb, ideal for AAV) to "
            "skin keratinocytes (topical AAV or ex vivo skin grafts); "
            "(2) AAV9 for neurodegeneration component (intrathecal); "
            "(3) CRISPR correction in patient keratinocyte stem cells for "
            "autologous skin grafting; "
            "(4) Strict UV protection remains cornerstone.  "
            "Delivery: topical AAV to skin + AAV9 IT for CNS."
        ),
        "clinical_programs": (
            "No gene therapy trials as of 2026.  "
            "Rigorous UV avoidance is sole management.  "
            "Preclinical AAV-XPA skin-directed gene therapy."
        ),
        "conditions": ["xeroderma_pigmentosum", "XP_A", "dna_repair",
                        "UV_sensitivity", "skin_cancer_predisposition",
                        "neurodegeneration"],
    },

    "XPC": {
        "gene_id": 7508,
        "chrom": "chr3",
        "start": 14_145_655,
        "end": 14_179_261,
        "strand": "-",
        "refseq": "NC_000003.12",
        "cytoband": "3p25.1",
        "exon_count": 16,
        "role": (
            "Xeroderma pigmentosum group C -- global genome NER (GG-NER) "
            "damage sensor, forms XPC-RAD23B-CETN2 complex.  Deficiency "
            "causes XP-C: high skin cancer risk but usually NO "
            "neurodegeneration (GG-NER specific, transcription-coupled NER "
            "intact).  Most common XP group in Europe/North Africa."
        ),
        "disease": "Xeroderma pigmentosum group C (XP-C)",
        "omim_disease": 278720,
        "omim_gene": 613208,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1643_1644delTG (p.Val548Alafs)",
                "rsid": "rs121913257",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Most common XP-C mutation in North Africa.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement to skin (cDNA ~2.8 kb fits AAV); "
            "(2) Ex vivo CRISPR correction in keratinocyte stem cells; "
            "(3) Topical gene therapy application.  "
            "Delivery: intradermal AAV or gene-corrected skin grafts."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["xeroderma_pigmentosum", "XP_C", "dna_repair",
                        "skin_cancer_predisposition"],
    },

    "ERCC6": {
        "gene_id": 2074,
        "chrom": "chr10",
        "start": 49_453_767,
        "end": 49_539_115,
        "strand": "-",
        "refseq": "NC_000010.11",
        "cytoband": "10q11.23",
        "exon_count": 21,
        "role": (
            "Excision repair cross-complementation group 6 (CSB) -- "
            "transcription-coupled NER (TC-NER) factor, SWI/SNF ATPase "
            "that remodels chromatin at stalled RNA Pol II.  Mutations cause "
            "Cockayne syndrome type B: dwarfism, microcephaly, progressive "
            "neurodegeneration (leukodystrophy, cerebellar atrophy), retinal "
            "degeneration, photosensitivity, progeroid features.  NO cancer "
            "predisposition (unlike XP)."
        ),
        "disease": "Cockayne syndrome type B (CS-B)",
        "omim_disease": 133540,
        "omim_gene": 609413,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various truncating mutations (no hotspot)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "~70% of Cockayne syndrome; most are ERCC6.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~4.5 kb, tight AAV fit); "
            "(2) AAV9 for neurological component; "
            "(3) CRISPR correction in patient fibroblasts (research).  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["cockayne_syndrome", "CS_B", "dna_repair",
                        "photosensitivity", "neurodegeneration", "progeroid"],
    },

    "ERCC2": {
        "gene_id": 2068,
        "chrom": "chr19",
        "start": 45_349_072,
        "end": 45_371_522,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19q13.32",
        "exon_count": 23,
        "role": (
            "Excision repair cross-complementation group 2 (XPD) -- TFIIH "
            "helicase subunit, essential for both NER and transcription "
            "initiation.  Different mutations cause three distinct diseases: "
            "XP-D (cancer-prone), trichothiodystrophy (TTD, brittle hair, "
            "ichthyosis, NO cancer), or XP/CS complex."
        ),
        "disease": "XP-D / Trichothiodystrophy / XP-CS complex",
        "omim_disease": 278730,
        "omim_gene": 126340,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg722Trp (R722W, TTD-specific)",
                "rsid": "rs121913019",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "C-terminal; causes TTD (NOT XP/cancer).",
            },
            {
                "name": "p.Arg683Trp (R683W, XP-specific)",
                "rsid": "rs121913020",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Helicase domain; causes classic XP-D with cancer.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.3 kb fits AAV); "
            "(2) Base editing for specific genotype-phenotype mutations; "
            "(3) UV protection for XP-D.  "
            "Delivery: AAV to skin + AAV9 for CNS."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["xeroderma_pigmentosum", "XP_D", "trichothiodystrophy",
                        "TTD", "dna_repair", "photosensitivity"],
    },

    "ERCC8": {
        "gene_id": 1161,
        "chrom": "chr5",
        "start": 60_870_677,
        "end": 60_944_099,
        "strand": "+",
        "refseq": "NC_000005.10",
        "cytoband": "5q12.1",
        "exon_count": 12,
        "role": (
            "Excision repair cross-complementation group 8 (CSA) -- WD40 "
            "repeat E3 ubiquitin ligase adaptor in TC-NER.  Mutations cause "
            "Cockayne syndrome type A: same phenotype as CS-B but ~30% of "
            "Cockayne syndrome."
        ),
        "disease": "Cockayne syndrome type A (CS-A)",
        "omim_disease": 216400,
        "omim_gene": 609412,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various truncating and missense",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "~30% of Cockayne syndrome.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.2 kb, ideal for AAV); "
            "(2) AAV9 for CNS component."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["cockayne_syndrome", "CS_A", "dna_repair",
                        "photosensitivity", "neurodegeneration"],
    },
}
