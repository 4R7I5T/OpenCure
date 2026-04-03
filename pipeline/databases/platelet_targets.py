"""
Platelet / megakaryocyte disorder gene-therapy targets -- GRCh38 (hg38).

Covers inherited bleeding disorders arising from platelet defects
(qualitative or quantitative).

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


PLATELET_TARGETS = {

    "ITGA2B": {
        "gene_id": 3674,
        "chrom": "chr17",
        "start": 44_343_998,
        "end": 44_361_651,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.31",
        "exon_count": 30,
        "role": (
            "Integrin alpha-IIb (GPIIb) -- platelet fibrinogen receptor "
            "alpha subunit, forms GPIIb/IIIa complex with ITGB3.  Biallelic "
            "mutations cause Glanzmann thrombasthenia (GT): absent platelet "
            "aggregation despite normal platelet count.  Severe mucocutaneous "
            "bleeding (epistaxis, menorrhagia, GI bleed)."
        ),
        "disease": "Glanzmann thrombasthenia (type I)",
        "omim_disease": 273800,
        "omim_gene": 607759,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "IVS15-1G>A (splice, Iraqi-Jewish founder)",
                "rsid": None,
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": "Common in Middle Eastern populations.",
            },
        ],
        "strategy": (
            "(1) Ex vivo lentiviral HSC gene therapy expressing ITGA2B in "
            "megakaryocytes (cDNA ~3.0 kb fits AAV/lenti); "
            "(2) AAV-mediated expression with megakaryocyte-specific "
            "promoter (GP1BA or PF4 promoter); "
            "(3) Platelet transfusion for acute bleeding.  "
            "Delivery: lentiviral HSC or AAV."
        ),
        "clinical_programs": (
            "Emicizumab (bypassing agent) off-label use reported.  "
            "No gene therapy trials as of 2026.  "
            "Preclinical lentiviral GT gene therapy (Paris group)."
        ),
        "conditions": ["glanzmann_thrombasthenia", "GT", "platelet_disorder",
                        "bleeding_disorder"],
    },

    "ITGB3": {
        "gene_id": 3690,
        "chrom": "chr17",
        "start": 47_253_827,
        "end": 47_313_743,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.32",
        "exon_count": 15,
        "role": (
            "Integrin beta-3 (GPIIIa) -- GPIIb/IIIa beta subunit.  "
            "Mutations cause Glanzmann thrombasthenia type II."
        ),
        "disease": "Glanzmann thrombasthenia (type II)",
        "omim_disease": 273800,
        "omim_gene": 173470,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various missense and truncating",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Less common than ITGA2B mutations.",
            },
        ],
        "strategy": (
            "(1) Ex vivo lentiviral HSC gene therapy (cDNA ~2.4 kb); "
            "(2) AAV-mediated expression with megakaryocyte promoter."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["glanzmann_thrombasthenia", "GT", "platelet_disorder",
                        "bleeding_disorder"],
    },

    "GP1BA": {
        "gene_id": 2811,
        "chrom": "chr17",
        "start": 4_932_494,
        "end": 4_935_544,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17p13.2",
        "exon_count": 2,
        "role": (
            "Glycoprotein Ib alpha -- von Willebrand factor receptor on "
            "platelets, essential for platelet adhesion to damaged vessel "
            "wall.  Biallelic mutations cause Bernard-Soulier syndrome (BSS): "
            "giant platelets, thrombocytopenia, severe bleeding."
        ),
        "disease": "Bernard-Soulier syndrome",
        "omim_disease": 231200,
        "omim_gene": 606672,
        "inheritance": "AR (classical); AD (some monoallelic variants)",
        "key_variants": [
            {
                "name": "Various truncating and missense",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Rare (~1:1,000,000); consanguineous families.",
            },
        ],
        "strategy": (
            "(1) Ex vivo lentiviral HSC gene therapy (cDNA ~1.9 kb); "
            "(2) AAV expression in megakaryocytes."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["bernard_soulier_syndrome", "BSS", "platelet_disorder",
                        "giant_platelets", "bleeding_disorder"],
    },

    "NBEAL2": {
        "gene_id": 23218,
        "chrom": "chr3",
        "start": 47_046_285,
        "end": 47_108_997,
        "strand": "+",
        "refseq": "NC_000003.12",
        "cytoband": "3p21.31",
        "exon_count": 54,
        "role": (
            "Neurobeachin-like 2 -- required for alpha-granule biogenesis "
            "in megakaryocytes.  Mutations cause gray platelet syndrome (GPS): "
            "large agranular platelets, mild thrombocytopenia, myelofibrosis.  "
            "Named for the gray appearance on Wright stain due to absent "
            "alpha-granules."
        ),
        "disease": "Gray platelet syndrome",
        "omim_disease": 139090,
        "omim_gene": 615285,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various truncating mutations",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Rare; ~100 cases reported worldwide.",
            },
        ],
        "strategy": (
            "(1) Gene large for AAV (~8.2 kb cDNA); dual-AAV or lentiviral; "
            "(2) Ex vivo CRISPR correction in HSCs."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["gray_platelet_syndrome", "GPS", "platelet_disorder",
                        "myelofibrosis"],
    },

    "HPS1": {
        "gene_id": 3257,
        "chrom": "chr10",
        "start": 98_379_652,
        "end": 98_404_960,
        "strand": "+",
        "refseq": "NC_000010.11",
        "cytoband": "10q24.2",
        "exon_count": 20,
        "role": (
            "Hermansky-Pudlak syndrome 1 -- BLOC-3 complex component required "
            "for lysosome-related organelle (LRO) biogenesis in melanocytes, "
            "platelets (dense granules), and lung (lamellar bodies).  HPS1: "
            "oculocutaneous albinism + platelet storage pool deficiency + "
            "pulmonary fibrosis (fatal, main cause of death in 4th-5th decade)."
        ),
        "disease": "Hermansky-Pudlak syndrome type 1",
        "omim_disease": 203300,
        "omim_gene": 604982,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1470_1486dup16 (Puerto Rican founder)",
                "rsid": "rs387906407",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "NW Puerto Rico founder; carrier rate ~1/21.  "
                    "~400+ HPS1 patients in PR."
                ),
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.1 kb, fits AAV); "
            "(2) Nintedanib (anti-fibrotic) for pulmonary fibrosis component "
            "(Hermansky-Pudlak Syndrome Network trial); "
            "(3) Ex vivo HSC gene therapy for platelet/melanocyte correction.  "
            "Delivery: AAV to lung (inhaled) + systemic."
        ),
        "clinical_programs": (
            "Nintedanib Phase 2/3 for HPS pulmonary fibrosis.  "
            "No gene therapy trials as of 2026."
        ),
        "conditions": ["hermansky_pudlak_syndrome", "HPS", "platelet_disorder",
                        "albinism", "pulmonary_fibrosis", "storage_pool_deficiency"],
    },
}
