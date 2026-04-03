"""
Porphyria gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers acute hepatic porphyrias and erythropoietic porphyrias.
Givosiran (Alnylam) is FDA-approved RNAi therapy for AIP.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


PORPHYRIA_TARGETS = {

    "HMBS": {
        "gene_id": 3145,
        "chrom": "chr11",
        "start": 119_084_875,
        "end": 119_093_270,
        "strand": "+",
        "refseq": "NC_000011.10",
        "cytoband": "11q23.3",
        "exon_count": 15,
        "role": (
            "Hydroxymethylbilane synthase (porphobilinogen deaminase) -- "
            "third enzyme in heme biosynthesis.  Haploinsufficiency causes "
            "acute intermittent porphyria (AIP): episodic neurovisceral "
            "attacks (severe abdominal pain, neuropsychiatric symptoms, "
            "autonomic dysfunction) triggered by cytochrome P450 inducers "
            "(drugs, hormones, fasting)."
        ),
        "disease": "Acute intermittent porphyria (AIP)",
        "omim_disease": 176000,
        "omim_gene": 609806,
        "inheritance": "AD (low penetrance ~10-20%)",
        "key_variants": [
            {
                "name": "p.Arg167Trp (R167W)",
                "rsid": "rs121918009",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common in Swedish population; ~50% enzyme activity.",
            },
            {
                "name": "p.Arg173Trp (R173W)",
                "rsid": "rs121918010",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Active site; severe enzyme deficiency.",
            },
        ],
        "strategy": (
            "(1) AAV5 liver-directed gene replacement (cDNA ~1.0 kb, ideal "
            "for AAV) -- Alnylam/Ultragenyx approach; "
            "(2) ALAS1 silencing (givosiran/Givlaari approach) prevents "
            "ALA/PBG accumulation without fixing underlying defect; "
            "(3) mRNA/LNP hepatic delivery for enzyme replacement; "
            "(4) CRISPRa to upregulate wild-type HMBS allele.  "
            "Delivery: AAV5 IV to liver."
        ),
        "clinical_programs": (
            "APPROVED: Givosiran (Givlaari, Alnylam) -- RNAi targeting ALAS1 "
            "in liver, FDA-approved 2019 for AIP.  Prevents attacks by "
            "reducing ALA/PBG (treats downstream, not root cause).  "
            "Ultragenyx UX701 (AAV5-HMBS) Phase 1/2 for AIP -- addresses "
            "root cause.  Gene therapy could be curative."
        ),
        "conditions": ["acute_intermittent_porphyria", "AIP", "porphyria",
                        "heme_biosynthesis", "neurovisceral_attack"],
    },

    "FECH": {
        "gene_id": 2235,
        "chrom": "chr18",
        "start": 57_569_003,
        "end": 57_609_857,
        "strand": "+",
        "refseq": "NC_000018.10",
        "cytoband": "18q21.31",
        "exon_count": 11,
        "role": (
            "Ferrochelatase -- terminal enzyme in heme biosynthesis (inserts "
            "Fe2+ into protoporphyrin IX to form heme).  Deficiency causes "
            "erythropoietic protoporphyria (EPP): severe phototoxicity from "
            "protoporphyrin IX accumulation in skin/erythrocytes.  Painful "
            "burning/stinging on sun exposure, hepatic protoporphyrin "
            "accumulation can cause liver failure."
        ),
        "disease": "Erythropoietic protoporphyria (EPP)",
        "omim_disease": 177000,
        "omim_gene": 612386,
        "inheritance": "AD (with low-expression IVS3-48C allele in trans)",
        "key_variants": [
            {
                "name": "IVS3-48T>C (low-expression allele)",
                "rsid": "rs2272783",
                "consequence": "splice_regulatory",
                "clinical_significance": "modifier",
                "notes": (
                    "Common hypomorphic allele (~10% Europeans).  EPP requires "
                    "a LOF mutation in trans with this low-expression allele."
                ),
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.3 kb, ideal for AAV) in "
            "erythroid precursors; "
            "(2) Ex vivo CRISPR correction in HSCs; "
            "(3) Afamelanotide (alpha-MSH analogue) approved for EPP "
            "(increases melanin, not gene therapy).  "
            "Delivery: lentiviral HSC or AAV8 liver for hepatic form."
        ),
        "clinical_programs": (
            "Afamelanotide (Scenesse) approved in EU/US for EPP (symptomatic).  "
            "No gene therapy trials as of 2026."
        ),
        "conditions": ["erythropoietic_protoporphyria", "EPP", "porphyria",
                        "photosensitivity", "heme_biosynthesis"],
    },

    "UROS": {
        "gene_id": 7390,
        "chrom": "chr10",
        "start": 125_497_415,
        "end": 125_532_990,
        "strand": "-",
        "refseq": "NC_000010.11",
        "cytoband": "10q26.2",
        "exon_count": 10,
        "role": (
            "Uroporphyrinogen III synthase -- fourth enzyme in heme pathway.  "
            "Deficiency causes congenital erythropoietic porphyria (CEP, "
            "Gunther disease): severe photomutilation, hemolytic anemia, "
            "red fluorescent teeth (erythrodontia), red urine from birth."
        ),
        "disease": "Congenital erythropoietic porphyria (CEP, Gunther disease)",
        "omim_disease": 263700,
        "omim_gene": 606938,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Cys73Arg (C73R)",
                "rsid": "rs104893753",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common (~33% of alleles); severe form.",
            },
        ],
        "strategy": (
            "(1) Ex vivo lentiviral HSC gene therapy (cDNA ~0.8 kb, ideal); "
            "(2) HSCT is curative (corrects erythroid defect); "
            "(3) AAV-directed erythroid expression.  "
            "Delivery: lentiviral HSC."
        ),
        "clinical_programs": (
            "HSCT curative for severe CEP.  "
            "No gene therapy trials as of 2026."
        ),
        "conditions": ["congenital_erythropoietic_porphyria", "CEP",
                        "gunther_disease", "porphyria", "photosensitivity"],
    },

    "PPOX": {
        "gene_id": 5498,
        "chrom": "chr1",
        "start": 160_098_038,
        "end": 160_103_312,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q23.3",
        "exon_count": 13,
        "role": (
            "Protoporphyrinogen oxidase -- sixth enzyme in heme biosynthesis.  "
            "Deficiency causes variegate porphyria (VP): neurovisceral attacks "
            "(like AIP) PLUS photocutaneous lesions.  Very common in South "
            "African Afrikaners (R59W founder, ~1:300)."
        ),
        "disease": "Variegate porphyria (VP)",
        "omim_disease": 176200,
        "omim_gene": 600923,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Arg59Trp (R59W, South African founder)",
                "rsid": "rs121908009",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Afrikaner founder (traced to 1688 settler); ~1:300 in SA.",
            },
        ],
        "strategy": (
            "(1) Givosiran (ALAS1 silencing) used off-label for VP attacks; "
            "(2) AAV liver gene replacement (cDNA ~1.5 kb, fits AAV); "
            "(3) CRISPRa to upregulate wild-type allele."
        ),
        "clinical_programs": "Givosiran off-label. No specific VP gene therapy trials.",
        "conditions": ["variegate_porphyria", "VP", "porphyria",
                        "neurovisceral_attack", "photosensitivity"],
    },

    "ALAS2": {
        "gene_id": 212,
        "chrom": "chrX",
        "start": 55_054_572,
        "end": 55_076_521,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xp11.21",
        "exon_count": 12,
        "role": (
            "5-Aminolevulinate synthase 2 -- erythroid-specific first enzyme "
            "of heme biosynthesis.  Loss-of-function causes X-linked "
            "sideroblastic anemia (XLSA): ring sideroblasts, iron overload.  "
            "Gain-of-function causes X-linked erythropoietic protoporphyria "
            "(XLEPP)."
        ),
        "disease": "X-linked sideroblastic anemia / XLEPP",
        "omim_disease": 300751,
        "omim_gene": 301300,
        "inheritance": "XL",
        "key_variants": [
            {
                "name": "Various missense (pyridoxine-responsive vs non-responsive)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Some respond to pyridoxine (vitamin B6) supplementation.",
            },
        ],
        "strategy": (
            "(1) Pyridoxine (B6) is effective for responsive genotypes; "
            "(2) AAV gene replacement (cDNA ~1.7 kb, fits AAV) in erythroid "
            "precursors for non-responsive forms.  "
            "Delivery: lentiviral HSC."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["sideroblastic_anemia", "XLSA", "XLEPP", "porphyria",
                        "iron_overload"],
    },
}
