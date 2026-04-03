"""
Congenital disorders of glycosylation (CDG) targets -- GRCh38 (hg38).

CDGs are multisystem metabolic disorders caused by defects in protein or
lipid glycosylation.  PMM2-CDG (CDG-Ia) is by far the most common,
accounting for ~80% of all diagnosed CDG.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


GLYCOSYLATION_TARGETS = {

    "PMM2": {
        "gene_id": 5373,
        "chrom": "chr16",
        "start": 8_891_670,
        "end": 8_943_013,
        "strand": "+",
        "refseq": "NC_000016.10",
        "cytoband": "16p13.2",
        "exon_count": 8,
        "role": (
            "Phosphomannomutase 2 -- cytosolic enzyme converting mannose-6-"
            "phosphate to mannose-1-phosphate, essential for dolichol-linked "
            "oligosaccharide (LLO) biosynthesis and N-glycosylation.  "
            "Deficiency causes PMM2-CDG (CDG-Ia), the most common CDG (~80% "
            "of cases): cerebellar hypoplasia, strabismus, inverted nipples, "
            "lipodystrophy, coagulopathy, hepatic dysfunction, intellectual "
            "disability.  ~1:20,000 incidence."
        ),
        "disease": "PMM2-CDG (Congenital disorder of glycosylation type Ia)",
        "omim_disease": 212065,
        "omim_gene": 601785,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg141His (R141H)",
                "rsid": "rs28936415",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common PMM2 mutation (~40% of European alleles).  "
                    "Near-null activity.  NEVER found homozygous (lethal) -- "
                    "always compound heterozygous with milder allele."
                ),
            },
            {
                "name": "p.Phe119Leu (F119L)",
                "rsid": "rs28936416",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Second most common; mild, residual activity ~20%.",
            },
            {
                "name": "p.Asp65Tyr (D65Y)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Active site adjacent; compound het with R141H = mild-moderate.",
            },
        ],
        "strategy": (
            "(1) AAV liver-directed gene replacement (cDNA ~0.7 kb, ideal "
            "for AAV) -- restoring hepatic PMM2 corrects coagulopathy and "
            "hepatic dysfunction; "
            "(2) AAV9 IV for systemic correction (CNS + liver); "
            "(3) Base editing for R141H (single nucleotide change); "
            "(4) Mannose-1-phosphate supplementation (bypasses PMM2 block "
            "but has delivery challenges due to phosphatase degradation); "
            "(5) Pharmacological chaperones to stabilize residual enzyme.  "
            "Delivery: AAV8 liver-directed + AAV9 IV for CNS."
        ),
        "clinical_programs": (
            "Magic Medicines AAV9-PMM2 gene therapy Phase 1 planned.  "
            "Acetazolamide showed improvement in cerebellar symptoms in "
            "small trial (Martina et al., JIMD 2023).  "
            "No approved therapies.  D-mannose supplementation does NOT "
            "work (wrong pathway)."
        ),
        "conditions": ["PMM2_CDG", "CDG_Ia", "congenital_disorder_of_glycosylation",
                        "glycosylation", "cerebellar_hypoplasia"],
    },

    "ALG1": {
        "gene_id": 56052,
        "chrom": "chr16",
        "start": 5_069_770,
        "end": 5_104_063,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16p13.3",
        "exon_count": 14,
        "role": (
            "ALG1 mannosyltransferase -- ER enzyme adding first mannose to "
            "dolichol-PP-GlcNAc2 in LLO biosynthesis.  Deficiency causes "
            "ALG1-CDG (CDG-Ik): severe epileptic encephalopathy, microcephaly, "
            "cardiomyopathy, coagulopathy."
        ),
        "disease": "ALG1-CDG (CDG-Ik)",
        "omim_disease": 608540,
        "omim_gene": 605907,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Ser258Leu (S258L)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Active site; severe neurological phenotype.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.3 kb, ideal for AAV); "
            "(2) Base editing for specific mutations."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["ALG1_CDG", "congenital_disorder_of_glycosylation",
                        "glycosylation", "epileptic_encephalopathy"],
    },

    "ALG6": {
        "gene_id": 29929,
        "chrom": "chr1",
        "start": 63_373_832,
        "end": 63_438_437,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1p31.3",
        "exon_count": 14,
        "role": (
            "ALG6 glucosyltransferase -- ER enzyme adding first glucose to "
            "LLO.  Deficiency causes ALG6-CDG (CDG-Ic): moderate ID, "
            "seizures, protein-losing enteropathy.  Second most common CDG."
        ),
        "disease": "ALG6-CDG (CDG-Ic, second most common CDG)",
        "omim_disease": 603147,
        "omim_gene": 604566,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.998C>T (p.Ala333Val, A333V)",
                "rsid": "rs121908413",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common ALG6 mutation (~50% of alleles).",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.5 kb, fits AAV); "
            "(2) Base editing for A333V."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["ALG6_CDG", "CDG_Ic", "congenital_disorder_of_glycosylation",
                        "glycosylation"],
    },

    "DPAGT1": {
        "gene_id": 1798,
        "chrom": "chr11",
        "start": 118_976_665,
        "end": 118_987_067,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11q23.3",
        "exon_count": 9,
        "role": (
            "Dolichyl-phosphate N-acetylglucosaminephosphotransferase 1 -- "
            "first enzyme in LLO biosynthesis (GlcNAc-1-P transfer to "
            "dolichol-P).  Deficiency causes DPAGT1-CDG: severe ID, seizures, "
            "microcephaly.  Also causes limb-girdle congenital myasthenic "
            "syndrome (CMS with tubular aggregates)."
        ),
        "disease": "DPAGT1-CDG / Congenital myasthenic syndrome",
        "omim_disease": 608093,
        "omim_gene": 191350,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various missense",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Hypomorphic alleles; complete LOF likely lethal.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.2 kb, ideal for AAV); "
            "(2) 3,4-DAP for myasthenic component."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["DPAGT1_CDG", "congenital_disorder_of_glycosylation",
                        "glycosylation", "congenital_myasthenic_syndrome"],
    },

    "SRD5A3": {
        "gene_id": 79644,
        "chrom": "chr4",
        "start": 56_178_449,
        "end": 56_206_178,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4q12",
        "exon_count": 6,
        "role": (
            "Steroid 5-alpha-reductase type 3 -- polyprenol reductase "
            "converting polyprenol to dolichol (essential for LLO synthesis).  "
            "Deficiency causes SRD5A3-CDG: cerebellar ataxia, ophthalmologic "
            "features (optic atrophy, coloboma), ID, coagulopathy."
        ),
        "disease": "SRD5A3-CDG",
        "omim_disease": 612379,
        "omim_gene": 611715,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Trp29Ter (W29X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Founder in Middle Eastern populations.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.0 kb, ideal for AAV); "
            "(2) Dolichol supplementation (theoretical bypass)."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["SRD5A3_CDG", "congenital_disorder_of_glycosylation",
                        "glycosylation", "cerebellar_ataxia"],
    },
}
