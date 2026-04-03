"""
Organic acidemia / additional metabolic gene-therapy targets -- GRCh38 (hg38).

Covers methylmalonic acidemia (MMA), glutaric aciduria type I, alkaptonuria,
and VLCAD deficiency -- inborn errors of metabolism with active gene therapy
programs or clear CRISPR therapeutic strategies.

Note: Propionic acidemia (PCCA/PCCB), isovaleric acidemia (IVD), MCAD
(ACADM), CBS (homocystinuria), CPT2, and HADHA are already in the
rare_disease pipeline.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - LogicBio (now Ultragenyx), Selecta Biosciences MMA programs
"""


ORGANIC_ACIDEMIA_TARGETS = {

    "MUT": {
        "gene_id": 4594,
        "chrom": "chr6",
        "start": 49_398_691,
        "end": 49_431_282,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6p12.3",
        "exon_count": 13,
        "role": (
            "Methylmalonyl-CoA mutase -- mitochondrial enzyme converting "
            "L-methylmalonyl-CoA to succinyl-CoA (requires adenosylcobalamin "
            "B12 cofactor).  Deficiency causes isolated methylmalonic acidemia "
            "(MMA, mut0/mut- subtypes): life-threatening metabolic crises "
            "(hyperammonemia, metabolic acidosis, ketosis), progressive renal "
            "failure, basal ganglia necrosis.  Incidence ~1:50,000."
        ),
        "disease": "Methylmalonic acidemia (mut0/mut-)",
        "omim_disease": 251000,
        "omim_gene": 609058,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Asn219Tyr (N219Y)",
                "rsid": "rs121918260",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "mut- subtype; some residual activity; B12-responsive in some.",
            },
            {
                "name": "Various null mutations (mut0)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "No residual activity; severe neonatal presentation.",
            },
        ],
        "strategy": (
            "(1) AAV liver-directed gene replacement (cDNA ~2.3 kb fits AAV) "
            "-- hepatic MUT expression clears toxic methylmalonic acid; "
            "(2) LogicBio/Ultragenyx LB-001 AAV gene therapy for MMA; "
            "(3) mRNA/LNP hepatic delivery (Moderna approach for PA/MMA); "
            "(4) CRISPR insertion at albumin safe harbor locus.  "
            "Delivery: AAV8/AAV5 IV to liver."
        ),
        "clinical_programs": (
            "Ultragenyx (formerly LogicBio) LB-001 (AAVhu37-MUT) Phase 1/2 "
            "for MMA (NCT04581785).  "
            "Selecta Biosciences SEL-302 (AAV-MUT + ImmTOR immune "
            "tolerance) Phase 1.  "
            "Moderna mRNA-3705 (MUT mRNA in LNP) Phase 1/2 "
            "(NCT04899310).  "
            "Liver or combined liver-kidney transplant for severe mut0."
        ),
        "conditions": ["methylmalonic_acidemia", "MMA", "organic_acidemia",
                        "metabolic", "B12_metabolism"],
    },

    "MMAA": {
        "gene_id": 166785,
        "chrom": "chr4",
        "start": 145_578_937,
        "end": 145_606_741,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4q31.21",
        "exon_count": 7,
        "role": (
            "Methylmalonic aciduria type A protein -- mitochondrial GTPase "
            "required for adenosylcobalamin (AdoCbl) synthesis and MUT "
            "activation.  Deficiency causes cblA-type MMA: B12-responsive "
            "methylmalonic acidemia (milder than mut0/mut-)."
        ),
        "disease": "MMA cblA type (B12-responsive)",
        "omim_disease": 251100,
        "omim_gene": 607481,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.433C>T (p.Arg145Ter)",
                "rsid": "rs121918261",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Most common cblA mutation.",
            },
        ],
        "strategy": (
            "(1) High-dose hydroxycobalamin (IM B12) is first-line and "
            "effective for most cblA patients; "
            "(2) AAV gene replacement (cDNA ~1.2 kb, ideal for AAV) for "
            "B12-non-responsive patients."
        ),
        "clinical_programs": "B12 supplementation is standard of care. No gene therapy trials.",
        "conditions": ["methylmalonic_acidemia", "MMA", "cblA",
                        "organic_acidemia", "B12_metabolism"],
    },

    "GCDH": {
        "gene_id": 2639,
        "chrom": "chr19",
        "start": 12_890_581,
        "end": 12_898_471,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19p13.13",
        "exon_count": 11,
        "role": (
            "Glutaryl-CoA dehydrogenase -- mitochondrial enzyme in lysine, "
            "hydroxylysine, and tryptophan catabolism.  Deficiency causes "
            "glutaric aciduria type I (GA1): accumulation of glutaric acid "
            "and 3-hydroxyglutaric acid.  Encephalopathic crises cause "
            "irreversible striatal necrosis (dystonic cerebral palsy).  "
            "Newborn screening and early treatment prevent crises."
        ),
        "disease": "Glutaric aciduria type I (GA1)",
        "omim_disease": 231670,
        "omim_gene": 608801,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg402Trp (R402W)",
                "rsid": "rs121434369",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "European founder; high residual activity but still at risk for crisis.",
            },
            {
                "name": "p.Ala421Val (A421V)",
                "rsid": "rs121434370",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Old Order Amish founder; carrier rate ~1/10.",
            },
        ],
        "strategy": (
            "(1) AAV liver-directed gene replacement (cDNA ~1.3 kb, ideal "
            "for AAV) -- hepatic expression metabolises glutaryl-CoA; "
            "(2) Lysine-restricted diet + carnitine supplementation is "
            "effective if started neonatally (NBS critical); "
            "(3) Base editing for R402W or A421V."
        ),
        "clinical_programs": (
            "No gene therapy trials as of 2026.  "
            "Newborn screening + dietary management prevents crises "
            "in most patients if identified early."
        ),
        "conditions": ["glutaric_aciduria_type_I", "GA1", "organic_acidemia",
                        "metabolic", "striatal_necrosis"],
    },

    "HGD": {
        "gene_id": 3081,
        "chrom": "chr3",
        "start": 120_301_766,
        "end": 120_356_397,
        "strand": "+",
        "refseq": "NC_000003.12",
        "cytoband": "3q13.33",
        "exon_count": 14,
        "role": (
            "Homogentisate 1,2-dioxygenase -- enzyme in tyrosine catabolism "
            "converting homogentisic acid (HGA) to maleylacetoacetate.  "
            "Deficiency causes alkaptonuria (AKU): lifelong HGA accumulation "
            "-> ochronosis (black pigmentation of cartilage/connective tissue), "
            "severe early-onset osteoarthritis, cardiac valve calcification.  "
            "One of Garrod's original 4 inborn errors of metabolism (1902)."
        ),
        "disease": "Alkaptonuria (AKU)",
        "omim_disease": 203500,
        "omim_gene": 607474,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Met368Val (M368V)",
                "rsid": "rs121917987",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Common in European populations.",
            },
        ],
        "strategy": (
            "(1) AAV liver-directed gene replacement (cDNA ~1.3 kb, ideal "
            "for AAV) to restore hepatic HGA catabolism; "
            "(2) Nitisinone (NTBC) blocks HPD upstream, preventing HGA "
            "production (same drug used for tyrosinemia I); "
            "(3) Base editing for specific mutations.  "
            "Delivery: AAV8 IV to liver."
        ),
        "clinical_programs": (
            "Nitisinone approved in EU (2020) for AKU after SONIA 2 trial.  "
            "No gene therapy trials."
        ),
        "conditions": ["alkaptonuria", "AKU", "organic_acidemia",
                        "ochronosis", "metabolic"],
    },

    "ACADVL": {
        "gene_id": 37,
        "chrom": "chr17",
        "start": 7_217_125,
        "end": 7_222_863,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17p13.1",
        "exon_count": 20,
        "role": (
            "Very-long-chain acyl-CoA dehydrogenase -- mitochondrial enzyme "
            "for first step of VLCFA beta-oxidation (C14-C20).  Deficiency "
            "causes VLCAD deficiency: cardiomyopathy/sudden death in neonates, "
            "hypoketotic hypoglycemia, rhabdomyolysis.  On newborn screening "
            "panels worldwide.  Incidence ~1:40,000-1:120,000."
        ),
        "disease": "VLCAD deficiency",
        "omim_disease": 201475,
        "omim_gene": 609575,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Val283Ala (V283A)",
                "rsid": "rs77931234",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mutation; residual enzyme activity ~10-15%.  "
                    "Milder phenotype (myopathic/hypoglycemic)."
                ),
            },
        ],
        "strategy": (
            "(1) AAV liver/muscle gene replacement (cDNA ~2.0 kb, fits AAV); "
            "(2) Medium-chain triglyceride (MCT) diet bypasses VLCAD block; "
            "(3) Triheptanoin (Dojolvi) FDA-approved 2020 for LC-FAOD.  "
            "Delivery: AAV8 to liver + AAV9 muscle."
        ),
        "clinical_programs": (
            "Triheptanoin (Dojolvi, Ultragenyx) FDA-approved 2020 for "
            "long-chain fatty acid oxidation disorders.  "
            "No gene therapy trials."
        ),
        "conditions": ["VLCAD_deficiency", "fatty_acid_oxidation_disorder",
                        "organic_acidemia", "metabolic", "rhabdomyolysis"],
    },
}
