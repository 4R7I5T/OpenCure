"""
Prion disease and genetic sleep disorder targets -- GRCh38 (hg38).

Covers genetic prion diseases (FFI, GSS, genetic CJD) and narcolepsy.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


PRION_SLEEP_TARGETS = {

    "PRNP": {
        "gene_id": 5621,
        "chrom": "chr20",
        "start": 4_666_796,
        "end": 4_682_233,
        "strand": "+",
        "refseq": "NC_000020.11",
        "cytoband": "20p13",
        "exon_count": 2,
        "role": (
            "Prion protein (PrPC) -- GPI-anchored glycoprotein on neuronal "
            "surface, normal function includes copper binding, myelin "
            "maintenance, and neuroprotection.  Misfolding to PrPSc causes "
            "transmissible spongiform encephalopathies.  Mutations cause "
            "genetic prion diseases: fatal familial insomnia (FFI), "
            "Gerstmann-Straussler-Scheinker (GSS), and genetic Creutzfeldt-"
            "Jakob disease (gCJD).  ALL are fatal with no treatment."
        ),
        "disease": "Fatal familial insomnia / GSS / genetic CJD",
        "omim_disease": 600072,
        "omim_gene": 176640,
        "inheritance": "AD (with codon 129 polymorphism modifier)",
        "key_variants": [
            {
                "name": "p.Asp178Asn + codon 129 Met (D178N-129M = FFI)",
                "rsid": "rs74315403",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "D178N with Met at codon 129 in cis = FFI (progressive "
                    "insomnia, dysautonomia, dementia).  D178N with Val at "
                    "129 in cis = gCJD.  Same mutation, different disease "
                    "based on 129 polymorphism."
                ),
            },
            {
                "name": "p.Pro102Leu (P102L = GSS)",
                "rsid": "rs74315400",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common GSS mutation; cerebellar ataxia + dementia.",
            },
            {
                "name": "p.Glu200Lys (E200K = gCJD)",
                "rsid": "rs28933385",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common gCJD mutation worldwide.  Libyan Jewish and "
                    "Slovak founder populations.  Onset ~60 years."
                ),
            },
        ],
        "strategy": (
            "(1) ASO-mediated PRNP reduction -- Ionis ION717 approach; "
            "PRNP knockdown is safe (PrP-null mice are largely normal) and "
            "prevents/delays prion disease in animal models; "
            "(2) CRISPRi to silence PRNP expression in CNS; "
            "(3) Allele-specific silencing of mutant PRNP allele; "
            "(4) CRISPR knockout of PRNP -- radical but potentially curative "
            "for presymptomatic carriers (heterozygous knockout sufficient).  "
            "Delivery: intrathecal ASO or AAV9 IV/IT."
        ),
        "clinical_programs": (
            "Ionis ION717 -- anti-PRNP ASO, Phase 1/2 for genetic prion "
            "disease (NCT06153173).  First disease-modifying trial for "
            "prion disease.  Sonia Vallabh and Eric Minikel (Broad Institute) "
            "are patient-scientists driving this program.  "
            "No CRISPR trials as of 2026.  "
            "All prion diseases are currently 100% fatal."
        ),
        "conditions": ["fatal_familial_insomnia", "FFI", "GSS",
                        "genetic_CJD", "prion_disease", "sleep_disorder",
                        "spongiform_encephalopathy"],
    },

    "HCRT": {
        "gene_id": 3060,
        "chrom": "chr17",
        "start": 42_187_284,
        "end": 42_188_703,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.2",
        "exon_count": 2,
        "role": (
            "Hypocretin/orexin -- neuropeptide produced by ~70,000 "
            "hypothalamic neurons, essential for sleep-wake stability.  "
            "Narcolepsy type 1 is caused by selective autoimmune destruction "
            "of orexin-producing neurons (associated with HLA-DQB1*06:02).  "
            "Extremely rare genetic narcolepsy from HCRT mutations exists.  "
            "Loss of orexin signaling causes cataplexy, excessive daytime "
            "sleepiness, sleep paralysis, and hypnagogic hallucinations."
        ),
        "disease": "Narcolepsy type 1 (orexin deficiency)",
        "omim_disease": 161400,
        "omim_gene": 602358,
        "inheritance": "complex/autoimmune (sporadic); AR (extremely rare genetic)",
        "key_variants": [
            {
                "name": "Autoimmune destruction of orexin neurons",
                "rsid": None,
                "consequence": "cell_loss",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Not a classical monogenic disease.  >98% of narcolepsy "
                    "type 1 is autoimmune.  HLA-DQB1*06:02 is present in "
                    "~98% of patients (vs ~25% of general population).  "
                    "Gene therapy targets orexin restoration, not gene correction."
                ),
            },
        ],
        "strategy": (
            "(1) AAV-mediated HCRT/orexin gene replacement to surviving "
            "hypothalamic neurons or ectopic expression sites (cDNA ~0.4 kb, "
            "ideal for AAV); "
            "(2) AAV-HCRT intracerebroventricular showed rescue in "
            "orexin-knockout mice; "
            "(3) Cell replacement: transplant of iPSC-derived orexin neurons; "
            "(4) Orexin receptor agonists as pharmacological approach "
            "(Takeda TAK-994 approach, but discontinued).  "
            "Delivery: AAV stereotactic to lateral hypothalamus or ICV."
        ),
        "clinical_programs": (
            "Takeda TAK-994 (oral orexin 2 receptor agonist) Phase 2 -- "
            "discontinued due to hepatotoxicity.  "
            "Takeda TAK-861 (next-gen OX2R agonist) Phase 2 ongoing.  "
            "AAV-HCRT preclinical (Bhatt lab, Stanford).  "
            "No gene therapy trials as of 2026."
        ),
        "conditions": ["narcolepsy", "narcolepsy_type_1", "orexin_deficiency",
                        "sleep_disorder", "cataplexy"],
    },
}
