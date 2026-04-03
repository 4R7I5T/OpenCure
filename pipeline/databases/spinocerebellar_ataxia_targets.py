"""
Spinocerebellar ataxia (SCA) gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers autosomal dominant spinocerebellar ataxias caused by coding or
non-coding repeat expansions or point mutations.  SCA3/Machado-Joseph
disease is the most common SCA worldwide.

Most SCAs are polyglutamine (polyQ) diseases caused by CAG repeat expansion.
Therapeutic strategy is gene silencing (not replacement) because disease
is caused by toxic gain-of-function protein.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - ClinicalTrials.gov
  - Ionis, Novartis, Vico Therapeutics SCA ASO programs
"""


SPINOCEREBELLAR_ATAXIA_TARGETS = {

    "ATXN1": {
        "gene_id": 6310,
        "chrom": "chr6",
        "start": 16_299_112,
        "end": 16_761_712,
        "strand": "-",
        "refseq": "NC_000006.12",
        "cytoband": "6p22.3",
        "exon_count": 9,
        "role": (
            "Ataxin-1 -- nuclear protein involved in transcriptional regulation "
            "via interaction with Capicua (CIC).  CAG repeat expansion causes "
            "SCA1: progressive cerebellar ataxia, dysarthria, bulbar dysfunction.  "
            "PolyQ-expanded ataxin-1 forms toxic nuclear inclusions and "
            "dysregulates CIC-dependent transcription."
        ),
        "disease": "Spinocerebellar ataxia type 1 (SCA1)",
        "omim_disease": 164400,
        "omim_gene": 601556,
        "inheritance": "AD (CAG repeat expansion)",
        "key_variants": [
            {
                "name": "CAG repeat expansion (>39 repeats, with CAT interruption loss)",
                "rsid": None,
                "consequence": "repeat_expansion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Normal: 6-44 CAGs with 1-3 CAT interruptions.  "
                    "Pathogenic: >39 uninterrupted CAGs.  Anticipation "
                    "(earlier onset with larger repeats).  ~6% of all SCA."
                ),
            },
        ],
        "strategy": (
            "(1) ASO-mediated ATXN1 reduction -- allele-nonselective reduction "
            "is tolerable (ATXN1-null mice have mild phenotype); "
            "(2) CRISPRi to silence ATXN1 expression in Purkinje cells; "
            "(3) CRISPR excision of expanded CAG tract; "
            "(4) Modulation of CIC (downstream effector) as therapeutic bypass.  "
            "Delivery: intrathecal ASO or AAV9 with cerebellar-tropic capsid."
        ),
        "clinical_programs": (
            "No clinical trials as of 2026.  "
            "Ionis ATXN1 ASO preclinical."
        ),
        "conditions": ["SCA1", "spinocerebellar_ataxia", "cerebellar_ataxia",
                        "polyglutamine_disease", "trinucleotide_repeat"],
    },

    "ATXN3": {
        "gene_id": 4287,
        "chrom": "chr14",
        "start": 92_058_552,
        "end": 92_106_621,
        "strand": "-",
        "refseq": "NC_000014.9",
        "cytoband": "14q32.12",
        "exon_count": 11,
        "role": (
            "Ataxin-3 (Josephin) -- deubiquitinase involved in protein quality "
            "control.  CAG repeat expansion causes SCA3/Machado-Joseph disease "
            "(MJD), the most common SCA worldwide (~21% of all dominant ataxias).  "
            "PolyQ-expanded ataxin-3 is neurotoxic despite retaining catalytic "
            "activity.  High prevalence in Azores/Portuguese populations."
        ),
        "disease": "SCA3 / Machado-Joseph disease (MJD)",
        "omim_disease": 109150,
        "omim_gene": 607047,
        "inheritance": "AD (CAG repeat expansion)",
        "key_variants": [
            {
                "name": "CAG repeat expansion (>55 repeats)",
                "rsid": None,
                "consequence": "repeat_expansion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Normal: 12-44 CAGs.  Pathogenic: 55-87 CAGs.  "
                    "Intermediate 45-54: reduced penetrance.  "
                    "Azorean founder mutation; ~1:4000 in Azores.  "
                    "Most common SCA worldwide."
                ),
            },
        ],
        "strategy": (
            "(1) ASO-mediated ATXN3 reduction -- allele-nonselective knockdown "
            "shown safe in primates; "
            "(2) Ionis/Biomarin ION541 anti-ATXN3 ASO in development; "
            "(3) Vico Therapeutics VO659 anti-sense oligonucleotide Phase 1/2; "
            "(4) CRISPRi to silence ATXN3; "
            "(5) CRISPR CAG repeat excision.  "
            "Delivery: intrathecal ASO (preferred) or AAV9."
        ),
        "clinical_programs": (
            "Vico Therapeutics VO659 (anti-ATXN3 ASO) Phase 1/2 "
            "(NCT05822908).  "
            "Ionis anti-ATXN3 ASO preclinical.  "
            "No CRISPR trials as of 2026."
        ),
        "conditions": ["SCA3", "machado_joseph_disease", "MJD",
                        "spinocerebellar_ataxia", "cerebellar_ataxia",
                        "polyglutamine_disease"],
    },

    "CACNA1A": {
        "gene_id": 773,
        "chrom": "chr19",
        "start": 13_206_442,
        "end": 13_506_960,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19p13.13",
        "exon_count": 47,
        "role": (
            "Cav2.1 alpha-1A subunit -- P/Q-type voltage-gated calcium channel, "
            "principal presynaptic calcium channel at CNS synapses.  Different "
            "mutations cause different diseases: small CAG expansions (20-33) = "
            "SCA6 (pure cerebellar ataxia); missense GOF = familial hemiplegic "
            "migraine type 1 (FHM1); LOF = episodic ataxia type 2 (EA2)."
        ),
        "disease": "SCA6 / Episodic ataxia type 2 / Familial hemiplegic migraine",
        "omim_disease": 183086,
        "omim_gene": 601011,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "CAG expansion 20-33 (SCA6)",
                "rsid": None,
                "consequence": "repeat_expansion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Normal: 4-18.  SCA6: 20-33 (smallest pathogenic "
                    "polyQ expansion).  Pure cerebellar ataxia, late onset "
                    "(~50 years), very slow progression."
                ),
            },
            {
                "name": "p.Arg192Gln (R192Q, FHM1)",
                "rsid": "rs121908217",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "GOF; enhanced Cav2.1 current; hemiplegic migraine.",
            },
            {
                "name": "Truncating mutations (EA2)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Haploinsufficiency -> episodic ataxia with nystagmus.",
            },
        ],
        "strategy": (
            "(1) Acetazolamide is effective for EA2 (reduces attack frequency); "
            "(2) For SCA6: CRISPRi to reduce polyQ-CACNA1A expression; "
            "(3) Gene too large for AAV (~7.5 kb cDNA); "
            "(4) Allele-specific silencing for GOF FHM1 mutations.  "
            "Delivery: intrathecal ASO or AAV9."
        ),
        "clinical_programs": (
            "Acetazolamide for EA2.  "
            "4-aminopyridine (fampridine) off-label for EA2.  "
            "No gene therapy trials."
        ),
        "conditions": ["SCA6", "episodic_ataxia", "EA2", "FHM1",
                        "hemiplegic_migraine", "spinocerebellar_ataxia",
                        "cerebellar_ataxia"],
    },

    "ATXN7": {
        "gene_id": 6314,
        "chrom": "chr3",
        "start": 63_865_653,
        "end": 64_002_710,
        "strand": "+",
        "refseq": "NC_000003.12",
        "cytoband": "3p14.1",
        "exon_count": 13,
        "role": (
            "Ataxin-7 -- component of STAGA transcriptional coactivator "
            "complex (HAT module).  CAG repeat expansion causes SCA7: "
            "cerebellar ataxia with progressive retinal degeneration "
            "(cone-rod dystrophy), the only SCA with visual loss.  "
            "Severe infantile form with >100 repeats."
        ),
        "disease": "SCA7 (spinocerebellar ataxia with retinal degeneration)",
        "omim_disease": 164500,
        "omim_gene": 607640,
        "inheritance": "AD (CAG expansion)",
        "key_variants": [
            {
                "name": "CAG expansion (>37 repeats)",
                "rsid": None,
                "consequence": "repeat_expansion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Normal: 4-35.  Pathogenic: >37.  Extreme anticipation "
                    "(parental repeat can expand massively, especially paternal)."
                ),
            },
        ],
        "strategy": (
            "(1) ASO-mediated ATXN7 reduction (intrathecal for CNS, "
            "intravitreal for retina); "
            "(2) CRISPRi to silence ATXN7; "
            "(3) Subretinal AAV for retinal degeneration component.  "
            "Delivery: intrathecal ASO + intravitreal."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["SCA7", "spinocerebellar_ataxia", "cerebellar_ataxia",
                        "retinal_degeneration", "polyglutamine_disease"],
    },

    "KCNC3": {
        "gene_id": 3748,
        "chrom": "chr19",
        "start": 50_228_496,
        "end": 50_260_963,
        "strand": "-",
        "refseq": "NC_000019.10",
        "cytoband": "19q13.33",
        "exon_count": 4,
        "role": (
            "Kv3.3 potassium channel -- fast-repolarizing channel highly "
            "expressed in cerebellar Purkinje cells and deep nuclei.  "
            "Mutations cause SCA13: childhood-onset cerebellar ataxia with "
            "intellectual disability."
        ),
        "disease": "SCA13",
        "omim_disease": 605259,
        "omim_gene": 176264,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Arg420His (R420H)",
                "rsid": "rs121918339",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Pore domain; dominant-negative; childhood-onset.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant-negative allele; "
            "(2) AAV gene replacement (cDNA ~2.2 kb fits AAV) with "
            "suppression-and-replacement; "
            "(3) Potassium channel modulators."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["SCA13", "spinocerebellar_ataxia", "cerebellar_ataxia",
                        "channelopathy"],
    },

    "ITPR1": {
        "gene_id": 3708,
        "chrom": "chr3",
        "start": 4_494_082,
        "end": 4_838_750,
        "strand": "+",
        "refseq": "NC_000003.12",
        "cytoband": "3p26.1",
        "exon_count": 58,
        "role": (
            "Inositol 1,4,5-trisphosphate receptor type 1 -- ER calcium "
            "release channel, the predominant IP3R in Purkinje cells.  "
            "Heterozygous mutations cause SCA15/SCA29 (slowly progressive "
            "cerebellar ataxia).  Homozygous LOF causes Gillespie syndrome "
            "(ataxia + aniridia + ID)."
        ),
        "disease": "SCA15/SCA29 / Gillespie syndrome",
        "omim_disease": 606658,
        "omim_gene": 147265,
        "inheritance": "AD (SCA15/29); AR (Gillespie)",
        "key_variants": [
            {
                "name": "Heterozygous deletions and missense",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Haploinsufficiency in SCA15; slow progression.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type allele; "
            "(2) Gene too large for AAV (~8.4 kb cDNA); dual-AAV approach."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["SCA15", "SCA29", "gillespie_syndrome",
                        "spinocerebellar_ataxia", "cerebellar_ataxia"],
    },
}
