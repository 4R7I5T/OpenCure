"""
Skeletal muscle channelopathy gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers non-dystrophic myotonias, periodic paralyses, and malignant
hyperthermia susceptibility.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


SKELETAL_CHANNELOPATHY_TARGETS = {

    "RYR1": {
        "gene_id": 6261,
        "chrom": "chr19",
        "start": 38_433_691,
        "end": 38_587_564,
        "strand": "-",
        "refseq": "NC_000019.10",
        "cytoband": "19q13.2",
        "exon_count": 106,
        "role": (
            "Ryanodine receptor 1 -- skeletal muscle sarcoplasmic reticulum "
            "calcium release channel.  Gain-of-function mutations cause "
            "malignant hyperthermia susceptibility (MHS1) and central core "
            "disease (CCD).  Loss-of-function causes multi-minicore disease "
            "and congenital fiber-type disproportion.  Most common inherited "
            "skeletal muscle disorder."
        ),
        "disease": "Malignant hyperthermia / Central core disease / Multi-minicore",
        "omim_disease": 145600,
        "omim_gene": 180901,
        "inheritance": "AD (MHS, CCD); AR (multi-minicore)",
        "key_variants": [
            {
                "name": "p.Arg614Cys (R614C)",
                "rsid": "rs121918592",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Classic MH mutation; N-terminal hotspot.",
            },
            {
                "name": "p.Arg2163His (R2163H)",
                "rsid": "rs121918594",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Central hotspot; MH/CCD.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of GOF mutant allele for MH "
            "(preventive gene therapy -- alternative to lifelong volatile "
            "anesthetic avoidance); "
            "(2) Gene too large for AAV (~15 kb cDNA); dual-AAV or "
            "exon skipping for AR forms; "
            "(3) Dantrolene remains acute treatment for MH crisis.  "
            "Delivery: AAV9 IV or intramuscular."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "Dantrolene (IV) for MH crisis.  "
            "Avoidance of triggering anesthetics (volatile agents, "
            "succinylcholine) is standard of care."
        ),
        "conditions": ["malignant_hyperthermia", "MH", "central_core_disease",
                        "multi_minicore", "skeletal_channelopathy", "myopathy"],
    },

    "CACNA1S": {
        "gene_id": 779,
        "chrom": "chr1",
        "start": 201_008_638,
        "end": 201_081_700,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q32.1",
        "exon_count": 44,
        "role": (
            "Voltage-gated calcium channel Cav1.1 alpha-1S subunit -- "
            "skeletal muscle L-type calcium channel/voltage sensor for "
            "excitation-contraction coupling.  Mutations cause hypokalemic "
            "periodic paralysis type 1 (HypoPP1) and MH susceptibility type 5."
        ),
        "disease": "Hypokalemic periodic paralysis type 1 / MHS5",
        "omim_disease": 170400,
        "omim_gene": 114208,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Arg528His (R528H)",
                "rsid": "rs121912709",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common HypoPP1 mutation (~50%); S4 voltage sensor.",
            },
            {
                "name": "p.Arg1239His (R1239H)",
                "rsid": "rs121912710",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Second most common; S4 of domain IV.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant allele (gating pore "
            "current mechanism); "
            "(2) Base editing for R528H or R1239H (S4 arginine mutations); "
            "(3) Acetazolamide/dichlorphenamide effective for most patients.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": (
            "Dichlorphenamide (Keveyis) FDA-approved for HypoPP.  "
            "No gene therapy trials."
        ),
        "conditions": ["hypokalemic_periodic_paralysis", "HypoPP",
                        "periodic_paralysis", "skeletal_channelopathy", "MHS"],
    },

    "SCN4A": {
        "gene_id": 6329,
        "chrom": "chr17",
        "start": 63_934_553,
        "end": 63_970_298,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q23.3",
        "exon_count": 24,
        "role": (
            "Skeletal muscle sodium channel Nav1.4 -- voltage-gated sodium "
            "channel responsible for action potential initiation in skeletal "
            "muscle.  GOF mutations cause hyperkalemic periodic paralysis "
            "(HyperPP), paramyotonia congenita (PMC), and sodium channel "
            "myotonia.  LOF causes congenital myasthenic syndrome."
        ),
        "disease": "HyperPP / Paramyotonia congenita / Sodium channel myotonia",
        "omim_disease": 170500,
        "omim_gene": 603967,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Thr704Met (T704M)",
                "rsid": "rs121918333",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common HyperPP mutation (~50%); enhanced slow inactivation.",
            },
            {
                "name": "p.Met1592Val (M1592V)",
                "rsid": "rs121918334",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Second common HyperPP (~30%).",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of GOF mutant allele; "
            "(2) Base editing for T704M or M1592V; "
            "(3) Mexiletine (sodium channel blocker) effective for most.  "
            "Delivery: AAV9 IV to skeletal muscle."
        ),
        "clinical_programs": (
            "Mexiletine off-label for myotonia/HyperPP.  "
            "Dichlorphenamide approved for periodic paralysis.  "
            "No gene therapy trials."
        ),
        "conditions": ["hyperkalemic_periodic_paralysis", "HyperPP",
                        "paramyotonia_congenita", "sodium_channel_myotonia",
                        "skeletal_channelopathy", "periodic_paralysis"],
    },

    "CLCN1": {
        "gene_id": 1180,
        "chrom": "chr7",
        "start": 143_300_699,
        "end": 143_325_330,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q34",
        "exon_count": 23,
        "role": (
            "Chloride voltage-gated channel 1 -- skeletal muscle chloride "
            "channel responsible for ~80% of resting membrane conductance.  "
            "Mutations cause myotonia congenita: Thomsen disease (AD, milder) "
            "or Becker myotonia (AR, more common and severe).  Muscle "
            "stiffness on initiating movement, warm-up phenomenon."
        ),
        "disease": "Myotonia congenita (Thomsen/Becker)",
        "omim_disease": 160800,
        "omim_gene": 118425,
        "inheritance": "AD (Thomsen) / AR (Becker)",
        "key_variants": [
            {
                "name": "p.Phe413Cys (F413C, Thomsen founder)",
                "rsid": "rs121908466",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Thomsen's original family (Danish founder).",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement for AR Becker form (cDNA ~2.9 kb, "
            "fits AAV); "
            "(2) Allele-specific silencing for AD Thomsen; "
            "(3) Mexiletine effective first-line therapy for most patients.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": "Mexiletine SOC. No gene therapy trials.",
        "conditions": ["myotonia_congenita", "thomsen_disease", "becker_myotonia",
                        "skeletal_channelopathy", "non_dystrophic_myotonia"],
    },

    "KCNJ2": {
        "gene_id": 3759,
        "chrom": "chr17",
        "start": 70_176_964,
        "end": 70_187_169,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q24.3",
        "exon_count": 3,
        "role": (
            "Inward rectifier potassium channel Kir2.1 -- sets resting "
            "membrane potential in muscle and heart.  Mutations cause "
            "Andersen-Tawil syndrome (ATS/LQT7): periodic paralysis + "
            "cardiac arrhythmias (bidirectional VT, prolonged QTc) + "
            "dysmorphic features (micrognathia, clinodactyly)."
        ),
        "disease": "Andersen-Tawil syndrome (ATS type 1 / LQT7)",
        "omim_disease": 170390,
        "omim_gene": 600681,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Arg218Trp (R218W)",
                "rsid": "rs121909250",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Selectivity filter; dominant-negative.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant-negative allele; "
            "(2) AAV gene replacement (cDNA ~1.3 kb fits AAV) with "
            "suppression-and-replacement; "
            "(3) Flecainide and acetazolamide for symptom control."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["andersen_tawil_syndrome", "ATS", "LQT7",
                        "periodic_paralysis", "skeletal_channelopathy",
                        "cardiac_arrhythmia"],
    },
}
