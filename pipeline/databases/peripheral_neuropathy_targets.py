"""
Peripheral neuropathy gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers Charcot-Marie-Tooth disease (CMT) subtypes and hereditary sensory
and autonomic neuropathies (HSAN).  CMT is the most common inherited
neuromuscular disorder (~1:2500 prevalence).

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - ClinicalTrials.gov
  - DTx Pharma/Novartis, Ionis, Sarepta CMT programs
"""


PERIPHERAL_NEUROPATHY_TARGETS = {

    # ===================================================================
    # CMT1A -- PMP22 (most common CMT, ~50%)
    # ===================================================================
    "PMP22": {
        "gene_id": 5376,
        "chrom": "chr17",
        "start": 15_132_888,
        "end": 15_168_639,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17p12",
        "exon_count": 5,
        "role": (
            "Peripheral myelin protein 22 -- integral membrane glycoprotein "
            "in compact myelin of Schwann cells.  1.4 Mb duplication on 17p12 "
            "containing PMP22 causes CMT1A (demyelinating, ~50% of all CMT).  "
            "PMP22 deletion causes HNPP.  Point mutations cause CMT1E "
            "(severe).  PMP22 is exquisitely dosage-sensitive."
        ),
        "disease": "CMT1A (duplication) / CMT1E (point mutation) / HNPP (deletion)",
        "omim_disease": 118220,
        "omim_gene": 601097,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "17p12 duplication (1.4 Mb)",
                "rsid": None,
                "consequence": "duplication",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Accounts for ~70-80% of CMT1.  NAHR between flanking "
                    "CMT1A-REP repeats.  3 copies of PMP22 -> overexpression "
                    "-> Schwann cell dysfunction -> demyelination."
                ),
            },
            {
                "name": "p.Leu16Pro (Trembler-J)",
                "rsid": "rs104894613",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Severe CMT1E; dominant-negative misfolded PMP22.",
            },
        ],
        "strategy": (
            "(1) ASO or siRNA-mediated PMP22 reduction to normalize dosage "
            "in CMT1A duplication patients (reducing 3-copy -> 2-copy "
            "equivalent) -- DTx Pharma/Novartis approach; "
            "(2) CRISPRi to fine-tune PMP22 expression to optimal level; "
            "(3) CRISPR excision of duplicated PMP22 copy; "
            "(4) Allele-specific silencing of mutant allele for CMT1E.  "
            "Delivery: subcutaneous siRNA (conjugated to Schwann cell ligand) "
            "or intrathecal ASO."
        ),
        "clinical_programs": (
            "DTx Pharma (acquired by Novartis) DTXP-003 -- siRNA targeting "
            "PMP22 for CMT1A, preclinical/IND-enabling.  "
            "Ionis anti-PMP22 ASO preclinical (showed benefit in CMT1A rats).  "
            "No approved gene therapies as of 2026.  "
            "PXT3003 (Pharnext, baclofen+naltrexone+sorbitol) Phase 3 failed."
        ),
        "conditions": ["CMT1A", "CMT1E", "HNPP", "charcot_marie_tooth",
                        "peripheral_neuropathy", "demyelinating_neuropathy"],
    },

    # ===================================================================
    # CMT2A -- MFN2 (most common axonal CMT)
    # ===================================================================
    "MFN2": {
        "gene_id": 9927,
        "chrom": "chr1",
        "start": 11_980_181,
        "end": 12_013_503,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1p36.22",
        "exon_count": 19,
        "role": (
            "Mitofusin 2 -- outer mitochondrial membrane GTPase essential "
            "for mitochondrial fusion, ER-mitochondria tethering, and axonal "
            "mitochondrial transport.  Mutations cause CMT2A (axonal, most "
            "common axonal CMT ~20%), optic atrophy plus, and hereditary "
            "motor-sensory neuropathy type VI."
        ),
        "disease": "CMT2A (axonal Charcot-Marie-Tooth)",
        "omim_disease": 609260,
        "omim_gene": 608507,
        "inheritance": "AD (most); AR (rare)",
        "key_variants": [
            {
                "name": "p.Arg94Gln (R94Q)",
                "rsid": "rs119103249",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "GTPase domain; most common MFN2 mutation; early-onset severe.",
            },
            {
                "name": "p.Arg280His (R280H)",
                "rsid": "rs119103250",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "GTPase domain; moderate severity.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.3 kb fits AAV) with neuron-"
            "specific promoter targeting DRG sensory and motor neurons; "
            "(2) Allele-specific silencing of dominant-negative mutant allele; "
            "(3) Base editing for recurrent mutations (R94Q).  "
            "Delivery: AAV9 intrathecal or IV."
        ),
        "clinical_programs": "No gene therapy or CRISPR trials as of 2026.",
        "conditions": ["CMT2A", "charcot_marie_tooth", "axonal_neuropathy",
                        "peripheral_neuropathy", "mitochondrial_dynamics"],
    },

    # ===================================================================
    # CMTX1 -- GJB1 (Connexin 32)
    # ===================================================================
    "GJB1": {
        "gene_id": 2705,
        "chrom": "chrX",
        "start": 71_223_251,
        "end": 71_233_297,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xq13.1",
        "exon_count": 2,
        "role": (
            "Connexin 32 -- gap junction protein in Schwann cells (non-compact "
            "myelin) and hepatocytes.  Mutations cause CMTX1, second most "
            "common form of CMT (~10%).  Males more severely affected; "
            "females often have milder symptoms or are asymptomatic."
        ),
        "disease": "CMTX1 (X-linked Charcot-Marie-Tooth)",
        "omim_disease": 302800,
        "omim_gene": 304040,
        "inheritance": "XL (semi-dominant)",
        "key_variants": [
            {
                "name": ">450 pathogenic mutations (no hotspot)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Distributed throughout coding region; missense predominate.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~0.8 kb, ideal for AAV) with "
            "Schwann cell-specific promoter (MPZ/P0 promoter); "
            "(2) AAV-GJB1 restored connexin 32 expression in Gjb1-null mice "
            "and improved nerve conduction (Sargiannidou et al.).  "
            "Delivery: intrathecal or intranerve AAV9."
        ),
        "clinical_programs": (
            "No clinical trials as of 2026.  "
            "Strong preclinical AAV-GJB1 data in Gjb1-null mice."
        ),
        "conditions": ["CMTX1", "charcot_marie_tooth", "peripheral_neuropathy",
                        "X_linked_CMT"],
    },

    # ===================================================================
    # CMT1B -- MPZ (Myelin Protein Zero)
    # ===================================================================
    "MPZ": {
        "gene_id": 4359,
        "chrom": "chr1",
        "start": 161_274_744,
        "end": 161_280_643,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q23.3",
        "exon_count": 6,
        "role": (
            "Myelin protein zero (P0) -- most abundant peripheral nerve myelin "
            "protein (~50% of PNS myelin).  Immunoglobulin superfamily; "
            "mediates myelin compaction via homophilic adhesion.  Mutations "
            "cause CMT1B (demyelinating) or CMT2I/J (axonal, late-onset)."
        ),
        "disease": "CMT1B / CMT2I / CMT2J / Dejerine-Sottas (DSS)",
        "omim_disease": 118200,
        "omim_gene": 159440,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Ser63Cys (S63C)",
                "rsid": "rs104894617",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Ig domain; unfolded protein response in Schwann cells.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant mutant allele; "
            "(2) AAV gene replacement (cDNA ~0.7 kb, ideal for AAV); "
            "(3) Base editing for specific missense mutations.  "
            "Delivery: AAV with Schwann cell promoter."
        ),
        "clinical_programs": "No gene therapy trials as of 2026.",
        "conditions": ["CMT1B", "charcot_marie_tooth", "peripheral_neuropathy",
                        "dejerine_sottas"],
    },

    # ===================================================================
    # CMT4C -- SH3TC2 (most common AR CMT)
    # ===================================================================
    "SH3TC2": {
        "gene_id": 79628,
        "chrom": "chr5",
        "start": 149_007_280,
        "end": 149_049_388,
        "strand": "+",
        "refseq": "NC_000005.10",
        "cytoband": "5q32",
        "exon_count": 17,
        "role": (
            "SH3 domain and tetratricopeptide repeats 2 -- Schwann cell "
            "protein involved in endosomal recycling and myelination.  "
            "Biallelic mutations cause CMT4C, the most common autosomal "
            "recessive demyelinating CMT.  Early-onset scoliosis is distinctive."
        ),
        "disease": "CMT4C (autosomal recessive demyelinating CMT)",
        "omim_disease": 601596,
        "omim_gene": 608206,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg954Ter (R954X)",
                "rsid": "rs121918360",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Most common in European Roma/Romani populations.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~3.8 kb fits AAV); "
            "(2) CRISPR correction for R954X.  "
            "Delivery: AAV9 intrathecal."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["CMT4C", "charcot_marie_tooth", "peripheral_neuropathy",
                        "autosomal_recessive_CMT"],
    },

    # ===================================================================
    # HSAN1 -- SPTLC1
    # ===================================================================
    "SPTLC1": {
        "gene_id": 10558,
        "chrom": "chr9",
        "start": 92_059_100,
        "end": 92_145_085,
        "strand": "-",
        "refseq": "NC_000009.12",
        "cytoband": "9q22.31",
        "exon_count": 15,
        "role": (
            "Serine palmitoyltransferase long chain base subunit 1 -- first "
            "enzyme in sphingolipid biosynthesis.  Gain-of-function mutations "
            "cause promiscuous use of alanine/glycine instead of serine, "
            "producing neurotoxic deoxy-sphingolipids.  Causes HSAN1: "
            "progressive sensory neuropathy with lancinating pain, ulcers, "
            "and osteomyelitis."
        ),
        "disease": "Hereditary sensory neuropathy type 1 (HSAN1)",
        "omim_disease": 162400,
        "omim_gene": 605712,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Cys133Trp (C133W)",
                "rsid": "rs104894283",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common; alters substrate specificity.",
            },
            {
                "name": "p.Cys133Tyr (C133Y)",
                "rsid": "rs104894284",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Same codon; slightly milder phenotype.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant SPTLC1 (GOF, not LOF); "
            "(2) Substrate competition: L-serine supplementation reduces "
            "deoxy-sphingolipid levels (clinical trials); "
            "(3) Base editing to correct C133W/Y.  "
            "Delivery: AAV DRG-targeted or oral L-serine."
        ),
        "clinical_programs": (
            "L-serine supplementation Phase 1/2 (NCT01733407) -- showed "
            "biomarker improvement.  No CRISPR trials."
        ),
        "conditions": ["HSAN1", "hereditary_sensory_neuropathy",
                        "peripheral_neuropathy", "sphingolipid_disorder"],
    },
}
