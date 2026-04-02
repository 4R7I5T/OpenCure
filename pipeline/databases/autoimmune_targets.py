"""
Autoinflammatory and autoimmune disease gene targets — GRCh38 coordinates.

Comprehensive panel covering monogenic autoinflammatory syndromes (SAIDs),
with emphasis on NOD2-associated diseases (Yao syndrome, Blau syndrome,
Crohn's disease) and the broader periodic fever / inflammasome spectrum.

Modeled after the Blueprint Genetics 47-gene autoinflammatory panel with
additional targets relevant to gene-therapy and diagnostic pipelines.

Sources: OMIM, NCBI Gene, ClinVar, Infevers database, Nomani et al. 2024
(PMC11449693), Yao 2019 (PMID: 31084224).

All coordinates are GRCh38/hg38.
"""

# ---------------------------------------------------------------------------
# Primary: NOD2 and associated diseases
# ---------------------------------------------------------------------------
AUTOIMMUNE_TARGETS = {
    "NOD2": {
        "chrom": "chr16",
        "start": 50693606,
        "end": 50733075,
        "role": "Cytosolic innate immune receptor — MDP/peptidoglycan sensing, "
                "NF-kB and p38 MAPK activation, autophagy",
        "strategy": "Diagnostic: genotype all coding exons + IVS8+158 intronic "
                     "region. Therapeutic: base-editing for gain-of-function "
                     "variants, CRISPRi for overexpression phenotypes",
        "omim_gene": 605956,
        "aliases": ["CARD15", "NLRC2", "IBD1", "BLAU"],
        "inheritance": "complex",
        "pathway": "NF-kB, p38_MAPK, autophagy, innate_immunity",
        "associated_diseases": [
            "yao_syndrome",
            "blau_syndrome",
            "crohn_disease",
            "early_onset_sarcoidosis",
        ],
        "conditions": ["yao_syndrome", "autoinflammatory", "IBD",
                        "blau_syndrome"],
        "protein_domains": {
            "CARD1":  {"start_aa": 1,   "end_aa": 110,  "function": "RIPK2 signaling"},
            "CARD2":  {"start_aa": 111, "end_aa": 194,  "function": "RIPK2 signaling"},
            "NBD":    {"start_aa": 195, "end_aa": 425,  "function": "ATP binding, oligomerization"},
            "HD1":    {"start_aa": 426, "end_aa": 485,  "function": "autoinhibition"},
            "WHD":    {"start_aa": 486, "end_aa": 602,  "function": "signal transduction"},
            "HD2":    {"start_aa": 603, "end_aa": 743,  "function": "LRR linkage"},
            "LRR":    {"start_aa": 744, "end_aa": 1020, "function": "MDP/PAMP recognition"},
        },
        # Must capture intronic regions — IVS8+158 is 158bp into intron 8
        "extended_capture_regions": [
            {"chrom": "chr16", "start": 50722705, "end": 50723301,
             "reason": "Intron 8 — contains IVS8+158 (rs5743289), the "
                       "primary Yao syndrome variant"},
        ],
    },

    # ---------------------------------------------------------------------------
    # Inflammasome / periodic fever genes
    # ---------------------------------------------------------------------------
    "MEFV": {
        "chrom": "chr16",
        "start": 3242027,
        "end": 3256633,
        "role": "Pyrin — inflammasome regulation, cytoskeletal innate immunity",
        "strategy": "Diagnostic: screen exon 2, 3, 5, 10 hotspots. "
                     "Therapeutic: base-editing to correct M694V and other "
                     "gain-of-function variants in exon 10",
        "omim_gene": 608107,
        "inheritance": "AR",
        "pathway": "pyrin_inflammasome, IL-1beta",
        "associated_diseases": ["familial_mediterranean_fever"],
        "conditions": ["FMF", "autoinflammatory", "periodic_fever"],
    },
    "NLRP3": {
        "chrom": "chr1",
        "start": 247416077,
        "end": 247448817,
        "role": "Cryopyrin — NLRP3 inflammasome sensor, IL-1beta processing",
        "strategy": "Diagnostic: screen NACHT domain exons 3-6. "
                     "Therapeutic: CRISPRi silencing or base-editing of "
                     "gain-of-function variants causing constitutive activation",
        "omim_gene": 606416,
        "inheritance": "AD",
        "pathway": "NLRP3_inflammasome, IL-1beta, caspase-1",
        "associated_diseases": ["FCAS", "MWS", "NOMID_CINCA"],
        "conditions": ["CAPS", "autoinflammatory", "periodic_fever"],
    },
    "NLRC4": {
        "chrom": "chr2",
        "start": 32224449,
        "end": 32265743,
        "role": "NLRC4 inflammasome — bacterial flagellin/T3SS sensing",
        "strategy": "Diagnostic: screen NACHT domain. Therapeutic: CRISPRi "
                     "silencing of gain-of-function variants",
        "omim_gene": 606831,
        "inheritance": "AD",
        "pathway": "NLRC4_inflammasome, IL-1beta, IL-18",
        "associated_diseases": ["NLRC4-MAS", "FCAS4"],
        "conditions": ["autoinflammatory", "macrophage_activation"],
    },
    "NLRP12": {
        "chrom": "chr19",
        "start": 53793583,
        "end": 53824402,
        "role": "NLRP12 — negative regulator of NF-kB and inflammasome",
        "strategy": "Diagnostic: full coding sequence. Modifier gene for "
                     "Yao syndrome (digenic NOD2+NLRP12 combinations in 22%)",
        "omim_gene": 609648,
        "inheritance": "AD",
        "pathway": "NF-kB_regulation, inflammasome",
        "associated_diseases": ["FCAS2", "NLRP12-AID"],
        "conditions": ["autoinflammatory", "periodic_fever",
                        "yao_syndrome_modifier"],
    },
    "MVK": {
        "chrom": "chr12",
        "start": 109573272,
        "end": 109598125,
        "role": "Mevalonate kinase — cholesterol/isoprenoid biosynthesis, "
                "pyrin inflammasome regulation",
        "strategy": "Diagnostic: screen for V377I and I268T. Therapeutic: "
                     "base-editing to correct loss-of-function variants",
        "omim_gene": 251170,
        "inheritance": "AR",
        "pathway": "mevalonate, pyrin_inflammasome",
        "associated_diseases": ["HIDS", "mevalonic_aciduria"],
        "conditions": ["autoinflammatory", "periodic_fever", "MKD"],
    },

    # ---------------------------------------------------------------------------
    # TNF receptor / NF-kB pathway
    # ---------------------------------------------------------------------------
    "TNFRSF1A": {
        "chrom": "chr12",
        "start": 6328771,
        "end": 6342076,
        "role": "TNF receptor 1 — NF-kB and apoptosis signaling",
        "strategy": "Diagnostic: screen extracellular cysteine-rich domains. "
                     "Therapeutic: base-editing to correct misfolding mutations",
        "omim_gene": 191190,
        "inheritance": "AD",
        "pathway": "TNF_signaling, NF-kB, apoptosis",
        "associated_diseases": ["TRAPS"],
        "conditions": ["TRAPS", "autoinflammatory", "periodic_fever"],
    },
    "TNFAIP3": {
        "chrom": "chr6",
        "start": 137866349,
        "end": 137883312,
        "role": "A20 ubiquitin-editing enzyme — NF-kB negative regulation",
        "strategy": "Diagnostic: screen for truncating/loss-of-function variants. "
                     "Therapeutic: CRISPRa to restore expression from remaining "
                     "allele in haploinsufficiency",
        "omim_gene": 191163,
        "inheritance": "AD",
        "pathway": "NF-kB_regulation, ubiquitin",
        "associated_diseases": ["HA20"],
        "conditions": ["autoinflammatory", "behcet_like", "autoimmune"],
    },
    "OTULIN": {
        "chrom": "chr5",
        "start": 14664718,
        "end": 14716525,
        "role": "Linear deubiquitinase — LUBAC/NF-kB regulation",
        "strategy": "Diagnostic: full coding screen. Therapeutic: "
                     "base-editing to correct loss-of-function",
        "omim_gene": 615712,
        "inheritance": "AR",
        "pathway": "NF-kB_regulation, linear_ubiquitin",
        "associated_diseases": ["ORAS"],
        "conditions": ["autoinflammatory"],
    },

    # ---------------------------------------------------------------------------
    # Proteasome / interferon pathway
    # ---------------------------------------------------------------------------
    "PSMB8": {
        "chrom": "chr6",
        "start": 32840717,
        "end": 32844679,
        "role": "Immunoproteasome beta subunit — protein degradation, "
                "type I interferon regulation",
        "strategy": "Diagnostic: screen for T75M and other loss-of-function. "
                     "Therapeutic: base-editing or HDR correction",
        "omim_gene": 177046,
        "inheritance": "AR",
        "pathway": "proteasome, type_I_interferon",
        "associated_diseases": ["CANDLE", "PRAAS"],
        "conditions": ["autoinflammatory", "interferonopathy"],
    },
    "STING1": {
        "chrom": "chr5",
        "start": 139475533,
        "end": 139482758,
        "role": "STING — cGAS-STING innate DNA sensing, type I IFN induction",
        "strategy": "Diagnostic: screen for gain-of-function in exons 5-7. "
                     "Therapeutic: CRISPRi silencing or base-editing of "
                     "constitutively active variants",
        "omim_gene": 612374,
        "aliases": ["TMEM173"],
        "inheritance": "AD",
        "pathway": "cGAS-STING, type_I_interferon",
        "associated_diseases": ["SAVI"],
        "conditions": ["autoinflammatory", "interferonopathy", "vasculopathy"],
    },
    "TREX1": {
        "chrom": "chr3",
        "start": 48465830,
        "end": 48467645,
        "role": "3'-5' DNA exonuclease — cytosolic DNA degradation, "
                "prevents cGAS-STING autoactivation",
        "strategy": "Diagnostic: screen single exon. Therapeutic: "
                     "base-editing to correct loss-of-function",
        "omim_gene": 606609,
        "inheritance": "AR/AD",
        "pathway": "cGAS-STING, type_I_interferon, DNA_metabolism",
        "associated_diseases": ["AGS1", "SLE", "FCL"],
        "conditions": ["interferonopathy", "autoimmune", "SLE"],
    },

    # ---------------------------------------------------------------------------
    # IL-1 family / IL-36 pathway
    # ---------------------------------------------------------------------------
    "IL1RN": {
        "chrom": "chr2",
        "start": 113099360,
        "end": 113134014,
        "role": "IL-1 receptor antagonist — blocks IL-1alpha/beta signaling",
        "strategy": "Diagnostic: screen for deletions and truncating variants. "
                     "Therapeutic: gene replacement via HDR at safe harbor",
        "omim_gene": 147679,
        "inheritance": "AR",
        "pathway": "IL-1_signaling",
        "associated_diseases": ["DIRA"],
        "conditions": ["autoinflammatory", "neonatal_onset"],
    },
    "IL36RN": {
        "chrom": "chr2",
        "start": 113058638,
        "end": 113064744,
        "role": "IL-36 receptor antagonist — blocks IL-36 inflammatory signaling",
        "strategy": "Diagnostic: screen for L27P and other loss-of-function. "
                     "Therapeutic: base-editing to correct pathogenic variants",
        "omim_gene": 605507,
        "inheritance": "AR",
        "pathway": "IL-36_signaling",
        "associated_diseases": ["DITRA"],
        "conditions": ["autoinflammatory", "pustular_psoriasis"],
    },

    # ---------------------------------------------------------------------------
    # Ubiquitin / somatic
    # ---------------------------------------------------------------------------
    "UBA1": {
        "chrom": "chrX",
        "start": 47190847,
        "end": 47215128,
        "role": "E1 ubiquitin-activating enzyme — ubiquitin-proteasome "
                "pathway initiation",
        "strategy": "Diagnostic: screen Met41 codon for somatic variants "
                     "(p.Met41Thr/Val/Leu). Note: somatic mosaicism "
                     "requires high-depth sequencing (>100x)",
        "omim_gene": 314370,
        "inheritance": "somatic_X-linked",
        "pathway": "ubiquitin-proteasome, myeloid_inflammation",
        "associated_diseases": ["VEXAS"],
        "conditions": ["autoinflammatory", "VEXAS", "myelodysplastic"],
    },

    # ---------------------------------------------------------------------------
    # Other monogenic autoinflammatory
    # ---------------------------------------------------------------------------
    "ADA2": {
        "chrom": "chr22",
        "start": 17178790,
        "end": 17221848,
        "role": "Adenosine deaminase 2 — extracellular adenosine metabolism, "
                "endothelial integrity",
        "strategy": "Diagnostic: full coding screen. Therapeutic: "
                     "HDR correction or gene replacement",
        "omim_gene": 607575,
        "aliases": ["CECR1"],
        "inheritance": "AR",
        "pathway": "adenosine_metabolism, endothelial_homeostasis",
        "associated_diseases": ["DADA2"],
        "conditions": ["autoinflammatory", "vasculitis", "stroke"],
    },
    "PLCG2": {
        "chrom": "chr16",
        "start": 81739041,
        "end": 81962685,
        "role": "Phospholipase C gamma-2 — B-cell, NK-cell, and mast cell "
                "signaling downstream of immune receptors",
        "strategy": "Diagnostic: screen for genomic deletions (PLAID) and "
                     "gain-of-function missense (APLAID). Therapeutic: "
                     "CRISPRi for gain-of-function",
        "omim_gene": 600220,
        "inheritance": "AD",
        "pathway": "PLC_signaling, B-cell_activation",
        "associated_diseases": ["PLAID", "APLAID"],
        "conditions": ["autoinflammatory", "immune_dysregulation"],
    },
    "PSTPIP1": {
        "chrom": "chr15",
        "start": 76994679,
        "end": 77037474,
        "role": "PEST-domain phosphatase interacting protein 1 — pyrin "
                "interaction, actin remodeling",
        "strategy": "Diagnostic: screen for A230T, E250Q. Therapeutic: "
                     "base-editing to correct gain-of-function",
        "omim_gene": 606347,
        "inheritance": "AD",
        "pathway": "pyrin_inflammasome, actin_remodeling",
        "associated_diseases": ["PAPA"],
        "conditions": ["autoinflammatory", "pyoderma_gangrenosum"],
    },
    "LPIN2": {
        "chrom": "chr18",
        "start": 2916994,
        "end": 3013144,
        "role": "Lipin-2 — phospholipid metabolism, inflammasome regulation",
        "strategy": "Diagnostic: full coding screen. Therapeutic: "
                     "HDR correction of loss-of-function variants",
        "omim_gene": 605519,
        "inheritance": "AR",
        "pathway": "phospholipid_metabolism, inflammasome",
        "associated_diseases": ["majeed_syndrome"],
        "conditions": ["autoinflammatory", "osteomyelitis"],
    },
}


# ---------------------------------------------------------------------------
# Digenic / modifier gene combinations relevant to Yao syndrome
# (22% of NOD2-positive YAOS patients carry variants in these genes)
# ---------------------------------------------------------------------------
YAO_MODIFIER_GENES = ["MEFV", "NLRP3", "NLRP12", "TNFRSF1A"]

YAO_DIGENIC_PATTERNS = {
    "NOD2_MEFV": {
        "frequency": "most common digenic combination",
        "example": "IVS8+158 + MEFV E148Q",
        "clinical_note": "May present with mixed FMF/YAOS features",
    },
    "NOD2_NLRP12": {
        "frequency": "common",
        "example": "IVS8+158/R702W + NLRP12 F402L",
        "clinical_note": "Cold-triggered flares may be more prominent",
    },
    "NOD2_NLRP3": {
        "frequency": "common",
        "example": "IVS8+158 + NLRP3 Q705K",
        "clinical_note": "Urticaria-like rash may be more prominent",
    },
    "NOD2_TNFRSF1A": {
        "frequency": "less common",
        "example": "Various combinations",
        "clinical_note": "Prolonged fever episodes",
    },
}


# ---------------------------------------------------------------------------
# Disease-to-gene mapping for diagnostic routing
# ---------------------------------------------------------------------------
AUTOINFLAMMATORY_DISEASES = {
    "yao_syndrome": {
        "primary_gene": "NOD2",
        "modifier_genes": ["MEFV", "NLRP3", "NLRP12", "TNFRSF1A"],
        "omim": 617321,
        "inheritance": "complex/GTD",
        "key_pathway": "NF-kB, p38_MAPK, IL-6",
        "treatment_targets": ["IL-1_inhibitor", "IL-6_inhibitor",
                              "JAK_inhibitor", "glucocorticoid"],
    },
    "FMF": {
        "primary_gene": "MEFV",
        "omim": 249100,
        "inheritance": "AR",
        "key_pathway": "pyrin_inflammasome",
        "treatment_targets": ["colchicine", "IL-1_inhibitor"],
    },
    "CAPS": {
        "primary_gene": "NLRP3",
        "omim": 120100,
        "inheritance": "AD",
        "key_pathway": "NLRP3_inflammasome",
        "treatment_targets": ["IL-1_inhibitor"],
    },
    "TRAPS": {
        "primary_gene": "TNFRSF1A",
        "omim": 142680,
        "inheritance": "AD",
        "key_pathway": "TNF_signaling",
        "treatment_targets": ["IL-1_inhibitor", "anti-TNF"],
    },
    "blau_syndrome": {
        "primary_gene": "NOD2",
        "omim": 186580,
        "inheritance": "AD",
        "key_pathway": "NF-kB_constitutive",
        "treatment_targets": ["anti-TNF", "IL-1_inhibitor"],
    },
    "VEXAS": {
        "primary_gene": "UBA1",
        "omim": 301054,
        "inheritance": "somatic_X-linked",
        "key_pathway": "ubiquitin-proteasome",
        "treatment_targets": ["JAK_inhibitor", "azacitidine"],
    },
    "DADA2": {
        "primary_gene": "ADA2",
        "omim": 615688,
        "inheritance": "AR",
        "key_pathway": "adenosine_metabolism",
        "treatment_targets": ["anti-TNF"],
    },
    "HA20": {
        "primary_gene": "TNFAIP3",
        "omim": 616744,
        "inheritance": "AD",
        "key_pathway": "NF-kB_regulation",
        "treatment_targets": ["anti-TNF", "JAK_inhibitor"],
    },
    "SAVI": {
        "primary_gene": "STING1",
        "omim": 615934,
        "inheritance": "AD",
        "key_pathway": "cGAS-STING_interferon",
        "treatment_targets": ["JAK_inhibitor"],
    },
    "MKD": {
        "primary_gene": "MVK",
        "omim": 260920,
        "inheritance": "AR",
        "key_pathway": "mevalonate",
        "treatment_targets": ["IL-1_inhibitor"],
    },
}
