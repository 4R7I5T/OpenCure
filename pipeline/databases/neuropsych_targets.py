"""
Neuropsychiatric gene targets for circuit-level modulation — GRCh38 coordinates.

These targets are organized by functional layer for a multi-construct
neuromodulation strategy combining CRISPRi/CRISPRa with chemogenetic
(DREADD) control.

Layers:
  1. Neurotransmitter receptor/transporter modulation
  2. Neurohormone receptor sensitivity
  3. Transcription factor and neurotrophic modulation
  4. On-demand chemogenetic control (DREADD)

All interventions require informed patient consent and IRB/ethics approval.

Sources: OMIM, GeneCards, UniProt, Allen Brain Atlas.
"""

# ---------------------------------------------------------------------------
# Layer 1: Neurotransmitter receptor / transporter modulation
# ---------------------------------------------------------------------------
NEUROPSYCH_TARGETS = {
    "CHRM1": {
        "chrom": "chr11",
        "start": 62432780,
        "end": 62460950,
        "role": "Muscarinic M1 receptor — cholinergic neurotransmission",
        "strategy": "CRISPRi silencing to reduce M1 receptor expression in "
                     "DLPFC and hippocampus",
        "layer": "neurotransmitter_receptor",
        "modality": "CRISPRi",
        "target_regions": ["DLPFC", "hippocampus"],
        "conditions": ["cognitive_impairment", "schizophrenia"],
        "key_variants": ["rs2067477", "rs2075748"],
    },
    "COMT": {
        "chrom": "chr22",
        "start": 19929262,
        "end": 19969975,
        "role": "Catechol-O-methyltransferase — prefrontal dopamine degradation",
        "strategy": "CRISPRa activation to increase COMT expression and "
                     "normalize dopamine turnover in DLPFC",
        "layer": "neurotransmitter_degradation",
        "modality": "CRISPRa",
        "target_regions": ["DLPFC"],
        "conditions": ["schizophrenia", "ADHD", "cognitive_impairment"],
        "key_variants": ["rs4680_Val158Met"],
    },
    "CACNA1C": {
        "chrom": "chr12",
        "start": 2161515,
        "end": 2322801,
        "role": "L-type voltage-gated calcium channel alpha-1C subunit — "
                "neuronal excitability and synaptic plasticity",
        "strategy": "CRISPRi silencing to attenuate gain-of-function variants "
                     "linked to psychiatric risk",
        "layer": "ion_channel",
        "modality": "CRISPRi",
        "target_regions": ["hippocampus", "cortex"],
        "conditions": ["bipolar_disorder", "schizophrenia",
                       "major_depressive_disorder"],
        "key_variants": ["rs1006737", "rs4765905"],
    },

    # ---------------------------------------------------------------------------
    # Layer 2: Neurohormone receptor sensitivity
    # ---------------------------------------------------------------------------
    "OXTR": {
        "chrom": "chr3",
        "start": 8775433,
        "end": 8788544,
        "role": "Oxytocin receptor — social cognition and bonding",
        "strategy": "CRISPRa activation to increase oxytocin receptor density "
                     "in hypothalamus and amygdala",
        "layer": "neurohormone_receptor",
        "modality": "CRISPRa",
        "target_regions": ["hypothalamus", "amygdala"],
        "conditions": ["autism_spectrum", "social_anxiety",
                       "attachment_disorders"],
        "key_variants": ["rs53576", "rs2254298"],
    },
    "SLC6A4": {
        "chrom": "chr17",
        "start": 30194319,
        "end": 30236002,
        "role": "Serotonin transporter (SERT) — serotonin reuptake",
        "strategy": "CRISPRi silencing to reduce SERT expression in dorsal "
                     "raphe, increasing synaptic serotonin (analogous to SSRI)",
        "layer": "neurotransmitter_transporter",
        "modality": "CRISPRi",
        "target_regions": ["dorsal_raphe"],
        "conditions": ["major_depressive_disorder", "anxiety", "OCD"],
        "key_variants": ["5-HTTLPR_short", "rs25531"],
    },
    "OPRM1": {
        "chrom": "chr6",
        "start": 154039662,
        "end": 154106498,
        "role": "Mu-opioid receptor — pain and reward processing",
        "strategy": "CRISPRa activation of mu-opioid receptor to restore "
                     "endogenous opioid sensitivity",
        "layer": "opioid_receptor",
        "modality": "CRISPRa",
        "target_regions": ["periaqueductal_gray", "nucleus_accumbens"],
        "conditions": ["opioid_use_disorder", "chronic_pain",
                       "anhedonia"],
        "key_variants": ["rs1799971_A118G"],
    },

    # ---------------------------------------------------------------------------
    # Layer 3: Transcription factor and neurotrophic modulation
    # ---------------------------------------------------------------------------
    "CREB1": {
        "chrom": "chr2",
        "start": 207529963,
        "end": 207605989,
        "role": "cAMP response element-binding protein — activity-dependent "
                "transcription",
        "strategy": "CRISPRi silencing to modulate CREB-dependent transcription "
                     "in hippocampus CA1",
        "layer": "transcription_factor",
        "modality": "CRISPRi",
        "target_regions": ["hippocampus_CA1"],
        "conditions": ["anxiety", "fear_conditioning", "addiction"],
        "key_variants": ["rs2253206", "rs2551922"],
    },
    "BDNF": {
        "chrom": "chr11",
        "start": 27654893,
        "end": 27722058,
        "role": "Brain-derived neurotrophic factor — synaptic plasticity "
                "and neuronal survival",
        "strategy": "CRISPRa activation to enhance BDNF expression in "
                     "treatment-resistant depression",
        "layer": "neurotrophic_factor",
        "modality": "CRISPRa",
        "target_regions": ["hippocampus", "prefrontal_cortex"],
        "conditions": ["treatment_resistant_depression",
                       "cognitive_decline", "PTSD"],
        "key_variants": ["rs6265_Val66Met"],
    },
    "CHAT": {
        "chrom": "chr10",
        "start": 123239119,
        "end": 123282787,
        "role": "Choline acetyltransferase — acetylcholine synthesis",
        "strategy": "CRISPRa activation to boost cholinergic tone in "
                     "basal forebrain circuits",
        "layer": "cholinergic_modulation",
        "modality": "CRISPRa",
        "target_regions": ["basal_forebrain", "hippocampus"],
        "conditions": ["cognitive_decline", "alzheimers_prodromal"],
        "key_variants": ["rs3810950", "rs868750"],
    },
}


# ---------------------------------------------------------------------------
# Layer 4: DREADD chemogenetic control constructs
# ---------------------------------------------------------------------------
DREADD_CONSTRUCTS = {
    "DREADD_inhibitory": {
        "type": "DREADD_insertion",
        "receptor": "hM4Di",
        "reporter": "IRES-mCherry",
        "promoter": "CamKII",
        "backbone": "AAV9",
        "target_regions": ["DLPFC", "hippocampus"],
        "activation_ligand": "deschloroclozapine (DCZ)",
        "purpose": "On-demand silencing of targeted excitatory neurons via "
                   "Gi-coupled DREADD. Provides reversible external control "
                   "switch independent of CRISPR constructs.",
        "mechanism": "DCZ activates hM4Di -> Gi signaling -> reduced neuronal "
                     "firing via GIRK channel opening and cAMP reduction",
        "crispr_role": "Cas9 HDR inserts DREADD transgene at AAVS1 safe harbor",
    },
    "DREADD_excitatory": {
        "type": "DREADD_insertion",
        "receptor": "hM3Dq",
        "reporter": "IRES-GFP",
        "promoter": "OXT",
        "backbone": "AAV9",
        "target_regions": ["paraventricular_nucleus"],
        "activation_ligand": "deschloroclozapine (DCZ)",
        "purpose": "On-demand activation of oxytocin neurons to enhance "
                   "pro-social signaling on demand.",
        "mechanism": "DCZ activates hM3Dq -> Gq signaling -> increased "
                     "neuronal firing via PLC/IP3/DAG cascade",
        "crispr_role": "Cas9 HDR inserts DREADD transgene at ROSA26 safe harbor",
    },
}


# ---------------------------------------------------------------------------
# Safe harbor loci for transgene insertion (DREADD / reporter cassettes)
# ---------------------------------------------------------------------------
SAFE_HARBOR_LOCI = {
    "AAVS1": {
        "chrom": "chr19",
        "start": 55115756,
        "end": 55115856,
        "role": "PPP1R12C intron 1 — validated safe harbor for transgene "
                "insertion with stable expression",
        "strategy": "HDR-mediated DREADD transgene insertion",
    },
    "ROSA26_human": {
        "chrom": "chr3",
        "start": 9396280,
        "end": 9396380,
        "role": "Human ROSA26 equivalent — alternative safe harbor for "
                "transgene insertion",
        "strategy": "HDR-mediated DREADD transgene insertion",
    },
}


# ---------------------------------------------------------------------------
# CRISPRi/CRISPRa promoter target regions (promoter + first exon)
# Used by guide design to know where to target dCas9-KRAB / dCas9-VPR
# ---------------------------------------------------------------------------
CRISPRI_TARGET_REGIONS = {
    gene: {
        "chrom": info["chrom"],
        "start": max(1, info["start"] - 1500),
        "end": info["start"] + 500,
        "modality": info["modality"],
        "purpose": info["strategy"],
    }
    for gene, info in NEUROPSYCH_TARGETS.items()
}


# Tissue-specific promoters appropriate for neuropsychiatric targets
NEURO_PROMOTERS = {
    "CamKII":  {"size_bp": 1300, "specificity": "excitatory_neurons"},
    "hSyn1":   {"size_bp": 470,  "specificity": "pan_neuronal"},
    "TPH2":    {"size_bp": 1800, "specificity": "serotonergic_neurons"},
    "OXT":     {"size_bp": 900,  "specificity": "oxytocin_neurons"},
    "TH":      {"size_bp": 1100, "specificity": "dopaminergic_neurons"},
    "GAD67":   {"size_bp": 1400, "specificity": "GABAergic_neurons"},
    "ChAT":    {"size_bp": 1200, "specificity": "cholinergic_neurons"},
}
