"""
Curated cancer gene targets — GRCh38 coordinates.

Sources: COSMIC Cancer Gene Census, OncoKB, ClinVar.
Each entry includes the gene's role, common mutation types in tumors,
and the therapeutic CRISPR strategy (knockout of oncogene, restoration of
tumor suppressor, etc.).
"""

CANCER_TARGETS = {
    # --- Tumor suppressors (goal: restore / reactivate) ---
    "TP53": {
        "chrom": "chr17",
        "start": 7668402,
        "end": 7687550,
        "role": "Tumor suppressor — cell-cycle arrest and apoptosis",
        "strategy": "CRISPRa to reactivate silenced wild-type allele, or "
                     "base-editing to correct hotspot missense mutations",
        "common_cancers": ["breast", "lung", "colorectal", "ovarian", "pancreatic"],
        "mutation_types": ["missense", "nonsense", "frameshift"],
    },
    "BRCA1": {
        "chrom": "chr17",
        "start": 43044295,
        "end": 43170245,
        "role": "DNA double-strand break repair via homologous recombination",
        "strategy": "HDR-mediated correction of pathogenic variants",
        "common_cancers": ["breast", "ovarian"],
        "mutation_types": ["frameshift", "nonsense", "splice_site"],
    },
    "BRCA2": {
        "chrom": "chr13",
        "start": 32315086,
        "end": 32400266,
        "role": "DNA repair — homologous recombination",
        "strategy": "HDR-mediated correction of pathogenic variants",
        "common_cancers": ["breast", "ovarian", "prostate", "pancreatic"],
        "mutation_types": ["frameshift", "nonsense"],
    },
    "PTEN": {
        "chrom": "chr10",
        "start": 87863113,
        "end": 87971930,
        "role": "PI3K/AKT pathway negative regulator",
        "strategy": "CRISPRa to restore expression in PTEN-deleted tumors",
        "common_cancers": ["glioblastoma", "prostate", "endometrial"],
        "mutation_types": ["deletion", "nonsense", "missense"],
    },
    "RB1": {
        "chrom": "chr13",
        "start": 48303747,
        "end": 48481890,
        "role": "Retinoblastoma protein — cell-cycle G1/S checkpoint",
        "strategy": "CRISPRa to reactivate or base-editing to correct mutations",
        "common_cancers": ["retinoblastoma", "lung_small_cell", "bladder"],
        "mutation_types": ["nonsense", "deletion", "splice_site"],
    },
    "APC": {
        "chrom": "chr5",
        "start": 112707498,
        "end": 112846239,
        "role": "Wnt signaling negative regulator",
        "strategy": "Base-editing to correct truncating mutations",
        "common_cancers": ["colorectal", "desmoid"],
        "mutation_types": ["nonsense", "frameshift"],
    },
    "VHL": {
        "chrom": "chr3",
        "start": 10141778,
        "end": 10153670,
        "role": "HIF degradation — oxygen sensing pathway",
        "strategy": "HDR correction of pathogenic VHL mutations",
        "common_cancers": ["renal_clear_cell"],
        "mutation_types": ["missense", "deletion", "nonsense"],
    },

    # --- Oncogenes (goal: disrupt / silence) ---
    "KRAS": {
        "chrom": "chr12",
        "start": 25204789,
        "end": 25250936,
        "role": "RAS/MAPK signaling — cell proliferation",
        "strategy": "CRISPRi silencing of mutant allele or allele-specific "
                     "Cas13 knockdown of mutant mRNA",
        "common_cancers": ["pancreatic", "lung", "colorectal"],
        "mutation_types": ["missense_hotspot"],
    },
    "EGFR": {
        "chrom": "chr7",
        "start": 55019017,
        "end": 55211628,
        "role": "Receptor tyrosine kinase — growth signaling",
        "strategy": "CRISPRi silencing or allele-specific knockout of "
                     "activating mutations (L858R, exon 19 del, T790M)",
        "common_cancers": ["lung_NSCLC", "glioblastoma"],
        "mutation_types": ["missense_hotspot", "in_frame_deletion"],
    },
    "BRAF": {
        "chrom": "chr7",
        "start": 140719327,
        "end": 140924929,
        "role": "RAF/MAPK signaling kinase",
        "strategy": "Allele-specific CRISPRi of V600E mutant allele",
        "common_cancers": ["melanoma", "colorectal", "thyroid"],
        "mutation_types": ["missense_hotspot"],
    },
    "PIK3CA": {
        "chrom": "chr3",
        "start": 179148114,
        "end": 179240093,
        "role": "PI3K catalytic subunit — growth/survival signaling",
        "strategy": "Allele-specific knockout of activating mutations",
        "common_cancers": ["breast", "endometrial", "colorectal"],
        "mutation_types": ["missense_hotspot"],
    },
    "MYC": {
        "chrom": "chr8",
        "start": 127735434,
        "end": 127742951,
        "role": "Transcription factor — proliferation and metabolism",
        "strategy": "CRISPRi silencing of amplified MYC",
        "common_cancers": ["burkitt_lymphoma", "breast", "lung"],
        "mutation_types": ["amplification", "translocation"],
    },
    "HER2": {
        "chrom": "chr17",
        "start": 39687914,
        "end": 39730426,
        "role": "ERBB2 receptor tyrosine kinase — growth signaling",
        "strategy": "CRISPRi silencing in HER2-amplified tumors",
        "common_cancers": ["breast", "gastric"],
        "mutation_types": ["amplification", "missense"],
    },
    "BCL2": {
        "chrom": "chr18",
        "start": 63123346,
        "end": 63320128,
        "role": "Anti-apoptotic protein",
        "strategy": "CRISPRi to reduce overexpression and restore apoptosis",
        "common_cancers": ["follicular_lymphoma", "CLL"],
        "mutation_types": ["translocation", "overexpression"],
    },
    "ALK": {
        "chrom": "chr2",
        "start": 29192774,
        "end": 29921566,
        "role": "Receptor tyrosine kinase — frequently fused in cancer",
        "strategy": "CRISPR disruption of fusion breakpoint",
        "common_cancers": ["lung_NSCLC", "neuroblastoma", "ALCL"],
        "mutation_types": ["fusion", "missense"],
    },

    # --- DNA repair / genome stability ---
    "ATM": {
        "chrom": "chr11",
        "start": 108222484,
        "end": 108369102,
        "role": "DNA damage response kinase",
        "strategy": "HDR correction of loss-of-function mutations",
        "common_cancers": ["breast", "pancreatic", "CLL"],
        "mutation_types": ["nonsense", "frameshift", "splice_site"],
    },
    "PALB2": {
        "chrom": "chr16",
        "start": 23603160,
        "end": 23641310,
        "role": "BRCA2 partner — homologous recombination repair",
        "strategy": "HDR correction of truncating mutations",
        "common_cancers": ["breast", "pancreatic"],
        "mutation_types": ["frameshift", "nonsense"],
    },

    # --- DNA mismatch repair / Lynch syndrome ---
    "MLH1": {
        "chrom": "chr3",
        "start": 36993332,
        "end": 37050918,
        "role": "MutL homolog 1 — DNA mismatch repair",
        "strategy": "CRISPRa to reactivate epigenetically silenced MLH1 "
                     "(promoter hypermethylation is the most common cause); "
                     "HDR correction of germline truncating mutations",
        "common_cancers": ["colorectal", "endometrial", "ovarian", "gastric"],
        "mutation_types": ["nonsense", "frameshift", "splice_site",
                           "promoter_methylation"],
    },
    "MSH2": {
        "chrom": "chr2",
        "start": 47403067,
        "end": 47709830,
        "role": "MutS homolog 2 — DNA mismatch repair",
        "strategy": "HDR correction of truncating mutations; gene replacement "
                     "feasible (cDNA ~2.8 kb)",
        "common_cancers": ["colorectal", "endometrial", "ovarian",
                           "urinary_tract"],
        "mutation_types": ["nonsense", "frameshift", "large_deletion"],
    },
    "MSH6": {
        "chrom": "chr2",
        "start": 47695530,
        "end": 47721031,
        "role": "MutS homolog 6 — DNA mismatch repair (MSH2-MSH6 heterodimer)",
        "strategy": "HDR correction of truncating mutations",
        "common_cancers": ["colorectal", "endometrial"],
        "mutation_types": ["missense", "nonsense", "frameshift"],
    },
    "PMS2": {
        "chrom": "chr7",
        "start": 5970925,
        "end": 6009106,
        "role": "PMS1 homolog 2 — DNA mismatch repair (MLH1-PMS2 heterodimer)",
        "strategy": "HDR correction; note pseudogene PMS2CL complicates "
                     "sequencing — long-read or gene-specific PCR required",
        "common_cancers": ["colorectal", "endometrial"],
        "mutation_types": ["nonsense", "frameshift", "splice_site"],
    },
}
