"""
Hematological genetic disease targets for CRISPR/gene-therapy — GRCh38 coordinates.

Covers sickle cell disease, beta-thalassemia, hemophilia A/B, Fanconi anemia,
Diamond-Blackfan anemia, and pyruvate kinase deficiency.

All coordinates are GRCh38/hg38. Gene strand and Ensembl IDs are included for
pipeline cross-referencing.  Variant positions are reported on the **plus strand**
of the reference genome unless otherwise noted.

Sources:
  - Ensembl REST API (GRCh38, accessed April 2026)
  - ClinVar / dbSNP (NCBI)
  - Frangoul et al., NEJM 2021 (Casgevy / exa-cel)
  - Musallam et al., Nature Medicine 2022 (BCL11A enhancer editing)
  - Pipe et al., NEJM 2023 (fidanacogene elaparvovec)
  - Ozelo et al., NEJM 2022 (etranacogene dezaparvovec)
  - Base-editing of HBG1/HBG2, NEJM April 2026 (risto-cel / BEAM-101)

All interventions require informed patient consent and IRB/ethics approval.
"""

# ===========================================================================
#  BETA-GLOBIN GENE CLUSTER  (chr11p15.4, minus strand)
#
#  5' ──LCR──HBE1──HBG2──HBG1──HBD──HBB── 3'  (gene order, 5'→3' on mRNA)
#  Higher genomic coords ◄────────────────► Lower genomic coords  (minus strand)
# ===========================================================================

BETA_GLOBIN_CLUSTER = {
    "description": "Human beta-globin gene cluster on chr11p15.4 (minus strand)",
    "chrom": "chr11",
    "strand": "-1",
    "cluster_region": {
        "start": 5225464,   # HBB 3' end
        "end": 5313839,     # approximate LCR 5' boundary (HS5)
    },

    # --- Locus Control Region (LCR) ---
    # HS positions are approximate, derived from HBE1 canonical TSS (5269945)
    # and published relative distances (-6, -11, -15, -18, -22 kb).
    "LCR": {
        "description": "Beta-globin Locus Control Region — five DNase I "
                       "hypersensitive sites (HS1-HS5) upstream of HBE1",
        "NCBI_gene_id": 109580095,
        "approximate_region": {"start": 5276000, "end": 5313839},
        "hypersensitive_sites": {
            "HS1": {"approx_center": 5275945, "relative": "-6 kb from HBE1 TSS",
                    "function": "Weak enhancer activity"},
            "HS2": {"approx_center": 5280645, "relative": "-10.7 kb from HBE1 TSS",
                    "function": "Strong enhancer, NF-E2/AP-1 binding"},
            "HS3": {"approx_center": 5284645, "relative": "-14.7 kb from HBE1 TSS",
                    "function": "Enhancer, GATA-1 binding"},
            "HS4": {"approx_center": 5288245, "relative": "-18.3 kb from HBE1 TSS",
                    "function": "Insulator/enhancer, CTCF binding"},
            "HS5": {"approx_center": 5291945, "relative": "-22 kb from HBE1 TSS",
                    "function": "Insulator, CTCF-dependent boundary"},
        },
    },

    # --- Individual genes (canonical transcript coordinates) ---
    "genes": {
        "HBE1": {
            "name": "Hemoglobin subunit epsilon 1",
            "ensembl_gene": "ENSG00000213931",
            "ensembl_transcript": "ENST00000396895",
            "chrom": "chr11", "strand": -1,
            "tx_start": 5268345, "tx_end": 5269945,
            "expression_stage": "Embryonic (yolk sac)",
        },
        "HBG2": {
            "name": "Hemoglobin subunit gamma 2 (G-gamma, Gly-136)",
            "ensembl_gene": "ENSG00000196565",
            "ensembl_transcript": "ENST00000336906",
            "chrom": "chr11", "strand": -1,
            "tx_start": 5253188, "tx_end": 5254781,
            "expression_stage": "Fetal",
        },
        "HBG1": {
            "name": "Hemoglobin subunit gamma 1 (A-gamma, Ala-136)",
            "ensembl_gene": "ENSG00000213934",
            "ensembl_transcript": "ENST00000330597",
            "chrom": "chr11", "strand": -1,
            "tx_start": 5248269, "tx_end": 5249857,
            "expression_stage": "Fetal",
        },
        "HBD": {
            "name": "Hemoglobin subunit delta",
            "ensembl_gene": "ENSG00000223609",
            "ensembl_transcript": "ENST00000650601",
            "chrom": "chr11", "strand": -1,
            "tx_start": 5232838, "tx_end": 5234483,
            "expression_stage": "Adult (minor, ~3%)",
        },
        "HBB": {
            "name": "Hemoglobin subunit beta",
            "ensembl_gene": "ENSG00000244734",
            "ensembl_transcript": "ENST00000335295",
            "chrom": "chr11", "strand": -1,
            "tx_start": 5225464, "tx_end": 5227071,
            "exons": [
                {"exon": 1, "start": 5226930, "end": 5227071},
                {"exon": 2, "start": 5226577, "end": 5226799},
                {"exon": 3, "start": 5225464, "end": 5225726},
            ],
            "expression_stage": "Adult (major, ~97%)",
        },
    },
}


# ===========================================================================
#  1.  SICKLE CELL DISEASE  (SCD)
# ===========================================================================

SICKLE_CELL_DISEASE = {
    "disease": "Sickle Cell Disease",
    "inheritance": "Autosomal recessive",
    "omim": 603903,

    # ----- Primary causative variant -----
    "causative_variant": {
        "gene": "HBB",
        "name": "HbS / E6V (Glu6Val)",
        "hgvs_c": "NM_000518.5:c.20A>T",
        "hgvs_p": "p.Glu7Val",  # includes Met; mature protein = Glu6Val
        "rsid": "rs334",
        "chrom": "chr11",
        "pos_grch38": 5227002,
        "ref": "T", "alt": "A",  # plus-strand; gene on minus strand → coding A>T
        "allele_string": "T/A",
        "clinical_significance": ["pathogenic"],
        "consequence": "missense_variant",
        "protein_effect": "Glu→Val at position 6 of mature beta-globin causes "
                          "HbS polymerization under deoxygenation",
    },

    # ----- Therapeutic targets & strategies -----
    "therapeutic_targets": {

        # (A) BCL11A enhancer disruption  — Casgevy / exa-cel (APPROVED)
        "BCL11A_enhancer": {
            "gene": "BCL11A",
            "ensembl_gene": "ENSG00000119866",
            "gene_coords": {"chrom": "chr2", "start": 60450520, "end": 60554467,
                            "strand": -1},
            "target_region": "+58 DHS erythroid enhancer (intron 2)",
            "enhancer_coords_grch38": {
                "broad_region": {"start": 60494251, "end": 60495546},
                "functional_core": "~115 bp within the broad region",
                "GATA1_motif": "CTG(n9)GATA composite element within core",
            },
            "other_DHS": {
                "+55": "Additional DHS site, synergistic with +58",
                "+62": "Most distal erythroid DHS",
            },
            "strategy": "CRISPR-Cas9 knockout of GATA1 binding site in the "
                         "+58 erythroid enhancer → loss of BCL11A expression "
                         "in erythroid cells → de-repression of HBG1/HBG2 → "
                         "HbF reactivation (pancellular) → anti-sickling",
            "therapy_name": "Casgevy (exagamglogene autotemcel / exa-cel)",
            "developer": "Vertex / CRISPR Therapeutics",
            "clinical_status": "FDA-approved Dec 2023 (SCD), EMA-approved Feb 2024",
            "trial_id": "NCT03745287",
            "mechanism_detail": "sgRNA directs Cas9 to the GATA1 motif in the "
                                "+58 enhancer; indels disrupt GATA1/TAL1 binding "
                                "→ BCL11A silenced specifically in erythroid "
                                "lineage (not in B-cells/neurons where BCL11A "
                                "has other roles)",
            "pediatric_expansion": "Regulatory submissions for ages 5-11 expected "
                                   "H1 2026",
        },

        # (B) HBG1/HBG2 promoter base editing  — BEAM-101 / risto-cel
        "HBG_promoter_base_editing": {
            "genes": ["HBG1", "HBG2"],
            "target_region": "Promoter BCL11A binding site (-198 to -195 region)",
            "key_HPFH_variants": {
                "-198_T>C": {
                    "gene": "HBG1",
                    "rsid": "rs35710727",
                    "pos_grch38": 5250055,  # HBG1 promoter
                    "ref": "A", "alt": "G",  # plus strand
                    "gene_strand_change": "T>C (British HPFH)",
                    "effect": "Disrupts BCL11A repressor binding → HbF induction",
                    "clinical_significance": ["pathogenic"],  # causes HPFH
                },
                "-196_C>T": {
                    "gene": "HBG1",
                    "rsid": "rs35983258",
                    "pos_grch38": 5250053,
                    "ref": "G", "alt": "A",
                    "gene_strand_change": "C>T",
                    "effect": "Disrupts BCL11A binding → HPFH",
                    "clinical_significance": ["pathogenic"],
                },
                "-195_C>G": {
                    "gene": "HBG1",
                    "rsid": "rs35321913",
                    "pos_grch38": 5250052,
                    "ref": "G", "alt": "C",
                    "gene_strand_change": "C>G",
                    "effect": "Disrupts BCL11A binding → HPFH",
                    "clinical_significance": ["pathogenic"],
                },
                "-175_T>C": {
                    "gene": "HBG1",
                    "rsid": "rs34153902",
                    "pos_grch38": 5250032,
                    "ref": "A", "alt": "G",  # plus strand
                    "gene_strand_change": "T>C (creates de novo TAL1 binding site)",
                    "effect": "Creates TAL1/E-box motif → LCR recruitment → "
                              "HbF activation (Black HPFH)",
                },
            },
            "strategy": "Adenine base editing (ABE8e) to install A>G edits at "
                         "BCL11A binding sites in HBG1/2 promoters — recreates "
                         "natural HPFH mutations without double-strand breaks",
            "therapy_name": "ristoglogene autogetemcel (risto-cel, formerly BEAM-101)",
            "developer": "Beam Therapeutics",
            "clinical_status": "Phase 1/2 (BEACON trial, NCT05456880); "
                               "NEJM publication April 2026; FDA Orphan Drug "
                               "Designation granted",
            "trial_id": "NCT05456880",
            "clinical_results": "17 patients dosed as of Feb 2025; HbF >60%, "
                                "HbS <40% in all evaluable patients; 26 patients "
                                "dosed as of June 2025",
        },

        # (C) Direct SCD correction  — BEAM-102
        "HbS_direct_correction": {
            "gene": "HBB",
            "target": "rs334 (E6V) direct correction to HbG-Makassar (E6V→E6A)",
            "strategy": "Base editing to convert HbS codon to HbG-Makassar — "
                         "a benign non-sickling hemoglobin variant",
            "therapy_name": "BEAM-102",
            "developer": "Beam Therapeutics",
            "clinical_status": "Preclinical (as of 2025)",
        },

        # (D) Other therapeutic targets for HbF induction
        "KLF1": {
            "gene": "KLF1",
            "ensembl_gene": "ENSG00000105610",
            "chrom": "chr19", "start": 12884422, "end": 12887201, "strand": -1,
            "role": "KLF transcription factor 1 (erythroid Kruppel-like factor) "
                    "— activates adult beta-globin and ZBTB7A expression; "
                    "haploinsufficiency causes HPFH",
            "strategy": "CRISPRi/partial knockdown to reduce KLF1 → lower "
                         "BCL11A and ZBTB7A → HbF de-repression. Caution: "
                         "complete KO causes severe anemia",
            "clinical_status": "Research/preclinical",
        },
        "ZBTB7A_LRF": {
            "gene": "ZBTB7A (LRF)",
            "ensembl_gene": "ENSG00000178951",
            "chrom": "chr19", "start": 4043303, "end": 4067636, "strand": -1,
            "role": "BCL11A-independent repressor of fetal hemoglobin; binds "
                    "directly to HBG1/2 promoters at a site distinct from BCL11A",
            "strategy": "CRISPR disruption of ZBTB7A binding site in HBG "
                         "promoters, or direct ZBTB7A knockdown in erythroid "
                         "cells. Shows ~28% HbF induction (comparable to BCL11A "
                         "disruption at ~26%)",
            "clinical_status": "Research/preclinical",
        },
    },
}


# ===========================================================================
#  2.  BETA-THALASSEMIA
# ===========================================================================

BETA_THALASSEMIA = {
    "disease": "Beta-Thalassemia",
    "inheritance": "Autosomal recessive",
    "omim": 613985,
    "gene": "HBB",
    "chrom": "chr11",
    "gene_start": 5225464,
    "gene_end": 5229395,

    # ---- Common pathogenic variants (API-verified GRCh38 positions) ----
    "pathogenic_variants": {

        "IVS_I_110_G>A": {
            "name": "IVS-I-110 (G>A) — beta-plus-thalassemia",
            "hgvs": "NM_000518.5:c.93-21G>A",
            "rsid": "rs35004220",
            "pos_grch38": 5226820,
            "ref": "C", "alt": "T",  # plus strand (gene on minus)
            "consequence": "Creates aberrant splice acceptor site in intron 1 "
                           "→ abnormal mRNA with premature stop",
            "clinical_significance": ["pathogenic"],
            "prevalence": "Most common in Mediterranean populations (40-41% in Turkey)",
            "severity": "beta+",
        },

        "IVS_I_5_G>C": {
            "name": "IVS-I-5 (G>C) — beta-plus-thalassemia",
            "hgvs": "NM_000518.5:c.92+5G>C",
            "rsid": "rs33915217",
            "pos_grch38": 5226925,
            "ref": "C", "alt": "A",  # plus strand; multiple alts: A/G/T
            "allele_string": "C/A/G/T",
            "consequence": "Splice donor site disruption in intron 1",
            "clinical_significance": ["pathogenic"],
            "prevalence": "Common in Asian Indian populations (~17.5%)",
            "severity": "beta+",
        },

        "IVS_I_6_T>C": {
            "name": "IVS-I-6 (T>C) — beta-plus-thalassemia",
            "hgvs": "NM_000518.5:c.92+6T>C",
            "rsid": "rs35724775",
            "pos_grch38": 5226924,
            "ref": "A", "alt": "G",  # plus strand
            "allele_string": "A/G/T",
            "consequence": "Splice donor site disruption in intron 1",
            "clinical_significance": ["pathogenic"],
            "severity": "beta+",
        },

        "IVS_II_654_C>T": {
            "name": "IVS-II-654 (C>T) — beta-plus-thalassemia",
            "hgvs": "NM_000518.5:c.316-197C>T",
            "rsid": "rs34690599",
            "pos_grch38": 5225832,
            "ref": "G", "alt": "A",  # plus strand
            "allele_string": "G/A/C",
            "consequence": "Creates cryptic splice site in intron 2",
            "clinical_significance": ["pathogenic"],
            "prevalence": "Common in Chinese/Southeast Asian populations",
            "severity": "beta+",
        },

        "Codon_39_C>T": {
            "name": "Codon 39 (C>T) nonsense — beta-zero-thalassemia",
            "hgvs": "NM_000518.5:c.118C>T",
            "rsid": "rs11549407",
            "pos_grch38": 5226774,
            "ref": "G", "alt": "A",  # plus strand
            "allele_string": "G/A/C/T",
            "consequence": "Premature stop codon (Gln→Stop) in exon 2 → "
                           "nonsense-mediated mRNA decay → no beta-globin",
            "clinical_significance": ["pathogenic"],
            "prevalence": "Common in Mediterranean (Sardinia, Middle East)",
            "severity": "beta0",
        },
    },

    # ---- Therapeutic strategies for beta-thalassemia ----
    "therapeutic_strategies": {
        "Casgevy_BCL11A": {
            "description": "Same BCL11A +58 enhancer disruption as SCD — "
                           "FDA-approved for transfusion-dependent beta-thalassemia "
                           "(TDT) in Jan 2024. Reactivates HbF to compensate for "
                           "deficient/absent adult beta-globin",
            "clinical_status": "FDA-approved (TDT), EMA-approved",
            "trial_id": "NCT03655678",
        },
        "betibeglogene_autotemcel": {
            "description": "Lentiviral gene addition therapy (Zynteglo) — adds "
                           "functional beta-globin (betaA-T87Q) gene via "
                           "lentiviral vector into patient HSCs",
            "therapy_name": "Zynteglo (betibeglogene autotemcel)",
            "developer": "bluebird bio",
            "clinical_status": "FDA-approved Aug 2022 for TDT",
        },
        "base_editing_correction": {
            "description": "Adenine base editing to directly correct splice-site "
                           "mutations (e.g., IVS-I-110 G>A → restore normal "
                           "splicing). Also applicable to nonsense mutations via "
                           "ABE or CBE approaches",
            "clinical_status": "Preclinical/early clinical for specific variants",
        },
        "CRISPR_splice_correction": {
            "description": "CRISPR/Cas9 disruption of aberrant regulatory elements "
                           "created by splice mutations (e.g., destroy the cryptic "
                           "splice site created by IVS-I-110) to restore normal "
                           "splicing patterns",
            "clinical_status": "Preclinical",
        },
    },
}


# ===========================================================================
#  3.  HEMOPHILIA A
# ===========================================================================

HEMOPHILIA_A = {
    "disease": "Hemophilia A",
    "inheritance": "X-linked recessive",
    "omim": 306700,
    "gene": "F8",
    "ensembl_gene": "ENSG00000185010",
    "ensembl_transcript": "ENST00000360256",
    "chrom": "chrX",
    "strand": -1,
    "gene_start": 154835788,
    "gene_end": 155026940,
    "tx_start": 154835792,  # canonical transcript
    "tx_end": 155022723,
    "num_exons": 26,
    "gene_size_kb": 186,  # one of the largest human genes

    # ---- Exon structure (minus strand, genomic coords descend with exon number) ----
    "exons": [
        {"exon": 1,  "start": 155022410, "end": 155022723},
        {"exon": 2,  "start": 154999479, "end": 154999600},
        {"exon": 14, "start": 154928571, "end": 154931676},  # largest exon (3.1 kb)
        {"exon": 22, "start": 154896077, "end": 154896232},
        {"exon": 23, "start": 154863083, "end": 154863227},
        {"exon": 26, "start": 154835792, "end": 154837752},
    ],

    # ---- Common inversions ----
    "structural_variants": {

        "intron_22_inversion": {
            "description": "Inversion of ~600 kb segment due to recombination "
                           "between int22h-1 (within intron 22) and one of two "
                           "extragenic homologous copies (int22h-2 or int22h-3) "
                           "located ~400-500 kb telomeric",
            "intron_22_coords": {
                "start": 154863228,  # exon 23 end + 1
                "end": 154896076,    # exon 22 start - 1
                "size_kb": 32.8,
            },
            "int22h1_region": "Within intron 22 (~9.5 kb repeat element)",
            "frequency": "~45% of severe hemophilia A cases",
            "consequence": "Splits F8 into two fragments with opposite "
                           "transcription directions → no functional FVIII",
            "severity": "Severe (<1% FVIII activity)",
        },

        "intron_1_inversion": {
            "description": "Inversion due to recombination between int1h-1 "
                           "(within intron 1) and int1h-2 (telomeric homolog)",
            "intron_1_coords": {
                "start": 154999601,  # exon 2 end + 1
                "end": 155022409,    # exon 1 start - 1
                "size_kb": 22.8,
            },
            "frequency": "~2-5% of severe hemophilia A cases",
            "severity": "Severe (<1% FVIII activity)",
        },
    },

    # ---- Therapeutic strategies ----
    "therapeutic_strategies": {
        "AAV_gene_addition": {
            "therapy_name": "Roctavian (valoctocogene roxaparvovec)",
            "developer": "BioMarin",
            "mechanism": "AAV5 vector delivering B-domain-deleted F8 cDNA "
                         "driven by liver-selective promoter",
            "clinical_status": "FDA-approved June 2023; EMA conditionally "
                               "approved Aug 2022",
            "trial_id": "NCT03370913 (GENEr8-1)",
            "efficacy": "Reduced treated bleeds from ~5/year to <1/year; "
                        "some decline in FVIII levels over 3-4 years",
        },
        "CRISPR_inversion_correction": {
            "mechanism": "CRISPR/Cas9 dual-cut to re-invert the F8 intron 22 "
                         "or intron 1 inversion in patient iPSCs/HSPCs, then "
                         "differentiate into functional cells",
            "clinical_status": "Preclinical (iPSC proof-of-concept demonstrated)",
            "challenges": "Large inversion size, delivery to hepatocytes in vivo, "
                          "pre-existing Cas9 immunity",
        },
        "antithrombin_knockdown": {
            "mechanism": "In vivo CRISPR/Cas9 (LNP delivery) to knock down "
                         "SERPINC1 (antithrombin) in liver → rebalances "
                         "hemostasis without replacing FVIII",
            "reference": "Antithrombin gene editing via LNP (Science Advances, 2022)",
            "clinical_status": "Preclinical",
        },
    },
}


# ===========================================================================
#  4.  HEMOPHILIA B
# ===========================================================================

HEMOPHILIA_B = {
    "disease": "Hemophilia B",
    "inheritance": "X-linked recessive",
    "omim": 306900,
    "gene": "F9",
    "ensembl_gene": "ENSG00000101981",
    "ensembl_transcript": "ENST00000218099",
    "chrom": "chrX",
    "strand": 1,
    "gene_start": 139530739,
    "gene_end": 139563459,
    "num_exons": 8,

    "exons": [
        {"exon": 1, "start": 139530739, "end": 139530852},
        {"exon": 2, "start": 139537010, "end": 139537173},
        {"exon": 3, "start": 139537362, "end": 139537386},
        {"exon": 4, "start": 139541076, "end": 139541189},
        {"exon": 5, "start": 139548363, "end": 139548491},
        {"exon": 6, "start": 139551062, "end": 139551264},
        {"exon": 7, "start": 139560741, "end": 139560855},
        {"exon": 8, "start": 139561524, "end": 139563459},
    ],

    # ---- Key variants ----
    "key_variants": {

        "Padua_R338L": {
            "name": "FIX-Padua (R338L) — gain-of-function variant",
            "hgvs": "NM_000133.4:c.1150C>T and c.1151G>T",
            "rsid": "rs137852283",
            "pos_grch38": 139561836,  # c.1151G>T position
            "exon": 8,
            "ref": "G", "alt": "T",  # plus strand
            "protein_change": "p.Arg384Leu (preprotein) / p.Arg338Leu (mature FIX)",
            "consequence": "Gain-of-function: 8-fold increased FIX coagulation "
                           "activity; used in gene therapy vector transgenes",
            "clinical_significance": ["pathogenic"],  # causes thrombophilia
            "in_gene_therapy": True,
            "note": "Both FDA-approved hemophilia B gene therapies deliver the "
                    "FIX-Padua variant transgene (R338L) for enhanced activity",
        },
    },

    # ---- Therapeutic strategies ----
    "therapeutic_strategies": {

        "etranacogene_dezaparvovec": {
            "therapy_name": "Hemgenix (etranacogene dezaparvovec)",
            "developer": "CSL Behring / uniQure",
            "mechanism": "AAV5 vector with codon-optimized FIX-Padua (R338L) "
                         "transgene under liver-specific promoter",
            "clinical_status": "FDA-approved Nov 2022",
            "trial_id": "NCT03569891 (HOPE-B)",
        },
        "fidanacogene_elaparvovec": {
            "therapy_name": "Beqvez (fidanacogene elaparvovec)",
            "developer": "Pfizer",
            "mechanism": "AAV vector with high-activity FIX-Padua (R338L) variant",
            "clinical_status": "FDA-approved April 2024",
            "trial_id": "NCT03861273 (BENEGENE-2)",
        },
        "CRISPR_F9_insertion": {
            "therapy_name": "REGV131-LNP1265",
            "mechanism": "In vivo CRISPR/Cas9 via LNP to insert functional F9 "
                         "sequence at a safe-harbor locus in hepatocytes",
            "clinical_status": "First-in-human Phase 1 (enrolled 2025)",
        },
        "CRISPR_Padua_knockin": {
            "mechanism": "CRISPR-mediated in situ introduction of the F9-Padua "
                         "R338L variant via HDR in patient iPSCs for autologous "
                         "cell therapy",
            "clinical_status": "Preclinical (proof-of-concept in iPSCs)",
        },
    },
}


# ===========================================================================
#  5.  FANCONI ANEMIA  (Complementation Group A)
# ===========================================================================

FANCONI_ANEMIA = {
    "disease": "Fanconi Anemia, Complementation Group A",
    "inheritance": "Autosomal recessive",
    "omim": 227650,
    "gene": "FANCA",
    "ensembl_gene": "ENSG00000187741",
    "ensembl_transcript": "ENST00000389301",
    "chrom": "chr16",
    "strand": -1,
    "gene_start": 89726683,
    "gene_end": 89816977,
    "gene_size_kb": 90,
    "num_pathogenic_variants_in_clinvar": 475,  # approximate

    "key_details": {
        "frequency": "FANCA mutations account for 60-70% of all Fanconi anemia cases",
        "mutation_spectrum": "Large deletions, splice-site mutations, nonsense, "
                            "and missense mutations throughout the gene",
        "pathway": "FANCA is part of the FA core complex required for DNA "
                   "interstrand crosslink (ICL) repair; deficiency causes "
                   "genomic instability and bone marrow failure",
    },

    "therapeutic_strategies": {
        "lentiviral_gene_addition": {
            "description": "Autologous HSPC transduction with lentiviral vector "
                           "carrying functional FANCA cDNA, followed by "
                           "re-infusion (reduced-intensity conditioning due to "
                           "DNA repair sensitivity)",
            "developer": "Rocket Pharmaceuticals",
            "clinical_status": "Phase 1/2 ongoing (NCT04248439); 9+ patients "
                               "treated; 49% survival of BM colony-forming units "
                               "at 12-36 months after ICL challenge",
            "trial_id": "NCT04248439",
        },
        "base_editing": {
            "description": "ABE8e adenine base editing to correct specific FANCA "
                           "nonsense mutations (e.g., converting stop codons to "
                           "missense) in patient HSPCs — proof of concept "
                           "demonstrated for two prevalent FANCA mutations",
            "clinical_status": "Preclinical (Nature Communications, 2022)",
            "advantage": "NHEJ-based correction is enhanced in FA cells due to "
                         "impaired HDR, making base editing particularly suitable",
        },
        "CRISPR_correction": {
            "description": "CRISPR/Cas9 NHEJ-based correction of FA mutations "
                           "in patient cells — leverages the fact that FA cells "
                           "have impaired HDR but functional NHEJ",
            "clinical_status": "Preclinical",
        },
    },
}


# ===========================================================================
#  6.  DIAMOND-BLACKFAN ANEMIA  (RPS19)
# ===========================================================================

DIAMOND_BLACKFAN_ANEMIA = {
    "disease": "Diamond-Blackfan Anemia",
    "inheritance": "Autosomal dominant",
    "omim": 105650,
    "gene": "RPS19",
    "ensembl_gene": "ENSG00000105372",
    "ensembl_transcript": "ENST00000598742",
    "chrom": "chr19",
    "strand": 1,
    "gene_start": 41860255,
    "gene_end": 41872925,

    "key_details": {
        "frequency": "RPS19 mutations account for ~25% of DBA cases",
        "other_genes": "RPL5, RPL11, RPS26, RPS24, RPL35A, and others "
                       "(all ribosomal protein genes)",
        "pathway": "Ribosomal protein haploinsufficiency → impaired ribosome "
                   "biogenesis → p53 activation → erythroid progenitor apoptosis",
    },

    "therapeutic_strategies": {
        "lentiviral_gene_addition": {
            "description": "LV carrying functional RPS19 under regulated promoter "
                           "(SJEFS-S19 vector) — restores ribosome assembly and "
                           "erythropoiesis in patient HSPCs",
            "developer": "St. Jude Children's Research Hospital",
            "clinical_status": "Preclinical (2024-2025 publications showing "
                               "efficacy in humanized mouse models)",
        },
        "GATA1_therapy": {
            "description": "Regulated GATA1 expression as a universal gene therapy "
                           "for DBA — bypasses the need to target specific RP genes; "
                           "up to 21-fold stimulation of erythroid production",
            "reference": "Cell Stem Cell, Nov 2024",
            "clinical_status": "Preclinical",
        },
        "CRISPR_correction": {
            "description": "CRISPR/Cas9 correction of specific RPS19 mutations "
                           "in iPSC-derived HSPCs via HDR or base editing",
            "clinical_status": "Preclinical (disease models established 2022-2023)",
        },
    },
}


# ===========================================================================
#  7.  PYRUVATE KINASE DEFICIENCY
# ===========================================================================

PYRUVATE_KINASE_DEFICIENCY = {
    "disease": "Pyruvate Kinase Deficiency",
    "inheritance": "Autosomal recessive",
    "omim": 266200,
    "gene": "PKLR",
    "ensembl_gene": "ENSG00000143627",
    "ensembl_transcript": "ENST00000342741",
    "chrom": "chr1",
    "strand": -1,
    "gene_start": 155289293,
    "gene_end": 155301438,

    "key_details": {
        "pathway": "PKLR encodes the liver/red cell isoform of pyruvate kinase; "
                   "deficiency impairs glycolysis in RBCs → chronic hemolytic "
                   "anemia, reticulocytosis, splenomegaly",
        "mutation_spectrum": ">300 pathogenic variants described (missense, "
                            "splice, nonsense, deletions)",
    },

    "therapeutic_strategies": {
        "small_molecule_activator": {
            "therapy_name": "Mitapivat (Pyrukynd)",
            "developer": "Agios Pharmaceuticals",
            "mechanism": "Allosteric activator of pyruvate kinase — stabilizes "
                         "the R-state of the enzyme, enhancing residual activity",
            "clinical_status": "FDA-approved Feb 2022 (adults); ACTIVATE-Kids "
                               "Phase 3 met primary endpoint Feb 2025 (ages 1-17)",
        },
        "lentiviral_gene_addition": {
            "description": "Autologous HSPC transduction with LV carrying "
                           "corrected PKLR cDNA, myeloablative conditioning, "
                           "re-infusion",
            "clinical_status": "Phase 1 completed (2 adults, 24-month follow-up "
                               "showing normal Hb and transfusion-free); "
                               "Phase 2 planned",
        },
        "CRISPR_correction": {
            "description": "Variant-specific CRISPR/Cas9 or base editing "
                           "correction in patient HSPCs",
            "clinical_status": "Preclinical",
        },
    },
}


# ===========================================================================
#  SUMMARY TABLE:  All genes for pipeline BED-file generation
# ===========================================================================

HEMATOLOGICAL_GENE_COORDS = {
    # gene: (chrom, start, end, strand, ensembl_id)
    "HBB":    ("chr11",   5225464,   5229395, -1, "ENSG00000244734"),
    "HBG1":   ("chr11",   5248269,   5249857, -1, "ENSG00000213934"),
    "HBG2":   ("chr11",   5253188,   5254781, -1, "ENSG00000196565"),
    "HBD":    ("chr11",   5232838,   5234483, -1, "ENSG00000223609"),
    "HBE1":   ("chr11",   5268345,   5269945, -1, "ENSG00000213931"),
    "BCL11A": ("chr2",  60450520, 60554467, -1, "ENSG00000119866"),
    "KLF1":   ("chr19", 12884422, 12887201, -1, "ENSG00000105610"),
    "ZBTB7A": ("chr19",  4043303,  4067636, -1, "ENSG00000178951"),
    "F8":     ("chrX", 154835788, 155026940, -1, "ENSG00000185010"),
    "F9":     ("chrX", 139530739, 139563459,  1, "ENSG00000101981"),
    "FANCA":  ("chr16", 89726683, 89816977, -1, "ENSG00000187741"),
    "RPS19":  ("chr19", 41860255, 41872925,  1, "ENSG00000105372"),
    "PKLR":   ("chr1", 155289293, 155301438, -1, "ENSG00000143627"),
}


# ===========================================================================
#  KEY VARIANT POSITIONS  (for pipeline VCF filtering / annotation)
# ===========================================================================

HEMATOLOGICAL_KEY_VARIANTS = [
    # (chrom, pos, rsid, ref, alt, gene, variant_name, disease)
    ("chr11",   5227002, "rs334",       "T", "A",  "HBB",  "E6V (HbS)",             "Sickle Cell Disease"),
    ("chr11",   5226820, "rs35004220",  "C", "T",  "HBB",  "IVS-I-110 G>A",         "Beta-Thalassemia"),
    ("chr11",   5226925, "rs33915217",  "C", "A",  "HBB",  "IVS-I-5 G>C",           "Beta-Thalassemia"),
    ("chr11",   5226924, "rs35724775",  "A", "G",  "HBB",  "IVS-I-6 T>C",           "Beta-Thalassemia"),
    ("chr11",   5225832, "rs34690599",  "G", "A",  "HBB",  "IVS-II-654 C>T",        "Beta-Thalassemia"),
    ("chr11",   5226774, "rs11549407",  "G", "A",  "HBB",  "Codon 39 (C>T) nonsense", "Beta-Thalassemia"),
    ("chr11",   5250055, "rs35710727",  "A", "G",  "HBG1", "-198 T>C (HPFH)",       "HPFH / Therapeutic"),
    ("chr11",   5250053, "rs35983258",  "G", "A",  "HBG1", "-196 C>T (HPFH)",       "HPFH / Therapeutic"),
    ("chr11",   5250052, "rs35321913",  "G", "C",  "HBG1", "-195 C>G (HPFH)",       "HPFH / Therapeutic"),
    ("chr11",   5250032, "rs34153902",  "A", "G",  "HBG1", "-175 T>C (HPFH)",       "HPFH / Therapeutic"),
    ("chrX", 139561836, "rs137852283",  "G", "T",  "F9",   "Padua R338L (gain-of-fn)", "Hemophilia B / Gene Therapy"),
]


# ===========================================================================
#  CRISPR THERAPEUTIC REGIONS  (for BED-file target definition)
# ===========================================================================

CRISPR_THERAPEUTIC_REGIONS = [
    # (chrom, start, end, name, strategy, clinical_status)
    ("chr2",  60494251, 60495546,
     "BCL11A_+58_enhancer", "Cas9 knockout (Casgevy)", "FDA-approved 2023"),
    ("chr11",  5250020,  5250070,
     "HBG1_promoter_BCL11A_binding", "Adenine base editing (BEAM-101)", "Phase 1/2"),
    ("chr11",  5255020,  5255070,
     "HBG2_promoter_BCL11A_binding", "Adenine base editing (BEAM-101)", "Phase 1/2"),
    ("chr11",  5226990,  5227010,
     "HBB_E6V_region", "Base editing to HbG-Makassar (BEAM-102)", "Preclinical"),
    ("chrX", 154863228, 154896076,
     "F8_intron22", "CRISPR inversion correction", "Preclinical"),
    ("chrX", 139561524, 139563459,
     "F9_exon8_Padua_region", "HDR knock-in of R338L", "Preclinical"),
]


# ===========================================================================
#  APPROVED / LATE-STAGE GENE THERAPIES  (quick-reference)
# ===========================================================================

APPROVED_GENE_THERAPIES = [
    {
        "therapy": "Casgevy (exagamglogene autotemcel)",
        "developer": "Vertex / CRISPR Therapeutics",
        "modality": "Ex vivo CRISPR-Cas9 (BCL11A enhancer knockout)",
        "diseases": ["Sickle Cell Disease", "Transfusion-Dependent Beta-Thalassemia"],
        "approval": "FDA Dec 2023 (SCD), Jan 2024 (TDT); EMA Feb 2024",
        "target_gene": "BCL11A (+58 enhancer)",
    },
    {
        "therapy": "Zynteglo (betibeglogene autotemcel)",
        "developer": "bluebird bio",
        "modality": "Ex vivo lentiviral gene addition (beta-A-T87Q-globin)",
        "diseases": ["Transfusion-Dependent Beta-Thalassemia"],
        "approval": "FDA Aug 2022",
        "target_gene": "HBB transgene",
    },
    {
        "therapy": "Lyfgenia (lovotibeglogene autotemcel)",
        "developer": "bluebird bio",
        "modality": "Ex vivo lentiviral gene addition (anti-sickling beta-globin)",
        "diseases": ["Sickle Cell Disease"],
        "approval": "FDA Dec 2023",
        "target_gene": "HBB transgene",
    },
    {
        "therapy": "Roctavian (valoctocogene roxaparvovec)",
        "developer": "BioMarin",
        "modality": "In vivo AAV5 gene addition (B-domain-deleted F8)",
        "diseases": ["Severe Hemophilia A"],
        "approval": "FDA June 2023; EMA conditional Aug 2022",
        "target_gene": "F8 transgene",
    },
    {
        "therapy": "Hemgenix (etranacogene dezaparvovec)",
        "developer": "CSL Behring / uniQure",
        "modality": "In vivo AAV5 gene addition (FIX-Padua R338L)",
        "diseases": ["Hemophilia B"],
        "approval": "FDA Nov 2022",
        "target_gene": "F9-Padua transgene",
    },
    {
        "therapy": "Beqvez (fidanacogene elaparvovec)",
        "developer": "Pfizer",
        "modality": "In vivo AAV gene addition (FIX-Padua R338L)",
        "diseases": ["Hemophilia B"],
        "approval": "FDA April 2024",
        "target_gene": "F9-Padua transgene",
    },
]

LATE_STAGE_GENE_THERAPIES = [
    {
        "therapy": "risto-cel (ristoglogene autogetemcel, BEAM-101)",
        "developer": "Beam Therapeutics",
        "modality": "Ex vivo adenine base editing (HBG1/2 promoter)",
        "diseases": ["Sickle Cell Disease"],
        "status": "Phase 1/2 (BEACON, NCT05456880); NEJM pub April 2026",
        "target_gene": "HBG1/HBG2 promoters",
    },
    {
        "therapy": "FANCA lentiviral gene therapy",
        "developer": "Rocket Pharmaceuticals",
        "modality": "Ex vivo lentiviral gene addition (FANCA)",
        "diseases": ["Fanconi Anemia Group A"],
        "status": "Phase 1/2 (NCT04248439)",
        "target_gene": "FANCA transgene",
    },
    {
        "therapy": "PKLR lentiviral gene therapy",
        "developer": "Various academic sponsors",
        "modality": "Ex vivo lentiviral gene addition (PKLR)",
        "diseases": ["Pyruvate Kinase Deficiency"],
        "status": "Phase 1 completed; Phase 2 planned",
        "target_gene": "PKLR transgene",
    },
]


# ===========================================================================
# UNIFIED PIPELINE-COMPATIBLE TARGET DICT
# Matches the {gene: {chrom, start, end, role, strategy, ...}} format
# used by the rest of the OpenCure pipeline.
# ===========================================================================

ALL_HEMATOLOGICAL_TARGETS = {
    "HBB": {
        "chrom": "chr11", "start": 5225464, "end": 5229395,
        "role": "Beta-globin — adult hemoglobin HbA component",
        "strategy": "Base-editing to correct E6V (rs334) for SCD; "
                     "HDR correction for beta-thalassemia splice/nonsense variants",
        "conditions": ["sickle_cell_disease", "beta_thalassemia"],
    },
    "BCL11A": {
        "chrom": "chr2", "start": 60450520, "end": 60554467,
        "role": "BCL11A — fetal hemoglobin repressor; erythroid enhancer "
                "is the Casgevy/exa-cel target",
        "strategy": "CRISPR disruption of +58 erythroid enhancer to "
                     "derepress HbF (Casgevy mechanism, FDA-approved)",
        "conditions": ["sickle_cell_disease", "beta_thalassemia"],
    },
    "HBG1": {
        "chrom": "chr11", "start": 5248269, "end": 5249857,
        "role": "Gamma-globin 1 (A-gamma) — fetal hemoglobin component",
        "strategy": "Base-editing of promoter HPFH sites to reactivate "
                     "HbF expression (risto-cel/BEAM-101 approach)",
        "conditions": ["sickle_cell_disease", "beta_thalassemia"],
    },
    "HBG2": {
        "chrom": "chr11", "start": 5253188, "end": 5254781,
        "role": "Gamma-globin 2 (G-gamma) — fetal hemoglobin component",
        "strategy": "Base-editing of promoter HPFH sites to reactivate HbF",
        "conditions": ["sickle_cell_disease", "beta_thalassemia"],
    },
    "KLF1": {
        "chrom": "chr19", "start": 12884422, "end": 12887201,
        "role": "Erythroid Kruppel-like factor — activates adult globin, "
                "represses fetal globin",
        "strategy": "CRISPRi silencing to shift globin balance toward HbF",
        "conditions": ["sickle_cell_disease", "beta_thalassemia"],
    },
    "ZBTB7A": {
        "chrom": "chr19", "start": 4043303, "end": 4067636,
        "role": "LRF/ZBTB7A — fetal hemoglobin repressor at HBG promoter",
        "strategy": "CRISPR disruption of ZBTB7A binding site in HBG "
                     "promoters to derepress HbF",
        "conditions": ["sickle_cell_disease", "beta_thalassemia"],
    },
    "F8": {
        "chrom": "chrX", "start": 154835788, "end": 155026940,
        "role": "Coagulation factor VIII — intrinsic coagulation cascade",
        "strategy": "CRISPR inversion correction for intron 22/1 inversions; "
                     "AAV gene addition (Roctavian, FDA-approved)",
        "conditions": ["hemophilia_A"],
    },
    "F9": {
        "chrom": "chrX", "start": 139530739, "end": 139563459,
        "role": "Coagulation factor IX — intrinsic coagulation cascade",
        "strategy": "AAV gene addition with FIX-Padua (R338L) transgene "
                     "(Hemgenix/Beqvez, FDA-approved); CRISPR insertion "
                     "at albumin safe harbor",
        "conditions": ["hemophilia_B"],
    },
    "FANCA": {
        "chrom": "chr16", "start": 89726683, "end": 89816977,
        "role": "Fanconi anemia complementation group A — DNA repair",
        "strategy": "Ex vivo lentiviral gene addition in HSCs; "
                     "ABE base-editing to correct pathogenic variants",
        "conditions": ["fanconi_anemia"],
    },
    "RPS19": {
        "chrom": "chr19", "start": 41860255, "end": 41872925,
        "role": "Ribosomal protein S19 — ribosome biogenesis",
        "strategy": "Lentiviral gene addition; GATA1-based universal "
                     "gene therapy approach",
        "conditions": ["diamond_blackfan_anemia"],
    },
    "PKLR": {
        "chrom": "chr1", "start": 155289293, "end": 155301438,
        "role": "Pyruvate kinase (liver and red cell) — glycolysis",
        "strategy": "Lentiviral gene addition in HSCs; mitapivat "
                     "(small molecule activator) approved",
        "conditions": ["pyruvate_kinase_deficiency"],
    },
}
