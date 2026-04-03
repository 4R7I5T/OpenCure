"""
Lipid disorder gene-therapy / CRISPR targets -- GRCh38 (hg38) coordinates.

Covers familial hypercholesterolemia, Lp(a) excess, sitosterolemia,
and related monogenic dyslipidemias with active gene therapy or gene
editing programs.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - Verve Therapeutics (PCSK9/ANGPTL3 base editing programs)
  - ClinicalTrials.gov
"""


LIPID_TARGETS = {

    # -------------------------------------------------------------------
    # 1.  LDLR -- LDL Receptor
    # -------------------------------------------------------------------
    "LDLR": {
        "gene_id": 3949,
        "chrom": "chr19",
        "start": 11_089_362,
        "end": 11_133_820,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19p13.2",
        "exon_count": 18,
        "role": (
            "LDL receptor -- hepatocyte surface receptor for LDL-cholesterol "
            "clearance.  Mutations cause familial hypercholesterolemia (FH).  "
            "Heterozygous FH (1:250) -> LDL-C ~200-350 mg/dL, premature ASCVD.  "
            "Homozygous FH (1:300,000) -> LDL-C >500 mg/dL, childhood MI."
        ),
        "disease": "Familial hypercholesterolemia (FH)",
        "omim_disease": 143890,
        "omim_gene": 606945,
        "inheritance": "AD (semi-dominant)",
        "key_variants": [
            {
                "name": ">2000 pathogenic variants described",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Extreme allelic heterogeneity.  Classes: null (class 1), "
                    "transport-deficient (class 2, e.g., FH-Lebanese), "
                    "binding-deficient (class 3), internalization-deficient "
                    "(class 4), recycling-deficient (class 5)."
                ),
            },
            {
                "name": "c.681C>G (FH-Quebec, class 2)",
                "rsid": "rs28942078",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "French Canadian founder; ER retention.",
            },
        ],
        "strategy": (
            "(1) AAV liver gene replacement for homozygous FH (cDNA ~2.6 kb "
            "fits AAV) -- most direct approach for receptor-negative patients; "
            "(2) CRISPRa to upregulate wild-type LDLR allele in heterozygous FH; "
            "(3) In vivo base editing of PCSK9 as alternative (reduces PCSK9 -> "
            "upregulates LDLR indirectly -- Verve approach); "
            "(4) Base editing for specific LDLR mutations.  "
            "Delivery: AAV8 IV to liver or LNP."
        ),
        "clinical_programs": (
            "Regeneron AAV-LDLR (Phase 1/2 for HoFH, NCT02651675).  "
            "Verve VERVE-101 (PCSK9 base editing) as indirect LDLR upregulation.  "
            "Standard of care: statins, ezetimibe, PCSK9 inhibitors "
            "(evolocumab, alirocumab), lomitapide (HoFH)."
        ),
        "conditions": ["familial_hypercholesterolemia", "FH", "HoFH",
                        "lipid_disorder", "hypercholesterolemia", "ASCVD"],
    },

    # -------------------------------------------------------------------
    # 2.  PCSK9 -- Proprotein Convertase Subtilisin/Kexin Type 9
    # -------------------------------------------------------------------
    "PCSK9": {
        "gene_id": 255738,
        "chrom": "chr1",
        "start": 55_039_475,
        "end": 55_064_852,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p32.3",
        "exon_count": 12,
        "role": (
            "PCSK9 -- serine protease that promotes LDLR degradation.  "
            "Gain-of-function mutations cause autosomal dominant FH (FH3).  "
            "Loss-of-function is protective: natural knockout individuals "
            "(e.g., PCSK9-Q152H homozygotes) have LDL-C ~15 mg/dL with no "
            "adverse health effects, validating PCSK9 as a safe target."
        ),
        "disease": "FH3 (gain-of-function) / cardioprotective target",
        "omim_disease": 603776,
        "omim_gene": 607786,
        "inheritance": "AD (GOF mutations)",
        "key_variants": [
            {
                "name": "p.Asp374Tyr (D374Y, GOF)",
                "rsid": "rs137852912",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most severe GOF mutation; LDL-C >300 mg/dL.",
            },
            {
                "name": "p.Tyr142Ter (Y142X, LOF, protective)",
                "rsid": "rs67608943",
                "consequence": "nonsense",
                "clinical_significance": "protective",
                "notes": "African American; 88% reduction in CHD risk.",
            },
            {
                "name": "p.Cys679Ter (C679X, LOF, protective)",
                "rsid": "rs28362286",
                "consequence": "nonsense",
                "clinical_significance": "protective",
                "notes": "Complete LOF; LDL-C ~15 mg/dL when homozygous.",
            },
        ],
        "strategy": (
            "(1) In vivo base editing to introduce LOF mutation in hepatocytes "
            "-- ONE-TIME treatment to permanently lower LDL-C (Verve VERVE-101 "
            "approach); "
            "(2) CRISPRi to silence PCSK9 expression; "
            "(3) CRISPR knockout of PCSK9 in liver (more aggressive but proven "
            "safe by natural LOF carriers); "
            "(4) Allele-specific disruption of GOF allele for FH3 patients.  "
            "Delivery: LNP IV to liver (Verve approach) or AAV8."
        ),
        "clinical_programs": (
            "Verve Therapeutics VERVE-101 -- adenine base editing to disrupt "
            "PCSK9 in liver, Phase 1b (heart-1 trial, NCT05398029).  "
            "FIRST IN VIVO BASE EDITING TRIAL IN HUMANS.  "
            "Initial results showed ~55% PCSK9 reduction, ~50% LDL-C reduction.  "
            "Verve VERVE-102 (improved LNP) Phase 1.  "
            "PCSK9 antibodies (evolocumab, alirocumab) FDA-approved."
        ),
        "conditions": ["familial_hypercholesterolemia", "FH", "lipid_disorder",
                        "PCSK9_GOF", "ASCVD", "hypercholesterolemia",
                        "cardiovascular_prevention"],
    },

    # -------------------------------------------------------------------
    # 3.  APOB -- Apolipoprotein B
    # -------------------------------------------------------------------
    "APOB": {
        "gene_id": 338,
        "chrom": "chr2",
        "start": 21_001_429,
        "end": 21_044_073,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2p24.1",
        "exon_count": 29,
        "role": (
            "Apolipoprotein B-100 -- structural protein of LDL particles, "
            "ligand for LDLR binding.  Gain-of-function mutations in the "
            "LDLR-binding domain cause familial defective ApoB (FDB, FH-like).  "
            "APOB is the sole structural protein of LDL and VLDL."
        ),
        "disease": "Familial defective ApoB-100 (FDB)",
        "omim_disease": 144010,
        "omim_gene": 107730,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Arg3527Gln (R3527Q)",
                "rsid": "rs5742904",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common FDB mutation; reduces LDLR binding ~95%.  "
                    "Prevalence ~1:500 in Central Europeans."
                ),
            },
        ],
        "strategy": (
            "(1) Base editing to correct R3527Q (single nucleotide change); "
            "(2) Verve VERVE-201 approach: base editing of ANGPTL3 as "
            "alternative LDL-lowering strategy; "
            "(3) ASO-mediated APOB reduction (mipomersen, approved).  "
            "Delivery: LNP to liver."
        ),
        "clinical_programs": (
            "Mipomersen (anti-APOB ASO) FDA-approved for HoFH.  "
            "No CRISPR trials targeting APOB directly as of 2026.  "
            "Verve exploring APOB base editing preclinically."
        ),
        "conditions": ["familial_defective_apoB", "FDB", "lipid_disorder",
                        "hypercholesterolemia", "FH"],
    },

    # -------------------------------------------------------------------
    # 4.  ANGPTL3 -- Angiopoietin-Like 3
    # -------------------------------------------------------------------
    "ANGPTL3": {
        "gene_id": 27329,
        "chrom": "chr1",
        "start": 62_597_520,
        "end": 62_606_313,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p31.3",
        "exon_count": 7,
        "role": (
            "ANGPTL3 -- inhibitor of lipoprotein lipase (LPL) and endothelial "
            "lipase.  Natural loss-of-function causes familial combined "
            "hypolipidemia (FHBL2): low LDL-C, low HDL-C, low triglycerides, "
            "with no adverse effects and ~40% reduced ASCVD risk.  "
            "Ideal CRISPR knockout target."
        ),
        "disease": "Therapeutic target for pan-lipid lowering",
        "omim_disease": 605019,
        "omim_gene": 604774,
        "inheritance": "AR (LOF = protective hypolipidemia)",
        "key_variants": [
            {
                "name": "p.Ser17Ter (S17X, LOF)",
                "rsid": "rs121917817",
                "consequence": "nonsense",
                "clinical_significance": "protective",
                "notes": "Complete LOF in Campodimele, Italy; very low lipids.",
            },
        ],
        "strategy": (
            "(1) In vivo base editing to create LOF in hepatocytes -- Verve "
            "VERVE-201 approach; "
            "(2) CRISPR knockout of ANGPTL3 in liver; "
            "(3) CRISPRi to silence ANGPTL3.  "
            "Delivery: LNP IV to liver."
        ),
        "clinical_programs": (
            "Verve Therapeutics VERVE-201 -- base editing of ANGPTL3, "
            "Phase 1 planned (IND filed 2025).  "
            "Evinacumab (anti-ANGPTL3 antibody) FDA-approved for HoFH.  "
            "Vupanorsen (anti-ANGPTL3 ASO) discontinued."
        ),
        "conditions": ["lipid_disorder", "pan_lipid_lowering", "HoFH",
                        "triglycerides", "ASCVD", "cardiovascular_prevention"],
    },

    # -------------------------------------------------------------------
    # 5.  LPA -- Lipoprotein(a)
    # -------------------------------------------------------------------
    "LPA": {
        "gene_id": 4018,
        "chrom": "chr6",
        "start": 160_531_482,
        "end": 160_664_275,
        "strand": "-",
        "refseq": "NC_000006.12",
        "cytoband": "6q25.3",
        "exon_count": 40,
        "role": (
            "Lipoprotein(a) -- LDL-like particle with apolipoprotein(a) "
            "covalently bound.  Elevated Lp(a) (>50 mg/dL, ~20% of population) "
            "is a causal, independent risk factor for ASCVD, aortic stenosis, "
            "and stroke.  Genetically determined -- no lifestyle modification "
            "lowers Lp(a) substantially."
        ),
        "disease": "Elevated Lp(a) / ASCVD risk",
        "omim_gene": 152200,
        "inheritance": "Quantitative trait (KIV-2 repeat number)",
        "key_variants": [
            {
                "name": "KIV-2 copy number variation",
                "rsid": None,
                "consequence": "copy_number_variation",
                "clinical_significance": "risk_determinant",
                "notes": (
                    "Inverse correlation: fewer KIV-2 repeats -> smaller "
                    "apo(a) isoform -> higher Lp(a) levels -> higher ASCVD risk.  "
                    "Genetically determined, ~90% heritability."
                ),
            },
            {
                "name": "rs10455872 (intron)",
                "rsid": "rs10455872",
                "consequence": "regulatory",
                "clinical_significance": "risk_factor",
                "notes": "Strong predictor of elevated Lp(a); aortic stenosis risk.",
            },
        ],
        "strategy": (
            "(1) CRISPRi or CRISPR knockout of LPA in liver to eliminate "
            "Lp(a) production (no known physiological requirement); "
            "(2) ASO-mediated LPA silencing (Novartis pelacarsen approach); "
            "(3) siRNA-mediated LPA reduction (Silence SLN360, Amgen olpasiran).  "
            "Delivery: LNP to liver or ASO subcutaneous."
        ),
        "clinical_programs": (
            "Novartis pelacarsen (anti-LPA ASO) Phase 3 Lp(a)HORIZON "
            "(NCT04023552) -- ASCVD outcomes trial.  "
            "Amgen olpasiran (siRNA) Phase 2 OCEAN(a) showed >95% Lp(a) "
            "reduction.  Phase 3 ongoing.  "
            "Silence Therapeutics SLN360 Phase 2.  "
            "No CRISPR trials for LPA as of 2026."
        ),
        "conditions": ["elevated_Lpa", "lipid_disorder", "ASCVD",
                        "aortic_stenosis", "cardiovascular_prevention"],
    },

    # -------------------------------------------------------------------
    # 6.  ABCG5 / ABCG8 -- Sitosterolemia
    # -------------------------------------------------------------------
    "ABCG5": {
        "gene_id": 64240,
        "chrom": "chr2",
        "start": 43_816_293,
        "end": 43_842_302,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2p21",
        "exon_count": 13,
        "role": (
            "Sterolin-1 (ABCG5) -- half-transporter forming heterodimer with "
            "ABCG8 to pump plant sterols and cholesterol out of enterocytes "
            "and hepatocytes.  Deficiency causes sitosterolemia: plant sterol "
            "accumulation -> xanthomas, premature ASCVD."
        ),
        "disease": "Sitosterolemia",
        "omim_disease": 210250,
        "omim_gene": 605459,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various truncating mutations",
                "rsid": None,
                "consequence": "nonsense/frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Rare -- estimated 1:200,000.",
            },
        ],
        "strategy": (
            "(1) AAV liver gene replacement (cDNA ~2.0 kb fits AAV); "
            "(2) Base editing for specific mutations.  "
            "Delivery: AAV8 liver-directed."
        ),
        "clinical_programs": (
            "Ezetimibe standard of care (blocks NPC1L1 sterol absorption).  "
            "No gene therapy trials."
        ),
        "conditions": ["sitosterolemia", "lipid_disorder", "plant_sterol",
                        "xanthomas"],
    },

    "ABCG8": {
        "gene_id": 64241,
        "chrom": "chr2",
        "start": 43_838_864,
        "end": 43_878_466,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2p21",
        "exon_count": 13,
        "role": (
            "Sterolin-2 (ABCG8) -- ABCG5 heterodimer partner for sterol "
            "efflux.  Mutations also cause sitosterolemia."
        ),
        "disease": "Sitosterolemia",
        "omim_disease": 210250,
        "omim_gene": 605460,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gly574Arg (G574R)",
                "rsid": "rs137852988",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts heterodimer function.",
            },
        ],
        "strategy": (
            "(1) AAV liver gene replacement (cDNA ~2.0 kb); "
            "(2) Dual ABCG5/ABCG8 replacement if both deficient."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["sitosterolemia", "lipid_disorder", "plant_sterol"],
    },
}
