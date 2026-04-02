"""
Mental-health pharmacogenomic & gene-therapy targets — GRCh38 coordinates.

These are legitimate therapeutic targets for consented patients.
Use cases:
  - Pharmacogenomics: predict drug response from patient genotype
  - Gene therapy: correct variants causing treatment-resistant conditions
  - Therapeutic modulation: up/down-regulate to restore healthy function

All interventions require informed patient consent and IRB/ethics approval.
"""

MENTAL_HEALTH_TARGETS = {
    # --- Depression / Mood disorders ---
    "SLC6A4": {
        "chrom": "chr17",
        "start": 30194319,
        "end": 30236002,
        "role": "Serotonin transporter (5-HTT) — SSRI drug target",
        "strategy": "Pharmacogenomic: 5-HTTLPR genotyping predicts SSRI "
                     "response. Gene therapy: CRISPRi to reduce SERT in "
                     "treatment-resistant depression (analogous to SSRI mechanism)",
        "conditions": ["major_depressive_disorder", "anxiety", "OCD"],
        "key_variants": ["5-HTTLPR_short", "rs25531"],
        "therapeutic_class": "antidepressant_response",
    },
    "BDNF": {
        "chrom": "chr11",
        "start": 27654893,
        "end": 27722058,
        "role": "Brain-derived neurotrophic factor — synaptic plasticity",
        "strategy": "CRISPRa to enhance BDNF expression in treatment-resistant "
                     "depression. Val66Met (rs6265) genotyping predicts "
                     "antidepressant response and cognitive outcomes",
        "conditions": ["treatment_resistant_depression", "cognitive_decline"],
        "key_variants": ["rs6265_Val66Met"],
        "therapeutic_class": "neuroplasticity",
    },
    "HTR2A": {
        "chrom": "chr13",
        "start": 46831546,
        "end": 46897076,
        "role": "Serotonin receptor 2A — target of atypical antipsychotics",
        "strategy": "Pharmacogenomic: rs6313/rs6311 predict antipsychotic "
                     "response and side effects. Gene therapy: modulate "
                     "expression for treatment-resistant schizophrenia",
        "conditions": ["schizophrenia", "depression", "psychosis"],
        "key_variants": ["rs6313", "rs6311"],
        "therapeutic_class": "antipsychotic_response",
    },
    "MTHFR": {
        "chrom": "chr1",
        "start": 11785723,
        "end": 11806103,
        "role": "Folate metabolism — methylation cycle critical for "
                "neurotransmitter synthesis",
        "strategy": "Pharmacogenomic: C677T and A1298C predict folate-related "
                     "depression and guide L-methylfolate supplementation. "
                     "Base-editing to correct C677T in severe cases",
        "conditions": ["depression", "folate_deficiency", "hyperhomocysteinemia"],
        "key_variants": ["rs1801133_C677T", "rs1801131_A1298C"],
        "therapeutic_class": "metabolic_psychiatry",
    },

    # --- Anxiety / Stress disorders ---
    "FKBP5": {
        "chrom": "chr6",
        "start": 35541362,
        "end": 35697008,
        "role": "FK506-binding protein 5 — HPA axis / glucocorticoid "
                "receptor sensitivity",
        "strategy": "Pharmacogenomic: rs1360780 predicts PTSD risk and "
                     "cortisol response. Epigenetic editing (dCas9-TET1) "
                     "to normalize methylation in trauma-exposed patients",
        "conditions": ["PTSD", "anxiety", "stress_disorders"],
        "key_variants": ["rs1360780", "rs9296158"],
        "therapeutic_class": "stress_response",
    },
    "NR3C1": {
        "chrom": "chr5",
        "start": 143277931,
        "end": 143435512,
        "role": "Glucocorticoid receptor — stress hormone signaling",
        "strategy": "Epigenetic editing to restore NR3C1 promoter methylation "
                     "patterns disrupted by early-life stress. CRISPRa to "
                     "enhance expression in cortisol-resistant depression",
        "conditions": ["depression", "PTSD", "anxiety"],
        "key_variants": ["rs6195", "exon_1F_methylation"],
        "therapeutic_class": "stress_response",
    },
    "GABRA2": {
        "chrom": "chr4",
        "start": 46232216,
        "end": 46394194,
        "role": "GABA-A receptor alpha-2 subunit — inhibitory neurotransmission",
        "strategy": "Pharmacogenomic: variants predict benzodiazepine response "
                     "and alcohol dependence risk. CRISPRa to enhance GABAergic "
                     "tone in treatment-resistant anxiety",
        "conditions": ["anxiety", "alcohol_use_disorder", "insomnia"],
        "key_variants": ["rs279858", "rs279871"],
        "therapeutic_class": "anxiolytic_response",
    },

    # --- Schizophrenia / Psychosis ---
    "COMT": {
        "chrom": "chr22",
        "start": 19929262,
        "end": 19969975,
        "role": "Catechol-O-methyltransferase — prefrontal dopamine metabolism",
        "strategy": "Pharmacogenomic: Val158Met (rs4680) predicts cognitive "
                     "function and antipsychotic response. Guides treatment "
                     "selection for schizophrenia patients",
        "conditions": ["schizophrenia", "ADHD", "cognitive_impairment"],
        "key_variants": ["rs4680_Val158Met"],
        "therapeutic_class": "dopamine_modulation",
    },
    "DRD2": {
        "chrom": "chr11",
        "start": 113409605,
        "end": 113475691,
        "role": "Dopamine D2 receptor — primary antipsychotic target",
        "strategy": "Pharmacogenomic: Taq1A (rs1800497) predicts antipsychotic "
                     "response and tardive dyskinesia risk",
        "conditions": ["schizophrenia", "bipolar_disorder"],
        "key_variants": ["rs1800497_Taq1A", "rs6277"],
        "therapeutic_class": "antipsychotic_response",
    },
    "DISC1": {
        "chrom": "chr1",
        "start": 231762561,
        "end": 232177019,
        "role": "Disrupted in Schizophrenia 1 — neurodevelopment and "
                "synaptic function",
        "strategy": "HDR correction of rare pathogenic variants in familial "
                     "schizophrenia cases",
        "conditions": ["schizophrenia", "schizoaffective_disorder"],
        "key_variants": ["rs821616", "Ser704Cys"],
        "therapeutic_class": "neurodevelopmental",
    },

    # --- Addiction ---
    "OPRM1": {
        "chrom": "chr6",
        "start": 154039662,
        "end": 154106498,
        "role": "Mu-opioid receptor — pain and reward processing",
        "strategy": "Pharmacogenomic: A118G (rs1799971) predicts opioid "
                     "dosing requirements and naltrexone response for "
                     "alcohol/opioid use disorder treatment",
        "conditions": ["opioid_use_disorder", "alcohol_use_disorder",
                       "chronic_pain"],
        "key_variants": ["rs1799971_A118G"],
        "therapeutic_class": "addiction_pharmacogenomics",
    },
    "ALDH2": {
        "chrom": "chr12",
        "start": 111766887,
        "end": 111817532,
        "role": "Aldehyde dehydrogenase 2 — alcohol metabolism",
        "strategy": "Pharmacogenomic: rs671 (Glu504Lys) affects alcohol "
                     "metabolism and disulfiram-like response. Gene therapy: "
                     "base-editing to restore ALDH2 activity in deficient patients",
        "conditions": ["alcohol_use_disorder", "acetaldehyde_toxicity"],
        "key_variants": ["rs671_Glu504Lys"],
        "therapeutic_class": "addiction_pharmacogenomics",
    },

    # --- ADHD ---
    "SLC6A3": {
        "chrom": "chr5",
        "start": 1392790,
        "end": 1445430,
        "role": "Dopamine transporter (DAT1) — stimulant drug target",
        "strategy": "Pharmacogenomic: 3'UTR VNTR predicts methylphenidate "
                     "response in ADHD",
        "conditions": ["ADHD"],
        "key_variants": ["DAT1_VNTR_10R/9R"],
        "therapeutic_class": "stimulant_response",
    },
}
