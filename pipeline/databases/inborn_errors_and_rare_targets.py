"""
Curated gene-therapy / CRISPR targets for inborn errors of metabolism and
rare monogenic diseases -- GRCh38 (hg38) coordinates.

Covers 16 disease categories (53 genes total):
  1.  Renal tubular disorders
  2.  Glycogen storage diseases
  3.  Peroxisomal disorders
  4.  Urea cycle defects
  5.  Organic acidemias
  6.  Fatty acid oxidation defects
  7.  Amino acid disorders
  8.  Carbohydrate metabolism disorders
  9.  Sphingolipidoses
  10. Immunoglobulin deficiencies
  11. Platelet / bleeding disorders
  12. Iron metabolism disorders
  13. Copper metabolism disorders
  14. Porphyrias
  15. Ciliopathies
  16. RASopathies

Coordinates are from NCBI Gene, GRCh38.p14 (RefSeq annotation RS_2025_08),
retrieved 2026-04-02.

IMPORTANT -- verify every coordinate against current NCBI / Ensembl releases
before production use.  Numbering can shift between patch levels.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14)
  - OMIM
  - ClinVar
  - ClinicalTrials.gov
  - CRISPR Medicine News
  - Innovative Genomics Institute clinical-trial tracker
  - Published literature through early 2026
"""


# ============================================================================
# Unified dictionary -- every gene keyed by HGNC symbol
# ============================================================================

ALL_INBORN_ERRORS_TARGETS = {

    # ========================================================================
    # 1.  RENAL TUBULAR DISORDERS
    # ========================================================================

    "SLC12A3": {
        "chrom": "chr16",
        "start": 56_865_207,
        "end": 56_915_850,
        "role": (
            "Thiazide-sensitive sodium-chloride cotransporter (NCC) in the "
            "distal convoluted tubule.  Loss-of-function causes Gitelman "
            "syndrome -- autosomal recessive salt-wasting tubulopathy with "
            "hypokalemic metabolic alkalosis, hypomagnesemia, and hypocalciuria."
        ),
        "strategy": (
            "Ex vivo base editing or prime editing of patient renal progenitor "
            "cells, or in vivo delivery of corrective constructs via kidney-"
            "tropic LNPs.  For common founder mutations (e.g. c.1180+1G>T in "
            "European cohorts), adenine base editing is feasible.  Alternatively, "
            "CRISPRa to upregulate residual SLC12A3 expression in hypomorphic "
            "alleles.  No gene therapy trials registered as of early 2026; "
            "standard of care is lifelong oral Mg/K supplementation."
        ),
        "conditions": [
            "gitelman_syndrome",
            "renal_tubular_disorder",
            "salt_wasting_tubulopathy",
            "hypokalemia",
        ],
    },

    "CLCNKB": {
        "chrom": "chr1",
        "start": 16_043_782,
        "end": 16_057_326,
        "role": (
            "Chloride voltage-gated channel Kb -- basolateral chloride channel "
            "in the thick ascending limb and distal convoluted tubule.  Biallelic "
            "loss-of-function causes Bartter syndrome type III (classic Bartter), "
            "a salt-wasting nephropathy with hypokalemic alkalosis, "
            "hypercalciuria, and nephrocalcinosis."
        ),
        "strategy": (
            "In vivo gene replacement via kidney-targeted AAV or LNP delivery "
            "of full-length CLCNKB cDNA (~2.1 kb, within AAV packaging limits).  "
            "For point mutations, adenine or cytosine base editing in renal "
            "tubular cells.  No clinical trials as of early 2026; current "
            "management is indomethacin + electrolyte replacement."
        ),
        "conditions": [
            "bartter_syndrome_type_III",
            "renal_tubular_disorder",
            "salt_wasting_tubulopathy",
            "hypokalemia",
        ],
    },

    "KCNJ1": {
        "chrom": "chr11",
        "start": 128_838_020,
        "end": 128_867_296,
        "role": (
            "Inward-rectifying potassium channel Kir1.1 (ROMK).  Apical K+ "
            "recycling channel in the thick ascending limb essential for NaCl "
            "reabsorption.  Biallelic loss-of-function causes Bartter syndrome "
            "type II (antenatal/neonatal Bartter) with severe polyhydramnios, "
            "premature birth, and neonatal salt wasting."
        ),
        "strategy": (
            "In vivo gene replacement with kidney-tropic AAV delivering ROMK "
            "cDNA (~1.2 kb).  Base editing for recurrent missense mutations.  "
            "Neonatal onset demands early intervention; preclinical kidney-"
            "targeted delivery remains a major hurdle.  No clinical trials "
            "registered.  Supportive care: indomethacin, K-sparing diuretics."
        ),
        "conditions": [
            "bartter_syndrome_type_II",
            "renal_tubular_disorder",
            "salt_wasting_tubulopathy",
            "neonatal_bartter",
        ],
    },

    # ========================================================================
    # 2.  GLYCOGEN STORAGE DISEASES
    # ========================================================================

    "AGL": {
        "chrom": "chr1",
        "start": 99_849_258,
        "end": 99_924_023,
        "role": (
            "Amylo-1,6-glucosidase / 4-alpha-glucanotransferase (glycogen "
            "debranching enzyme).  Biallelic loss-of-function causes GSD III "
            "(Cori/Forbes disease) -- hepatomyopathy with fasting hypoglycemia, "
            "hepatomegaly, progressive myopathy and cardiomyopathy."
        ),
        "strategy": (
            "Liver-directed dual-AAV gene replacement (AGL cDNA ~4.6 kb exceeds "
            "single AAV capacity; split-intein or trans-splicing dual-vector "
            "approach).  Ultragenyx UX053 dual-AAV program in preclinical "
            "development.  Alternative: in vivo base editing for common founder "
            "mutations (e.g. c.16C>T in North African cohort).  Cornstarch "
            "therapy is current standard of care."
        ),
        "conditions": [
            "GSD_III",
            "cori_disease",
            "forbes_disease",
            "glycogen_storage_disease",
            "hepatomyopathy",
        ],
    },

    "GBE1": {
        "chrom": "chr3",
        "start": 81_489_703,
        "end": 81_761_645,
        "role": (
            "1,4-alpha-glucan-branching enzyme (glycogen branching enzyme).  "
            "Deficiency causes GSD IV (Andersen disease) -- accumulation of "
            "abnormal amylopectin-like glycogen (polyglucosan bodies) leading to "
            "progressive hepatic cirrhosis, or neuromuscular/cardiac forms."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement delivering GBE1 cDNA (~2.2 kb).  "
            "The severe hepatic form often requires liver transplant in infancy.  "
            "Gene therapy could avert transplant if delivered early.  For the "
            "adult polyglucosan body disease form (APBD), CNS-directed AAV9 "
            "delivery may be needed.  No active clinical trials as of early 2026."
        ),
        "conditions": [
            "GSD_IV",
            "andersen_disease",
            "glycogen_storage_disease",
            "polyglucosan_body_disease",
            "hepatic_cirrhosis",
        ],
    },

    "PYGM": {
        "chrom": "chr11",
        "start": 64_746_389,
        "end": 64_760_715,
        "role": (
            "Muscle glycogen phosphorylase.  Deficiency causes GSD V (McArdle "
            "disease) -- exercise intolerance with myalgia, rhabdomyolysis risk, "
            "and 'second wind' phenomenon.  Most common glycogen storage "
            "myopathy (~1:100,000)."
        ),
        "strategy": (
            "Intramuscular or systemic AAV gene replacement (PYGM cDNA ~2.5 kb).  "
            "Spark/Genethon preclinical AAV programs.  The p.Arg50Ter (R50X) "
            "nonsense mutation accounts for ~50% of alleles in European patients "
            "and is amenable to adenine base editing (premature stop -> sense "
            "codon) or antisense-mediated readthrough.  Exercise avoidance and "
            "oral sucrose are current management."
        ),
        "conditions": [
            "GSD_V",
            "mcardle_disease",
            "glycogen_storage_disease",
            "exercise_intolerance",
            "rhabdomyolysis",
        ],
    },

    "PHKA2": {
        "chrom": "chrX",
        "start": 18_892_298,
        "end": 18_984_114,
        "role": (
            "Phosphorylase kinase regulatory subunit alpha 2 (liver isoform).  "
            "Mutations cause GSD IXa -- the most common hepatic phosphorylase "
            "kinase deficiency.  X-linked; affected boys have hepatomegaly, "
            "fasting ketotic hypoglycemia, growth retardation.  Often improves "
            "with age."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (PHKA2 cDNA ~3.5 kb).  The "
            "relatively benign course and spontaneous improvement in many "
            "patients makes the risk-benefit calculation for gene therapy "
            "nuanced.  Base editing for recurrent mutations in severe cases.  "
            "Dietary management (frequent meals, cornstarch) is standard."
        ),
        "conditions": [
            "GSD_IXa",
            "glycogen_storage_disease",
            "phosphorylase_kinase_deficiency",
            "hepatomegaly",
        ],
    },

    # ========================================================================
    # 3.  PEROXISOMAL DISORDERS
    # ========================================================================

    "PEX1": {
        "chrom": "chr7",
        "start": 92_487_025,
        "end": 92_528_520,
        "role": (
            "Peroxisomal biogenesis factor 1 -- AAA-ATPase required for "
            "peroxisomal matrix protein import.  Biallelic loss-of-function "
            "causes Zellweger spectrum disorder (ZSD), ranging from lethal "
            "neonatal Zellweger syndrome to milder infantile Refsum disease.  "
            "PEX1 is the most commonly mutated gene in ZSD (~65% of cases)."
        ),
        "strategy": (
            "Ex vivo lentiviral gene replacement in patient fibroblasts "
            "(proof-of-concept shown).  In vivo liver-directed AAV-PEX1 or "
            "systemic AAV9-PEX1 for multi-organ disease.  The common "
            "p.Gly843Asp hypomorphic allele (~30% of alleles) is amenable to "
            "adenine base editing.  No clinical trials as of early 2026; "
            "supportive care only."
        ),
        "conditions": [
            "zellweger_spectrum_disorder",
            "peroxisomal_biogenesis_disorder",
            "infantile_refsum",
            "neonatal_adrenoleukodystrophy",
        ],
    },

    "ABCD1": {
        "chrom": "chrX",
        "start": 153_724_856,
        "end": 153_744_755,
        "role": (
            "ATP-binding cassette transporter D1 -- peroxisomal membrane "
            "transporter for very-long-chain fatty acids (VLCFA).  Mutations "
            "cause X-linked adrenoleukodystrophy (X-ALD / Lorenzo's oil "
            "disease) -- VLCFA accumulation leads to cerebral demyelination "
            "(childhood cerebral ALD) or progressive myelopathy (AMN)."
        ),
        "strategy": (
            "Ex vivo lentiviral gene replacement in autologous CD34+ HSPCs "
            "followed by myeloablative conditioning and transplant.  Bluebird "
            "Bio's Skysona (elivaldogene autotemcel / Lenti-D) received FDA "
            "accelerated approval (2022) for early cerebral ALD in boys.  "
            "Next-generation approaches: base editing of common mutations, "
            "AAV9-mediated CNS gene therapy for AMN.  Lorenzo's oil (dietary "
            "VLCFA reduction) has limited efficacy."
        ),
        "conditions": [
            "x_linked_adrenoleukodystrophy",
            "cerebral_ALD",
            "adrenomyeloneuropathy",
            "lorenzos_oil_disease",
            "peroxisomal_disorder",
        ],
    },

    "PEX7": {
        "chrom": "chr6",
        "start": 136_822_592,
        "end": 136_913_934,
        "role": (
            "Peroxisomal biogenesis factor 7 -- PTS2 receptor that imports "
            "peroxisomal matrix proteins bearing a type-2 targeting signal.  "
            "Biallelic loss-of-function causes rhizomelic chondrodysplasia "
            "punctata type 1 (RCDP1) -- severe skeletal dysplasia with "
            "proximal limb shortening, cataracts, and intellectual disability."
        ),
        "strategy": (
            "Systemic AAV gene replacement delivering PEX7 cDNA.  Multi-organ "
            "involvement (bone, brain, lens) complicates delivery.  The common "
            "c.875T>A (p.Leu292Ter) mutation in European populations may be "
            "targetable by base editing or antisense readthrough.  No clinical "
            "trials; supportive/palliative care is current standard."
        ),
        "conditions": [
            "RCDP_type_1",
            "rhizomelic_chondrodysplasia_punctata",
            "peroxisomal_disorder",
            "skeletal_dysplasia",
        ],
    },

    # ========================================================================
    # 4.  UREA CYCLE DEFECTS
    # ========================================================================

    "CPS1": {
        "chrom": "chr2",
        "start": 210_477_685,
        "end": 210_679_107,
        "role": (
            "Carbamoyl-phosphate synthetase 1 -- mitochondrial enzyme catalyzing "
            "the first committed step of the urea cycle (ammonia + bicarbonate "
            "-> carbamoyl phosphate).  Deficiency causes CPS1 deficiency, "
            "the rarest but often most severe urea cycle defect with neonatal "
            "hyperammonemia."
        ),
        "strategy": (
            "Liver-directed dual-AAV gene replacement (CPS1 cDNA ~4.5 kb "
            "approaches single-AAV limits; dual-vector or helper-dependent "
            "adenovirus may be needed).  mRNA-LNP hepatic delivery of CPS1 "
            "mRNA is an alternative bridging strategy.  Nitrogen scavengers "
            "(sodium benzoate, phenylbutyrate) and protein restriction are "
            "current management; liver transplant is curative."
        ),
        "conditions": [
            "CPS1_deficiency",
            "urea_cycle_defect",
            "hyperammonemia",
            "neonatal_metabolic_crisis",
        ],
    },

    "ARG1": {
        "chrom": "chr6",
        "start": 131_573_226,
        "end": 131_584_329,
        "role": (
            "Arginase 1 -- cytosolic enzyme catalyzing the final step of the "
            "urea cycle (arginine -> urea + ornithine).  Deficiency causes "
            "arginase deficiency (hyperargininemia) -- progressive spastic "
            "diplegia, intellectual disability, and episodic hyperammonemia.  "
            "Distinct from other UCD by prominent neurological features."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (ARG1 cDNA ~1.0 kb, easily "
            "packaged in AAV8/AAV5).  Preclinical AAV-ARG1 studies in arg1-/- "
            "mice show correction of hyperargininemia and spasticity.  Base "
            "editing feasible for recurrent point mutations.  Pegzilarginase "
            "(enzyme replacement) is in clinical development (Aeglea).  Protein "
            "restriction + nitrogen scavengers are standard care."
        ),
        "conditions": [
            "arginase_deficiency",
            "hyperargininemia",
            "urea_cycle_defect",
            "spastic_diplegia",
        ],
    },

    "SLC25A13": {
        "chrom": "chr7",
        "start": 96_120_220,
        "end": 96_322_098,
        "role": (
            "Citrin -- mitochondrial aspartate-glutamate carrier AGC2.  "
            "Deficiency causes citrullinemia type II (CTLN2) in adults "
            "(recurrent hyperammonemia, neuropsychiatric symptoms, fatty liver) "
            "and neonatal intrahepatic cholestasis caused by citrin deficiency "
            "(NICCD).  Prevalent in East Asian populations."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (SLC25A13 cDNA ~2.0 kb).  "
            "The adult-onset CTLN2 phenotype is triggered by high-carbohydrate "
            "diet; dietary intervention (high-protein, low-carb, MCT "
            "supplementation) is first-line.  Liver transplant is curative for "
            "severe CTLN2.  No gene therapy trials registered; preclinical "
            "AAV-citrin studies in murine models show promise."
        ),
        "conditions": [
            "citrullinemia_type_II",
            "citrin_deficiency",
            "NICCD",
            "urea_cycle_defect",
            "neonatal_cholestasis",
        ],
    },

    # ========================================================================
    # 5.  ORGANIC ACIDEMIAS
    # ========================================================================

    "MMUT": {
        "chrom": "chr6",
        "start": 49_430_360,
        "end": 49_463_253,
        "role": (
            "Methylmalonyl-CoA mutase (formerly MUT) -- mitochondrial enzyme "
            "converting methylmalonyl-CoA to succinyl-CoA (vitamin B12-dependent).  "
            "Deficiency causes methylmalonic acidemia (MMA) -- recurrent metabolic "
            "crises with ketoacidosis, hyperammonemia, and progressive renal and "
            "neurological damage."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement: LogicBio LB-001 (AAVhu37-MMUT) "
            "reached Phase 1/2 clinical trial.  mRNA-LNP approaches (Moderna "
            "mRNA-3705) also in clinical trials for repeated hepatic expression.  "
            "Dual-function correction: liver gene therapy reduces toxic "
            "metabolites systemically even though kidney/brain also affected.  "
            "B12-responsive forms may respond to cobalamin therapy."
        ),
        "conditions": [
            "methylmalonic_acidemia",
            "organic_acidemia",
            "MMA",
            "metabolic_crisis",
            "renal_impairment",
        ],
    },

    "PCCA": {
        "chrom": "chr13",
        "start": 100_089_093,
        "end": 100_530_435,
        "role": (
            "Propionyl-CoA carboxylase subunit alpha.  Together with PCCB forms "
            "the mitochondrial biotin-dependent enzyme converting propionyl-CoA "
            "to D-methylmalonyl-CoA.  Deficiency of either subunit causes "
            "propionic acidemia (PA) -- life-threatening metabolic crises with "
            "ketoacidosis, hyperammonemia, and cardiomyopathy."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (PCCA cDNA ~2.2 kb).  "
            "Preclinical studies in Pcca-/- mice with AAV8-PCCA show long-term "
            "metabolic correction.  Homology Medicines exploring in vivo gene "
            "editing approaches.  mRNA-LNP hepatic delivery is an alternative.  "
            "Protein restriction, carnitine supplementation, and emergency "
            "protocols for metabolic crises are standard care.  Liver transplant "
            "reduces crisis frequency but does not prevent cardiomyopathy."
        ),
        "conditions": [
            "propionic_acidemia",
            "organic_acidemia",
            "metabolic_crisis",
            "cardiomyopathy",
        ],
    },

    "PCCB": {
        "chrom": "chr3",
        "start": 136_250_340,
        "end": 136_330_169,
        "role": (
            "Propionyl-CoA carboxylase subunit beta.  Biallelic mutations cause "
            "propionic acidemia indistinguishable from PCCA-mutant PA.  PCCB "
            "mutations account for ~75% of PA cases worldwide."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (PCCB cDNA ~1.6 kb).  Same "
            "therapeutic rationale as PCCA.  The c.1218_1231del14ins12 mutation "
            "is common in European PA patients.  Base editing or prime editing "
            "could correct specific recurrent variants.  Clinical management "
            "identical to PCCA deficiency."
        ),
        "conditions": [
            "propionic_acidemia",
            "organic_acidemia",
            "metabolic_crisis",
            "cardiomyopathy",
        ],
    },

    "IVD": {
        "chrom": "chr15",
        "start": 40_405_795,
        "end": 40_435_947,
        "role": (
            "Isovaleryl-CoA dehydrogenase -- mitochondrial FAD-dependent enzyme "
            "in leucine catabolism.  Deficiency causes isovaleric acidemia (IVA) "
            "-- acute neonatal presentation with 'sweaty feet' odor, metabolic "
            "acidosis, and pancytopenia, or chronic intermittent form."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (IVD cDNA ~1.2 kb).  "
            "Relatively small cDNA is easily packaged.  Leucine-restricted diet "
            "and glycine/carnitine supplementation are effective in most patients.  "
            "Gene therapy may benefit severe neonatal cases unresponsive to "
            "dietary management.  No clinical trials as of early 2026."
        ),
        "conditions": [
            "isovaleric_acidemia",
            "organic_acidemia",
            "leucine_metabolism",
            "metabolic_crisis",
        ],
    },

    # ========================================================================
    # 6.  FATTY ACID OXIDATION DEFECTS
    # ========================================================================

    "ACADM": {
        "chrom": "chr1",
        "start": 75_724_709,
        "end": 75_763_679,
        "role": (
            "Medium-chain acyl-CoA dehydrogenase.  Deficiency causes MCAD "
            "deficiency -- the most common fatty acid oxidation disorder "
            "(~1:10,000 by newborn screening).  Fasting intolerance with "
            "hypoketotic hypoglycemia, Reye-like episodes; potentially fatal "
            "if undiagnosed."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (ACADM cDNA ~1.3 kb).  The "
            "c.985A>G (p.Lys329Glu, K329E) mutation accounts for ~80% of "
            "alleles in Northern European patients and is an ideal target for "
            "adenine base editing (G->A reversion at the DNA level).  Newborn "
            "screening and fasting avoidance are highly effective; gene therapy "
            "would benefit non-compliant or diagnostically missed patients."
        ),
        "conditions": [
            "MCAD_deficiency",
            "fatty_acid_oxidation_defect",
            "hypoketotic_hypoglycemia",
            "newborn_screening_target",
        ],
    },

    "ACADL": {
        "chrom": "chr2",
        "start": 210_187_923,
        "end": 210_225_447,
        "role": (
            "Long-chain acyl-CoA dehydrogenase.  LCAD deficiency is rare and "
            "historically controversial (many early cases reclassified as VLCAD "
            "deficiency).  True ACADL deficiency presents with fasting-induced "
            "hypoketotic hypoglycemia and cardiomyopathy."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (ACADL cDNA ~1.3 kb).  Given "
            "rarity and diagnostic overlap with VLCAD deficiency, genetic "
            "confirmation is essential before targeting.  Dietary management "
            "(fasting avoidance, MCT supplementation) is standard.  No gene "
            "therapy clinical trials."
        ),
        "conditions": [
            "LCAD_deficiency",
            "fatty_acid_oxidation_defect",
            "hypoketotic_hypoglycemia",
            "cardiomyopathy",
        ],
    },

    "HADHA": {
        "chrom": "chr2",
        "start": 26_190_635,
        "end": 26_244_632,
        "role": (
            "Hydroxyacyl-CoA dehydrogenase trifunctional multienzyme complex "
            "subunit alpha.  With HADHB, forms the mitochondrial trifunctional "
            "protein (MTP).  HADHA mutations cause LCHAD deficiency or complete "
            "MTP deficiency -- fasting intolerance, cardiomyopathy, "
            "rhabdomyolysis, peripheral neuropathy, and retinopathy."
        ),
        "strategy": (
            "Liver- and muscle-directed AAV gene replacement (HADHA cDNA ~2.3 kb).  "
            "Multi-organ involvement complicates gene therapy.  The common "
            "c.1528G>C (p.Glu510Gln, E510Q) mutation (~90% of LCHAD alleles) "
            "is a prime target for cytosine base editing.  MCT-based diet and "
            "fasting avoidance are current standard; triheptanoin (Dojolvi) is "
            "FDA-approved for LC-FAO disorders."
        ),
        "conditions": [
            "LCHAD_deficiency",
            "MTP_deficiency",
            "fatty_acid_oxidation_defect",
            "cardiomyopathy",
            "retinopathy",
        ],
    },

    "CPT2": {
        "chrom": "chr1",
        "start": 53_196_824,
        "end": 53_214_197,
        "role": (
            "Carnitine palmitoyltransferase 2 -- inner mitochondrial membrane "
            "enzyme completing long-chain fatty acid import into the matrix.  "
            "Three clinical forms: lethal neonatal (multi-organ), severe "
            "infantile (hepatocardiomuscular), and myopathic adult form "
            "(exercise-induced rhabdomyolysis -- most common)."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (CPT2 cDNA ~2.0 kb) for "
            "severe hepatic forms.  Muscle-directed AAV for the myopathic form.  "
            "The p.Ser113Leu (c.338C>T) thermolabile variant accounts for ~60% "
            "of alleles in the adult myopathic form and is amenable to cytosine "
            "base editing.  Bezafibrate (PPARa agonist) upregulates residual "
            "CPT2 activity and is used off-label.  Triheptanoin approved for "
            "LC-FAO disorders."
        ),
        "conditions": [
            "CPT2_deficiency",
            "fatty_acid_oxidation_defect",
            "rhabdomyolysis",
            "exercise_intolerance",
        ],
    },

    # ========================================================================
    # 7.  AMINO ACID DISORDERS
    # ========================================================================

    "ASL": {
        "chrom": "chr7",
        "start": 66_075_819,
        "end": 66_093_576,
        "role": (
            "Argininosuccinate lyase -- urea cycle enzyme cleaving "
            "argininosuccinate to arginine + fumarate.  Deficiency causes "
            "argininosuccinic aciduria (ASA) -- hyperammonemia, liver disease, "
            "and neurocognitive impairment even with good metabolic control "
            "(due to cell-autonomous NO/nitric oxide deficiency)."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (ASL cDNA ~1.4 kb).  "
            "Preclinical AAV-ASL in Asl-/- mice restores ureagenesis and "
            "normalizes argininosuccinate.  However, brain-autonomous NO "
            "deficiency may require additional CNS-directed therapy.  Arginine "
            "supplementation, nitrogen scavengers, and protein restriction "
            "are standard.  Liver transplant corrects urea cycle but not "
            "neurocognitive trajectory."
        ),
        "conditions": [
            "argininosuccinic_aciduria",
            "amino_acid_disorder",
            "urea_cycle_defect",
            "hyperammonemia",
        ],
    },

    "SLC7A7": {
        "chrom": "chr14",
        "start": 22_773_222,
        "end": 22_819_791,
        "role": (
            "y+L amino acid transporter 1 -- basolateral cationic amino acid "
            "transporter in intestine and kidney.  Deficiency causes lysinuric "
            "protein intolerance (LPI) -- impaired absorption/reabsorption of "
            "lysine, arginine, ornithine leading to hyperammonemia, protein "
            "aversion, osteoporosis, and pulmonary alveolar proteinosis."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement to restore hepatic arginine/"
            "ornithine transport and urea cycle function.  Intestinal correction "
            "would require tissue-specific promoters or separate targeting.  "
            "Citrulline supplementation bypasses the transport defect for urea "
            "cycle support.  No gene therapy trials.  PAP (pulmonary alveolar "
            "proteinosis) is the major cause of mortality and may require "
            "macrophage-directed gene therapy."
        ),
        "conditions": [
            "lysinuric_protein_intolerance",
            "amino_acid_disorder",
            "hyperammonemia",
            "pulmonary_alveolar_proteinosis",
        ],
    },

    "SLC3A1": {
        "chrom": "chr2",
        "start": 44_275_480,
        "end": 44_322_437,
        "role": (
            "Neutral and basic amino acid transport protein rBAT -- heavy "
            "subunit of the heterodimeric amino acid transporter b0,+AT/rBAT "
            "in the proximal tubule and small intestine.  Mutations cause "
            "cystinuria type A -- impaired renal reabsorption of cystine "
            "leading to recurrent cystine kidney stones."
        ),
        "strategy": (
            "Kidney-directed gene therapy to restore proximal tubule cystine "
            "reabsorption.  Kidney-tropic AAV serotypes (AAV2, AAV9) or LNP "
            "delivery of SLC3A1 mRNA.  Technically challenging due to need for "
            "apical membrane expression in proximal tubular cells.  Current "
            "management: hyperhydration, urine alkalinization, tiopronin/"
            "D-penicillamine.  No clinical trials for gene therapy."
        ),
        "conditions": [
            "cystinuria",
            "amino_acid_disorder",
            "kidney_stones",
            "nephrolithiasis",
        ],
    },

    # ========================================================================
    # 8.  CARBOHYDRATE METABOLISM DISORDERS
    # ========================================================================

    "GALK1": {
        "chrom": "chr17",
        "start": 75_751_469,
        "end": 75_765_192,
        "role": (
            "Galactokinase -- first enzyme in the Leloir pathway of galactose "
            "metabolism (galactose -> galactose-1-phosphate).  Deficiency causes "
            "galactokinase deficiency (galactosemia type II) -- galactitol "
            "accumulation leads to bilateral cataracts; milder than classic "
            "galactosemia (GALT deficiency)."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (GALK1 cDNA ~1.2 kb).  "
            "Galactose-restricted diet prevents cataracts and is highly "
            "effective, making gene therapy risk-benefit less favorable.  "
            "Gene therapy may be considered for patients with poor dietary "
            "compliance or late diagnosis.  No clinical trials."
        ),
        "conditions": [
            "galactokinase_deficiency",
            "galactosemia_type_II",
            "carbohydrate_metabolism",
            "cataracts",
        ],
    },

    "SLC2A1": {
        "chrom": "chr1",
        "start": 42_925_353,
        "end": 42_958_868,
        "role": (
            "Glucose transporter 1 (GLUT1) -- primary glucose transporter at "
            "the blood-brain barrier.  Haploinsufficiency causes GLUT1 "
            "deficiency syndrome (De Vivo disease) -- epilepsy, developmental "
            "delay, microcephaly, and movement disorder due to cerebral energy "
            "failure from impaired glucose transport into the brain."
        ),
        "strategy": (
            "CNS-directed AAV9 gene replacement delivering SLC2A1 under a "
            "neuron/endothelial promoter (intrathecal or intravenous with "
            "BBB-crossing AAV).  CRISPRa to upregulate the remaining functional "
            "allele is an alternative for haploinsufficiency.  Ketogenic diet "
            "bypasses the glucose transport defect by providing ketone bodies "
            "as alternative brain fuel -- first-line standard of care."
        ),
        "conditions": [
            "GLUT1_deficiency",
            "de_vivo_disease",
            "carbohydrate_metabolism",
            "epilepsy",
            "movement_disorder",
        ],
    },

    "ALDOB": {
        "chrom": "chr9",
        "start": 101_420_560,
        "end": 101_435_774,
        "role": (
            "Aldolase B (fructose-bisphosphate aldolase B) -- liver enzyme "
            "catalyzing fructose-1-phosphate cleavage.  Deficiency causes "
            "hereditary fructose intolerance (HFI) -- fructose/sucrose "
            "ingestion triggers hypoglycemia, vomiting, liver failure, and "
            "renal tubular dysfunction."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (ALDOB cDNA ~1.1 kb).  "
            "Preclinical AAV8-ALDOB in Aldob-/- mice shows complete rescue.  "
            "The p.Ala149Pro (A149P) mutation accounts for ~65% of European "
            "alleles and is targetable by base editing.  Strict fructose/"
            "sucrose/sorbitol avoidance is curative but socially burdensome; "
            "gene therapy could liberate diet."
        ),
        "conditions": [
            "hereditary_fructose_intolerance",
            "carbohydrate_metabolism",
            "liver_disease",
            "hypoglycemia",
        ],
    },

    # ========================================================================
    # 9.  SPHINGOLIPIDOSES
    # ========================================================================

    "SMPD1": {
        "chrom": "chr11",
        "start": 6_390_474,
        "end": 6_394_996,
        "role": (
            "Acid sphingomyelinase (ASM) -- lysosomal enzyme hydrolyzing "
            "sphingomyelin to ceramide + phosphocholine.  Deficiency causes "
            "Niemann-Pick disease type A (severe infantile neurovisceral, "
            "fatal by ~3 years) and type B (chronic visceral, hepatosplenomegaly "
            "and pulmonary disease without primary CNS involvement)."
        ),
        "strategy": (
            "Enzyme replacement therapy: olipudase alfa (Xenpozyme) received "
            "FDA approval (2022) for non-CNS manifestations of acid "
            "sphingomyelinase deficiency (ASMD type B).  For type A (CNS), "
            "intrathecal AAV-SMPD1 or AAV9-mediated CNS gene therapy is needed.  "
            "Ex vivo HSCT-based lentiviral gene therapy under investigation.  "
            "Liver-directed AAV for type B as alternative to lifelong ERT."
        ),
        "conditions": [
            "niemann_pick_A",
            "niemann_pick_B",
            "ASMD",
            "sphingolipidosis",
            "lysosomal_storage_disease",
        ],
    },

    "ASAH1": {
        "chrom": "chr8",
        "start": 18_055_992,
        "end": 18_084_961,
        "role": (
            "Acid ceramidase -- lysosomal enzyme hydrolyzing ceramide to "
            "sphingosine + fatty acid.  Deficiency causes Farber disease "
            "(Farber lipogranulomatosis) -- painful joint swelling, "
            "subcutaneous nodules, hoarse voice, and progressive "
            "neurodegeneration.  Also causes spinal muscular atrophy with "
            "progressive myoclonic epilepsy (SMA-PME)."
        ),
        "strategy": (
            "Systemic AAV gene replacement (ASAH1 cDNA ~1.2 kb) for visceral "
            "and CNS manifestations.  Ex vivo lentiviral gene therapy in HSPCs "
            "for hematopoietic/macrophage correction (analogous to Fabry/Gaucher "
            "approaches).  Recombinant acid ceramidase enzyme replacement is in "
            "preclinical development.  HSCT has shown partial benefit in mild "
            "Farber disease.  No clinical trials for gene therapy."
        ),
        "conditions": [
            "farber_disease",
            "farber_lipogranulomatosis",
            "SMA_PME",
            "sphingolipidosis",
            "lysosomal_storage_disease",
        ],
    },

    "HEXB": {
        "chrom": "chr5",
        "start": 74_640_023,
        "end": 74_721_288,
        "role": (
            "Hexosaminidase subunit beta -- with HEXA forms hexosaminidase A "
            "(alpha-beta), and self-associates to form hexosaminidase B "
            "(beta-beta).  HEXB deficiency causes Sandhoff disease -- GM2 "
            "ganglioside accumulation with progressive neurodegeneration "
            "clinically similar to Tay-Sachs but with additional visceral "
            "involvement."
        ),
        "strategy": (
            "CNS-directed AAV gene therapy: intrathecal/intracerebroventricular "
            "AAV9-HEXB.  Preclinical success in Sandhoff mice and cats.  "
            "Axovant/Sio Gene Therapies AXO-AAV-GM2 (bicistronic HEXA+HEXB) "
            "Phase 1/2 trial for GM2 gangliosidoses (NCT04669535) -- treats both "
            "Tay-Sachs and Sandhoff.  Substrate reduction therapy (miglustat) "
            "has limited CNS penetration."
        ),
        "conditions": [
            "sandhoff_disease",
            "GM2_gangliosidosis",
            "sphingolipidosis",
            "lysosomal_storage_disease",
            "neurodegeneration",
        ],
    },

    # ========================================================================
    # 10.  IMMUNOGLOBULIN DEFICIENCIES
    # ========================================================================

    "BTK": {
        "chrom": "chrX",
        "start": 101_349_450,
        "end": 101_390_796,
        "role": (
            "Bruton tyrosine kinase -- cytoplasmic tyrosine kinase essential "
            "for B-cell receptor signaling and B-cell maturation.  Loss-of-"
            "function causes X-linked agammaglobulinemia (XLA / Bruton "
            "agammaglobulinemia) -- absent mature B cells and profound "
            "hypogammaglobulinemia with recurrent bacterial infections."
        ),
        "strategy": (
            "Ex vivo lentiviral gene replacement in autologous CD34+ HSPCs.  "
            "BTK cDNA (~2.0 kb) under a B-cell-specific promoter restores "
            "B-cell development.  Preclinical lenti-BTK in Btk-/- mice "
            "rescues B-cell compartment.  CRISPR-mediated knock-in at safe "
            "harbor locus (AAVS1) is an alternative.  Lifelong IVIG "
            "replacement is current standard; gene therapy could be curative."
        ),
        "conditions": [
            "x_linked_agammaglobulinemia",
            "bruton_agammaglobulinemia",
            "XLA",
            "immunoglobulin_deficiency",
            "primary_immunodeficiency",
        ],
    },

    "PIK3CD": {
        "chrom": "chr1",
        "start": 9_627_258,
        "end": 9_729_114,
        "role": (
            "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit "
            "delta.  Gain-of-function mutations cause activated PI3K delta "
            "syndrome (APDS / PASLI) -- combined immunodeficiency with "
            "lymphoproliferation, recurrent sinopulmonary infections, "
            "herpesvirus susceptibility, and lymphoma predisposition."
        ),
        "strategy": (
            "APDS is a gain-of-function disorder; CRISPRi to downregulate the "
            "mutant allele, or allele-specific CRISPR disruption of the "
            "activating mutation (e.g. E1021K).  Leniolisib (Joenja), a "
            "selective PI3Kdelta inhibitor, received FDA approval (2023) -- "
            "first targeted oral therapy for APDS.  Gene therapy approaches "
            "are less advanced given the availability of leniolisib."
        ),
        "conditions": [
            "APDS",
            "activated_PI3K_delta_syndrome",
            "PASLI",
            "immunoglobulin_deficiency",
            "primary_immunodeficiency",
            "lymphoproliferation",
        ],
    },

    "IGHM": {
        "chrom": "chr14",
        "start": 105_851_966,
        "end": 105_856_217,
        "role": (
            "Immunoglobulin heavy constant mu -- encodes the heavy chain of IgM, "
            "the first immunoglobulin expressed during B-cell development.  "
            "Biallelic loss-of-function causes autosomal recessive "
            "agammaglobulinemia -- absent B cells and pan-hypogammaglobulinemia, "
            "clinically resembling XLA."
        ),
        "strategy": (
            "Ex vivo gene replacement in CD34+ HSPCs with lentiviral IGHM "
            "construct.  Technically challenging because IGHM must be expressed "
            "in the correct genomic context for V(D)J recombination and class "
            "switching.  Knock-in at the endogenous IGHM locus via CRISPR-HDR "
            "may be preferable.  IVIG replacement is standard; HSCT is curative."
        ),
        "conditions": [
            "agammaglobulinemia",
            "immunoglobulin_deficiency",
            "primary_immunodeficiency",
            "absent_B_cells",
        ],
    },

    # ========================================================================
    # 11.  PLATELET / BLEEDING DISORDERS
    # ========================================================================

    "ITGA2B": {
        "chrom": "chr17",
        "start": 44_372_181,
        "end": 44_389_649,
        "role": (
            "Integrin alpha-IIb -- with ITGB3 forms the platelet fibrinogen "
            "receptor GPIIb/IIIa (integrin alphaIIb-beta3).  Biallelic loss "
            "causes Glanzmann thrombasthenia (GT) -- severe mucocutaneous "
            "bleeding due to absent platelet aggregation."
        ),
        "strategy": (
            "Ex vivo lentiviral gene replacement in autologous CD34+ HSPCs "
            "using a megakaryocyte-specific promoter (e.g. GPIba or ITGA2B "
            "promoter).  Preclinical lenti-ITGA2B in GT patient iPSC-derived "
            "megakaryocytes restores fibrinogen binding.  In vivo approaches: "
            "liver-directed AAV to produce secreted thrombopoietic factors.  "
            "Platelet transfusion and antifibrinolytics are current management; "
            "rFVIIa (NovoSeven) for severe bleeds.  Gene therapy trials in "
            "planning."
        ),
        "conditions": [
            "glanzmann_thrombasthenia",
            "platelet_disorder",
            "bleeding_disorder",
            "mucocutaneous_bleeding",
        ],
    },

    "ITGB3": {
        "chrom": "chr17",
        "start": 47_253_827,
        "end": 47_313_743,
        "role": (
            "Integrin beta-3 -- with ITGA2B forms GPIIb/IIIa.  Biallelic "
            "loss also causes Glanzmann thrombasthenia.  ITGB3 mutations "
            "account for ~20-25% of GT cases."
        ),
        "strategy": (
            "Same ex vivo HSPC lentiviral gene replacement approach as ITGA2B.  "
            "ITGB3 cDNA (~2.4 kb) under megakaryocyte-specific promoter.  "
            "Both subunits must be co-expressed for functional receptor; "
            "correction of either defective gene restores the complex.  "
            "Eptifibatide/abciximab (GPIIb/IIIa inhibitors, used in cardiology) "
            "confirm the therapeutic target."
        ),
        "conditions": [
            "glanzmann_thrombasthenia",
            "platelet_disorder",
            "bleeding_disorder",
        ],
    },

    "GP1BA": {
        "chrom": "chr17",
        "start": 4_932_277,
        "end": 4_935_023,
        "role": (
            "Glycoprotein Ib platelet subunit alpha -- major subunit of the "
            "GPIb-IX-V complex, the platelet von Willebrand factor (vWF) "
            "receptor.  Biallelic loss causes Bernard-Soulier syndrome (BSS) "
            "-- macrothrombocytopenia with severe bleeding due to defective "
            "platelet adhesion to subendothelium."
        ),
        "strategy": (
            "Ex vivo lentiviral gene replacement in HSPCs with megakaryocyte-"
            "specific promoter driving GP1BA.  Challenges include restoring "
            "proper stoichiometry of the multi-subunit GPIb-IX-V complex.  "
            "DDAVP and platelet transfusions for bleeding episodes.  "
            "Recombinant vWF may help.  No gene therapy trials registered."
        ),
        "conditions": [
            "bernard_soulier_syndrome",
            "macrothrombocytopenia",
            "platelet_disorder",
            "bleeding_disorder",
        ],
    },

    "FERMT3": {
        "chrom": "chr11",
        "start": 64_205_920,
        "end": 64_223_891,
        "role": (
            "Fermitin family member 3 (kindlin-3) -- integrin activator "
            "essential for inside-out signaling in leukocytes, platelets, and "
            "osteoclasts.  Deficiency causes leukocyte adhesion deficiency "
            "type III (LAD-III) -- recurrent severe infections (like LAD-I), "
            "Glanzmann-like bleeding, and osteopetrosis."
        ),
        "strategy": (
            "Ex vivo lentiviral gene replacement in autologous CD34+ HSPCs "
            "(FERMT3 cDNA ~2.0 kb).  Corrects both leukocyte adhesion and "
            "platelet aggregation defects simultaneously.  HSCT is curative "
            "but carries transplant-related morbidity.  Gene therapy could "
            "avoid GVHD risk.  No clinical trials as of early 2026."
        ),
        "conditions": [
            "LAD_III",
            "leukocyte_adhesion_deficiency",
            "platelet_disorder",
            "bleeding_disorder",
            "osteopetrosis",
            "immunodeficiency",
        ],
    },

    # ========================================================================
    # 12.  IRON METABOLISM DISORDERS
    # ========================================================================

    "HFE": {
        "chrom": "chr6",
        "start": 26_087_429,
        "end": 26_098_343,
        "role": (
            "Homeostatic iron regulator -- MHC class I-like molecule that "
            "modulates hepcidin expression in response to iron status.  The "
            "p.Cys282Tyr (C282Y) homozygous genotype causes hereditary "
            "hemochromatosis type 1 -- progressive iron overload affecting "
            "liver, heart, pancreas, and joints.  Most common genetic disorder "
            "in Northern Europeans (~1:200 homozygous)."
        ),
        "strategy": (
            "Gene therapy is not the primary approach given effectiveness of "
            "phlebotomy.  However, for patients intolerant of phlebotomy: "
            "CRISPRa upregulation of HAMP (hepcidin) to reduce intestinal iron "
            "absorption, or liver-directed AAV-HFE gene replacement.  Base "
            "editing to correct C282Y (rs1800562, G>A) is technically "
            "straightforward.  Phlebotomy remains first-line; iron chelation "
            "(deferasirox) for those who cannot tolerate phlebotomy."
        ),
        "conditions": [
            "hereditary_hemochromatosis",
            "iron_overload",
            "hemochromatosis_type_1",
            "HFE_hemochromatosis",
        ],
    },

    "TFR2": {
        "chrom": "chr7",
        "start": 100_620_420,
        "end": 100_641_552,
        "role": (
            "Transferrin receptor 2 -- hepatic iron sensor that signals to "
            "hepcidin via the BMP-SMAD pathway.  Biallelic loss-of-function "
            "causes hemochromatosis type 3 -- iron overload similar to HFE "
            "hemochromatosis but typically with earlier onset and more "
            "severe phenotype."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (TFR2 cDNA ~2.5 kb) to "
            "restore hepcidin regulation.  Base editing for specific mutations.  "
            "CRISPRa activation of HAMP as a downstream bypass strategy.  "
            "Phlebotomy and iron chelation are effective standard treatments.  "
            "Gene therapy could benefit patients with severe juvenile-onset "
            "phenotype.  No clinical trials."
        ),
        "conditions": [
            "hemochromatosis_type_3",
            "iron_overload",
            "hereditary_hemochromatosis",
        ],
    },

    "SLC40A1": {
        "chrom": "chr2",
        "start": 189_560_590,
        "end": 189_580_786,
        "role": (
            "Ferroportin (SLC40A1) -- the only known cellular iron exporter, "
            "expressed on enterocytes, macrophages, and hepatocytes.  "
            "Loss-of-function mutations cause ferroportin disease (hemochromatosis "
            "type 4A, autosomal dominant) -- macrophage iron loading with "
            "hyperferritinemia.  Gain-of-function (hepcidin-resistant) mutations "
            "cause type 4B -- hepatocyte iron loading resembling classic HH."
        ),
        "strategy": (
            "Autosomal dominant -- allele-specific CRISPR disruption of the "
            "mutant allele for gain-of-function type 4B.  For loss-of-function "
            "type 4A, CRISPRa or gene replacement to boost wild-type allele "
            "expression.  Phlebotomy is effective for type 4B; type 4A patients "
            "tolerate phlebotomy poorly (anemia).  Low-dose iron chelation "
            "may be used.  No gene therapy trials."
        ),
        "conditions": [
            "ferroportin_disease",
            "hemochromatosis_type_4",
            "iron_overload",
            "hereditary_hemochromatosis",
        ],
    },

    "HAMP": {
        "chrom": "chr19",
        "start": 35_282_528,
        "end": 35_285_143,
        "role": (
            "Hepcidin antimicrobial peptide -- master regulator of systemic "
            "iron homeostasis.  Secreted by hepatocytes; binds and degrades "
            "ferroportin to limit iron absorption and macrophage iron release.  "
            "Biallelic loss-of-function causes hemochromatosis type 2B (juvenile "
            "hemochromatosis) -- severe early-onset iron overload with "
            "cardiomyopathy and hypogonadism, often fatal by age 30 without "
            "treatment."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (HAMP cDNA is tiny, ~0.3 kb, "
            "easily packaged with liver-specific promoter).  mRNA-LNP delivery "
            "of HAMP mRNA for repeated dosing.  Synthetic hepcidin mimetics "
            "(e.g. rusfertide/PTG-300) are in clinical trials and may obviate "
            "gene therapy.  Aggressive phlebotomy + iron chelation are standard "
            "but compliance is challenging in young patients."
        ),
        "conditions": [
            "hemochromatosis_type_2B",
            "juvenile_hemochromatosis",
            "iron_overload",
            "cardiomyopathy",
        ],
    },

    # ========================================================================
    # 13.  COPPER METABOLISM DISORDERS
    # ========================================================================

    "ATP7A": {
        "chrom": "chrX",
        "start": 77_910_693,
        "end": 78_050_395,
        "role": (
            "Copper-transporting ATPase 1 -- P-type ATPase responsible for "
            "intestinal copper absorption and delivery of copper to secreted "
            "cuproenzymes (dopamine beta-hydroxylase, lysyl oxidase, tyrosinase, "
            "cytochrome c oxidase).  Loss-of-function causes Menkes disease -- "
            "X-linked lethal copper deficiency with progressive "
            "neurodegeneration, connective tissue abnormalities, kinky hair, "
            "and death by age 3 in classic form."
        ),
        "strategy": (
            "CNS-directed AAV9 gene replacement (ATP7A cDNA ~4.5 kb, at the "
            "limit of AAV packaging).  Intracerebroventricular AAV9-ATP7A in "
            "Menkes mouse (mottled-brindled) shows neurological rescue.  NIH "
            "clinical trial of subcutaneous copper histidinate shows benefit if "
            "started neonatally (before neurodegeneration).  Combined early "
            "copper replacement + gene therapy may be optimal.  Dual-AAV "
            "strategies being explored for the large cDNA."
        ),
        "conditions": [
            "menkes_disease",
            "copper_deficiency",
            "copper_metabolism",
            "neurodegeneration",
            "connective_tissue_disorder",
        ],
    },

    # ========================================================================
    # 14.  PORPHYRIAS
    # ========================================================================

    "HMBS": {
        "chrom": "chr11",
        "start": 119_084_881,
        "end": 119_093_549,
        "role": (
            "Hydroxymethylbilane synthase (porphobilinogen deaminase) -- third "
            "enzyme in the heme biosynthesis pathway.  Heterozygous loss-of-"
            "function causes acute intermittent porphyria (AIP) -- recurrent "
            "neurovisceral attacks (abdominal pain, neuropathy, psychiatric "
            "symptoms) triggered by drugs, hormones, or fasting."
        ),
        "strategy": (
            "Liver-directed AAV5 gene replacement: Alnylam/Ultragenyx AAV5-HMBS "
            "in Phase 1/2 (ILLUMINATE-1 trial).  Small interfering RNA: "
            "givosiran (Givlaari), FDA-approved (2019), targets ALAS1 mRNA to "
            "reduce toxic ALA/PBG accumulation -- first approved RNAi therapy "
            "for porphyria.  Gene therapy aims to avoid lifelong givosiran "
            "injections.  Hemin infusions for acute attacks."
        ),
        "conditions": [
            "acute_intermittent_porphyria",
            "AIP",
            "porphyria",
            "neurovisceral_attacks",
        ],
    },

    "PPOX": {
        "chrom": "chr1",
        "start": 161_165_728,
        "end": 161_178_013,
        "role": (
            "Protoporphyrinogen oxidase -- penultimate enzyme in heme "
            "biosynthesis.  Heterozygous loss-of-function causes variegate "
            "porphyria (VP) -- both acute neurovisceral attacks (like AIP) "
            "AND cutaneous photosensitivity (skin fragility, blistering).  "
            "Prevalent in South African Afrikaners due to founder effect "
            "(p.Arg59Trp, ~1:300)."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (PPOX cDNA ~1.6 kb).  "
            "Givosiran (ALAS1 siRNA) also reduces attacks in VP by decreasing "
            "upstream porphyrin precursors, though not formally approved for VP.  "
            "Avoidance of triggers (porphyrinogenic drugs, fasting, alcohol) "
            "is critical.  Sunscreen for cutaneous symptoms.  Gene therapy "
            "could prevent both acute and cutaneous manifestations."
        ),
        "conditions": [
            "variegate_porphyria",
            "porphyria",
            "neurovisceral_attacks",
            "photosensitivity",
        ],
    },

    "UROD": {
        "chrom": "chr1",
        "start": 45_012_254,
        "end": 45_015_575,
        "role": (
            "Uroporphyrinogen decarboxylase -- fifth enzyme in heme "
            "biosynthesis.  Heterozygous mutations plus environmental triggers "
            "(iron, alcohol, hepatitis C, estrogens) cause porphyria cutanea "
            "tarda (PCT) -- the most common porphyria.  Cutaneous "
            "photosensitivity with skin blistering, hypertrichosis, and "
            "hyperpigmentation.  Homozygous UROD mutations cause hepato-"
            "erythropoietic porphyria (HEP, severe childhood form)."
        ),
        "strategy": (
            "PCT is primarily managed by removing triggers: phlebotomy to "
            "reduce iron, hydroxychloroquine, hepatitis C treatment.  Gene "
            "therapy mainly relevant for HEP (homozygous UROD deficiency), "
            "which is refractory to standard PCT treatments.  Liver-directed "
            "AAV-UROD (cDNA ~1.1 kb) or base editing for the severe HEP form."
        ),
        "conditions": [
            "porphyria_cutanea_tarda",
            "hepatoerythropoietic_porphyria",
            "porphyria",
            "photosensitivity",
        ],
    },

    "CPOX": {
        "chrom": "chr3",
        "start": 98_570_488,
        "end": 98_593_611,
        "role": (
            "Coproporphyrinogen oxidase -- sixth enzyme in heme biosynthesis.  "
            "Heterozygous loss-of-function causes hereditary coproporphyria "
            "(HCP) -- acute neurovisceral porphyria attacks and sometimes "
            "cutaneous photosensitivity.  Less common than AIP or VP."
        ),
        "strategy": (
            "Liver-directed AAV gene replacement (CPOX cDNA ~1.1 kb).  Same "
            "therapeutic logic as other acute hepatic porphyrias.  Givosiran "
            "(ALAS1 siRNA) may reduce attacks as with AIP.  Trigger avoidance "
            "is primary management.  Hemin infusions for acute attacks.  No "
            "gene therapy trials specific to HCP."
        ),
        "conditions": [
            "hereditary_coproporphyria",
            "HCP",
            "porphyria",
            "neurovisceral_attacks",
        ],
    },

    # ========================================================================
    # 15.  CILIOPATHIES
    # ========================================================================

    "BBS1": {
        "chrom": "chr11",
        "start": 66_510_635,
        "end": 66_533_598,
        "role": (
            "Bardet-Biedl syndrome 1 protein -- component of the BBSome, an "
            "octameric complex required for ciliary membrane protein trafficking.  "
            "Biallelic mutations cause Bardet-Biedl syndrome (BBS) -- retinal "
            "dystrophy, obesity, polydactyly, renal anomalies, intellectual "
            "disability, and hypogonadism."
        ),
        "strategy": (
            "Subretinal AAV gene replacement (AAV-BBS1) for retinal dystrophy "
            "component -- analogous to Luxturna (AAV-RPE65).  Systemic gene "
            "therapy for obesity/renal/cognitive features is much more complex.  "
            "The p.Met390Arg (M390R) mutation accounts for ~80% of BBS1 alleles "
            "and is targetable by base editing.  Setmelanotide (MC4R agonist) "
            "is FDA-approved (2020) for BBS-associated obesity."
        ),
        "conditions": [
            "bardet_biedl_syndrome",
            "ciliopathy",
            "retinal_dystrophy",
            "obesity",
            "renal_anomaly",
        ],
    },

    "IFT80": {
        "chrom": "chr3",
        "start": 160_256_986,
        "end": 160_399_225,
        "role": (
            "Intraflagellar transport protein 80 -- component of the IFT-B "
            "complex required for anterograde ciliary transport.  Biallelic "
            "mutations cause Jeune asphyxiating thoracic dystrophy (ATD / "
            "short-rib thoracic dysplasia type 2) -- narrow thorax with "
            "respiratory insufficiency, short limbs, renal cystic disease, "
            "and retinal degeneration."
        ),
        "strategy": (
            "Systemic AAV gene replacement is technically challenging due to "
            "multi-organ skeletal/renal/retinal involvement.  For the retinal "
            "component, subretinal AAV-IFT80 may be feasible.  Surgical "
            "thoracic expansion (lateral thoracic expansion / VEPTR) for "
            "respiratory compromise.  Renal management for progressive CKD.  "
            "No gene therapy trials.  Prognosis depends on severity of "
            "thoracic restriction."
        ),
        "conditions": [
            "jeune_syndrome",
            "asphyxiating_thoracic_dystrophy",
            "short_rib_thoracic_dysplasia",
            "ciliopathy",
            "renal_cystic_disease",
        ],
    },

    "NPHP1": {
        "chrom": "chr2",
        "start": 110_123_348,
        "end": 110_205_013,
        "role": (
            "Nephrocystin-1 -- ciliary/centrosomal protein involved in cell "
            "signaling and cell-cell/cell-matrix adhesion.  Biallelic loss "
            "(most commonly homozygous ~290 kb deletion) causes "
            "nephronophthisis type 1 -- autosomal recessive cystic kidney "
            "disease progressing to end-stage renal disease (ESRD) by ~13 "
            "years.  Most common genetic cause of ESRD in children."
        ),
        "strategy": (
            "In vivo kidney-directed AAV gene replacement (NPHP1 cDNA ~2.2 kb).  "
            "Major challenges: targeting sufficient renal tubular cells for "
            "functional rescue, and the common cause being a large deletion "
            "(not amenable to base editing).  Full gene replacement or large "
            "fragment knock-in via CRISPR-mediated insertion required.  Renal "
            "transplant is current definitive treatment.  Gene therapy could "
            "delay/prevent ESRD if delivered early."
        ),
        "conditions": [
            "nephronophthisis",
            "ciliopathy",
            "cystic_kidney_disease",
            "end_stage_renal_disease",
        ],
    },

    # ========================================================================
    # 16.  RASOPATHIES
    # ========================================================================

    "PTPN11": {
        "chrom": "chr12",
        "start": 112_418_947,
        "end": 112_509_918,
        "role": (
            "Tyrosine-protein phosphatase non-receptor type 11 (SHP-2).  "
            "Gain-of-function mutations cause Noonan syndrome type 1 (~50% of "
            "all Noonan syndrome) -- short stature, congenital heart defects "
            "(pulmonary stenosis, HCM), characteristic facial features, "
            "bleeding diathesis, and variable intellectual disability.  Also "
            "causes LEOPARD syndrome and juvenile myelomonocytic leukemia."
        ),
        "strategy": (
            "Gain-of-function disease -- allele-specific CRISPR disruption of "
            "the mutant allele, or CRISPRi to reduce expression of the "
            "activating allele.  SHP-2 inhibitors (e.g. TNO155, RMC-4630) are "
            "in oncology trials and may have application in RASopathies.  MEK "
            "inhibitors (trametinib) used off-label for severe HCM in Noonan.  "
            "Growth hormone approved for Noonan short stature.  Gene therapy "
            "approaches are early-stage."
        ),
        "conditions": [
            "noonan_syndrome",
            "rasopathy",
            "congenital_heart_defect",
            "LEOPARD_syndrome",
            "JMML",
        ],
    },

    "RAF1": {
        "chrom": "chr3",
        "start": 12_583_601,
        "end": 12_664_117,
        "role": (
            "RAF proto-oncogene serine/threonine-protein kinase (C-RAF).  "
            "Gain-of-function mutations cause Noonan syndrome type 5 (~5% of "
            "Noonan) -- notable for high prevalence of hypertrophic "
            "cardiomyopathy (HCM, ~80% of RAF1-Noonan patients) which can be "
            "severe and life-threatening."
        ),
        "strategy": (
            "Gain-of-function -- allele-specific CRISPR disruption or CRISPRi.  "
            "RAF inhibitors and MEK inhibitors (trametinib, selumetinib) are "
            "the pharmacological approach for severe HCM.  Cardiac-specific "
            "AAV9 delivery of CRISPRi components could attenuate mutant RAF1 "
            "in cardiomyocytes.  No gene therapy trials specific to RAF1-Noonan."
        ),
        "conditions": [
            "noonan_syndrome",
            "rasopathy",
            "hypertrophic_cardiomyopathy",
        ],
    },

    "SOS1": {
        "chrom": "chr2",
        "start": 38_981_549,
        "end": 39_124_868,
        "role": (
            "Son of Sevenless homolog 1 -- RAS guanine nucleotide exchange "
            "factor (GEF) that activates RAS by promoting GDP-to-GTP exchange.  "
            "Gain-of-function mutations cause Noonan syndrome type 4 (~10% of "
            "Noonan) -- typically milder phenotype with prominent ectodermal "
            "features (keratosis pilaris, sparse eyebrows, curly hair)."
        ),
        "strategy": (
            "Gain-of-function -- allele-specific CRISPR disruption or CRISPRi.  "
            "SOS1 inhibitors (e.g. BI-3406) are in early oncology development.  "
            "The generally milder phenotype of SOS1-Noonan makes gene therapy "
            "risk-benefit less compelling.  Growth hormone for short stature, "
            "cardiac surveillance, developmental support are standard management."
        ),
        "conditions": [
            "noonan_syndrome",
            "rasopathy",
            "ectodermal_features",
        ],
    },

    "HRAS": {
        "chrom": "chr11",
        "start": 532_242,
        "end": 535_576,
        "role": (
            "HRas proto-oncogene GTPase -- small GTPase in the RAS-MAPK "
            "signaling pathway.  Gain-of-function mutations (predominantly "
            "p.Gly12Ser) cause Costello syndrome -- intellectual disability, "
            "distinctive facial features, loose skin, cardiomyopathy, "
            "predisposition to embryonal rhabdomyosarcoma and bladder carcinoma."
        ),
        "strategy": (
            "Gain-of-function -- allele-specific CRISPR disruption of the "
            "mutant HRAS allele.  The recurrent p.Gly12Ser mutation is a "
            "defined target for allele-specific guide RNA design.  MEK "
            "inhibitors (trametinib) used compassionately for severe HCM.  "
            "Tumor surveillance is essential (annual screening).  Farnesyl-"
            "transferase inhibitors (tipifarnib, lonafarnib) trialed for "
            "RASopathies with limited success.  Gene therapy approaches are "
            "preclinical."
        ),
        "conditions": [
            "costello_syndrome",
            "rasopathy",
            "cardiomyopathy",
            "cancer_predisposition",
            "rhabdomyosarcoma",
        ],
    },
}
