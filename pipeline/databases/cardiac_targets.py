"""
Inherited cardiac and cardiovascular disease gene targets — GRCh38 coordinates.

Covers cardiomyopathies (HCM, DCM, ACM), channelopathies (LQTS, Brugada,
CPVT), and connective tissue aortopathies (Marfan, Loeys-Dietz, vEDS).

All coordinates are GRCh38/hg38 from NCBI Gene / Ensembl.

Sources:
  - ClinVar / ClinGen Cardiomyopathy Expert Panel
  - Tenaya Therapeutics TN-201 (MYBPC3), TN-401 (PKP2) clinical programs
  - Verve Therapeutics VERVE-101/102 (PCSK9 base editing)
  - Intellia NTLA-2001 (TTR knockout for ATTR-CM)

All interventions require informed patient consent and IRB/ethics approval.
"""

ALL_CARDIAC_TARGETS = {
    # =========================================================================
    # HYPERTROPHIC CARDIOMYOPATHY (HCM)
    # =========================================================================
    "MYH7": {
        "chrom": "chr14",
        "start": 23412739,
        "end": 23435717,
        "role": "Beta-myosin heavy chain — sarcomeric motor protein, "
                "primary force generator in cardiomyocytes",
        "strategy": "Allele-specific silencing of dominant-negative mutant "
                     "allele; base editing for recurrent missense (R403Q, "
                     "R719W, R663H); suppression-and-replacement approach",
        "conditions": ["HCM", "hypertrophic_cardiomyopathy", "cardiomyopathy"],
        "inheritance": "AD",
        "key_variants": ["rs121913625_R403Q", "rs121913626_R719W",
                         "rs121913627_R663H"],
        "clinical_programs": "Preclinical CRISPR; mavacamten (CAMZYOS) "
                             "approved as myosin inhibitor",
    },
    "MYBPC3": {
        "chrom": "chr11",
        "start": 47331394,
        "end": 47352694,
        "role": "Myosin-binding protein C3 — sarcomeric structural protein, "
                "most common HCM gene (haploinsufficiency)",
        "strategy": "AAV9 gene replacement (cDNA ~3.8 kb fits AAV); "
                     "Tenaya TN-201 Phase 1/2 (NCT05836259); exon skipping "
                     "for specific truncating variants",
        "conditions": ["HCM", "hypertrophic_cardiomyopathy", "cardiomyopathy"],
        "inheritance": "AD",
        "key_variants": ["rs397516037_W792fs", "rs121912508_R943X",
                         "rs115345946_R502W"],
        "clinical_programs": "TN-201 Phase 1/2; Lexeo LX2006 preclinical",
    },

    # =========================================================================
    # DILATED CARDIOMYOPATHY (DCM)
    # =========================================================================
    "TTN": {
        "chrom": "chr2",
        "start": 178525989,
        "end": 178807423,
        "role": "Titin — largest human protein, sarcomeric spring/scaffold; "
                "TTN truncating variants cause ~25% of familial DCM",
        "strategy": "Exon skipping for truncating variants (analogous to DMD "
                     "approach); gene too large (~100 kb cDNA) for AAV "
                     "replacement; base/prime editing for specific variants",
        "conditions": ["DCM", "dilated_cardiomyopathy", "cardiomyopathy"],
        "inheritance": "AD",
        "key_variants": ["TTNtv_in_A-band_exons"],
        "clinical_programs": "Early preclinical",
    },
    "LMNA": {
        "chrom": "chr1",
        "start": 156104878,
        "end": 156134876,
        "role": "Lamin A/C — nuclear lamina structural protein; mutations "
                "cause DCM with conduction disease and sudden death",
        "strategy": "Allele-specific silencing of dominant-negative allele; "
                     "gene replacement (cDNA ~2.0 kb fits AAV); CRISPRi "
                     "for mutant allele transcriptional repression",
        "conditions": ["DCM", "dilated_cardiomyopathy", "cardiomyopathy",
                        "laminopathy"],
        "inheritance": "AD",
        "key_variants": ["rs57045855_R453W", "rs121912485_R225X",
                         "rs28928900_R249Q"],
        "clinical_programs": "Preclinical; managed with ICD + HF therapy",
    },

    # =========================================================================
    # ARRHYTHMOGENIC CARDIOMYOPATHY (ACM)
    # =========================================================================
    "PKP2": {
        "chrom": "chr12",
        "start": 32829663,
        "end": 32920279,
        "role": "Plakophilin-2 — desmosomal cell junction protein; most "
                "common ACM gene (haploinsufficiency)",
        "strategy": "AAV9 gene replacement (cDNA ~2.6 kb fits AAV); "
                     "Tenaya TN-401 and Rocket RP-A601 in development",
        "conditions": ["ACM", "arrhythmogenic_cardiomyopathy", "ARVC",
                        "cardiomyopathy"],
        "inheritance": "AD",
        "key_variants": ["rs121434420_C796R", "rs121434422_Q49fs"],
        "clinical_programs": "TN-401 Phase 1 planned; RP-A601 IND-enabling",
    },
    "DSP": {
        "chrom": "chr6",
        "start": 7541868,
        "end": 7586947,
        "role": "Desmoplakin — desmosomal plaque protein linking desmosomes "
                "to intermediate filaments",
        "strategy": "Dual AAV split-intein (cDNA ~8.6 kb exceeds AAV limit); "
                     "exon skipping for specific truncating variants; "
                     "base/prime editing for point mutations",
        "conditions": ["ACM", "arrhythmogenic_cardiomyopathy",
                        "cardiomyopathy"],
        "inheritance": "AD",
        "key_variants": ["rs121912461_R1267X"],
        "clinical_programs": "Preclinical",
    },
    "DSG2": {
        "chrom": "chr18",
        "start": 31494742,
        "end": 31527610,
        "role": "Desmoglein-2 — desmosomal cadherin, cell-cell adhesion",
        "strategy": "AAV gene replacement (cDNA ~3.4 kb fits AAV)",
        "conditions": ["ACM", "arrhythmogenic_cardiomyopathy",
                        "cardiomyopathy"],
        "inheritance": "AD/AR",
        "key_variants": ["rs121908002_splice", "rs121908003_F531C"],
        "clinical_programs": "Preclinical",
    },

    # =========================================================================
    # CHANNELOPATHIES — LONG QT SYNDROME
    # =========================================================================
    "KCNQ1": {
        "chrom": "chr11",
        "start": 2444990,
        "end": 2849109,
        "role": "Potassium voltage-gated channel (IKs) — cardiac "
                "repolarization; LQT1 (most common LQTS type, ~35%)",
        "strategy": "AAV gene replacement (cDNA ~2.0 kb fits AAV); CRISPRa "
                     "to upregulate wild-type allele; suppression-and-"
                     "replacement for dominant-negative mutations",
        "conditions": ["LQTS", "LQT1", "long_QT_syndrome", "channelopathy",
                        "arrhythmia"],
        "inheritance": "AD",
        "key_variants": ["rs120074175_A341V"],
        "clinical_programs": "Preclinical; beta-blockers standard of care",
    },
    "KCNH2": {
        "chrom": "chr7",
        "start": 150642049,
        "end": 150675403,
        "role": "Potassium voltage-gated channel (IKr/hERG) — cardiac "
                "repolarization; LQT2 (~30% of LQTS)",
        "strategy": "AAV gene replacement (cDNA ~3.5 kb fits AAV); "
                     "suppression-and-replacement for dominant-negative; "
                     "pharmacological chaperones as alternative",
        "conditions": ["LQTS", "LQT2", "long_QT_syndrome", "channelopathy",
                        "arrhythmia"],
        "inheritance": "AD",
        "key_variants": ["rs121912775_A561T", "rs36210422_G628S"],
        "clinical_programs": "Preclinical iPSC-CM correction studies",
    },
    "SCN5A": {
        "chrom": "chr3",
        "start": 38548062,
        "end": 38649667,
        "role": "Cardiac sodium channel Nav1.5 — depolarization/conduction; "
                "LQT3 (gain-of-function), Brugada (loss-of-function), DCM",
        "strategy": "Allele-specific disruption of gain-of-function (LQT3); "
                     "gene replacement borderline AAV (~6 kb cDNA); "
                     "base editing for specific mutations; mexiletine for LQT3",
        "conditions": ["LQTS", "LQT3", "brugada_syndrome", "channelopathy",
                        "arrhythmia", "DCM"],
        "inheritance": "AD",
        "key_variants": ["rs121918852_delKPQ", "rs121918860_G717R",
                         "rs121918851_E1784K"],
        "clinical_programs": "Preclinical",
    },

    # =========================================================================
    # CATECHOLAMINERGIC POLYMORPHIC VT (CPVT)
    # =========================================================================
    "RYR2": {
        "chrom": "chr1",
        "start": 237042823,
        "end": 237833982,
        "role": "Ryanodine receptor 2 — SR calcium release channel; "
                "gain-of-function mutations cause exercise-triggered VT",
        "strategy": "Allele-specific silencing of gain-of-function allele; "
                     "gene too large (~15 kb cDNA) for AAV; CRISPRi for "
                     "mutant allele; base editing for hotspot mutations",
        "conditions": ["CPVT", "catecholaminergic_polymorphic_VT",
                        "channelopathy", "arrhythmia"],
        "inheritance": "AD",
        "key_variants": ["rs121918589_R4159Q", "rs121918587_S2246L"],
        "clinical_programs": "Preclinical; flecainide + beta-blockers SOC",
    },
    "CASQ2": {
        "chrom": "chr1",
        "start": 115700422,
        "end": 115768963,
        "role": "Calsequestrin 2 — SR calcium buffer; loss-of-function "
                "causes autosomal recessive CPVT",
        "strategy": "AAV9 gene replacement (cDNA ~1.2 kb, ideal AAV target); "
                     "preclinical correction in CASQ2-KO mice normalizes "
                     "arrhythmia",
        "conditions": ["CPVT", "catecholaminergic_polymorphic_VT",
                        "channelopathy", "arrhythmia"],
        "inheritance": "AR",
        "key_variants": ["rs121918598_R33X", "rs121918600_D338G"],
        "clinical_programs": "Preclinical with promising animal data",
    },

    # =========================================================================
    # AORTOPATHIES / CONNECTIVE TISSUE
    # =========================================================================
    "FBN1": {
        "chrom": "chr15",
        "start": 48408306,
        "end": 48645709,
        "role": "Fibrillin-1 — extracellular matrix glycoprotein; mutations "
                "cause Marfan syndrome (aortic aneurysm, lens subluxation)",
        "strategy": "Allele-specific silencing of dominant-negative allele; "
                     "dual AAV gene replacement (cDNA ~8.6 kb exceeds AAV); "
                     "TGF-beta pathway modulation (losartan); >1800 variants "
                     "identified, extreme allelic heterogeneity",
        "conditions": ["marfan_syndrome", "aortopathy",
                        "connective_tissue_disorder"],
        "inheritance": "AD",
        "key_variants": ["various_Cys_substitutions_in_cbEGF_domains"],
        "clinical_programs": "Preclinical; losartan/ARBs standard of care",
    },
    "COL3A1": {
        "chrom": "chr2",
        "start": 188974367,
        "end": 189012746,
        "role": "Type III collagen — vascular/organ structural collagen; "
                "mutations cause vascular Ehlers-Danlos (vEDS, life-"
                "threatening arterial rupture)",
        "strategy": "Allele-specific silencing of dominant-negative allele "
                     "(most mutations are Gly->X in triple helix domain); "
                     "suppression-and-replacement approach; AAV possible "
                     "(cDNA ~4.4 kb)",
        "conditions": ["vascular_ehlers_danlos", "vEDS",
                        "connective_tissue_disorder", "aortopathy"],
        "inheritance": "AD",
        "key_variants": ["Gly_substitutions_in_Gly-X-Y_repeat"],
        "clinical_programs": "Preclinical; celiprolol only specific pharma",
    },
    "TGFBR1": {
        "chrom": "chr9",
        "start": 99104038,
        "end": 99154192,
        "role": "TGF-beta receptor type 1 — TGF-beta signaling; mutations "
                "cause Loeys-Dietz syndrome (aggressive aortic aneurysm)",
        "strategy": "AAV gene replacement (cDNA ~1.5 kb, ideal AAV target); "
                     "allele-specific silencing; TGF-beta pathway modulation "
                     "(losartan/ARBs)",
        "conditions": ["loeys_dietz_syndrome", "aortopathy",
                        "connective_tissue_disorder"],
        "inheritance": "AD",
        "key_variants": ["rs121434558_R487W"],
        "clinical_programs": "Preclinical; losartan + surgery SOC",
    },
    "TGFBR2": {
        "chrom": "chr3",
        "start": 30606381,
        "end": 30694142,
        "role": "TGF-beta receptor type 2 — TGF-beta signaling; Loeys-Dietz "
                "syndrome type 2",
        "strategy": "AAV gene replacement (cDNA ~1.7 kb fits AAV); "
                     "allele-specific silencing",
        "conditions": ["loeys_dietz_syndrome", "aortopathy",
                        "connective_tissue_disorder"],
        "inheritance": "AD",
        "key_variants": ["rs121434541_R528H", "rs121434540_R528C"],
        "clinical_programs": "Preclinical",
    },
}
