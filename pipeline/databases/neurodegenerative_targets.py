"""
Neurodegenerative disease gene-therapy / CRISPR targets -- GRCh38 (hg38) coordinates.

Covers Alzheimer's disease, Parkinson's disease, frontotemporal dementia (FTD),
and amyotrophic lateral sclerosis (ALS) targets not already present in the
hereditary disease database (which covers HTT/Huntington, SMN1/SMA, SOD1/ALS).

Categories:
  1. Alzheimer's Disease (APP, PSEN1, PSEN2, APOE, TREM2)
  2. Parkinson's Disease (SNCA, LRRK2, GBA1, PRKN, PINK1, PARK7)
  3. Frontotemporal Dementia (MAPT, GRN, C9orf72)
  4. ALS — additional targets (FUS, TARDBP, ATXN2)

Coordinates are from NCBI Gene, GRCh38.p14 (RefSeq annotation RS_2025_08).
Variant positions from NCBI dbSNP / ClinVar.

IMPORTANT -- verify every coordinate against current NCBI / Ensembl releases
before production use.  Numbering can shift between patch levels.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14)
  - OMIM
  - ClinVar / dbSNP
  - ClinicalTrials.gov
  - Alzheimer's Disease Neuroimaging Initiative (ADNI)
  - Michael J. Fox Foundation Parkinson's gene therapy tracker
  - Published literature through early 2026
"""


# ============================================================================
# 1. ALZHEIMER'S DISEASE
# ============================================================================

ALZHEIMERS_TARGETS = {

    # -------------------------------------------------------------------
    # 1a.  APP -- Amyloid Precursor Protein
    # -------------------------------------------------------------------
    "APP": {
        "gene_id": 351,
        "chrom": "chr21",
        "start": 25_880_550,
        "end": 26_170_723,
        "strand": "+",
        "refseq": "NC_000021.9",
        "cytoband": "21q21.3",
        "exon_count": 18,
        "role": (
            "Amyloid precursor protein -- type I transmembrane glycoprotein.  "
            "Sequential cleavage by beta-secretase (BACE1) and gamma-secretase "
            "produces amyloid-beta (Abeta) peptides.  Accumulation of Abeta42 "
            "triggers amyloid plaques, neuroinflammation, and neuronal death.  "
            "Autosomal dominant mutations (>50 known) cluster near the "
            "secretase cleavage sites and increase Abeta42/40 ratio."
        ),
        "disease": "Early-onset familial Alzheimer's disease (AD1)",
        "omim_disease": 104300,
        "omim_gene": 104760,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Val717Ile (London mutation)",
                "hgvs_coding": "NM_000484.4:c.2149G>A",
                "rsid": "rs63750264",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Gamma-secretase site; increases Abeta42/40 ratio.",
            },
            {
                "name": "p.Lys670Asn/Met671Leu (Swedish mutation)",
                "rsid": "rs63750671",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Beta-secretase site; increases total Abeta production "
                    "6-8 fold.  Basis for many AD mouse models."
                ),
            },
            {
                "name": "p.Ala673Thr (Icelandic protective variant)",
                "rsid": "rs63750847",
                "consequence": "missense",
                "clinical_significance": "protective",
                "notes": (
                    "Reduces Abeta production ~40%.  Heterozygous carriers "
                    "have lower AD risk and slower cognitive decline.  "
                    "Prime editing to install A673T is a therapeutic strategy."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 50 pathogenic missense mutations in ClinVar, nearly all in "
            "exons 16-17 near secretase cleavage sites.  Duplications of the "
            "APP locus (including in Down syndrome / trisomy 21) also cause "
            "early-onset AD.  The A673T protective variant is an active "
            "therapeutic target."
        ),
        "strategy": (
            "(1) Prime editing to install the protective A673T variant, reducing "
            "Abeta production ~40% -- most promising single-gene approach; "
            "(2) CRISPRi to reduce APP expression in neurons (brain-penetrant "
            "AAV9 or intrathecal delivery); "
            "(3) Allele-specific silencing of mutant APP allele for familial AD; "
            "(4) Base editing to correct specific pathogenic missense mutations.  "
            "Delivery: AAV9 or AAVrh10 with neuron-specific promoter (SYN1)."
        ),
        "clinical_programs": (
            "No CRISPR clinical trials for APP as of 2026.  "
            "Lecanemab (anti-Abeta antibody) FDA-approved 2023.  "
            "Donanemab FDA-approved 2024.  ASO approaches (BIIB080/IONIS) "
            "in Phase 1 for tau."
        ),
        "conditions": ["alzheimers_disease", "early_onset_AD", "familial_AD",
                        "neurodegenerative", "dementia"],
    },

    # -------------------------------------------------------------------
    # 1b.  PSEN1 -- Presenilin 1
    # -------------------------------------------------------------------
    "PSEN1": {
        "gene_id": 5663,
        "chrom": "chr14",
        "start": 73_136_417,
        "end": 73_223_691,
        "strand": "-",
        "refseq": "NC_000014.9",
        "cytoband": "14q24.2",
        "exon_count": 12,
        "role": (
            "Presenilin 1 -- catalytic subunit of gamma-secretase complex.  "
            "Cleaves APP to generate Abeta peptides.  Most common cause of "
            "early-onset familial Alzheimer's disease (AD3).  Over 300 "
            "pathogenic mutations identified, nearly all missense, causing "
            "onset as early as age 25-30."
        ),
        "disease": "Early-onset familial Alzheimer's disease (AD3)",
        "omim_disease": 607822,
        "omim_gene": 104311,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Ala246Glu (A246E)",
                "rsid": "rs63749810",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "TM6 domain; onset ~55 years; well-characterized.",
            },
            {
                "name": "p.Glu280Ala (E280A, Paisa mutation)",
                "rsid": "rs63750082",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Founder mutation in Colombian kindred (~5000 carriers).  "
                    "Onset ~44 years.  Subject of Alzheimer's Prevention "
                    "Initiative (API) trials."
                ),
            },
            {
                "name": "p.Met146Leu (M146L)",
                "rsid": "rs63750231",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "TM2 domain; onset ~40 years; one of the earliest described.",
            },
        ],
        "mutational_landscape": (
            ">300 pathogenic variants, more than any other AD gene.  ~95% "
            "missense.  Mutations alter gamma-secretase cleavage to increase "
            "Abeta42/40 ratio.  Extreme allelic heterogeneity makes allele-"
            "specific approaches challenging."
        ),
        "strategy": (
            "(1) Allele-specific CRISPRi/Cas13 knockdown of mutant PSEN1 mRNA "
            "while preserving wild-type allele; "
            "(2) Base editing for recurrent mutations (e.g., E280A in the "
            "Colombian kindred); "
            "(3) HDR correction of specific pathogenic variants in patient "
            "iPSC-derived neurons for autologous transplant; "
            "(4) CRISPRa of wild-type allele.  "
            "Delivery: AAV9 intrathecal or AAVrh10 systemic with CNS tropism."
        ),
        "clinical_programs": (
            "No CRISPR trials for PSEN1 as of 2026.  API Colombia trial "
            "testing crenezumab in E280A carriers.  Banner Alzheimer's "
            "Prevention Initiative active."
        ),
        "conditions": ["alzheimers_disease", "early_onset_AD", "familial_AD",
                        "neurodegenerative", "dementia"],
    },

    # -------------------------------------------------------------------
    # 1c.  PSEN2 -- Presenilin 2
    # -------------------------------------------------------------------
    "PSEN2": {
        "gene_id": 5664,
        "chrom": "chr1",
        "start": 226_870_594,
        "end": 226_903_828,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1q42.13",
        "exon_count": 12,
        "role": (
            "Presenilin 2 -- gamma-secretase catalytic subunit homolog.  "
            "Mutations cause AD4 (early-onset familial Alzheimer's, rare).  "
            "Fewer pathogenic variants than PSEN1 but same mechanism."
        ),
        "disease": "Early-onset familial Alzheimer's disease (AD4)",
        "omim_disease": 606889,
        "omim_gene": 600759,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Asn141Ile (N141I, Volga German mutation)",
                "rsid": "rs63750215",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Founder mutation in Volga German families; onset ~52 years.",
            },
            {
                "name": "p.Met239Val (M239V)",
                "rsid": "rs63750399",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Italian kindred; onset ~50 years.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant PSEN2 allele; "
            "(2) Base editing for N141I founder mutation; "
            "(3) HDR correction in iPSC-derived neurons.  "
            "Same delivery considerations as PSEN1."
        ),
        "clinical_programs": "No CRISPR/gene therapy trials as of 2026.",
        "conditions": ["alzheimers_disease", "early_onset_AD", "familial_AD",
                        "neurodegenerative", "dementia"],
    },

    # -------------------------------------------------------------------
    # 1d.  APOE -- Apolipoprotein E
    # -------------------------------------------------------------------
    "APOE": {
        "gene_id": 348,
        "chrom": "chr19",
        "start": 44_905_754,
        "end": 44_909_393,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19q13.32",
        "exon_count": 4,
        "role": (
            "Apolipoprotein E -- lipid transport protein, major genetic risk "
            "factor for late-onset Alzheimer's.  APOE4 (rs429358-C, rs7412-C) "
            "increases AD risk 3-15x.  APOE2 is protective.  APOE mediates "
            "Abeta clearance, neuroinflammation, and BBB integrity.  "
            "APOE4 homozygotes: ~60% lifetime AD risk."
        ),
        "disease": "Late-onset Alzheimer's disease (major risk factor)",
        "omim_disease": 104310,
        "omim_gene": 107741,
        "inheritance": "complex/risk_factor",
        "key_variants": [
            {
                "name": "APOE4 (epsilon 4)",
                "rsid": "rs429358",
                "consequence": "missense (Cys112Arg)",
                "clinical_significance": "risk_factor",
                "notes": (
                    "APOE4 allele: Arg112/Arg158.  Homozygotes have ~12x AD "
                    "risk.  Present in ~25% of population, ~40% of AD patients."
                ),
            },
            {
                "name": "APOE2 (epsilon 2, protective)",
                "rsid": "rs7412",
                "consequence": "missense (Arg158Cys)",
                "clinical_significance": "protective",
                "notes": "APOE2 allele: Cys112/Cys158.  Reduces AD risk ~40%.",
            },
            {
                "name": "APOE Christchurch (R136S)",
                "rsid": "rs121918393",
                "consequence": "missense",
                "clinical_significance": "protective",
                "notes": (
                    "Identified in a PSEN1 E280A carrier who remained cognitively "
                    "normal until age 70+.  Reduces APOE-HSPG binding.  "
                    "Prime editing target for neuroprotection."
                ),
            },
        ],
        "strategy": (
            "(1) Base editing to convert APOE4 -> APOE3 (Arg112Cys single "
            "nucleotide change) or APOE4 -> APOE2 -- most impactful single "
            "edit for AD prevention; "
            "(2) Prime editing to install Christchurch variant (R136S) as "
            "neuroprotective modification; "
            "(3) CRISPRi to reduce APOE expression in astrocytes (lower "
            "APOE4 = less amyloid seeding); "
            "(4) AAV-delivered APOE2 overexpression (Lexeo LX1001 approach).  "
            "Delivery: AAV9 intracisternal; AAVrh10 IV."
        ),
        "clinical_programs": (
            "Lexeo Therapeutics LX1001 -- AAV-APOE2 intracisternal injection "
            "for APOE4 homozygotes, Phase 1/2 (NCT05400330).  "
            "Verve Therapeutics exploring APOE base editing (preclinical).  "
            "No CRISPR clinical trials as of 2026."
        ),
        "conditions": ["alzheimers_disease", "late_onset_AD",
                        "neurodegenerative", "dementia", "APOE4"],
    },

    # -------------------------------------------------------------------
    # 1e.  TREM2 -- Triggering Receptor Expressed on Myeloid Cells 2
    # -------------------------------------------------------------------
    "TREM2": {
        "gene_id": 54209,
        "chrom": "chr6",
        "start": 41_126_243,
        "end": 41_130_923,
        "strand": "-",
        "refseq": "NC_000006.12",
        "cytoband": "6p21.1",
        "exon_count": 5,
        "role": (
            "TREM2 -- microglial receptor, regulates phagocytosis of Abeta, "
            "synaptic pruning, and neuroinflammation.  R47H variant triples "
            "AD risk (comparable to one APOE4 allele).  Loss-of-function "
            "mutations cause Nasu-Hakola disease (polycystic lipomembranous "
            "osteodysplasia with sclerosing leukoencephalopathy)."
        ),
        "disease": "Alzheimer's disease risk / Nasu-Hakola disease",
        "omim_disease": 605514,
        "omim_gene": 605086,
        "inheritance": "complex (AD risk); AR (Nasu-Hakola)",
        "key_variants": [
            {
                "name": "p.Arg47His (R47H)",
                "rsid": "rs75932628",
                "consequence": "missense",
                "clinical_significance": "risk_factor",
                "notes": "~3x AD risk; impairs lipid ligand binding.",
            },
            {
                "name": "p.Arg62His (R62H)",
                "rsid": "rs143332484",
                "consequence": "missense",
                "clinical_significance": "risk_factor",
                "notes": "Modest AD risk increase; more common than R47H.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to boost TREM2 expression in microglia, enhancing "
            "Abeta clearance; "
            "(2) Base editing to correct R47H in patients carrying that variant; "
            "(3) AAV-mediated TREM2 overexpression (cDNA ~0.7 kb, ideal for AAV).  "
            "Delivery: AAV9 intrathecal; microglia-tropic AAV capsids in "
            "development (AAV.MG)."
        ),
        "clinical_programs": (
            "Vigil Neuroscience anti-TREM2 agonist antibody (VGL101) Phase 1.  "
            "No CRISPR trials as of 2026.  Alector AL002 TREM2 agonist "
            "discontinued after Phase 2 futility."
        ),
        "conditions": ["alzheimers_disease", "neurodegenerative", "dementia",
                        "nasu_hakola_disease", "neuroinflammation"],
    },
}


# ============================================================================
# 2. PARKINSON'S DISEASE
# ============================================================================

PARKINSONS_TARGETS = {

    # -------------------------------------------------------------------
    # 2a.  SNCA -- Alpha-synuclein
    # -------------------------------------------------------------------
    "SNCA": {
        "gene_id": 6622,
        "chrom": "chr4",
        "start": 89_724_099,
        "end": 89_838_315,
        "strand": "-",
        "refseq": "NC_000004.12",
        "cytoband": "4q22.1",
        "exon_count": 6,
        "role": (
            "Alpha-synuclein -- small presynaptic protein, key component of "
            "Lewy bodies.  Missense mutations (A53T, A30P, E46K) and gene "
            "multiplications cause autosomal dominant Parkinson's disease "
            "(PARK1/PARK4).  Alpha-synuclein aggregation is central to PD "
            "pathology and a major therapeutic target."
        ),
        "disease": "Parkinson's disease (PARK1/PARK4)",
        "omim_disease": 168601,
        "omim_gene": 163890,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Ala53Thr (A53T)",
                "rsid": "rs104893877",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "First PD gene identified (Contursi kindred); aggressive course.",
            },
            {
                "name": "p.Ala30Pro (A30P)",
                "rsid": "rs104893878",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "German family; alters membrane binding.",
            },
            {
                "name": "SNCA duplication/triplication",
                "rsid": None,
                "consequence": "copy_number_gain",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Gene dose correlates with severity: duplication = typical PD, "
                    "triplication = early-onset aggressive PD with dementia."
                ),
            },
        ],
        "strategy": (
            "(1) CRISPRi to reduce SNCA expression (lower alpha-synuclein = "
            "less aggregation) -- especially for duplications/triplications; "
            "(2) Allele-specific silencing of mutant allele for point mutations; "
            "(3) CRISPR-mediated deletion of SNCA enhancer elements to reduce "
            "expression by ~50% (mimicking heterozygous knockout, which is "
            "protective in mouse models); "
            "(4) ASO/siRNA approaches as complementary strategy.  "
            "Delivery: AAV9 or AAVrh10 to substantia nigra (stereotactic) "
            "or intrathecal."
        ),
        "clinical_programs": (
            "Alnylam/Regeneron ALN-APP (ASO for SNCA) Phase 1.  "
            "Prevail Therapeutics (Lilly) PR001 AAV-GBA1 for GBA-PD Phase 2.  "
            "No CRISPR trials for SNCA as of 2026."
        ),
        "conditions": ["parkinsons_disease", "parkinsonism", "PARK1",
                        "synucleinopathy", "neurodegenerative", "lewy_body"],
    },

    # -------------------------------------------------------------------
    # 2b.  LRRK2 -- Leucine-Rich Repeat Kinase 2
    # -------------------------------------------------------------------
    "LRRK2": {
        "gene_id": 120892,
        "chrom": "chr12",
        "start": 40_618_813,
        "end": 40_763_087,
        "strand": "+",
        "refseq": "NC_000012.12",
        "cytoband": "12q12",
        "exon_count": 51,
        "role": (
            "LRRK2 (dardarin) -- large multidomain kinase/GTPase.  G2019S "
            "gain-of-function mutation is the most common genetic cause of PD "
            "(~1-2% of sporadic PD, up to 40% in Ashkenazi Jewish and North "
            "African Berber populations).  Kinase hyperactivity drives "
            "neurodegeneration."
        ),
        "disease": "Parkinson's disease (PARK8)",
        "omim_disease": 607060,
        "omim_gene": 609007,
        "inheritance": "AD (incomplete penetrance ~30-75%)",
        "key_variants": [
            {
                "name": "p.Gly2019Ser (G2019S)",
                "rsid": "rs34637584",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common PD mutation worldwide.  Kinase domain; "
                    "increases kinase activity 2-3 fold.  Target of LRRK2 "
                    "kinase inhibitors (DNL151/BIIB122)."
                ),
            },
            {
                "name": "p.Arg1441Cys (R1441C)",
                "rsid": "rs33939927",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "ROC-GTPase domain; Basque founder mutation.",
            },
            {
                "name": "p.Arg1441Gly (R1441G)",
                "rsid": "rs33939927",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common LRRK2 mutation in Basque population.",
            },
        ],
        "strategy": (
            "(1) CRISPRi to reduce LRRK2 kinase activity by lowering "
            "expression (partial knockdown safe based on primate data); "
            "(2) Base editing to correct G2019S (single G>A change); "
            "(3) Allele-specific Cas13 knockdown of mutant mRNA; "
            "(4) CRISPR disruption of kinase domain in mutant allele.  "
            "Delivery: AAV9 stereotactic to striatum/substantia nigra."
        ),
        "clinical_programs": (
            "Biogen BIIB122 (LRRK2 kinase inhibitor) Phase 2b.  "
            "Denali DNL151 LRRK2 inhibitor Phase 1b.  "
            "No CRISPR trials for LRRK2 as of 2026.  "
            "Prevail PR004 AAV-based LRRK2 reduction preclinical."
        ),
        "conditions": ["parkinsons_disease", "parkinsonism", "PARK8",
                        "LRRK2_PD", "neurodegenerative"],
    },

    # -------------------------------------------------------------------
    # 2c.  GBA1 -- Glucocerebrosidase
    # -------------------------------------------------------------------
    "GBA1": {
        "gene_id": 2629,
        "chrom": "chr1",
        "start": 155_234_452,
        "end": 155_244_627,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1q22",
        "exon_count": 11,
        "role": (
            "Glucocerebrosidase -- lysosomal enzyme cleaving glucosylceramide.  "
            "Biallelic loss-of-function = Gaucher disease (most common LSD).  "
            "Heterozygous GBA1 variants are the most common genetic risk factor "
            "for Parkinson's disease (~5-10% of PD patients carry GBA1 variants).  "
            "GBA1 deficiency impairs lysosomal alpha-synuclein clearance."
        ),
        "disease": "Parkinson's disease (GBA-PD) / Gaucher disease",
        "omim_disease": 168600,
        "omim_gene": 606463,
        "inheritance": "AD risk factor (PD); AR (Gaucher)",
        "key_variants": [
            {
                "name": "p.Asn370Ser (N370S)",
                "rsid": "rs76763715",
                "consequence": "missense",
                "clinical_significance": "pathogenic/risk_factor",
                "notes": (
                    "Most common Gaucher mutation (Type 1).  PD risk ~5x.  "
                    "Residual enzyme activity; non-neuropathic Gaucher."
                ),
            },
            {
                "name": "p.Leu444Pro (L444P)",
                "rsid": "rs421016",
                "consequence": "missense",
                "clinical_significance": "pathogenic/risk_factor",
                "notes": (
                    "Severe Gaucher mutation (Type 2/3 when homozygous).  "
                    "Higher PD risk than N370S."
                ),
            },
            {
                "name": "p.Glu326Lys (E326K)",
                "rsid": "rs2230288",
                "consequence": "missense",
                "clinical_significance": "risk_factor",
                "notes": "Mild PD risk factor; does not cause Gaucher disease alone.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated GBA1 gene replacement (cDNA ~1.6 kb, ideal for "
            "AAV) -- Prevail PR001 approach; "
            "(2) Base editing to correct specific variants (N370S, L444P); "
            "(3) CRISPRa to upregulate wild-type GBA1 expression; "
            "(4) Small molecule chaperones (ambroxol) as adjunct.  "
            "Delivery: AAV9 intracisternal magna (Prevail approach) or "
            "AAV9 IV with CNS tropism."
        ),
        "clinical_programs": (
            "Prevail Therapeutics (Lilly) PR001 -- AAV9-GBA1 intracisternal, "
            "Phase 1/2 for GBA-PD (NCT04127578) and Gaucher Type 2 "
            "(NCT04411654).  Passage Bio PBGM01 for neuropathic Gaucher.  "
            "No CRISPR trials for GBA1 as of 2026."
        ),
        "conditions": ["parkinsons_disease", "GBA_PD", "gaucher_disease",
                        "neurodegenerative", "lysosomal_storage"],
    },

    # -------------------------------------------------------------------
    # 2d.  PRKN (PARK2) -- Parkin
    # -------------------------------------------------------------------
    "PRKN": {
        "gene_id": 5071,
        "chrom": "chr6",
        "start": 161_347_420,
        "end": 162_727_766,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6q26",
        "exon_count": 12,
        "role": (
            "Parkin -- E3 ubiquitin ligase, essential for mitophagy (PINK1-Parkin "
            "pathway).  Loss-of-function causes autosomal recessive juvenile "
            "Parkinson's disease (PARK2), the most common form of early-onset PD.  "
            "Parkin ubiquitinates damaged mitochondria for autophagic clearance."
        ),
        "disease": "Juvenile/early-onset Parkinson's disease (PARK2)",
        "omim_disease": 600116,
        "omim_gene": 602544,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Exon 3-4 deletion",
                "rsid": None,
                "consequence": "large_deletion",
                "clinical_significance": "pathogenic",
                "notes": "Most common PRKN mutation (~50% of PARK2 cases).",
            },
            {
                "name": "p.Arg275Trp (R275W)",
                "rsid": "rs34424986",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "RING1 domain; impairs E3 ligase activity.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated PRKN gene replacement (cDNA ~1.4 kb, ideal for "
            "AAV) under neuron-specific promoter; "
            "(2) CRISPRa to upregulate remaining functional allele in "
            "compound heterozygotes; "
            "(3) Gene therapy restores mitophagy in PRKN-null neurons.  "
            "Delivery: AAV9 or AAV2 stereotactic to substantia nigra."
        ),
        "clinical_programs": (
            "No CRISPR or gene therapy trials for PRKN as of 2026.  "
            "PARK2 patients respond well to L-DOPA with slow progression."
        ),
        "conditions": ["parkinsons_disease", "juvenile_PD", "PARK2",
                        "early_onset_PD", "neurodegenerative"],
    },

    # -------------------------------------------------------------------
    # 2e.  PINK1 -- PTEN-Induced Kinase 1
    # -------------------------------------------------------------------
    "PINK1": {
        "gene_id": 65018,
        "chrom": "chr1",
        "start": 20_633_403,
        "end": 20_651_539,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1p36.12",
        "exon_count": 8,
        "role": (
            "PINK1 -- mitochondrial serine/threonine kinase.  Senses damaged "
            "mitochondria and phosphorylates ubiquitin + Parkin to initiate "
            "mitophagy.  Loss-of-function causes autosomal recessive early-"
            "onset Parkinson's disease (PARK6)."
        ),
        "disease": "Early-onset Parkinson's disease (PARK6)",
        "omim_disease": 605909,
        "omim_gene": 608309,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gln456Ter (Q456X)",
                "rsid": "rs45478900",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Truncating mutation; kinase domain loss.",
            },
            {
                "name": "p.Gly309Asp (G309D)",
                "rsid": "rs45530340",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Kinase domain; impairs catalytic activity.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated PINK1 gene replacement (cDNA ~1.7 kb, ideal for "
            "AAV) to restore mitophagy; "
            "(2) CRISPRa to boost residual expression; "
            "(3) Base editing for specific missense variants.  "
            "Delivery: AAV9 or AAV2 stereotactic to substantia nigra."
        ),
        "clinical_programs": "No gene therapy or CRISPR trials as of 2026.",
        "conditions": ["parkinsons_disease", "early_onset_PD", "PARK6",
                        "neurodegenerative", "mitophagy_defect"],
    },

    # -------------------------------------------------------------------
    # 2f.  PARK7 (DJ-1) -- Parkinsonism-Associated Deglycase
    # -------------------------------------------------------------------
    "PARK7": {
        "gene_id": 11315,
        "chrom": "chr1",
        "start": 7_961_654,
        "end": 7_985_505,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p36.23",
        "exon_count": 8,
        "role": (
            "DJ-1 -- oxidative stress sensor and mitochondrial protector.  "
            "Loss-of-function causes rare autosomal recessive early-onset PD "
            "(PARK7).  DJ-1 protects against oxidative damage and regulates "
            "mitochondrial function."
        ),
        "disease": "Early-onset Parkinson's disease (PARK7)",
        "omim_disease": 606324,
        "omim_gene": 602533,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Leu166Pro (L166P)",
                "rsid": "rs28940891",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Destabilizes DJ-1 dimer; most characterized mutation.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~0.6 kb, excellent AAV target); "
            "(2) Base editing for L166P.  "
            "Delivery: AAV9 stereotactic or intrathecal."
        ),
        "clinical_programs": "No gene therapy or CRISPR trials as of 2026.",
        "conditions": ["parkinsons_disease", "early_onset_PD", "PARK7",
                        "neurodegenerative"],
    },
}


# ============================================================================
# 3. FRONTOTEMPORAL DEMENTIA (FTD)
# ============================================================================

FTD_TARGETS = {

    # -------------------------------------------------------------------
    # 3a.  MAPT -- Microtubule-Associated Protein Tau
    # -------------------------------------------------------------------
    "MAPT": {
        "gene_id": 4137,
        "chrom": "chr17",
        "start": 45_894_382,
        "end": 46_028_334,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.31",
        "exon_count": 16,
        "role": (
            "Tau -- microtubule-binding protein essential for axonal transport.  "
            "Mutations cause frontotemporal dementia with parkinsonism linked "
            "to chromosome 17 (FTDP-17).  Tau aggregation (neurofibrillary "
            "tangles) is central to tauopathies including AD, PSP, CBD."
        ),
        "disease": "Frontotemporal dementia (FTDP-17)",
        "omim_disease": 600274,
        "omim_gene": 157140,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Pro301Leu (P301L)",
                "rsid": "rs63751273",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common MAPT mutation for FTD.  Exon 10; increases "
                    "4R tau and aggregation propensity.  Basis for rTg4510 mouse model."
                ),
            },
            {
                "name": "IVS10+16 C>T (intronic)",
                "rsid": "rs63751011",
                "consequence": "splice_regulatory",
                "clinical_significance": "pathogenic",
                "notes": "Increases exon 10 inclusion -> excess 4R tau.",
            },
            {
                "name": "p.Val337Met (V337M)",
                "rsid": "rs63751140",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Seattle family A; disrupts microtubule binding.",
            },
        ],
        "strategy": (
            "(1) ASO-mediated MAPT reduction (BIIB080/IONIS-MAPTRx approach); "
            "(2) CRISPRi to reduce total tau expression; "
            "(3) Allele-specific silencing of mutant MAPT for familial FTD; "
            "(4) CRISPR-mediated skipping of exon 10 to reduce 4R/3R tau ratio.  "
            "Delivery: intrathecal ASO or AAV9 with neuronal promoter."
        ),
        "clinical_programs": (
            "Biogen BIIB080 (IONIS-MAPTRx) -- anti-tau ASO, Phase 1/2 for "
            "mild AD (NCT03186989).  Showed ~50% CSF tau reduction.  "
            "No CRISPR trials for MAPT as of 2026."
        ),
        "conditions": ["frontotemporal_dementia", "FTD", "FTDP_17",
                        "tauopathy", "neurodegenerative", "PSP", "CBD"],
    },

    # -------------------------------------------------------------------
    # 3b.  GRN -- Progranulin
    # -------------------------------------------------------------------
    "GRN": {
        "gene_id": 2896,
        "chrom": "chr17",
        "start": 44_345_086,
        "end": 44_353_106,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.31",
        "exon_count": 13,
        "role": (
            "Progranulin -- secreted glycoprotein, neurotrophic factor and "
            "lysosomal regulator.  Haploinsufficiency (heterozygous loss-of-"
            "function) causes FTD-GRN, a common genetic form of FTD.  "
            "Homozygous loss causes neuronal ceroid lipofuscinosis type 11 "
            "(CLN11).  Progranulin levels correlate inversely with TDP-43 "
            "pathology."
        ),
        "disease": "Frontotemporal dementia (FTD-GRN) / CLN11",
        "omim_disease": 607485,
        "omim_gene": 138945,
        "inheritance": "AD (FTD); AR (CLN11)",
        "key_variants": [
            {
                "name": "p.Arg493Ter (R493X)",
                "rsid": "rs63750723",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Common founder mutation in Flanders/Netherlands.",
            },
            {
                "name": "c.813_816del (IVS7-1)",
                "rsid": None,
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": "Causes exon 8 skipping and NMD.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated GRN gene replacement (cDNA ~1.8 kb, fits AAV) "
            "to restore progranulin levels -- lead approach; "
            "(2) CRISPRa to upregulate the intact GRN allele in heterozygous "
            "patients (boosting to normal levels may suffice); "
            "(3) CRISPR disruption of sortilin (SORT1) to reduce progranulin "
            "lysosomal degradation and raise plasma levels.  "
            "Delivery: AAV9 or AAVhu68 intracisternal."
        ),
        "clinical_programs": (
            "Prevail PR006 (AAV9-GRN) Phase 1/2 for FTD-GRN (NCT04408625).  "
            "Passage Bio PBFT02 (AAVhu68-GRN) Phase 1/2 (NCT04747431).  "
            "Alector AL001 anti-sortilin antibody to raise GRN levels Phase 3.  "
            "No CRISPR trials for GRN as of 2026."
        ),
        "conditions": ["frontotemporal_dementia", "FTD", "FTD_GRN",
                        "CLN11", "neurodegenerative", "TDP43_proteinopathy"],
    },

    # -------------------------------------------------------------------
    # 3c.  C9orf72 -- Chromosome 9 Open Reading Frame 72
    # -------------------------------------------------------------------
    "C9orf72": {
        "gene_id": 203228,
        "chrom": "chr9",
        "start": 27_546_544,
        "end": 27_573_866,
        "strand": "-",
        "refseq": "NC_000009.12",
        "cytoband": "9p21.2",
        "exon_count": 12,
        "role": (
            "C9orf72 -- GTPase regulating autophagy and vesicle trafficking.  "
            "GGGGCC hexanucleotide repeat expansion in intron 1 is the most "
            "common genetic cause of both ALS and FTD (~40% familial ALS, "
            "~25% familial FTD).  Pathogenic mechanisms: (1) loss of C9orf72 "
            "function, (2) toxic RNA foci from sense/antisense transcripts, "
            "(3) dipeptide repeat (DPR) proteins from RAN translation."
        ),
        "disease": "ALS/FTD (C9orf72-related)",
        "omim_disease": 105550,
        "omim_gene": 614260,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "GGGGCC hexanucleotide repeat expansion",
                "rsid": None,
                "consequence": "repeat_expansion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Normal: 2-30 repeats.  Pathogenic: >30 (typically "
                    "hundreds to thousands).  Cannot be detected by standard "
                    "short-read sequencing; requires repeat-primed PCR or "
                    "long-read sequencing."
                ),
            },
        ],
        "strategy": (
            "(1) CRISPR excision of the expanded repeat region; "
            "(2) CRISPRi to silence sense and antisense transcription from "
            "the expanded allele while preserving normal allele expression; "
            "(3) Cas13-mediated degradation of toxic repeat RNA; "
            "(4) ASO targeting sense strand (BIIB078 approach).  "
            "Delivery: AAV9 intrathecal or IV; ASO intrathecal.  "
            "Note: CRISPR excision is complicated by repeat size (can be "
            ">10 kb) and risk of chromosomal rearrangement."
        ),
        "clinical_programs": (
            "Biogen BIIB078 (tofersen-like ASO for C9orf72) Phase 1 "
            "(NCT04288856) -- discontinued due to lack of efficacy.  "
            "Wave Life Sciences WVE-004 (stereopure ASO) Phase 1/2 -- "
            "discontinued.  Ionis/Biogen redesigning ASO approach.  "
            "No CRISPR trials for C9orf72 as of 2026."
        ),
        "conditions": ["ALS", "amyotrophic_lateral_sclerosis",
                        "frontotemporal_dementia", "FTD", "C9orf72_ALS_FTD",
                        "neurodegenerative", "motor_neuron_disease"],
    },
}


# ============================================================================
# 4. ALS -- ADDITIONAL TARGETS (SOD1 is in hereditary_disease_targets.py)
# ============================================================================

ALS_ADDITIONAL_TARGETS = {

    # -------------------------------------------------------------------
    # 4a.  FUS -- Fused in Sarcoma
    # -------------------------------------------------------------------
    "FUS": {
        "gene_id": 2521,
        "chrom": "chr16",
        "start": 31_180_113,
        "end": 31_194_867,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16p11.2",
        "exon_count": 15,
        "role": (
            "FUS -- RNA-binding protein involved in transcription, splicing, "
            "and DNA repair.  Mutations (mostly in the C-terminal NLS) cause "
            "aggressive juvenile-onset ALS (ALS6).  FUS mislocalizes from "
            "nucleus to cytoplasm, forming toxic aggregates."
        ),
        "disease": "Amyotrophic lateral sclerosis (ALS6)",
        "omim_disease": 608030,
        "omim_gene": 137070,
        "inheritance": "AD (most); AR (some juvenile)",
        "key_variants": [
            {
                "name": "p.Pro525Leu (P525L)",
                "rsid": "rs121909668",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "NLS domain; aggressive juvenile ALS, onset <25 years.",
            },
            {
                "name": "p.Arg521Cys (R521C)",
                "rsid": "rs121909664",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common FUS mutation; NLS disruption.",
            },
        ],
        "strategy": (
            "(1) Allele-specific CRISPRi/Cas13 to silence mutant FUS while "
            "preserving wild-type (FUS is essential); "
            "(2) Base editing for specific NLS mutations (P525L, R521C); "
            "(3) CRISPR-enhanced nuclear import signal repair.  "
            "Delivery: AAV9 intrathecal."
        ),
        "clinical_programs": "No gene therapy or CRISPR trials for FUS as of 2026.",
        "conditions": ["ALS", "amyotrophic_lateral_sclerosis", "ALS6",
                        "juvenile_ALS", "neurodegenerative", "motor_neuron_disease"],
    },

    # -------------------------------------------------------------------
    # 4b.  TARDBP -- TDP-43
    # -------------------------------------------------------------------
    "TARDBP": {
        "gene_id": 23435,
        "chrom": "chr1",
        "start": 11_012_344,
        "end": 11_025_739,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p36.22",
        "exon_count": 6,
        "role": (
            "TDP-43 (TAR DNA-binding protein 43) -- RNA-binding protein, major "
            "component of ubiquitinated inclusions in ~97% of ALS and ~45% of "
            "FTD cases.  Mutations cause ALS10.  TDP-43 pathology (cytoplasmic "
            "aggregation + nuclear depletion) is the hallmark of most ALS."
        ),
        "disease": "Amyotrophic lateral sclerosis (ALS10)",
        "omim_disease": 612069,
        "omim_gene": 605078,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Met337Val (M337V)",
                "rsid": "rs80356726",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Glycine-rich domain; enhanced aggregation propensity.",
            },
            {
                "name": "p.Ala382Thr (A382T)",
                "rsid": "rs80356731",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Sardinian founder mutation; 3% of Sardinian ALS.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant TARDBP allele; "
            "(2) Base editing for specific missense mutations; "
            "(3) CRISPRi approach must be allele-specific as TDP-43 is "
            "essential for RNA processing.  "
            "Delivery: AAV9 intrathecal."
        ),
        "clinical_programs": "No gene therapy or CRISPR trials for TARDBP as of 2026.",
        "conditions": ["ALS", "amyotrophic_lateral_sclerosis", "ALS10",
                        "neurodegenerative", "motor_neuron_disease",
                        "TDP43_proteinopathy"],
    },

    # -------------------------------------------------------------------
    # 4c.  ATXN2 -- Ataxin-2
    # -------------------------------------------------------------------
    "ATXN2": {
        "gene_id": 6311,
        "chrom": "chr12",
        "start": 111_452_214,
        "end": 111_599_676,
        "strand": "-",
        "refseq": "NC_000012.12",
        "cytoband": "12q24.12",
        "exon_count": 25,
        "role": (
            "Ataxin-2 -- polyglutamine protein, RNA metabolism regulator.  "
            "Long expansions (>34 repeats) cause spinocerebellar ataxia type 2 "
            "(SCA2).  Intermediate expansions (27-33 repeats) are a significant "
            "ALS risk factor and modifier.  ATXN2 lowering rescues TDP-43 "
            "toxicity in animal models."
        ),
        "disease": "ALS risk modifier / Spinocerebellar ataxia type 2",
        "omim_disease": 183090,
        "omim_gene": 601517,
        "inheritance": "AD (SCA2); risk modifier (ALS)",
        "key_variants": [
            {
                "name": "CAG repeat expansion (intermediate 27-33)",
                "rsid": None,
                "consequence": "repeat_expansion",
                "clinical_significance": "risk_factor",
                "notes": (
                    "Normal: 22-23 repeats.  ALS risk: 27-33 repeats.  "
                    "SCA2: >34 repeats.  ATXN2 lowering reduces TDP-43 "
                    "aggregation in models."
                ),
            },
        ],
        "strategy": (
            "(1) ASO-mediated ATXN2 reduction (Biogen BIIB105 approach); "
            "(2) CRISPRi to lower ATXN2 expression -- reducing TDP-43 "
            "toxicity independent of specific ALS genotype; "
            "(3) CRISPR excision of expanded CAG tract.  "
            "Delivery: intrathecal ASO or AAV9 intrathecal."
        ),
        "clinical_programs": (
            "Biogen BIIB105 (ION541) -- anti-ATXN2 ASO for ALS, Phase 1 "
            "(NCT04494256).  Preclinical data showed dramatic rescue "
            "of TDP-43 proteinopathy in mice.  No CRISPR trials."
        ),
        "conditions": ["ALS", "amyotrophic_lateral_sclerosis", "SCA2",
                        "spinocerebellar_ataxia", "neurodegenerative",
                        "motor_neuron_disease"],
    },
}


# ============================================================================
# Combined export
# ============================================================================

ALL_NEURODEGENERATIVE_TARGETS = {}
ALL_NEURODEGENERATIVE_TARGETS.update(ALZHEIMERS_TARGETS)
ALL_NEURODEGENERATIVE_TARGETS.update(PARKINSONS_TARGETS)
ALL_NEURODEGENERATIVE_TARGETS.update(FTD_TARGETS)
ALL_NEURODEGENERATIVE_TARGETS.update(ALS_ADDITIONAL_TARGETS)
