"""
Expanded disease gene-therapy / CRISPR targets -- GRCh38 (hg38) coordinates.

Covers monogenic diabetes (MODY/neonatal), monogenic obesity, ophthalmology
(beyond RPE65/CEP290), neuronal ceroid lipofuscinoses (Batten disease),
bone marrow failure syndromes, immune checkpoint targets for cancer
immunotherapy, and craniofacial disorders.

Each gene entry includes verified NCBI GRCh38.p14 coordinates, key
pathogenic variants, CRISPR/gene-therapy strategies, and clinical-trial
or approved-therapy status as of early 2026.

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
  - CRISPR Medicine News
  - Innovative Genomics Institute clinical-trial tracker
  - Published literature through early 2026
"""


# ============================================================================
# 1. MONOGENIC DIABETES TARGETS
#    Maturity-Onset Diabetes of the Young (MODY) and neonatal diabetes
# ============================================================================

MONOGENIC_DIABETES_TARGETS = {

    # -------------------------------------------------------------------
    # 1a.  GCK -- Glucokinase  (MODY2)
    # -------------------------------------------------------------------
    "GCK": {
        "gene_id": 2645,
        "chrom": "chr7",
        "start": 44_143_213,
        "end": 44_189_439,
        "strand": "-",
        "refseq": "NC_000007.14",
        "cytoband": "7p13",
        "exon_count": 15,
        "role": (
            "Glucokinase -- hexokinase IV.  Key glucose sensor in pancreatic "
            "beta cells and hepatocytes.  Catalyses the first committed step "
            "of glycolysis (glucose -> glucose-6-phosphate).  In beta cells, "
            "GCK activity sets the threshold for glucose-stimulated insulin "
            "secretion (GSIS).  Heterozygous loss-of-function mutations shift "
            "the GSIS set-point upward, causing mild fasting hyperglycaemia "
            "(MODY2).  Homozygous/compound-heterozygous loss causes permanent "
            "neonatal diabetes mellitus (PNDM)."
        ),
        "disease": "MODY2 (GCK-MODY) / Permanent Neonatal Diabetes (homozygous)",
        "omim_disease": 125851,
        "omim_gene": 138079,
        "inheritance": "AD (MODY2); AR (PNDM)",
        "key_variants": [
            {
                "name": "p.Gly261Arg (G261R)",
                "hgvs_coding": "NM_000162.5:c.781G>A",
                "rsid": "rs193922262",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Catalytic domain; abolishes kinase activity.",
            },
            {
                "name": "p.Thr228Met (T228M)",
                "hgvs_coding": "NM_000162.5:c.683C>T",
                "rsid": "rs80356650",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Recurrent European mutation; reduces Vmax ~80%.",
            },
            {
                "name": "p.Glu256Lys (E256K)",
                "rsid": "rs80356668",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Glucose binding domain; common in Asian MODY2 cohorts.",
            },
        ],
        "mutational_landscape": (
            "Over 800 pathogenic GCK variants in ClinVar/HGMD spanning all "
            "10 coding exons.  ~65% are missense, ~15% nonsense/frameshift, "
            "~10% splice-site.  Most MODY2 patients are heterozygous with mild "
            "stable fasting hyperglycaemia (5.5-8 mmol/L) that rarely requires "
            "pharmacotherapy.  Homozygous null = severe PNDM requiring insulin."
        ),
        "strategy": (
            "MODY2 is generally mild and rarely requires treatment beyond diet.  "
            "Gene therapy is mainly relevant for homozygous PNDM: "
            "(1) AAV8 liver-directed gene replacement delivering GCK cDNA "
            "(~1.4 kb, easily fits AAV) under a hepatocyte-specific promoter "
            "to restore hepatic glucose sensing; "
            "(2) Base editing (ABE/CBE) of specific missense mutations in "
            "autologous iPSC-derived beta cells for transplantation; "
            "(3) CRISPRa to upregulate the intact allele in MODY2 heterozygotes "
            "if treatment is warranted.  "
            "No CRISPR clinical trials as of 2026; all preclinical."
        ),
        "clinical_programs": (
            "No gene therapy trials.  GCK activators (dorzagliatin/HMS5552) "
            "approved in China (2022) for T2DM; not standard for MODY2.  "
            "MODY2 management is typically observation only."
        ),
        "conditions": ["MODY2", "monogenic_diabetes", "neonatal_diabetes",
                        "glucokinase_deficiency"],
    },

    # -------------------------------------------------------------------
    # 1b.  HNF1A -- Hepatocyte Nuclear Factor 1-alpha  (MODY3)
    # -------------------------------------------------------------------
    "HNF1A": {
        "gene_id": 6927,
        "chrom": "chr12",
        "start": 120_978_543,
        "end": 121_002_512,
        "strand": "+",
        "refseq": "NC_000012.12",
        "cytoband": "12q24.31",
        "exon_count": 9,
        "role": (
            "Hepatocyte nuclear factor 1-alpha -- homeodomain-containing "
            "transcription factor expressed in liver, kidney, intestine, and "
            "pancreatic islets.  Regulates genes involved in glucose transport "
            "(GLUT2/SLC2A2), glucose metabolism, and insulin secretion.  "
            "Loss-of-function causes the most common form of MODY (MODY3), "
            "characterised by progressive beta-cell dysfunction and "
            "sensitivity to sulfonylureas."
        ),
        "disease": "MODY3 (HNF1A-MODY)",
        "omim_disease": 600496,
        "omim_gene": 142410,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Pro291fsinsC",
                "hgvs_coding": "NM_000545.8:c.872dupC",
                "rsid": "rs137853233",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Poly-C tract in exon 4; most common HNF1A mutation "
                    "worldwide (~8% of all MODY3).  Hotspot for replication "
                    "slippage."
                ),
            },
            {
                "name": "p.Arg271Trp (R271W)",
                "hgvs_coding": "NM_000545.8:c.811C>T",
                "rsid": "rs137853232",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Homeodomain; impairs DNA binding.",
            },
            {
                "name": "p.Gly292fs (G292fs)",
                "hgvs_coding": "NM_000545.8:c.874_875delinsA",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Adjacent to P291fs hotspot.",
            },
        ],
        "mutational_landscape": (
            "Over 400 pathogenic variants in ClinVar.  ~60% missense, ~25% "
            "truncating (nonsense/frameshift), ~10% splice-site.  Mutations "
            "cluster in the dimerisation domain (exons 1-2), DNA-binding "
            "homeodomain (exons 2-4), and transactivation domain (exons 5-10).  "
            "MODY3 patients respond well to low-dose sulfonylureas."
        ),
        "strategy": (
            "(1) Low-dose sulfonylureas remain first-line therapy and are "
            "highly effective for decades.  Gene therapy primarily of research "
            "interest: (2) AAV-mediated HNF1A gene replacement in beta cells "
            "(cDNA ~1.9 kb, fits AAV) via pancreatic duct delivery; "
            "(3) Base editing for recurrent missense variants (e.g., R271W "
            "is a C>T amenable to ABE correction); "
            "(4) CRISPRa to boost expression from the wild-type allele.  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Sulfonylureas (gliclazide, glipizide) "
            "are standard pharmacotherapy.  Pharmacogenomic testing for MODY3 "
            "diagnosis increasingly integrated into diabetes care."
        ),
        "conditions": ["MODY3", "monogenic_diabetes", "HNF1A_MODY"],
    },

    # -------------------------------------------------------------------
    # 1c.  HNF4A -- Hepatocyte Nuclear Factor 4-alpha  (MODY1)
    # -------------------------------------------------------------------
    "HNF4A": {
        "gene_id": 3172,
        "chrom": "chr20",
        "start": 44_355_699,
        "end": 44_434_596,
        "strand": "+",
        "refseq": "NC_000020.11",
        "cytoband": "20q13.12",
        "exon_count": 15,
        "role": (
            "Hepatocyte nuclear factor 4-alpha -- nuclear receptor "
            "transcription factor.  Master regulator of hepatocyte "
            "differentiation and lipid/glucose metabolism.  In beta cells, "
            "HNF4A regulates HNF1A expression, creating a transcriptional "
            "cascade essential for insulin secretion.  Loss-of-function "
            "causes MODY1 with progressive beta-cell failure; neonatal "
            "macrosomia and transient hyperinsulinaemic hypoglycaemia "
            "are common at birth."
        ),
        "disease": "MODY1 (HNF4A-MODY)",
        "omim_disease": 125850,
        "omim_gene": 600281,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Arg127Trp (R127W)",
                "hgvs_coding": "NM_175914.5:c.379C>T",
                "rsid": "rs137853336",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "DNA-binding domain; disrupts target gene regulation.",
            },
            {
                "name": "p.Arg154Ter (R154X)",
                "hgvs_coding": "NM_175914.5:c.460C>T",
                "rsid": "rs137853337",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Truncating; classic MODY1 presentation.",
            },
        ],
        "mutational_landscape": (
            "~200 pathogenic variants in ClinVar.  Second-rarest of the common "
            "MODY subtypes.  Mutations distributed across all functional "
            "domains.  MODY1 closely resembles MODY3 clinically but with "
            "more prominent lipid abnormalities."
        ),
        "strategy": (
            "(1) Sulfonylureas are first-line therapy (as for MODY3); "
            "(2) AAV-mediated HNF4A gene supplementation in beta cells "
            "(cDNA ~1.4 kb, fits AAV); "
            "(3) CRISPRa to upregulate the intact allele; "
            "(4) Base editing for specific missense variants.  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Standard care is sulfonylureas.  "
            "Neonatal hypoglycaemia management with diazoxide."
        ),
        "conditions": ["MODY1", "monogenic_diabetes", "HNF4A_MODY"],
    },

    # -------------------------------------------------------------------
    # 1d.  HNF1B -- Hepatocyte Nuclear Factor 1-beta  (MODY5)
    # -------------------------------------------------------------------
    "HNF1B": {
        "gene_id": 6928,
        "chrom": "chr17",
        "start": 37_686_431,
        "end": 37_745_059,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q12",
        "exon_count": 11,
        "role": (
            "Hepatocyte nuclear factor 1-beta -- transcription factor "
            "essential for kidney, pancreas, liver, and genital tract "
            "development.  Regulates genes involved in nephron "
            "segmentation, bile acid metabolism, and insulin secretion.  "
            "Loss-of-function causes MODY5 (renal cysts and diabetes "
            "syndrome, RCAD).  Whole-gene deletions at 17q12 are common "
            "(~50% of cases) and also cause 17q12 microdeletion syndrome "
            "with neurodevelopmental features."
        ),
        "disease": "MODY5 / Renal cysts and diabetes syndrome (RCAD)",
        "omim_disease": 137920,
        "omim_gene": 189907,
        "inheritance": "AD (often de novo)",
        "key_variants": [
            {
                "name": "17q12 whole-gene deletion (~1.4 Mb)",
                "consequence": "deletion",
                "clinical_significance": "pathogenic",
                "frequency": "~50% of MODY5 cases",
                "notes": (
                    "Recurrent microdeletion mediated by flanking segmental "
                    "duplications.  Causes diabetes + renal malformations + "
                    "possible neurodevelopmental features."
                ),
            },
            {
                "name": "p.Arg177Ter (R177X)",
                "hgvs_coding": "NM_000458.4:c.529C>T",
                "rsid": "rs137853234",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Truncating; classic MODY5 with renal cysts.",
            },
        ],
        "mutational_landscape": (
            "~50% of MODY5 is caused by 17q12 microdeletion; the remainder "
            "are point mutations (missense, nonsense, frameshift) throughout "
            "the gene.  Unlike other MODY subtypes, MODY5 does NOT respond "
            "to sulfonylureas and usually requires insulin."
        ),
        "strategy": (
            "(1) HNF1B-MODY typically requires insulin therapy; "
            "(2) For point mutations: base/prime editing of specific variants; "
            "(3) For 17q12 deletions: gene replacement via AAV (cDNA ~1.7 kb "
            "fits AAV) directed to kidney and pancreas -- challenging due to "
            "multi-organ involvement; "
            "(4) iPSC-derived beta-cell transplantation with corrected HNF1B.  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Standard care is insulin + management "
            "of renal complications.  Genetic testing increasingly used for "
            "early diagnosis of renal cystic disease."
        ),
        "conditions": ["MODY5", "monogenic_diabetes", "RCAD",
                        "renal_cysts_and_diabetes", "17q12_deletion"],
    },

    # -------------------------------------------------------------------
    # 1e.  KCNJ11 -- Kir6.2  (Neonatal Diabetes)
    # -------------------------------------------------------------------
    "KCNJ11": {
        "gene_id": 3767,
        "chrom": "chr11",
        "start": 17_385_248,
        "end": 17_389_346,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11p15.1",
        "exon_count": 4,
        "role": (
            "Potassium inwardly-rectifying channel subfamily J member 11 "
            "(Kir6.2) -- pore-forming subunit of the pancreatic beta-cell "
            "ATP-sensitive potassium channel (KATP).  Forms a hetero-octamer "
            "with 4 ABCC8/SUR1 subunits.  Channel closure in response to "
            "rising ATP/ADP ratio triggers membrane depolarisation, Ca2+ "
            "influx, and insulin granule exocytosis.  Gain-of-function "
            "mutations keep the channel open, preventing insulin secretion "
            "(neonatal diabetes).  Loss-of-function causes congenital "
            "hyperinsulinism."
        ),
        "disease": "Permanent/Transient Neonatal Diabetes Mellitus (PNDM/TNDM)",
        "omim_disease": 606176,
        "omim_gene": 600937,
        "inheritance": "AD (gain-of-function, often de novo)",
        "key_variants": [
            {
                "name": "p.Arg201His (R201H)",
                "hgvs_coding": "NM_000525.4:c.602G>A",
                "rsid": "rs80356611",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common KCNJ11 PNDM mutation (~30% of cases).  "
                    "Reduces ATP sensitivity of KATP channel.  Most patients "
                    "respond to high-dose sulfonylureas (glibenclamide)."
                ),
            },
            {
                "name": "p.Arg201Cys (R201C)",
                "hgvs_coding": "NM_000525.4:c.601C>T",
                "rsid": "rs80356610",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Same codon as R201H; similar phenotype.",
            },
            {
                "name": "p.Val59Met (V59M)",
                "hgvs_coding": "NM_000525.4:c.175G>A",
                "rsid": "rs80356617",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Associated with DEND syndrome (Developmental delay, "
                    "Epilepsy, Neonatal Diabetes) when severe."
                ),
            },
        ],
        "mutational_landscape": (
            "~50 pathogenic gain-of-function variants in ClinVar.  Almost all "
            "are missense mutations affecting ATP binding, pore gating, or "
            "subunit interactions.  ~90% of KCNJ11-PNDM patients can switch "
            "from insulin to oral sulfonylureas, making this one of the most "
            "successful examples of pharmacogenomics in diabetes."
        ),
        "strategy": (
            "Sulfonylureas (high-dose glibenclamide/glyburide) are first-line "
            "and highly effective for ~90% of patients.  Gene therapy of "
            "limited relevance given pharmacological success.  Potential "
            "approaches: (1) Allele-specific silencing of gain-of-function "
            "allele via CRISPRi for sulfonylurea-resistant cases; "
            "(2) Base editing to correct specific missense mutations (R201H "
            "is G>A, ideal ABE target); "
            "(3) Prime editing for DEND syndrome variants where neurological "
            "involvement is not addressed by sulfonylureas.  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "Sulfonylureas remain transformative treatment.  Multiple "
            "observational studies confirm long-term efficacy (>15 years).  "
            "Gene therapy for sulfonylurea-resistant DEND syndrome is a "
            "potential future application."
        ),
        "conditions": ["neonatal_diabetes", "PNDM", "TNDM", "DEND_syndrome",
                        "monogenic_diabetes", "KATP_channel"],
    },

    # -------------------------------------------------------------------
    # 1f.  ABCC8 -- SUR1  (Neonatal Diabetes)
    # -------------------------------------------------------------------
    "ABCC8": {
        "gene_id": 6833,
        "chrom": "chr11",
        "start": 17_392_498,
        "end": 17_476_845,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11p15.1",
        "exon_count": 38,
        "role": (
            "ATP-binding cassette subfamily C member 8 (sulfonylurea "
            "receptor 1, SUR1) -- regulatory subunit of the pancreatic "
            "beta-cell KATP channel.  Four SUR1 subunits surround the "
            "Kir6.2 pore tetramer.  SUR1 contains two nucleotide-binding "
            "domains (NBDs) that sense MgADP and confer sulfonylurea "
            "sensitivity.  Gain-of-function mutations cause neonatal "
            "diabetes; loss-of-function causes congenital hyperinsulinism "
            "(CHI) -- the most common genetic cause of persistent "
            "hypoglycaemia in infancy."
        ),
        "disease": "Neonatal Diabetes / Congenital Hyperinsulinism",
        "omim_disease": 606176,
        "omim_gene": 600509,
        "inheritance": "AD (GOF -> NDM) / AR (LOF -> CHI)",
        "key_variants": [
            {
                "name": "p.Arg1380Leu (R1380L)",
                "hgvs_coding": "NM_000352.6:c.4139G>T",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "NBD2 domain; gain-of-function causing PNDM.",
            },
            {
                "name": "c.3989-9G>A (splice)",
                "consequence": "splice_variant",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Deep intronic; common founder mutation in Ashkenazi "
                    "Jewish congenital hyperinsulinism (~90% of CHI alleles "
                    "in this population)."
                ),
            },
            {
                "name": "p.Val187Asp (V187D)",
                "hgvs_coding": "NM_000352.6:c.560T>A",
                "rsid": "rs137852699",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "TMD0 domain; loss-of-function causing diffuse CHI.  "
                    "Homozygous patients may require near-total pancreatectomy."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 450 pathogenic variants in ClinVar.  Loss-of-function "
            "mutations (CHI) are more common than gain-of-function (NDM).  "
            "CHI mutations distributed across all domains, with clusters in "
            "NBD1, NBD2, and TMD0.  Focal CHI may be treatable by limited "
            "pancreatectomy if paternal LOH at 11p15 leads to focal lesion."
        ),
        "strategy": (
            "For neonatal diabetes (GOF): sulfonylureas as for KCNJ11.  "
            "For congenital hyperinsulinism (LOF): "
            "(1) Gene replacement via AAV delivering ABCC8 cDNA (~4.7 kb, "
            "near AAV capacity limit; dual-AAV or lentiviral may be needed); "
            "(2) Base editing for specific missense mutations causing "
            "diffuse CHI; "
            "(3) Ex vivo editing of iPSC-derived beta cells with corrected "
            "ABCC8, followed by encapsulated transplantation; "
            "(4) For focal CHI: no gene therapy needed -- focal lesion "
            "resection is curative.  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "Diazoxide is first-line for diffuse CHI; octreotide for "
            "diazoxide-unresponsive cases.  Near-total pancreatectomy "
            "for medically refractory diffuse CHI.  Gene therapy for "
            "diffuse CHI is an active preclinical area."
        ),
        "conditions": ["neonatal_diabetes", "congenital_hyperinsulinism",
                        "monogenic_diabetes", "KATP_channel", "CHI"],
    },
}


# ============================================================================
# 2. OBESITY TARGETS
#    Monogenic obesity genes in the leptin-melanocortin pathway
# ============================================================================

OBESITY_TARGETS = {

    # -------------------------------------------------------------------
    # 2a.  MC4R -- Melanocortin 4 Receptor
    # -------------------------------------------------------------------
    "MC4R": {
        "gene_id": 4160,
        "chrom": "chr18",
        "start": 60_371_062,
        "end": 60_372_775,
        "strand": "-",
        "refseq": "NC_000018.10",
        "cytoband": "18q21.32",
        "exon_count": 1,
        "role": (
            "Melanocortin 4 receptor -- G-protein-coupled receptor in the "
            "hypothalamic paraventricular nucleus.  Central integrator of "
            "the leptin-melanocortin satiety pathway.  Activation by "
            "alpha-MSH (from POMC cleavage) suppresses appetite and "
            "increases energy expenditure.  MC4R is the single most common "
            "monogenic cause of severe obesity, accounting for 2-6% of "
            "severe early-onset obesity cases.  The gene is intronless."
        ),
        "disease": "MC4R-associated monogenic obesity",
        "omim_disease": 601665,
        "omim_gene": 155541,
        "inheritance": "AD (incomplete penetrance) / AR (more severe)",
        "key_variants": [
            {
                "name": "p.Thr11Ser (T11S)",
                "rsid": "rs2229616",
                "consequence": "missense",
                "clinical_significance": "protective (gain-of-function)",
                "notes": "Gain-of-function; associated with lower BMI.",
            },
            {
                "name": "p.Ile251Leu (I251L)",
                "rsid": "rs52820871",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Loss-of-function; common in severe early-onset obesity.",
            },
            {
                "name": "p.Arg165Trp (R165W)",
                "rsid": "rs121913563",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Transmembrane domain 4; impairs receptor trafficking "
                    "and signalling.  Homozygous carriers have severe obesity."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 200 pathogenic/likely-pathogenic MC4R variants in ClinVar.  "
            "Nearly all are missense (intronless gene).  Most cause "
            "loss-of-function through impaired receptor trafficking, ligand "
            "binding, or G-protein coupling.  Heterozygous carriers show "
            "variable obesity; homozygous/compound-het carriers have severe "
            "early-onset obesity with hyperphagia."
        ),
        "strategy": (
            "(1) Setmelanotide (IMCIVREE, approved 2020) -- MC4R pathway "
            "agonist effective for POMC/PCSK1/LEPR deficiency obesity but "
            "NOT for MC4R loss-of-function (receptor itself is defective); "
            "(2) Gene replacement via AAV targeting hypothalamus -- "
            "challenging due to CNS delivery requirements; intrathecal or "
            "intra-hypothalamic stereotactic injection of AAV9-MC4R; "
            "(3) CRISPRa to upregulate the wild-type allele in heterozygous "
            "patients; "
            "(4) Base/prime editing for specific missense variants, "
            "delivered via AAV9 to hypothalamus.  "
            "All gene therapy approaches are preclinical as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Setmelanotide (Rhythm Pharmaceuticals) "
            "approved for POMC/PCSK1/LEPR deficiency but not MC4R deficiency.  "
            "Anti-obesity GLP-1 receptor agonists (semaglutide, tirzepatide) "
            "show partial efficacy in MC4R heterozygotes."
        ),
        "conditions": ["monogenic_obesity", "MC4R_deficiency",
                        "severe_early_onset_obesity", "leptin_melanocortin_pathway"],
    },

    # -------------------------------------------------------------------
    # 2b.  LEP -- Leptin
    # -------------------------------------------------------------------
    "LEP": {
        "gene_id": 3952,
        "chrom": "chr7",
        "start": 128_241_278,
        "end": 128_257_629,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q32.1",
        "exon_count": 3,
        "role": (
            "Leptin -- adipokine hormone secreted by white adipose tissue "
            "in proportion to fat mass.  Signals to hypothalamic neurons "
            "(via LEPR) to suppress appetite and increase energy expenditure.  "
            "Complete leptin deficiency causes severe early-onset obesity, "
            "hyperphagia, hypogonadism, and immune dysfunction.  One of the "
            "first monogenic obesity genes identified (1997)."
        ),
        "disease": "Congenital leptin deficiency",
        "omim_disease": 614962,
        "omim_gene": 164160,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg105Trp (R105W, delta-G133)",
                "hgvs_coding": "NM_000230.3:c.313C>T",
                "rsid": "rs17151919",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Original mutation identified in Pakistani consanguineous "
                    "families (Montague et al., Nature 1997).  Protein is "
                    "synthesised but not secreted."
                ),
            },
            {
                "name": "c.398delG (frameshift)",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Truncating; no functional leptin produced.",
            },
        ],
        "mutational_landscape": (
            "Very rare (~60 cases reported worldwide).  Most families are "
            "consanguineous.  Founder mutations in Pakistani (delta-G133), "
            "Turkish, and Egyptian kindreds.  All pathogenic variants "
            "abolish or severely reduce circulating leptin."
        ),
        "strategy": (
            "(1) Recombinant metreleptin (MYALEPT, approved 2014 for "
            "lipodystrophy; used off-label for congenital leptin deficiency) "
            "is transformative -- normalises weight, hunger, puberty, and "
            "immunity; requires daily injections for life.  "
            "(2) AAV-mediated LEP gene replacement in adipose tissue "
            "(cDNA ~0.5 kb, very small -- easily fits AAV) could provide "
            "a one-time cure; preclinical studies in ob/ob mice show "
            "long-term weight normalisation; "
            "(3) Base editing for the R105W founder mutation (C>T, ideal "
            "ABE target).  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "Metreleptin (Amylin/AstraZeneca/Amryt) available as enzyme "
            "replacement therapy.  Gene therapy preclinical in mouse models."
        ),
        "approved_therapy": "Metreleptin (MYALEPT) -- recombinant leptin",
        "conditions": ["monogenic_obesity", "leptin_deficiency",
                        "congenital_leptin_deficiency",
                        "leptin_melanocortin_pathway"],
    },

    # -------------------------------------------------------------------
    # 2c.  LEPR -- Leptin Receptor
    # -------------------------------------------------------------------
    "LEPR": {
        "gene_id": 3953,
        "chrom": "chr1",
        "start": 65_420_652,
        "end": 65_641_559,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p31.3",
        "exon_count": 24,
        "role": (
            "Leptin receptor (OB-R) -- type I cytokine receptor.  The long "
            "isoform (OB-Rb) mediates leptin signalling in hypothalamic "
            "neurons via JAK2-STAT3, PI3K, and MAPK pathways.  "
            "Loss-of-function causes severe early-onset obesity clinically "
            "identical to leptin deficiency, but unresponsive to "
            "recombinant leptin.  Patients respond to setmelanotide "
            "(downstream MC4R pathway agonist)."
        ),
        "disease": "Congenital leptin receptor deficiency",
        "omim_disease": 614963,
        "omim_gene": 601007,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Pro316Thr (P316T, A>C splice variant region)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Extracellular domain; impairs leptin binding.",
            },
            {
                "name": "p.His684Pro (H684P)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Impairs receptor activation.  Homozygous carriers have "
                    "severe obesity from infancy."
                ),
            },
        ],
        "mutational_landscape": (
            "~100 pathogenic variants reported.  Rare (~3% of severe "
            "early-onset obesity cases).  Mutations throughout extracellular "
            "and intracellular domains.  Homozygous/compound-heterozygous "
            "loss-of-function required for obesity phenotype."
        ),
        "strategy": (
            "(1) Setmelanotide (IMCIVREE, Rhythm Pharmaceuticals) -- "
            "FDA/EMA-approved for LEPR deficiency obesity (2020/2021).  "
            "Bypasses the defective receptor to activate MC4R directly.  "
            "Produces ~10% body weight loss and dramatic hunger reduction.  "
            "(2) AAV-mediated LEPR gene replacement in hypothalamic neurons "
            "(cDNA ~3.5 kb for long isoform, fits AAV); preclinical in "
            "db/db mice.  "
            "(3) Base editing for specific missense mutations.  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "Setmelanotide APPROVED for LEPR deficiency.  Phase 3 data: "
            "~45% of patients achieved >=10% weight loss at 1 year.  "
            "Gene therapy preclinical."
        ),
        "approved_therapy": "Setmelanotide (IMCIVREE) -- MC4R agonist",
        "conditions": ["monogenic_obesity", "leptin_receptor_deficiency",
                        "leptin_melanocortin_pathway"],
    },

    # -------------------------------------------------------------------
    # 2d.  POMC -- Proopiomelanocortin
    # -------------------------------------------------------------------
    "POMC": {
        "gene_id": 5443,
        "chrom": "chr2",
        "start": 25_160_860,
        "end": 25_168_580,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2p23.3",
        "exon_count": 4,
        "role": (
            "Proopiomelanocortin -- polypeptide precursor cleaved into "
            "multiple bioactive peptides: alpha-MSH, beta-MSH, ACTH, "
            "beta-endorphin, and others.  In the hypothalamus, POMC-derived "
            "alpha-MSH activates MC4R to suppress appetite.  POMC deficiency "
            "causes severe early-onset obesity, adrenal insufficiency "
            "(no ACTH), and red hair / pale skin (no alpha-MSH for MC1R).  "
            "Extremely rare (<50 cases reported)."
        ),
        "disease": "POMC deficiency obesity",
        "omim_disease": 609734,
        "omim_gene": 176830,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.7013G>A (p.Arg236Gly, disrupts ACTH/alpha-MSH)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Disrupts processing of ACTH and alpha-MSH.  "
                    "Homozygous carriers have classic POMC deficiency triad."
                ),
            },
            {
                "name": "p.Glu110Ter (E110X)",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Truncating; no functional POMC peptides produced.",
            },
        ],
        "mutational_landscape": (
            "Very rare.  Fewer than 50 cases worldwide.  Null mutations "
            "produce the full triad (obesity + adrenal insufficiency + "
            "red hair).  Heterozygous carriers may have partial obesity "
            "phenotype."
        ),
        "strategy": (
            "(1) Setmelanotide (IMCIVREE) -- APPROVED for POMC deficiency "
            "obesity (2020).  First targeted pharmacotherapy for a genetic "
            "obesity syndrome.  Bypasses the missing alpha-MSH.  "
            "(2) AAV-mediated POMC gene replacement in hypothalamic "
            "arcuate nucleus (cDNA ~1.0 kb, very small -- fits AAV); "
            "preclinical proof-of-concept in Pomc-null mice.  "
            "(3) Hydrocortisone replacement still needed for ACTH deficiency "
            "(setmelanotide does not address adrenal axis).  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "Setmelanotide APPROVED for POMC deficiency.  Hydrocortisone "
            "replacement for adrenal insufficiency.  Gene therapy preclinical."
        ),
        "approved_therapy": "Setmelanotide (IMCIVREE) -- MC4R agonist",
        "conditions": ["monogenic_obesity", "POMC_deficiency",
                        "adrenal_insufficiency", "leptin_melanocortin_pathway"],
    },

    # -------------------------------------------------------------------
    # 2e.  PCSK1 -- Proprotein Convertase Subtilisin/Kexin Type 1
    # -------------------------------------------------------------------
    "PCSK1": {
        "gene_id": 5122,
        "chrom": "chr5",
        "start": 96_390_333,
        "end": 96_433_248,
        "strand": "-",
        "refseq": "NC_000005.10",
        "cytoband": "5q15",
        "exon_count": 15,
        "role": (
            "Proprotein convertase subtilisin/kexin type 1 (PC1/3) -- "
            "neuroendocrine serine protease that cleaves prohormone "
            "precursors including POMC, proinsulin, proglucagon, and "
            "proGnRH.  Essential for processing alpha-MSH from POMC "
            "and mature insulin from proinsulin.  Loss-of-function causes "
            "severe obesity, malabsorptive diarrhoea in infancy, and "
            "endocrine dysfunction (hyperproinsulinaemia, hypogonadism, "
            "adrenal insufficiency)."
        ),
        "disease": "PCSK1 deficiency (proprotein convertase 1/3 deficiency)",
        "omim_disease": 600955,
        "omim_gene": 162150,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Asn222Asp (N222D)",
                "rsid": "rs6232",
                "consequence": "missense",
                "clinical_significance": "pathogenic/risk",
                "notes": (
                    "Reduces PC1/3 catalytic activity by ~50%.  Risk allele "
                    "for common obesity in GWAS (OR ~1.34) when homozygous."
                ),
            },
            {
                "name": "p.Ser690Thr (S690T)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "C-terminal domain; impairs enzyme maturation.",
            },
        ],
        "mutational_landscape": (
            "Very rare biallelic loss-of-function (~20 cases).  Common "
            "hypomorphic variants (N222D, rs6232) confer ~1.3-fold obesity "
            "risk in the general population.  Biallelic null mutations cause "
            "neonatal malabsorptive diarrhoea + obesity + endocrine failure."
        ),
        "strategy": (
            "(1) Setmelanotide (IMCIVREE) -- APPROVED for PCSK1 deficiency "
            "obesity (2022 label expansion).  Effective because PCSK1 "
            "deficiency prevents POMC processing to alpha-MSH, and "
            "setmelanotide acts as an alpha-MSH mimetic at MC4R.  "
            "(2) AAV-mediated PCSK1 gene replacement in enteroendocrine "
            "cells and hypothalamus (cDNA ~2.3 kb, fits AAV); "
            "(3) Base editing for specific missense variants.  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "Setmelanotide APPROVED for PCSK1 deficiency obesity (2022).  "
            "Pancreatic enzyme and hormone replacement for GI and endocrine "
            "manifestations.  Gene therapy preclinical."
        ),
        "approved_therapy": "Setmelanotide (IMCIVREE) -- MC4R agonist",
        "conditions": ["monogenic_obesity", "PCSK1_deficiency",
                        "proprotein_convertase_deficiency",
                        "leptin_melanocortin_pathway"],
    },
}


# ============================================================================
# 3. OPHTHALMOLOGY TARGETS
#    Beyond RPE65 and CEP290 (already in hereditary_disease_targets.py)
# ============================================================================

OPHTHALMOLOGY_TARGETS = {

    # -------------------------------------------------------------------
    # 3a.  PAX6 -- Paired Box 6  (Aniridia)
    # -------------------------------------------------------------------
    "PAX6": {
        "gene_id": 5080,
        "chrom": "chr11",
        "start": 31_789_026,
        "end": 31_817_961,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11p13",
        "exon_count": 25,
        "role": (
            "Paired box 6 -- master transcription factor for eye "
            "development.  Contains paired domain and homeodomain for "
            "DNA binding.  Essential for lens, iris, cornea, and retina "
            "morphogenesis.  Haploinsufficiency causes aniridia (absent "
            "iris) with progressive corneal opacification (aniridia-related "
            "keratopathy, ARK), cataracts, glaucoma, and foveal hypoplasia.  "
            "Homozygous loss is lethal (anencephaly).  PAX6 is also "
            "critical for pancreatic islet and CNS development."
        ),
        "disease": "Aniridia / Aniridia-related keratopathy (ARK)",
        "omim_disease": 106210,
        "omim_gene": 607108,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "p.Arg240Ter (R240X)",
                "hgvs_coding": "NM_000280.4:c.718C>T",
                "rsid": "rs121907916",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Recurrent CpG mutation; classic aniridia.",
            },
            {
                "name": "11p13 deletion encompassing PAX6",
                "consequence": "deletion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Contiguous gene deletion including WT1 causes WAGR "
                    "syndrome (Wilms tumour, Aniridia, Genitourinary "
                    "anomalies, Retardation)."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 500 unique PAX6 mutations catalogued.  ~70% are "
            "truncating (nonsense, frameshift, splice-site), consistent "
            "with haploinsufficiency mechanism.  Missense mutations in "
            "paired domain cause milder/atypical phenotypes.  ~1/3 are "
            "de novo."
        ),
        "strategy": (
            "(1) AAV-mediated PAX6 gene supplementation to corneal limbal "
            "stem cells for aniridia-related keratopathy (cDNA ~1.3 kb, "
            "fits AAV); delivery via subconjunctival or intrastromal "
            "injection; MUST use tightly regulated promoter to avoid "
            "PAX6 overexpression toxicity; "
            "(2) CRISPRa to upregulate the wild-type allele in limbal "
            "stem cells -- promising approach for haploinsufficiency; "
            "(3) Ex vivo gene-corrected limbal stem cell transplantation; "
            "(4) Ataluren (read-through of premature stop codons) -- "
            "tested in aniridia Phase 2 trial with mixed results.  "
            "No CRISPR clinical trials for aniridia as of 2026."
        ),
        "clinical_programs": (
            "Ataluren Phase 2 for aniridia (NCT02647359, completed with "
            "inconclusive results).  Multiple AAV gene supplementation "
            "programs in preclinical development.  CRISPR epigenomic "
            "upregulation (CRISPRa) at PAX6 locus in preclinical studies "
            "(University College London)."
        ),
        "conditions": ["aniridia", "aniridia_related_keratopathy",
                        "WAGR_syndrome", "ophthalmology", "corneal_disease"],
    },

    # -------------------------------------------------------------------
    # 3b.  ABCA4 -- ATP-Binding Cassette A4  (Stargardt Disease)
    # -------------------------------------------------------------------
    "ABCA4": {
        "gene_id": 24,
        "chrom": "chr1",
        "start": 93_992_834,
        "end": 94_121_148,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1p22.1",
        "exon_count": 50,
        "role": (
            "ATP-binding cassette transporter A4 -- retina-specific "
            "flippase in photoreceptor outer segment disc membranes.  "
            "Transports N-retinylidene-phosphatidylethanolamine (N-ret-PE) "
            "from the disc lumen to the cytoplasmic leaflet, enabling "
            "all-trans-retinal recycling.  Loss-of-function causes toxic "
            "accumulation of lipofuscin (A2E) in RPE cells, leading to "
            "progressive macular degeneration.  ABCA4 mutations cause "
            "Stargardt disease (most common inherited macular dystrophy), "
            "cone-rod dystrophy, and some cases of autosomal recessive "
            "retinitis pigmentosa."
        ),
        "disease": "Stargardt disease / Cone-rod dystrophy / AR-RP",
        "omim_disease": 248200,
        "omim_gene": 601691,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gly1961Glu (G1961E)",
                "hgvs_coding": "NM_000350.3:c.5882G>A",
                "rsid": "rs1800553",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "frequency": "Most common ABCA4 allele (~10-15% in Europeans)",
                "notes": (
                    "NBD2 domain; mild/moderate reduction in transport "
                    "activity.  Considered a 'mild' allele; disease "
                    "severity depends on second allele."
                ),
            },
            {
                "name": "c.5461-10T>C (splice variant)",
                "rsid": "rs61750646",
                "consequence": "splice_variant",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Common pathogenic variant in Northern Europeans.  "
                    "Creates a cryptic splice site leading to exon skipping."
                ),
            },
            {
                "name": "p.Leu541Pro;Ala1038Val (complex allele)",
                "consequence": "missense_complex",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Two missense changes in cis; individually mild but "
                    "severe in combination.  Common in Central Europeans."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 1,200 pathogenic variants in ClinVar.  Extreme allelic "
            "heterogeneity -- most patients are compound heterozygotes.  "
            "Genotype-phenotype correlation driven by residual ABCA4 "
            "activity: two null alleles -> severe/early cone-rod dystrophy; "
            "one mild + one null -> classic Stargardt; two mild -> late-onset."
        ),
        "strategy": (
            "ABCA4 cDNA is ~6.8 kb, exceeding AAV packaging capacity.  "
            "Strategies: "
            "(1) Dual-AAV approaches (split-intein trans-splicing or "
            "overlapping) to deliver full-length ABCA4; multiple groups "
            "in preclinical/early clinical development; "
            "(2) Lentiviral delivery (SAR422459/Sanofi, Star-LHON study) -- "
            "Phase I/II completed with subretinal lentivirus; "
            "(3) Non-viral lipid nanoparticle delivery of ABCA4 mRNA -- "
            "preclinical (Stoke Therapeutics/Bayer); "
            "(4) Antisense oligonucleotide to correct the common c.5461-10T>C "
            "splice variant (Stoke/ProQR approach); "
            "(5) CRISPR base editing for specific missense variants.  "
            "Challenge: photoreceptors are post-mitotic, so editing must "
            "be efficient in non-dividing cells."
        ),
        "clinical_programs": {
            "SAR422459_Sanofi": {
                "type": "Lentiviral ABCA4 gene replacement",
                "route": "subretinal",
                "status": "Phase I/II completed; modest efficacy, safety acceptable.",
            },
            "dual_AAV_programs": {
                "type": "Dual-AAV split-intein ABCA4",
                "sponsors": "Applied Genetic Technologies (AGTC), Iveric Bio (now Astellas)",
                "status": "Preclinical/IND-enabling as of 2026.",
            },
        },
        "conditions": ["stargardt_disease", "cone_rod_dystrophy",
                        "macular_dystrophy", "ophthalmology",
                        "retinal_dystrophy"],
    },

    # -------------------------------------------------------------------
    # 3c.  RS1 -- Retinoschisin  (X-linked Retinoschisis)
    # -------------------------------------------------------------------
    "RS1": {
        "gene_id": 6247,
        "chrom": "chrX",
        "start": 18_639_688,
        "end": 18_672_108,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xp22.13",
        "exon_count": 7,
        "role": (
            "Retinoschisin -- secreted protein that forms homo-octameric "
            "complexes on the photoreceptor and bipolar cell surface.  "
            "Essential for retinal layer adhesion and synaptic "
            "architecture.  Loss-of-function causes X-linked retinoschisis "
            "(XLRS), characterised by foveal schisis (splitting of retinal "
            "layers), reduced visual acuity, and sometimes peripheral "
            "schisis with vitreous haemorrhage.  Affects ~1:5,000-25,000 "
            "males."
        ),
        "disease": "X-linked retinoschisis (XLRS)",
        "omim_disease": 312700,
        "omim_gene": 300839,
        "inheritance": "XLR",
        "key_variants": [
            {
                "name": "p.Glu72Lys (E72K)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Discoidin domain; impairs octamer assembly.",
            },
            {
                "name": "p.Arg102Trp (R102W)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Recurrent; one of the most common XLRS mutations.",
            },
            {
                "name": "p.Cys59Ser (C59S)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts disulfide bonding critical for octamerisation.",
            },
        ],
        "mutational_landscape": (
            "Over 250 unique RS1 mutations reported.  ~50% missense, "
            "~30% splice-site, ~20% truncating.  Missense mutations "
            "cluster in the discoidin domain (exons 4-6).  Highly variable "
            "expressivity even within families."
        ),
        "strategy": (
            "(1) AAV8-RS1 intravitreal gene replacement (cDNA ~0.7 kb, "
            "very small -- easily fits AAV).  Several clinical programs: "
            "Applied Genetic Technologies (AGTC) Phase I/II and NEI/NIH "
            "Phase I/II both delivered intravitreal AAV8-RS1; "
            "(2) Novel capsids (AAV.SPR, AAV44.9) for improved "
            "photoreceptor transduction from intravitreal route; "
            "(3) Base editing for recurrent missense mutations.  "
            "RS1 is an attractive gene therapy target due to small cDNA "
            "size and extracellular protein (bystander effect from "
            "transduced cells)."
        ),
        "clinical_programs": {
            "AGTC_AAV8_RS1": {
                "nct": "NCT02416622",
                "type": "AAV8-RS1 intravitreal",
                "status": (
                    "Phase I/II completed; well tolerated but limited "
                    "efficacy at initial doses.  Dose escalation and "
                    "capsid optimization ongoing."
                ),
            },
            "NEI_NIH_RS1": {
                "nct": "NCT02317887",
                "type": "AAV8-RS1 intravitreal",
                "status": "Phase I/II completed; safety established.",
            },
        },
        "conditions": ["X_linked_retinoschisis", "XLRS", "ophthalmology",
                        "retinal_dystrophy"],
    },

    # -------------------------------------------------------------------
    # 3d.  BEST1 -- Bestrophin 1  (Best Disease)
    # -------------------------------------------------------------------
    "BEST1": {
        "gene_id": 7439,
        "chrom": "chr11",
        "start": 61_949_821,
        "end": 61_965_515,
        "strand": "+",
        "refseq": "NC_000011.10",
        "cytoband": "11q12.3",
        "exon_count": 13,
        "role": (
            "Bestrophin 1 -- calcium-activated chloride channel in the "
            "basolateral membrane of retinal pigment epithelium (RPE).  "
            "Regulates transepithelial ion transport, volume regulation, "
            "and the light peak of the electro-oculogram (EOG).  "
            "Mutations cause Best vitelliform macular dystrophy (BVMD), "
            "adult-onset vitelliform macular dystrophy, autosomal recessive "
            "bestrophinopathy, and retinitis pigmentosa."
        ),
        "disease": "Best vitelliform macular dystrophy (BVMD)",
        "omim_disease": 153700,
        "omim_gene": 607854,
        "inheritance": "AD (BVMD); AR (ARB)",
        "key_variants": [
            {
                "name": "p.Tyr227Asn (Y227N)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Transmembrane domain; impairs chloride channel function.",
            },
            {
                "name": "p.Ala243Val (A243V)",
                "rsid": "rs121912646",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Recurrent European mutation; classic BVMD.",
            },
        ],
        "mutational_landscape": (
            "Over 250 pathogenic variants.  Almost exclusively missense "
            "for dominant BVMD.  Truncating mutations more common in "
            "autosomal recessive bestrophinopathy.  Highly variable "
            "expressivity -- some carriers retain near-normal vision."
        ),
        "strategy": (
            "(1) AAV-mediated BEST1 gene supplementation in RPE cells "
            "(cDNA ~1.7 kb, fits AAV); subretinal delivery to target RPE; "
            "suppression-and-replacement needed for dominant-negative "
            "mutations; "
            "(2) CRISPRi to silence dominant-negative mutant allele "
            "combined with allele-specific resistant transgene; "
            "(3) Base editing for specific recurrent missense variants.  "
            "No CRISPR clinical trials as of 2026.  Multiple preclinical "
            "AAV programs."
        ),
        "clinical_programs": (
            "Preclinical AAV gene supplementation (multiple academic groups).  "
            "No clinical trials as of 2026."
        ),
        "conditions": ["best_disease", "vitelliform_macular_dystrophy",
                        "bestrophinopathy", "ophthalmology",
                        "retinal_dystrophy"],
    },

    # -------------------------------------------------------------------
    # 3e.  CYP1B1 -- Cytochrome P450 1B1  (Primary Congenital Glaucoma)
    # -------------------------------------------------------------------
    "CYP1B1": {
        "gene_id": 1545,
        "chrom": "chr2",
        "start": 38_067_509,
        "end": 38_076_151,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2p22.2",
        "exon_count": 3,
        "role": (
            "Cytochrome P450 family 1 subfamily B member 1 -- "
            "monooxygenase expressed in fetal and adult trabecular "
            "meshwork, iris, and ciliary body.  Metabolises retinoic "
            "acid, oestradiol, and melatonin.  Essential for anterior "
            "segment development.  Loss-of-function causes primary "
            "congenital glaucoma (PCG) -- the most common childhood "
            "glaucoma (~1:10,000 births; higher in consanguineous "
            "populations)."
        ),
        "disease": "Primary congenital glaucoma (PCG / buphthalmos)",
        "omim_disease": 231300,
        "omim_gene": 601771,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gly61Glu (G61E)",
                "hgvs_coding": "NM_000104.4:c.182G>A",
                "rsid": "rs28936415",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Abolishes enzyme activity; common in Middle Eastern PCG.",
            },
            {
                "name": "p.Arg368His (R368H)",
                "hgvs_coding": "NM_000104.4:c.1103G>A",
                "rsid": "rs28936416",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Haem-binding domain; common worldwide.",
            },
            {
                "name": "c.1209_1218insGTGCCCATCC (4340delG)",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Brazilian/Portuguese founder mutation.",
            },
        ],
        "mutational_landscape": (
            "Over 200 pathogenic variants.  High allelic heterogeneity "
            "but strong founder effects in consanguineous populations.  "
            "Saudi Arabia: p.G61E founder (~80% of PCG alleles).  Brazil: "
            "4340delG (~50%).  Genotype-phenotype correlation is limited."
        ),
        "strategy": (
            "Current treatment is surgical (goniotomy/trabeculotomy).  "
            "Gene therapy rationale: (1) AAV-mediated CYP1B1 gene "
            "replacement in trabecular meshwork (cDNA ~1.6 kb, fits AAV); "
            "intracameral delivery; would need to be given in neonatal "
            "period during anterior segment development; "
            "(2) Base editing for founder mutations (G61E is G>A, ideal "
            "ABE target); "
            "(3) Major challenge: CYP1B1 function is most critical during "
            "fetal development, so postnatal gene replacement may have "
            "limited benefit.  "
            "No clinical trials as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Standard care is surgery "
            "(goniotomy, trabeculotomy, trabeculectomy).  "
            "Gene therapy is a challenging target due to developmental "
            "window constraints."
        ),
        "conditions": ["primary_congenital_glaucoma", "PCG", "buphthalmos",
                        "ophthalmology", "glaucoma"],
    },

    # -------------------------------------------------------------------
    # 3f.  CNGA3 -- Cyclic Nucleotide-Gated Channel Alpha 3 (Achromatopsia)
    # -------------------------------------------------------------------
    "CNGA3": {
        "gene_id": 1261,
        "chrom": "chr2",
        "start": 98_346_456,
        "end": 98_398_601,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2q11.2",
        "exon_count": 9,
        "role": (
            "Cyclic nucleotide-gated channel subunit alpha 3 -- alpha "
            "subunit of the cone photoreceptor CNG channel.  Forms "
            "heteromeric channels with CNGB3 that mediate the cone "
            "phototransduction dark current.  Loss-of-function causes "
            "complete achromatopsia (total colour blindness, severe "
            "photophobia, nystagmus, ~20/200 visual acuity).  CNGA3 "
            "accounts for ~25% of achromatopsia cases."
        ),
        "disease": "Achromatopsia (ACHM2)",
        "omim_disease": 216900,
        "omim_gene": 600053,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg283Trp (R283W)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Pore region; abolishes channel conductance.",
            },
            {
                "name": "p.Arg427Cys (R427C)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "cGMP binding domain; impairs ligand gating.",
            },
        ],
        "mutational_landscape": (
            "Over 100 pathogenic variants.  Mix of missense, nonsense, "
            "and splice-site.  Founder mutations in specific populations "
            "(e.g., Pingelapese islanders carry CNGB3 mutations, not CNGA3)."
        ),
        "strategy": (
            "(1) AAV-mediated CNGA3 gene replacement in cone photoreceptors "
            "(cDNA ~2.1 kb, fits AAV); subretinal or intravitreal delivery; "
            "(2) STZ-CNGA3 (MeiraGTx/Janssen) Phase I/II clinical trial "
            "using AAV8-CNGA3 subretinal injection; "
            "(3) Tubingen/RD-CURE rAAV.CNGA3 Phase I/II; "
            "(4) Base editing for specific missense variants.  "
            "Achromatopsia is a leading gene therapy target because "
            "cone structure is preserved despite dysfunction."
        ),
        "clinical_programs": {
            "RD_CURE_CNGA3": {
                "nct": "NCT02610582",
                "sponsor": "University of Tubingen / STZ eyetrial",
                "type": "rAAV8.CNGA3 subretinal",
                "status": (
                    "Phase I/II; dose-escalation completed with acceptable "
                    "safety.  Evidence of improved cone function in some "
                    "patients (pupillometry, colour vision testing)."
                ),
            },
        },
        "conditions": ["achromatopsia", "ACHM2", "cone_dysfunction",
                        "colour_blindness", "ophthalmology",
                        "retinal_dystrophy"],
    },

    # -------------------------------------------------------------------
    # 3g.  CNGB3 -- Cyclic Nucleotide-Gated Channel Beta 3 (Achromatopsia)
    # -------------------------------------------------------------------
    "CNGB3": {
        "gene_id": 54714,
        "chrom": "chr8",
        "start": 86_574_179,
        "end": 86_743_634,
        "strand": "-",
        "refseq": "NC_000008.11",
        "cytoband": "8q21.3",
        "exon_count": 19,
        "role": (
            "Cyclic nucleotide-gated channel subunit beta 3 -- beta "
            "subunit of the cone photoreceptor CNG channel.  Modulates "
            "channel gating, ion selectivity, and ligand sensitivity.  "
            "Loss-of-function causes complete achromatopsia.  CNGB3 is "
            "the most common cause of achromatopsia (~50% of cases), "
            "with the p.Thr383fs founder mutation highly prevalent."
        ),
        "disease": "Achromatopsia (ACHM1)",
        "omim_disease": 262300,
        "omim_gene": 605080,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Thr383fs (c.1148delC)",
                "hgvs_coding": "NM_019098.5:c.1148delC",
                "rsid": "rs121912533",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "frequency": "Most common achromatopsia allele worldwide",
                "notes": (
                    "Single nucleotide deletion; accounts for ~70% of CNGB3 "
                    "mutations in Europeans and nearly 100% in Pingelapese "
                    "islanders (famous founder effect, 4-10% carrier rate)."
                ),
            },
            {
                "name": "p.Arg403Gln (R403Q)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Second most common in Europeans.",
            },
        ],
        "mutational_landscape": (
            "Over 80 pathogenic variants.  p.Thr383fs dominates (~70% of "
            "alleles in most populations).  Pingelapese atoll in Micronesia "
            "has the highest achromatopsia prevalence in the world due to "
            "CNGB3 founder effect (typhoon Lengkieki, ~1775)."
        ),
        "strategy": (
            "(1) AAV-mediated CNGB3 gene replacement in cone photoreceptors "
            "(cDNA ~2.4 kb, fits AAV); AGTC (now Beacon Therapeutics) "
            "Phase I/II clinical trial with intravitreal/subretinal "
            "AAV-CNGB3; "
            "(2) MeiraGTx/Janssen Phase I/II with AAV-CNGB3; "
            "(3) Base editing for the Thr383fs allele would require prime "
            "editing (deletion correction).  "
            "CNGB3 achromatopsia is one of the most advanced inherited "
            "retinal disease gene therapy programs."
        ),
        "clinical_programs": {
            "AGTC_CNGB3": {
                "nct": "NCT02599922",
                "sponsor": "AGTC / Beacon Therapeutics",
                "type": "AAV-CNGB3 subretinal",
                "status": (
                    "Phase I/II; dose-escalation completed.  Evidence of "
                    "improved light sensitivity and cone function in "
                    "paediatric cohort."
                ),
            },
        },
        "conditions": ["achromatopsia", "ACHM1", "cone_dysfunction",
                        "ophthalmology", "retinal_dystrophy"],
    },
}


# ============================================================================
# 4. NCL TARGETS
#    Neuronal ceroid lipofuscinoses (Batten disease family)
# ============================================================================

NCL_TARGETS = {

    # -------------------------------------------------------------------
    # 4a.  CLN3 -- CLN3, Battenin  (Juvenile NCL / CLN3 disease)
    # -------------------------------------------------------------------
    "CLN3": {
        "gene_id": 1201,
        "chrom": "chr16",
        "start": 28_466_653,
        "end": 28_492_082,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16p12.1",
        "exon_count": 16,
        "role": (
            "CLN3 / battenin -- lysosomal/endosomal transmembrane protein.  "
            "Exact function incompletely understood; implicated in endosomal/ "
            "lysosomal trafficking, autophagy regulation, and lipid "
            "metabolism.  Loss-of-function causes juvenile neuronal ceroid "
            "lipofuscinosis (JNCL, CLN3 disease / classic Batten disease), "
            "the most common NCL form.  Onset age 4-7 years with rapid "
            "vision loss, seizures, cognitive decline, motor deterioration, "
            "and premature death (typically 15-30 years)."
        ),
        "disease": "CLN3 disease / Juvenile NCL / Classic Batten disease",
        "omim_disease": 204200,
        "omim_gene": 607042,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "1 kb deletion (c.461-280_677+382del, ~966 bp)",
                "consequence": "deletion",
                "clinical_significance": "pathogenic",
                "frequency": "~85% of CLN3 disease alleles worldwide",
                "notes": (
                    "Common founder deletion removing exons 7-8 and flanking "
                    "intronic sequence.  Abolishes protein function.  "
                    "Homozygous in ~75% of patients; compound heterozygous "
                    "with another variant in ~20%."
                ),
            },
        ],
        "mutational_landscape": (
            "~60 pathogenic variants, but the 1 kb exon 7-8 deletion "
            "dominates (~85% of alleles).  Remaining variants include "
            "missense, nonsense, splice-site, and small indels throughout "
            "the gene.  Incidence ~1:100,000 births."
        ),
        "strategy": (
            "(1) AAV9-CLN3 gene replacement via intrathecal/intracerebroventricular "
            "delivery (cDNA ~1.3 kb, fits AAV).  Must achieve broad CNS "
            "distribution given diffuse neurodegeneration; "
            "(2) Weill Cornell / Nationwide Children's AAV9-CLN3 preclinical "
            "program; "
            "(3) Gene-modified haematopoietic stem cell transplant "
            "(cross-correction via secreted CLN3-containing exosomes -- "
            "limited evidence); "
            "(4) CRISPR excision repair: for the common 1 kb deletion, "
            "HDR-mediated reinsertion of exons 7-8 or micro-gene replacement "
            "approach; prime editing theoretically possible but large "
            "deletion correction is challenging.  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "No approved therapy.  AAV9-CLN3 gene therapy in preclinical "
            "development (multiple academic centres).  Cysteamine and "
            "N-acetylcysteine trialled with limited efficacy.  "
            "Miglustat investigated but no significant benefit."
        ),
        "conditions": ["juvenile_NCL", "CLN3_disease", "batten_disease",
                        "neuronal_ceroid_lipofuscinosis", "lysosomal"],
    },

    # -------------------------------------------------------------------
    # 4b.  TPP1 (CLN2) -- Tripeptidyl Peptidase 1  (Late Infantile NCL)
    # -------------------------------------------------------------------
    "TPP1": {
        "gene_id": 1200,
        "chrom": "chr11",
        "start": 6_612_768,
        "end": 6_619_422,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11p15.4",
        "exon_count": 13,
        "alias": "CLN2",
        "role": (
            "Tripeptidyl peptidase 1 -- soluble lysosomal serine protease "
            "that sequentially removes tripeptides from the N-terminus "
            "of small proteins/peptides.  Loss-of-function causes late "
            "infantile NCL (CLN2 disease), with onset age 2-4 years "
            "featuring seizures, ataxia, language regression, visual "
            "loss, and progressive motor/cognitive decline.  Untreated, "
            "death typically by age 6-12."
        ),
        "disease": "CLN2 disease / Late infantile NCL",
        "omim_disease": 204500,
        "omim_gene": 607998,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.509-1G>C (splice acceptor, intron 5)",
                "rsid": "rs766983415",
                "consequence": "splice_variant",
                "clinical_significance": "pathogenic",
                "frequency": "Most common CLN2 allele (~38% in Europeans)",
                "notes": "Abolishes exon 6 inclusion; no functional TPP1.",
            },
            {
                "name": "p.Arg208Ter (R208X)",
                "hgvs_coding": "NM_000391.4:c.622C>T",
                "rsid": "rs118204427",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "frequency": "Second most common (~19% of alleles)",
                "notes": "Premature termination; no enzyme activity.",
            },
        ],
        "mutational_landscape": (
            "Over 120 pathogenic variants.  Two mutations (c.509-1G>C and "
            "R208X) account for ~55% of alleles.  Most are null alleles.  "
            "Rare missense mutations with residual enzyme activity have "
            "later onset and slower progression."
        ),
        "strategy": (
            "(1) CERLIPONASE ALFA (BRINEURA) -- APPROVED by FDA (2017) and "
            "EMA (2017).  Recombinant TPP1 enzyme replacement therapy "
            "delivered by intracerebroventricular infusion every 2 weeks.  "
            "First approved treatment for any NCL.  Significantly slows "
            "motor and language decline.  "
            "(2) AAV gene therapy: AAV9-TPP1 intrathecal or ICV delivery "
            "(cDNA ~1.7 kb, fits AAV).  Weill Cornell phase I/II trial "
            "(NCT01414985/NCT01161576) with AAVrh.10-CLN2; "
            "(3) Base editing for c.509-1G>C (G>C at splice site -- CBE "
            "could correct to G>T/G>A to restore splicing, or ABE if "
            "targeting the antisense strand); "
            "(4) Intrathecal AAV9 for broader CNS distribution than "
            "Brineura's focal ICV delivery."
        ),
        "clinical_programs": {
            "Brineura": {
                "generic": "cerliponase alfa",
                "type": "Recombinant TPP1 enzyme replacement",
                "route": "intracerebroventricular (ICV) infusion",
                "approval": "FDA 2017, EMA 2017",
                "status": (
                    "APPROVED -- first-in-class ERT for NCL.  Biweekly ICV "
                    "infusion.  Clinical trials showed 2.1-point reduction "
                    "in rate of CLN2 Clinical Rating Scale decline per year "
                    "vs. natural history (p<0.001)."
                ),
            },
            "Weill_Cornell_AAVrh10": {
                "nct": "NCT01414985",
                "type": "AAVrh.10-CLN2 gene replacement",
                "route": "intracranial (12 burr hole sites)",
                "status": (
                    "Phase I/II completed; evidence of reduced rate of "
                    "decline in some patients.  Newer AAV9-CLN2 programs "
                    "with intrathecal delivery in development."
                ),
            },
        },
        "approved_therapy": "Cerliponase alfa (BRINEURA) -- ERT via ICV",
        "conditions": ["late_infantile_NCL", "CLN2_disease", "batten_disease",
                        "neuronal_ceroid_lipofuscinosis", "lysosomal"],
    },

    # -------------------------------------------------------------------
    # 4c.  PPT1 (CLN1) -- Palmitoyl-Protein Thioesterase 1 (Infantile NCL)
    # -------------------------------------------------------------------
    "PPT1": {
        "gene_id": 5538,
        "chrom": "chr1",
        "start": 40_071_461,
        "end": 40_097_252,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1p34.2",
        "exon_count": 9,
        "alias": "CLN1",
        "role": (
            "Palmitoyl-protein thioesterase 1 -- soluble lysosomal enzyme "
            "that removes thioester-linked fatty acyl groups (palmitate) "
            "from cysteine residues of proteins destined for degradation.  "
            "Loss-of-function causes infantile NCL (CLN1 disease), the "
            "most severe NCL form with onset at 6-18 months, rapid "
            "neurodegeneration, blindness, seizures, and death by age "
            "8-13 years.  Later-onset forms (juvenile, adult) occur with "
            "partial enzyme deficiency."
        ),
        "disease": "CLN1 disease / Infantile NCL",
        "omim_disease": 256730,
        "omim_gene": 600722,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg122Trp (R122W)",
                "hgvs_coding": "NM_000310.4:c.364C>T",
                "rsid": "rs28940893",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Finnish founder mutation; accounts for ~45% of CLN1 "
                    "alleles in Finland (NCL incidence 1:20,000 in Finland)."
                ),
            },
            {
                "name": "p.Thr75Pro (T75P)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Residual enzyme activity; associated with later onset.",
            },
        ],
        "mutational_landscape": (
            "Over 70 pathogenic variants.  Finnish founder R122W most "
            "common.  Truncating mutations cause classic infantile onset; "
            "missense with residual activity cause later-onset forms.  "
            "PPT1 is a soluble enzyme amenable to cross-correction "
            "(enzyme secreted from corrected cells can be taken up by "
            "uncorrected neighbours)."
        ),
        "strategy": (
            "(1) AAV-mediated PPT1 gene replacement via intrathecal/ICV "
            "delivery (cDNA ~0.9 kb, very small -- fits AAV).  Soluble "
            "enzyme allows cross-correction, increasing therapeutic reach; "
            "(2) ERT via ICV (analogous to Brineura for CLN2) -- "
            "preclinical studies of recombinant PPT1; "
            "(3) HSC gene therapy: ex vivo lentiviral transduction of "
            "CD34+ cells to create a CNS depot of PPT1-producing "
            "microglia (analogous to Libmeldy for MLD); "
            "(4) Base editing for R122W Finnish founder mutation (C>T, "
            "ideal ABE target).  "
            "No approved therapy; active preclinical programs."
        ),
        "clinical_programs": (
            "No approved therapy.  AAV gene therapy programs in preclinical "
            "development (Nationwide Children's, Weill Cornell).  "
            "Cysteamine/N-acetylcysteine combination trialled historically "
            "with limited benefit.  PPT1 ERT in preclinical studies."
        ),
        "conditions": ["infantile_NCL", "CLN1_disease", "batten_disease",
                        "neuronal_ceroid_lipofuscinosis", "lysosomal"],
    },

    # -------------------------------------------------------------------
    # 4d.  MFSD8 (CLN7) -- Major Facilitator Superfamily Domain 7
    # -------------------------------------------------------------------
    "MFSD8": {
        "gene_id": 256471,
        "chrom": "chr4",
        "start": 127_917_732,
        "end": 127_965_963,
        "strand": "-",
        "refseq": "NC_000004.12",
        "cytoband": "4q28.2",
        "exon_count": 13,
        "alias": "CLN7",
        "role": (
            "Major facilitator superfamily domain-containing protein 8 -- "
            "lysosomal membrane transporter.  Predicted 12-transmembrane "
            "domain protein of the MFS superfamily.  Exact substrate "
            "unknown; may transport small molecules across the lysosomal "
            "membrane.  Loss-of-function causes variant late-infantile NCL "
            "(CLN7 disease), with onset age 2-6 years and clinical "
            "features similar to CLN2 disease."
        ),
        "disease": "CLN7 disease / Variant late-infantile NCL",
        "omim_disease": 610951,
        "omim_gene": 611124,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Thr294Lys (T294K)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Common in Turkish consanguineous families.",
            },
            {
                "name": "c.863+1G>C (splice donor)",
                "consequence": "splice_variant",
                "clinical_significance": "pathogenic",
                "notes": "Abolishes normal splicing of exon 9.",
            },
        ],
        "mutational_landscape": (
            "Over 40 pathogenic variants.  Most cases in consanguineous "
            "families (Turkey, Pakistan, India, Roma population).  "
            "Mix of missense, splice-site, and truncating mutations."
        ),
        "strategy": (
            "(1) AAV9-MFSD8 gene replacement via intrathecal/ICV delivery "
            "(cDNA ~1.5 kb, fits AAV); "
            "(2) In 2019, a single-patient IND was used to treat a child "
            "with CLN7 disease with ASO (milasen) customised for her "
            "specific splice mutation (NEJM 2019, Kim et al.) -- "
            "demonstrated the N-of-1 therapeutic paradigm; "
            "(3) Base/prime editing for specific variants.  "
            "No standard CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "Milasen (single-patient ASO) -- landmark personalised medicine "
            "case.  No commercial gene therapy programs.  AAV9-MFSD8 in "
            "preclinical studies."
        ),
        "conditions": ["variant_late_infantile_NCL", "CLN7_disease",
                        "batten_disease", "neuronal_ceroid_lipofuscinosis",
                        "lysosomal"],
    },

    # -------------------------------------------------------------------
    # 4e.  CLN5 -- CLN5, Intracellular Protein
    # -------------------------------------------------------------------
    "CLN5": {
        "gene_id": 1203,
        "chrom": "chr13",
        "start": 76_992_081,
        "end": 77_005_117,
        "strand": "+",
        "refseq": "NC_000013.11",
        "cytoband": "13q22.3",
        "exon_count": 5,
        "role": (
            "CLN5 -- soluble lysosomal glycoprotein.  Recently identified "
            "as a lysosomal bis(monoacylglycero)phosphate (BMP) synthase.  "
            "BMP is a lipid unique to lysosomes that regulates cholesterol "
            "and sphingolipid metabolism.  Loss-of-function causes Finnish "
            "variant late-infantile NCL (CLN5 disease), with onset age "
            "4-7 years and slower progression than CLN2."
        ),
        "disease": "CLN5 disease / Finnish variant late-infantile NCL",
        "omim_disease": 256731,
        "omim_gene": 608102,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1175delAT (p.Tyr392X)",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "frequency": "Finnish founder mutation (~90% of Finnish alleles)",
                "notes": (
                    "Frameshift leading to premature termination.  NCL "
                    "incidence ~1:20,000 in Finland."
                ),
            },
        ],
        "mutational_landscape": (
            "~40 pathogenic variants.  Strong Finnish founder effect.  "
            "Also reported in Colombia, Portugal, Netherlands, Pakistan.  "
            "CLN5 is a soluble enzyme, so cross-correction is possible."
        ),
        "strategy": (
            "(1) AAV9-CLN5 gene replacement via intrathecal/ICV delivery "
            "(cDNA ~1.2 kb, small -- fits AAV); soluble enzyme enables "
            "cross-correction; "
            "(2) ERT with recombinant CLN5 (preclinical); "
            "(3) HSC gene therapy for CNS enzyme depot.  "
            "No clinical trials as of 2026."
        ),
        "clinical_programs": (
            "No approved therapy.  AAV gene therapy in preclinical "
            "development.  Lincoln Laboratory / Batten Disease Support "
            "and Research Association funding preclinical studies."
        ),
        "conditions": ["CLN5_disease", "finnish_variant_NCL",
                        "batten_disease", "neuronal_ceroid_lipofuscinosis",
                        "lysosomal"],
    },

    # -------------------------------------------------------------------
    # 4f.  CLN6 -- CLN6 Transmembrane ER Protein
    # -------------------------------------------------------------------
    "CLN6": {
        "gene_id": 54982,
        "chrom": "chr15",
        "start": 68_206_992,
        "end": 68_257_211,
        "strand": "-",
        "refseq": "NC_000015.10",
        "cytoband": "15q23",
        "exon_count": 8,
        "role": (
            "CLN6 -- endoplasmic reticulum (ER) transmembrane protein.  "
            "Functions as a receptor/adaptor for lysosomal enzyme sorting "
            "at the ER, facilitating ER-to-Golgi transport of lysosomal "
            "enzymes.  Loss-of-function causes variant late-infantile NCL "
            "(CLN6 disease) and adult-onset NCL (Kufs disease type A).  "
            "Because CLN6 is an ER-resident protein (not lysosomal), "
            "cross-correction by secreted enzyme is NOT possible -- "
            "each neuron must be individually corrected."
        ),
        "disease": "CLN6 disease / Variant late-infantile NCL / Kufs type A",
        "omim_disease": 601780,
        "omim_gene": 606725,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg252His (R252H)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Recurrent; affects luminal loop.",
            },
        ],
        "mutational_landscape": (
            "Over 70 pathogenic variants.  Widely distributed across "
            "different populations (South Asian, European, Latin American, "
            "Middle Eastern).  Includes missense, nonsense, and splice-site."
        ),
        "strategy": (
            "(1) AAV9-CLN6 gene replacement via intrathecal delivery "
            "(cDNA ~0.9 kb, very small -- fits AAV).  Nationwide Children's "
            "Hospital / Amicus Therapeutics AT-GTX-501 Phase I/II "
            "intrathecal AAV9-CLN6 trial; "
            "(2) Challenges: ER-resident protein means no cross-correction; "
            "every affected neuron needs transduction; intrathecal AAV9 "
            "achieves broad but incomplete CNS coverage; "
            "(3) Gene editing approaches theoretically possible but face "
            "the same delivery challenges.  "
            "Note: CLN6 was the first NCL gene to enter AAV gene therapy "
            "clinical trials."
        ),
        "clinical_programs": {
            "AT_GTX_501_Amicus": {
                "nct": "NCT02725580",
                "sponsor": "Nationwide Children's / Amicus Therapeutics",
                "type": "scAAV9.CB.CLN6 intrathecal",
                "status": (
                    "Phase I/II; initial data showed slower disease "
                    "progression in some patients vs. natural history.  "
                    "Longer follow-up shows variable durability.  Amicus "
                    "continues development as AT-GTX-501."
                ),
            },
        },
        "conditions": ["CLN6_disease", "variant_late_infantile_NCL",
                        "kufs_disease", "batten_disease",
                        "neuronal_ceroid_lipofuscinosis", "lysosomal"],
    },

    # -------------------------------------------------------------------
    # 4g.  CLN8 -- CLN8 Transmembrane ER/ERGIC Protein
    # -------------------------------------------------------------------
    "CLN8": {
        "gene_id": 2055,
        "chrom": "chr8",
        "start": 1_753_059,
        "end": 1_786_570,
        "strand": "+",
        "refseq": "NC_000008.11",
        "cytoband": "8p23.3",
        "exon_count": 6,
        "role": (
            "CLN8 -- ER-to-ERGIC (ER-Golgi intermediate compartment) "
            "cargo receptor.  Works in partnership with CLN6 to sort "
            "lysosomal enzymes at the ER for Golgi trafficking.  "
            "Loss-of-function causes variant late-infantile NCL (CLN8 "
            "disease) and Northern epilepsy (progressive epilepsy with "
            "cognitive decline, a milder allelic disorder found in Finland).  "
            "Like CLN6, CLN8 is ER-resident and NOT amenable to "
            "cross-correction."
        ),
        "disease": "CLN8 disease / Northern epilepsy",
        "omim_disease": 600143,
        "omim_gene": 607837,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg24Gly (R24G)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Finnish founder mutation causing Northern epilepsy "
                    "(milder phenotype).  Carrier frequency ~1:135 in "
                    "Finland."
                ),
            },
            {
                "name": "p.Trp263Cys (W263C)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Turkish founder mutation; severe late-infantile NCL.",
            },
        ],
        "mutational_landscape": (
            "~30 pathogenic variants.  Finnish (R24G) and Turkish (W263C) "
            "founder mutations most common.  Genotype-phenotype correlation: "
            "R24G homozygotes have Northern epilepsy (slow progression); "
            "other mutations cause more typical NCL."
        ),
        "strategy": (
            "(1) AAV9-CLN8 gene replacement via intrathecal delivery "
            "(cDNA ~0.9 kb, very small -- fits AAV); Nationwide Children's "
            "Hospital / Amicus Therapeutics developing AT-GTX-502; "
            "(2) Same cross-correction limitation as CLN6 -- ER-resident "
            "protein requires transduction of individual neurons; "
            "(3) Base editing for R24G (G>A, ideal ABE target) could "
            "benefit Finnish Northern epilepsy patients specifically."
        ),
        "clinical_programs": {
            "AT_GTX_502_Amicus": {
                "nct": "NCT04737460",
                "sponsor": "Nationwide Children's / Amicus Therapeutics",
                "type": "scAAV9.CLN8 intrathecal",
                "status": (
                    "Phase I/II; enrolling.  Preclinical data in CLN8 mice "
                    "showed extended survival and preservation of motor "
                    "function."
                ),
            },
        },
        "conditions": ["CLN8_disease", "northern_epilepsy",
                        "batten_disease", "neuronal_ceroid_lipofuscinosis",
                        "lysosomal"],
    },
}


# ============================================================================
# 5. BONE MARROW FAILURE TARGETS
#    (ELANE and TERT/TERC already covered in other databases -- see notes)
# ============================================================================

BONE_MARROW_FAILURE_TARGETS = {

    # -------------------------------------------------------------------
    # 5a.  RUNX1 -- Runt-Related Transcription Factor 1
    # -------------------------------------------------------------------
    "RUNX1": {
        "gene_id": 861,
        "chrom": "chr21",
        "start": 34_787_801,
        "end": 35_049_302,
        "strand": "-",
        "refseq": "NC_000021.9",
        "cytoband": "21q22.12",
        "exon_count": 12,
        "role": (
            "Runt-related transcription factor 1 (AML1/CBFA2) -- master "
            "regulator of definitive haematopoiesis.  Forms a heterodimer "
            "with CBFB to activate target genes essential for HSC "
            "emergence, megakaryopoiesis, and myeloid/lymphoid "
            "differentiation.  Germline loss-of-function causes familial "
            "platelet disorder with predisposition to acute myeloid "
            "leukaemia (FPD/AML).  Somatic RUNX1 mutations are among "
            "the most common in AML and MDS."
        ),
        "disease": "Familial Platelet Disorder with AML predisposition (FPD/AML)",
        "omim_disease": 601399,
        "omim_gene": 151385,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Arg201Gln (R201Q)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Runt homology domain (RHD); impairs DNA binding.  "
                    "Common recurrent germline mutation."
                ),
            },
            {
                "name": "p.Arg166Gln (R166Q)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "RHD; reduces transcriptional activity.",
            },
            {
                "name": "Large intragenic deletions",
                "consequence": "deletion",
                "clinical_significance": "pathogenic",
                "notes": "~15% of FPD/AML families; may be missed by sequencing.",
            },
        ],
        "mutational_landscape": (
            "Over 80 germline pathogenic variants.  ~60% missense/nonsense "
            "in RHD (exons 3-5), ~15% frameshifts, ~15% deletions.  "
            "Lifetime AML risk ~35-40%.  MDS/AML arises through second-hit "
            "mutations (often the same RUNX1 allele or cooperating "
            "mutations in ASXL1, BCOR, PHF6, RAS pathway)."
        ),
        "strategy": (
            "(1) Allogeneic HSCT is currently the only curative option and "
            "is recommended before AML transformation, especially with "
            "matched sibling donor who does not carry the familial variant; "
            "(2) Ex vivo gene correction of autologous HSPCs: base editing "
            "for missense RHD mutations (e.g., R201Q is G>A, ideal ABE "
            "target) followed by autologous transplant to eliminate the "
            "inherited AML predisposition; "
            "(3) CRISPRi to silence the mutant allele if dominant-negative; "
            "(4) Gene addition of wild-type RUNX1 cDNA (~1.5 kb, fits AAV "
            "or lentiviral) -- but haploinsufficiency is the mechanism, so "
            "gene addition may be sufficient.  "
            "No CRISPR clinical trials as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Allogeneic HSCT for AML prevention "
            "in high-risk families.  RUNX1 germline testing increasingly "
            "incorporated into AML workup."
        ),
        "conditions": ["familial_platelet_disorder", "FPD_AML",
                        "bone_marrow_failure", "AML_predisposition",
                        "thrombocytopenia"],
    },

    # -------------------------------------------------------------------
    # 5b.  GATA2 -- GATA Binding Protein 2
    # -------------------------------------------------------------------
    "GATA2": {
        "gene_id": 2624,
        "chrom": "chr3",
        "start": 128_479_422,
        "end": 128_493_201,
        "strand": "-",
        "refseq": "NC_000003.12",
        "cytoband": "3q21.3",
        "exon_count": 8,
        "role": (
            "GATA binding protein 2 -- zinc-finger transcription factor "
            "essential for HSC maintenance, lymphatic vessel development, "
            "and myeloid/lymphoid differentiation.  GATA2 haploinsufficiency "
            "causes a complex immunodeficiency syndrome with monocytopenia, "
            "B/NK cell deficiency, pulmonary alveolar proteinosis, "
            "lymphoedema, and high risk of MDS/AML (~75% lifetime risk).  "
            "Also causes Emberger syndrome (lymphoedema + MDS) and "
            "MonoMAC syndrome."
        ),
        "disease": "GATA2 deficiency / MonoMAC / Emberger syndrome",
        "omim_disease": 614172,
        "omim_gene": 137295,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "p.Thr354Met (T354M)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Zinc finger 2 (ZF2); most common germline GATA2 "
                    "missense.  Impairs DNA binding."
                ),
            },
            {
                "name": "p.Arg396Gln (R396Q)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "ZF2; impairs transcriptional activation.",
            },
            {
                "name": "Intron 4 enhancer deletions/mutations",
                "consequence": "regulatory",
                "clinical_significance": "pathogenic",
                "notes": (
                    "A distal 9.5 kb intron-4 enhancer regulates GATA2 "
                    "expression in HSCs.  Deletions or point mutations in "
                    "this enhancer cause GATA2 deficiency without coding "
                    "sequence mutations."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 100 germline pathogenic variants.  ~40% in ZF2 (exon 5), "
            "~20% frameshifts, ~10% whole-gene or partial deletions, "
            "~10% regulatory (intron 4 enhancer).  MDS/AML risk ~75%.  "
            "Median age of MDS onset ~20 years."
        ),
        "strategy": (
            "(1) Allogeneic HSCT is the standard of care and only cure, "
            "recommended before MDS transformation.  Reduced-intensity "
            "conditioning preferred given infectious complications.  "
            "(2) Ex vivo gene correction of autologous HSPCs: CRISPR base "
            "editing for ZF2 missense mutations or gene addition of "
            "wild-type GATA2 (cDNA ~1.5 kb, fits AAV/lentiviral) under "
            "regulated promoter to avoid overexpression; "
            "(3) CRISPRa to upregulate the intact GATA2 allele -- "
            "requires careful dosage control; "
            "(4) For intron-4 enhancer mutations: CRISPR-mediated enhancer "
            "repair or replacement.  "
            "No gene therapy trials as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Allogeneic HSCT outcomes improving "
            "with reduced-intensity conditioning protocols.  NIH Natural "
            "History Study (NCT01905826) characterising disease progression."
        ),
        "conditions": ["GATA2_deficiency", "MonoMAC", "Emberger_syndrome",
                        "bone_marrow_failure", "MDS_predisposition",
                        "immunodeficiency"],
    },

    # -------------------------------------------------------------------
    # 5c.  MPL -- Thrombopoietin Receptor
    # -------------------------------------------------------------------
    "MPL": {
        "gene_id": 4352,
        "chrom": "chr1",
        "start": 43_337_818,
        "end": 43_354_466,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p34.2",
        "exon_count": 12,
        "role": (
            "MPL proto-oncogene / thrombopoietin receptor -- type I "
            "cytokine receptor that binds thrombopoietin (TPO).  "
            "Activates JAK2-STAT signalling.  Essential for HSC "
            "maintenance, megakaryocyte differentiation, and platelet "
            "production.  Biallelic loss-of-function causes congenital "
            "amegakaryocytic thrombocytopenia (CAMT), presenting with "
            "severe thrombocytopenia at birth that progresses to "
            "pancytopenia/aplastic anaemia due to HSC exhaustion.  "
            "Somatic gain-of-function (W515L/K) drives myeloproliferative "
            "neoplasms."
        ),
        "disease": "Congenital Amegakaryocytic Thrombocytopenia (CAMT)",
        "omim_disease": 604498,
        "omim_gene": 159530,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg102Pro (R102P)",
                "hgvs_coding": "NM_005373.3:c.305G>C",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "CAMT type II -- partial loss of signalling; later "
                    "progression to aplastic anaemia.  Transient response "
                    "to eltrombopag possible."
                ),
            },
            {
                "name": "c.79+2T>C (splice donor, intron 2)",
                "consequence": "splice_variant",
                "clinical_significance": "pathogenic",
                "notes": (
                    "CAMT type I -- complete loss of MPL; severe neonatal "
                    "thrombocytopenia with early pancytopenia."
                ),
            },
        ],
        "mutational_landscape": (
            "~30 pathogenic variants.  CAMT type I (null mutations) has "
            "more severe course with earlier aplastic anaemia (~2-5 years) "
            "vs. CAMT type II (hypomorphic, ~5-10 years).  Very rare "
            "(~1:1,000,000 births)."
        ),
        "strategy": (
            "(1) Allogeneic HSCT is the only cure -- recommended early "
            "before aplastic anaemia; "
            "(2) Ex vivo lentiviral gene addition of MPL to autologous "
            "HSPCs (cDNA ~1.9 kb, fits lentiviral/AAV): could avoid "
            "allogeneic transplant risks; "
            "(3) TPO receptor agonists (eltrombopag, romiplostim) may "
            "provide temporary benefit in CAMT type II with residual MPL "
            "function, but are ineffective in type I (no receptor); "
            "(4) Base editing for specific missense mutations (R102P is "
            "G>C -- CBE could correct to G>T but not ideal; prime editing "
            "needed for precise correction).  "
            "No gene therapy trials as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Allogeneic HSCT with matched donor "
            "is standard.  Eltrombopag used as bridge therapy.  "
            "Lentiviral MPL gene addition in autologous HSPCs in "
            "preclinical development (multiple academic groups)."
        ),
        "conditions": ["congenital_amegakaryocytic_thrombocytopenia",
                        "CAMT", "bone_marrow_failure",
                        "aplastic_anaemia", "thrombocytopenia"],
    },

    # -------------------------------------------------------------------
    # 5d.  SBDS -- Shwachman-Bodian-Diamond Syndrome Protein
    # -------------------------------------------------------------------
    "SBDS": {
        "gene_id": 51119,
        "chrom": "chr7",
        "start": 66_987_680,
        "end": 66_995_586,
        "strand": "-",
        "refseq": "NC_000007.14",
        "cytoband": "7q11.21",
        "exon_count": 5,
        "role": (
            "SBDS ribosome maturation factor -- essential for 60S ribosomal "
            "subunit maturation.  Works with EFL1 GTPase to release eIF6 "
            "from the pre-60S particle, enabling 80S ribosome assembly.  "
            "Loss-of-function causes Shwachman-Diamond syndrome (SDS), "
            "characterised by exocrine pancreatic insufficiency (lipomatous "
            "replacement of acinar cells), bone marrow failure (neutropenia "
            "most common, also anaemia and thrombocytopenia), skeletal "
            "abnormalities, and predisposition to MDS/AML.  SDS is the "
            "third most common inherited bone marrow failure syndrome "
            "after Fanconi anaemia and Diamond-Blackfan anaemia."
        ),
        "disease": "Shwachman-Diamond Syndrome (SDS)",
        "omim_disease": 260400,
        "omim_gene": 607444,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.183_184TA>CT (p.Lys62Ter, K62X)",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "frequency": "~75% of SDS alleles",
                "notes": (
                    "Result of gene conversion from the adjacent SBDSP1 "
                    "pseudogene.  Creates a premature stop codon.  Almost "
                    "always compound heterozygous with c.258+2T>C."
                ),
            },
            {
                "name": "c.258+2T>C (splice donor, intron 2)",
                "consequence": "splice_variant",
                "clinical_significance": "pathogenic",
                "frequency": "~60% of SDS alleles",
                "notes": (
                    "Also derived from SBDSP1 pseudogene conversion.  "
                    "K62X + 258+2T>C compound het accounts for ~60% of "
                    "SDS patients."
                ),
            },
        ],
        "mutational_landscape": (
            "~90% of SBDS mutations arise from gene conversion events "
            "with the adjacent SBDSP1 pseudogene.  Two pseudogene-derived "
            "mutations (K62X and 258+2T>C) account for ~75% and ~60% of "
            "alleles respectively.  Biallelic null is likely embryonic "
            "lethal -- all patients retain some residual SBDS function.  "
            "MDS/AML risk ~15-30% by age 30."
        ),
        "strategy": (
            "(1) Allogeneic HSCT for severe bone marrow failure or MDS.  "
            "Reduced-intensity conditioning recommended due to increased "
            "treatment-related toxicity in SDS; "
            "(2) Ex vivo lentiviral gene addition of SBDS to autologous "
            "HSPCs (cDNA ~0.8 kb, very small -- fits any vector); "
            "(3) Base editing for K62X (A>C at stop codon position -- "
            "CBE may not directly correct; ABE or prime editing for "
            "precise reversion); CRISPR disruption of the pseudogene "
            "SBDSP1 to prevent ongoing gene conversion; "
            "(4) Pancreatic enzyme supplementation (PERT) for exocrine "
            "insufficiency is standard supportive care.  "
            "No gene therapy trials as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  PERT for pancreatic insufficiency.  "
            "G-CSF for neutropenia.  Allogeneic HSCT for severe "
            "cytopenias or MDS.  Gene therapy in preclinical studies "
            "(Sick Kids Toronto, Boston Children's)."
        ),
        "conditions": ["shwachman_diamond_syndrome", "SDS",
                        "bone_marrow_failure", "MDS_predisposition",
                        "pancreatic_insufficiency", "ribosomopathy"],
    },

    # -------------------------------------------------------------------
    # Cross-references to genes covered in other databases:
    # -------------------------------------------------------------------
    # "ELANE": see additional_rare_disease_targets.py
    #     (congenital neutropenia / severe chronic neutropenia)
    # "TERT": see infectious_immune_aging_targets.py
    #     (dyskeratosis congenita / telomere biology disorders)
    # "TERC": see infectious_immune_aging_targets.py
    #     (dyskeratosis congenita / telomere biology disorders)
}


# ============================================================================
# 6. IMMUNE CHECKPOINT TARGETS
#    For CRISPR-based cancer immunotherapy (ex vivo T-cell editing)
# ============================================================================

IMMUNE_CHECKPOINT_TARGETS = {

    # -------------------------------------------------------------------
    # 6a.  PDCD1 -- Programmed Cell Death 1 (PD-1)
    # -------------------------------------------------------------------
    "PDCD1": {
        "gene_id": 5133,
        "chrom": "chr2",
        "start": 241_849_884,
        "end": 241_858_894,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q37.3",
        "exon_count": 6,
        "role": (
            "Programmed cell death 1 (PD-1) -- type I transmembrane "
            "inhibitory receptor on activated T cells, B cells, and "
            "myeloid cells.  Engagement by ligands PD-L1 (CD274) or "
            "PD-L2 (PDCD1LG2) recruits SHP-2 phosphatase, which "
            "dephosphorylates TCR signalling molecules and attenuates "
            "T-cell effector function.  Tumour cells exploit PD-1/PD-L1 "
            "axis for immune evasion.  CRISPR knockout of PD-1 in "
            "adoptive T cells enhances anti-tumour activity."
        ),
        "disease": "Cancer immunotherapy target (not a disease gene)",
        "omim_gene": 600244,
        "inheritance": "N/A (somatic editing for immunotherapy)",
        "key_variants": [],
        "strategy": (
            "(1) CRISPR-Cas9 knockout of PDCD1 in patient-derived tumour-"
            "infiltrating lymphocytes (TILs) or CAR-T cells to prevent "
            "checkpoint-mediated exhaustion.  First-in-human CRISPR trial "
            "used PDCD1-edited T cells (Lu et al., Nature Medicine 2020); "
            "(2) Simultaneous knockout of PDCD1 + TRAC (to create "
            "universal allogeneic CAR-T cells, e.g., CTX110 / Caribou); "
            "(3) Base editing to disrupt PD-1 without DSBs (reducing "
            "translocation risk); "
            "(4) Antibody therapies remain standard: pembrolizumab "
            "(Keytruda), nivolumab (Opdivo), cemiplimab (Libtayo).  "
            "CRISPR editing offers potential for durable T-cell-intrinsic "
            "checkpoint resistance."
        ),
        "clinical_programs": {
            "first_in_human_CRISPR": {
                "ref": "Lu et al., Nature Medicine 2020; NCT02793856",
                "sponsor": "Sichuan University / Carl June",
                "type": "CRISPR-Cas9 PD-1 knockout autologous T cells",
                "status": (
                    "Phase I completed (12 patients, NSCLC); safety "
                    "established.  Edited T cells persisted in vivo.  "
                    "Modest clinical activity (1 PR, 2 SD)."
                ),
            },
            "NYCE_Penn": {
                "nct": "NCT03399448",
                "sponsor": "University of Pennsylvania",
                "type": "CRISPR KO of PDCD1 + TRAC + TRBC in NY-ESO-1 TCR-T cells",
                "status": (
                    "Phase I; triple-edited T cells safely infused.  "
                    "Demonstrated feasibility of multiplex CRISPR editing "
                    "in clinical setting."
                ),
            },
        },
        "conditions": ["cancer_immunotherapy", "immune_checkpoint",
                        "PD1", "T_cell_engineering", "CAR_T"],
    },

    # -------------------------------------------------------------------
    # 6b.  CD274 -- PD-L1 (Programmed Death Ligand 1)
    # -------------------------------------------------------------------
    "CD274": {
        "gene_id": 29126,
        "chrom": "chr9",
        "start": 5_450_542,
        "end": 5_470_554,
        "strand": "+",
        "refseq": "NC_000009.12",
        "cytoband": "9p24.1",
        "exon_count": 7,
        "role": (
            "CD274 / PD-L1 / B7-H1 -- type I transmembrane protein and "
            "primary ligand for PD-1.  Expressed on tumour cells, APCs, "
            "and various tissue cells.  Tumour PD-L1 expression is a "
            "biomarker for anti-PD-1/PD-L1 antibody response.  Genomic "
            "amplification of 9p24.1 (CD274/PDCD1LG2/JAK2) is common "
            "in Hodgkin lymphoma and mediastinal B-cell lymphoma.  "
            "Anti-PD-L1 antibodies: atezolizumab (Tecentriq), "
            "durvalumab (Imfinzi), avelumab (Bavencio)."
        ),
        "disease": "Cancer immunotherapy target",
        "omim_gene": 605402,
        "inheritance": "N/A (somatic editing / tumour biology)",
        "key_variants": [
            {
                "name": "9p24.1 amplification",
                "consequence": "copy_number_gain",
                "clinical_significance": "predictive",
                "notes": (
                    "Amplification at 9p24.1 (CD274/PDCD1LG2/JAK2) found "
                    "in ~40% of classic Hodgkin lymphoma and ~60% of "
                    "primary mediastinal B-cell lymphoma.  Predicts high "
                    "PD-L1 expression and excellent response to "
                    "anti-PD-1 therapy."
                ),
            },
        ],
        "strategy": (
            "(1) CRISPR knockout of CD274 on tumour cells to study immune "
            "evasion mechanisms (research tool); "
            "(2) CRISPRa upregulation of PD-L1 in regulatory T cells or "
            "transplanted tissues to promote immune tolerance; "
            "(3) Structural variants/amplification at 9p24.1 as a "
            "biomarker for patient selection for anti-PD-1/PD-L1 "
            "immunotherapy; "
            "(4) CRISPR screens identify CD274 as essential for tumour "
            "immune evasion, validating therapeutic targeting.  "
            "Anti-PD-L1 antibodies approved across >15 cancer types."
        ),
        "clinical_programs": (
            "Anti-PD-L1 antibodies approved: atezolizumab, durvalumab, "
            "avelumab.  CRISPR editing of CD274 primarily in research "
            "setting.  No clinical trials of CD274-edited cell therapies."
        ),
        "conditions": ["cancer_immunotherapy", "immune_checkpoint",
                        "PDL1", "biomarker", "Hodgkin_lymphoma"],
    },

    # -------------------------------------------------------------------
    # 6c.  CTLA4 -- reference to existing database
    # -------------------------------------------------------------------
    # CTLA4 is covered in infectious_immune_aging_targets.py
    # as both an immune checkpoint and a haploinsufficiency disease gene.
    # Coordinates: chr2:203,867,788-203,873,960 (GRCh38)
    # See that database for full entry.

    # -------------------------------------------------------------------
    # 6d.  LAG3 -- Lymphocyte Activation Gene 3
    # -------------------------------------------------------------------
    "LAG3": {
        "gene_id": 3902,
        "chrom": "chr12",
        "start": 6_772_520,
        "end": 6_778_455,
        "strand": "+",
        "refseq": "NC_000012.12",
        "cytoband": "12p13.31",
        "exon_count": 8,
        "role": (
            "Lymphocyte activation gene 3 (LAG-3, CD223) -- inhibitory "
            "receptor on activated T cells and NK cells.  Binds MHC class "
            "II with higher affinity than CD4, as well as FGL1 (fibrinogen-"
            "like protein 1).  Promotes T-cell exhaustion in the tumour "
            "microenvironment, often co-expressed with PD-1 on "
            "dysfunctional TILs.  LAG-3 and PD-1 have synergistic "
            "inhibitory effects -- dual blockade produces superior "
            "anti-tumour responses."
        ),
        "disease": "Cancer immunotherapy target",
        "omim_gene": 153337,
        "inheritance": "N/A (immunotherapy target)",
        "key_variants": [],
        "strategy": (
            "(1) Anti-LAG-3 antibody relatlimab (Opdualag, combined with "
            "nivolumab) -- APPROVED by FDA (2022) for advanced melanoma.  "
            "First LAG-3-targeting therapy; "
            "(2) CRISPR knockout of LAG3 in CAR-T cells or TILs to "
            "prevent exhaustion -- preclinical studies show enhanced "
            "persistence and anti-tumour function; "
            "(3) Dual PDCD1 + LAG3 knockout in CAR-T cells for maximal "
            "checkpoint resistance; "
            "(4) Base editing to disrupt LAG3 without DSBs.  "
            "Multiple anti-LAG-3 antibodies in clinical development."
        ),
        "clinical_programs": {
            "Opdualag": {
                "generic": "relatlimab + nivolumab",
                "approval": "FDA 2022 for advanced melanoma",
                "status": (
                    "RELATIVITY-047 Phase III: relatlimab + nivolumab "
                    "doubled PFS vs. nivolumab alone (10.1 vs. 4.6 mo) "
                    "in previously untreated advanced melanoma."
                ),
            },
        },
        "approved_therapy": "Relatlimab + nivolumab (Opdualag) -- anti-LAG-3 combo",
        "conditions": ["cancer_immunotherapy", "immune_checkpoint",
                        "LAG3", "T_cell_exhaustion", "melanoma"],
    },

    # -------------------------------------------------------------------
    # 6e.  HAVCR2 -- TIM-3 (T-cell Immunoglobulin and Mucin Domain 3)
    # -------------------------------------------------------------------
    "HAVCR2": {
        "gene_id": 84868,
        "chrom": "chr5",
        "start": 157_085_832,
        "end": 157_109_044,
        "strand": "-",
        "refseq": "NC_000005.10",
        "cytoband": "5q33.3",
        "exon_count": 7,
        "role": (
            "Hepatitis A virus cellular receptor 2 / TIM-3 -- type I "
            "transmembrane protein of the immunoglobulin superfamily.  "
            "Inhibitory receptor on T cells (Th1, CD8+, Tregs), NK cells, "
            "and myeloid cells.  Ligands include galectin-9, CEACAM1, "
            "HMGB1, and phosphatidylserine.  Co-expressed with PD-1 "
            "on terminally exhausted T cells in the tumour microenvironment.  "
            "Germline loss-of-function variants cause subcutaneous "
            "panniculitis-like T-cell lymphoma (SPTCL) with hemophagocytic "
            "lymphohistiocytosis."
        ),
        "disease": "Cancer immunotherapy target / SPTCL-HLH (germline LOF)",
        "omim_gene": 606652,
        "inheritance": "N/A (immunotherapy target); AR for SPTCL-HLH",
        "key_variants": [
            {
                "name": "p.Tyr82Cys (Y82C)",
                "rsid": "rs147827860",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Germline variant in IgV domain; associated with "
                    "subcutaneous panniculitis-like T-cell lymphoma and "
                    "HLH.  Enriched in East Asian and Polynesian populations."
                ),
            },
        ],
        "strategy": (
            "(1) Anti-TIM-3 antibodies in clinical trials: sabatolimab "
            "(MBG453, Novartis), cobolimab (GSK), LY3321367 (Lilly); "
            "(2) CRISPR knockout of HAVCR2 in CAR-T cells to prevent "
            "terminal exhaustion -- triple knockout (PD-1 + LAG-3 + TIM-3) "
            "is a major research goal; "
            "(3) Anti-TIM-3 particularly promising in AML/MDS where "
            "TIM-3 is expressed on leukemic stem cells (non-T-cell role); "
            "(4) For germline SPTCL-HLH: HSCT or potentially gene "
            "correction of autologous T cells."
        ),
        "clinical_programs": {
            "sabatolimab": {
                "sponsor": "Novartis",
                "type": "Anti-TIM-3 monoclonal antibody",
                "status": (
                    "Phase III STIMULUS trials in MDS.  Phase II in AML "
                    "combination with decitabine.  Most advanced TIM-3 "
                    "clinical program."
                ),
            },
        },
        "conditions": ["cancer_immunotherapy", "immune_checkpoint",
                        "TIM3", "T_cell_exhaustion", "AML", "MDS",
                        "SPTCL_HLH"],
    },

    # -------------------------------------------------------------------
    # 6f.  TIGIT -- T-cell Immunoreceptor with Ig and ITIM Domains
    # -------------------------------------------------------------------
    "TIGIT": {
        "gene_id": 201633,
        "chrom": "chr3",
        "start": 114_294_028,
        "end": 114_310_288,
        "strand": "+",
        "refseq": "NC_000003.12",
        "cytoband": "3q13.31",
        "exon_count": 5,
        "role": (
            "TIGIT -- inhibitory receptor on T cells and NK cells.  "
            "Competes with the activating receptor DNAM-1 (CD226) for "
            "binding to nectin ligands (PVR/CD155, PVRL2/CD112) on "
            "tumour and antigen-presenting cells.  TIGIT engagement "
            "delivers inhibitory signals via cytoplasmic ITIM domain "
            "and disrupts DNAM-1 homodimerisation.  Highly expressed "
            "on tumour-infiltrating Tregs and exhausted CD8+ T cells.  "
            "Anti-TIGIT combinations with anti-PD-1 are a major area "
            "of clinical investigation."
        ),
        "disease": "Cancer immunotherapy target",
        "omim_gene": 612859,
        "inheritance": "N/A (immunotherapy target)",
        "key_variants": [],
        "strategy": (
            "(1) Anti-TIGIT antibodies: tiragolumab (Roche), vibostolimab "
            "(Merck), domvanalimab (Arcus/Gilead), ociperlimab (BeiGene); "
            "(2) CRISPR knockout of TIGIT in CAR-T/TIL products to "
            "enhance NK cell and T-cell anti-tumour function; "
            "(3) Combination strategies: anti-TIGIT + anti-PD-L1/PD-1 "
            "showed promising Phase II results (CITYSCAPE trial: "
            "tiragolumab + atezolizumab in PD-L1-high NSCLC); "
            "(4) Phase III results have been mixed (SKYSCRAPER-01 missed "
            "primary endpoint in all-comer NSCLC), suggesting biomarker "
            "selection will be critical.  "
            "CRISPR-edited TIGIT-knockout T cells in preclinical studies."
        ),
        "clinical_programs": {
            "tiragolumab_Roche": {
                "status": (
                    "Phase III SKYSCRAPER program.  SKYSCRAPER-01 in NSCLC "
                    "missed co-primary PFS endpoint.  SKYSCRAPER-02 in "
                    "ES-SCLC also negative.  Biomarker-selected subgroups "
                    "still under investigation."
                ),
            },
            "domvanalimab_Arcus": {
                "status": (
                    "Phase III ARC-7 in NSCLC; domvanalimab + zimberelimab "
                    "(anti-PD-1) showed improved PFS vs. zimberelimab alone "
                    "in PD-L1-high patients.  FDA Breakthrough Therapy "
                    "Designation granted."
                ),
            },
        },
        "conditions": ["cancer_immunotherapy", "immune_checkpoint",
                        "TIGIT", "T_cell_exhaustion", "NK_cell",
                        "NSCLC"],
    },
}


# ============================================================================
# 7. CRANIOFACIAL TARGETS
# ============================================================================

CRANIOFACIAL_TARGETS = {

    # -------------------------------------------------------------------
    # 7a.  TCOF1 -- Treacle Ribosome Biogenesis Factor 1
    # -------------------------------------------------------------------
    "TCOF1": {
        "gene_id": 6949,
        "chrom": "chr5",
        "start": 150_357_697,
        "end": 150_400_293,
        "strand": "+",
        "refseq": "NC_000005.10",
        "cytoband": "5q32-q33.1",
        "exon_count": 29,
        "role": (
            "Treacle ribosome biogenesis factor 1 -- nucleolar "
            "phosphoprotein essential for rRNA transcription and "
            "ribosome biogenesis.  Treacle recruits upstream binding "
            "factor (UBF) and RNA polymerase I to rDNA repeats and "
            "participates in rRNA methylation via NOP56/fibrillarin.  "
            "Haploinsufficiency reduces ribosome production in cranial "
            "neural crest cells (cNCCs), triggering p53-mediated "
            "apoptosis and causing Treacher Collins syndrome (TCS) -- "
            "mandibulofacial dysostosis with downslanting palpebral "
            "fissures, malar/mandibular hypoplasia, microtia, and "
            "conductive hearing loss.  Severity is extremely variable "
            "even within families."
        ),
        "disease": "Treacher Collins syndrome (TCS1)",
        "omim_disease": 154500,
        "omim_gene": 606847,
        "inheritance": "AD (haploinsufficiency; ~60% de novo)",
        "key_variants": [
            {
                "name": "c.4369_4373delAAGAA (p.Lys1457fs)",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Exon 24 deletion hotspot; one of the most recurrent "
                    "TCOF1 mutations."
                ),
            },
            {
                "name": "Over 200 unique mutations (mostly truncating)",
                "consequence": "various",
                "notes": (
                    "Most TCOF1 mutations are private (family-specific).  "
                    "~60% frameshift, ~16% nonsense, ~11% splice-site.  "
                    "All consistent with haploinsufficiency."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 200 pathogenic variants, almost all truncating.  No clear "
            "genotype-phenotype correlation.  Incidence ~1:50,000 births.  "
            "~60% de novo.  Intrafamilial variability suggests strong "
            "modifier effects and stochastic developmental noise."
        ),
        "strategy": (
            "TCS is a developmental condition -- craniofacial structures "
            "form during embryogenesis.  Gene therapy would need to be "
            "prenatal or very early postnatal: "
            "(1) In-utero gene supplementation to cranial neural crest "
            "cells (theoretical; extremely challenging); "
            "(2) p53 inhibition in cNCCs (pifithrin-alpha rescues "
            "craniofacial defects in Tcof1+/- mice) -- potential "
            "pharmacological approach during critical developmental window; "
            "(3) CRISPRa to upregulate the wild-type TCOF1 allele in "
            "neural crest progenitors; "
            "(4) Current management is reconstructive surgery, hearing "
            "aids, and MDT support.  "
            "Gene therapy for TCS is largely theoretical as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Standard care is multidisciplinary: "
            "craniofacial reconstructive surgery (mandibular distraction "
            "osteogenesis, malar augmentation, ear reconstruction), "
            "hearing aids / bone-anchored hearing devices, speech therapy.  "
            "p53 modulation as a potential prenatal intervention is at "
            "proof-of-concept stage in animal models."
        ),
        "conditions": ["treacher_collins_syndrome", "mandibulofacial_dysostosis",
                        "craniofacial", "ribosomopathy", "neural_crest"],
    },

    # -------------------------------------------------------------------
    # 7b.  FGFR2 -- Fibroblast Growth Factor Receptor 2
    # -------------------------------------------------------------------
    "FGFR2": {
        "gene_id": 2263,
        "chrom": "chr10",
        "start": 121_478_330,
        "end": 121_598_458,
        "strand": "-",
        "refseq": "NC_000010.11",
        "cytoband": "10q26.13",
        "exon_count": 26,
        "role": (
            "Fibroblast growth factor receptor 2 -- receptor tyrosine "
            "kinase.  Two major isoforms: FGFR2-IIIb (epithelial, "
            "keratinocyte growth factor receptor) and FGFR2-IIIc "
            "(mesenchymal).  Critical for craniofacial suture biology, "
            "limb development, and organogenesis.  Gain-of-function "
            "mutations cause craniosynostosis syndromes: Crouzon syndrome "
            "(premature coronal suture fusion, midfacial hypoplasia), "
            "Apert syndrome (craniosynostosis + syndactyly), Pfeiffer "
            "syndrome, Beare-Stevenson, Jackson-Weiss.  Mutations cluster "
            "in the IgII-IgIII linker and IgIII domain."
        ),
        "disease": "Crouzon / Apert / Pfeiffer / craniosynostosis syndromes",
        "omim_gene": 176943,
        "inheritance": "AD (gain-of-function; many de novo, especially Apert)",
        "key_variants": [
            {
                "name": "p.Ser252Trp (S252W)",
                "hgvs_coding": "NM_000141.5:c.755C>G",
                "rsid": "rs121909218",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Apert syndrome -- accounts for ~65% of cases.  "
                    "IgII-IgIII linker; increases FGF binding affinity "
                    "and broadens ligand specificity.  Almost exclusively "
                    "de novo (paternal age effect)."
                ),
            },
            {
                "name": "p.Pro253Arg (P253R)",
                "hgvs_coding": "NM_000141.5:c.758C>G",
                "rsid": "rs121909219",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Apert syndrome -- accounts for ~33% of cases.  "
                    "Adjacent to S252W; same domain and mechanism."
                ),
            },
            {
                "name": "p.Cys342Tyr/Arg/Ser/Phe/Trp",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "IgIII domain cysteine -- multiple substitutions at "
                    "this residue cause Crouzon or Pfeiffer syndrome.  "
                    "Loss of disulfide bond alters receptor conformation."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 100 pathogenic variants.  Two mutations (S252W and P253R) "
            "account for ~98% of Apert syndrome.  Crouzon mutations more "
            "scattered across IgIII domain.  All are gain-of-function "
            "(constitutive activation or enhanced ligand binding)."
        ),
        "strategy": (
            "Craniosynostosis syndromes are developmental -- suture fusion "
            "occurs in utero and early infancy.  Approaches: "
            "(1) Allele-specific silencing of gain-of-function FGFR2 allele "
            "using CRISPRi or ASO -- would need to be delivered to "
            "suture mesenchyme before or during fusion; "
            "(2) FGFR2 kinase inhibitors (infigratinib, pemigatinib) -- "
            "potentially repurposed from oncology; FGFR inhibition in "
            "Crouzon mouse models prevents suture fusion; "
            "(3) Local FGFR inhibitor delivery to sutures (preclinical); "
            "(4) For Apert S252W: base editing (C>G, would need CBE on "
            "antisense strand or prime editing) -- prenatal delivery "
            "challenge remains.  "
            "Current management is cranial vault surgery."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Infigratinib (Lilly) being "
            "investigated for achondroplasia (FGFR3) -- potential "
            "extrapolation to FGFR2.  Standard care is cranial vault "
            "remodelling surgery (strip craniectomy, fronto-orbital "
            "advancement), midface advancement (Le Fort III), "
            "syndactyly release for Apert."
        ),
        "conditions": ["crouzon_syndrome", "apert_syndrome",
                        "pfeiffer_syndrome", "craniosynostosis",
                        "craniofacial"],
    },

    # -------------------------------------------------------------------
    # 7c.  TWIST1 -- Twist Family bHLH Transcription Factor 1
    # -------------------------------------------------------------------
    "TWIST1": {
        "gene_id": 7291,
        "chrom": "chr7",
        "start": 19_113_047,
        "end": 19_117_636,
        "strand": "-",
        "refseq": "NC_000007.14",
        "cytoband": "7p21.1",
        "exon_count": 4,
        "role": (
            "TWIST1 -- basic helix-loop-helix (bHLH) transcription factor.  "
            "Master regulator of mesoderm development, neural crest "
            "migration, and cranial suture patency.  TWIST1 inhibits "
            "osteogenic differentiation in suture mesenchyme by "
            "antagonising RUNX2.  Haploinsufficiency causes Saethre-Chotzen "
            "syndrome (coronal craniosynostosis, facial asymmetry, ptosis, "
            "brachydactyly, cutaneous syndactyly of digits 2-3).  "
            "Complete loss of TWIST1 causes more severe bilateral coronal "
            "synostosis."
        ),
        "disease": "Saethre-Chotzen syndrome (SCS)",
        "omim_disease": 101400,
        "omim_gene": 601622,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "7p21 deletion encompassing TWIST1",
                "consequence": "deletion",
                "clinical_significance": "pathogenic",
                "frequency": "~15% of SCS cases",
                "notes": (
                    "Large deletions may include neighbouring genes; "
                    "associated with more severe phenotype and learning "
                    "difficulties."
                ),
            },
            {
                "name": "p.Arg118Cys (R118C)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "bHLH domain; impairs DNA binding and dimerisation.",
            },
            {
                "name": "Multiple truncating mutations throughout gene",
                "consequence": "nonsense/frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Consistent with haploinsufficiency.  No genotype-"
                    "phenotype correlation."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 120 pathogenic variants.  ~35% missense (mostly bHLH "
            "domain), ~35% truncating, ~15% large deletions, ~15% "
            "splice-site.  Second most common syndromic craniosynostosis "
            "after Crouzon/Apert."
        ),
        "strategy": (
            "Like other craniosynostosis genes, the developmental "
            "window limits gene therapy approaches: "
            "(1) CRISPRa to upregulate the wild-type TWIST1 allele in "
            "suture mesenchyme -- most conceptually straightforward for "
            "haploinsufficiency; "
            "(2) Local gene supplementation to patent sutures "
            "(AAV-TWIST1, cDNA ~0.6 kb, very small) to prevent "
            "premature fusion; "
            "(3) Pharmacological approaches to inhibit osteogenesis at "
            "sutures -- BMP antagonists, RUNX2 inhibitors; "
            "(4) Standard care remains surgery.  "
            "All gene therapy approaches are theoretical as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Standard care is cranial vault "
            "surgery.  Some patients require orbital surgery for ptosis "
            "and strabismus."
        ),
        "conditions": ["saethre_chotzen_syndrome", "craniosynostosis",
                        "craniofacial"],
    },

    # -------------------------------------------------------------------
    # 7d.  IRF6 -- Interferon Regulatory Factor 6
    # -------------------------------------------------------------------
    "IRF6": {
        "gene_id": 3664,
        "chrom": "chr1",
        "start": 209_785_617,
        "end": 209_806_142,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1q32.2",
        "exon_count": 9,
        "role": (
            "Interferon regulatory factor 6 -- transcription factor "
            "essential for keratinocyte differentiation, periderm "
            "formation, and palatal shelf fusion.  IRF6 regulates the "
            "switch between proliferation and differentiation in oral "
            "epithelium.  Loss-of-function (haploinsufficiency) causes "
            "Van der Woude syndrome (VWS) -- the most common syndromic "
            "form of cleft lip/palate, characterised by lower lip pits "
            "and cleft lip with or without cleft palate.  More severe "
            "loss causes popliteal pterygium syndrome (PPS) -- VWS "
            "features plus popliteal webbing, genital anomalies, and "
            "syndactyly."
        ),
        "disease": "Van der Woude syndrome (VWS) / Popliteal pterygium syndrome (PPS)",
        "omim_disease_VWS": 119300,
        "omim_disease_PPS": 263650,
        "omim_gene": 607199,
        "inheritance": "AD (variable expressivity, ~80% penetrance)",
        "key_variants": [
            {
                "name": "p.Arg84Cys (R84C)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "DNA-binding domain; most common VWS mutation.  "
                    "Impairs IRF6 transcriptional activity."
                ),
            },
            {
                "name": "p.Arg84His (R84H)",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Same codon as R84C; VWS phenotype.",
            },
            {
                "name": "p.Arg250Pro (R250P) and other protein-interaction domain variants",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Protein-interaction domain mutations tend to cause "
                    "the more severe PPS phenotype."
                ),
            },
        ],
        "mutational_landscape": (
            "Over 200 pathogenic variants.  ~70% cause VWS, ~20% cause PPS.  "
            "VWS mutations predominantly affect DNA-binding domain; PPS "
            "mutations cluster in protein-interaction domain.  VWS is the "
            "most common syndromic clefting disorder (~1:35,000 births; "
            "~2% of all cleft lip/palate)."
        ),
        "strategy": (
            "Cleft lip/palate and lip pits form during embryonic weeks "
            "6-12; gene therapy would need to be prenatal: "
            "(1) Standard care is surgical repair of cleft lip (3-6 months), "
            "cleft palate (9-14 months), and lip pits.  Results are "
            "excellent with modern techniques; "
            "(2) For PPS: additional surgery for popliteal webbing and "
            "genital anomalies; "
            "(3) Research interest in IRF6 modulation to improve wound "
            "healing and reduce scarring post-surgery; "
            "(4) CRISPRa to upregulate IRF6 could theoretically be used "
            "to promote oral epithelial differentiation in tissue "
            "engineering for palatal reconstruction; "
            "(5) Common variant rs642961 in IRF6 enhancer is a risk "
            "factor for non-syndromic cleft lip in the general population.  "
            "Gene therapy for VWS/PPS is theoretical as of 2026."
        ),
        "clinical_programs": (
            "No gene therapy trials.  Standard care is surgical.  "
            "IRF6 common variant rs642961 used in genetic counselling "
            "for non-syndromic CL/P risk assessment."
        ),
        "conditions": ["van_der_woude_syndrome", "popliteal_pterygium_syndrome",
                        "cleft_lip_palate", "craniofacial",
                        "orofacial_clefting"],
    },
}


# ============================================================================
# 8. UNIFIED DATABASE -- ALL EXPANDED TARGETS
# ============================================================================

ALL_EXPANDED_TARGETS = {}
ALL_EXPANDED_TARGETS.update(MONOGENIC_DIABETES_TARGETS)
ALL_EXPANDED_TARGETS.update(OBESITY_TARGETS)
ALL_EXPANDED_TARGETS.update(OPHTHALMOLOGY_TARGETS)
ALL_EXPANDED_TARGETS.update(NCL_TARGETS)
ALL_EXPANDED_TARGETS.update(BONE_MARROW_FAILURE_TARGETS)
ALL_EXPANDED_TARGETS.update(IMMUNE_CHECKPOINT_TARGETS)
ALL_EXPANDED_TARGETS.update(CRANIOFACIAL_TARGETS)


# ============================================================================
# QUICK-REFERENCE COORDINATE TABLE
# ============================================================================

EXPANDED_COORDINATE_TABLE = {
    # Monogenic Diabetes
    "GCK":    {"chrom": "chr7",  "start": 44_143_213,  "end": 44_189_439},
    "HNF1A":  {"chrom": "chr12", "start": 120_978_543, "end": 121_002_512},
    "HNF4A":  {"chrom": "chr20", "start": 44_355_699,  "end": 44_434_596},
    "HNF1B":  {"chrom": "chr17", "start": 37_686_431,  "end": 37_745_059},
    "KCNJ11": {"chrom": "chr11", "start": 17_385_248,  "end": 17_389_346},
    "ABCC8":  {"chrom": "chr11", "start": 17_392_498,  "end": 17_476_845},
    # Obesity
    "MC4R":   {"chrom": "chr18", "start": 60_371_062,  "end": 60_372_775},
    "LEP":    {"chrom": "chr7",  "start": 128_241_278, "end": 128_257_629},
    "LEPR":   {"chrom": "chr1",  "start": 65_420_652,  "end": 65_641_559},
    "POMC":   {"chrom": "chr2",  "start": 25_160_860,  "end": 25_168_580},
    "PCSK1":  {"chrom": "chr5",  "start": 96_390_333,  "end": 96_433_248},
    # Ophthalmology
    "PAX6":   {"chrom": "chr11", "start": 31_789_026,  "end": 31_817_961},
    "ABCA4":  {"chrom": "chr1",  "start": 93_992_834,  "end": 94_121_148},
    "RS1":    {"chrom": "chrX",  "start": 18_639_688,  "end": 18_672_108},
    "BEST1":  {"chrom": "chr11", "start": 61_949_821,  "end": 61_965_515},
    "CYP1B1": {"chrom": "chr2",  "start": 38_067_509,  "end": 38_076_151},
    "CNGA3":  {"chrom": "chr2",  "start": 98_346_456,  "end": 98_398_601},
    "CNGB3":  {"chrom": "chr8",  "start": 86_574_179,  "end": 86_743_634},
    # NCL / Batten Disease
    "CLN3":   {"chrom": "chr16", "start": 28_466_653,  "end": 28_492_082},
    "TPP1":   {"chrom": "chr11", "start": 6_612_768,   "end": 6_619_422},
    "PPT1":   {"chrom": "chr1",  "start": 40_071_461,  "end": 40_097_252},
    "MFSD8":  {"chrom": "chr4",  "start": 127_917_732, "end": 127_965_963},
    "CLN5":   {"chrom": "chr13", "start": 76_992_081,  "end": 77_005_117},
    "CLN6":   {"chrom": "chr15", "start": 68_206_992,  "end": 68_257_211},
    "CLN8":   {"chrom": "chr8",  "start": 1_753_059,   "end": 1_786_570},
    # Bone Marrow Failure
    "RUNX1":  {"chrom": "chr21", "start": 34_787_801,  "end": 35_049_302},
    "GATA2":  {"chrom": "chr3",  "start": 128_479_422, "end": 128_493_201},
    "MPL":    {"chrom": "chr1",  "start": 43_337_818,  "end": 43_354_466},
    "SBDS":   {"chrom": "chr7",  "start": 66_987_680,  "end": 66_995_586},
    # Immune Checkpoints
    "PDCD1":  {"chrom": "chr2",  "start": 241_849_884, "end": 241_858_894},
    "CD274":  {"chrom": "chr9",  "start": 5_450_542,   "end": 5_470_554},
    "LAG3":   {"chrom": "chr12", "start": 6_772_520,   "end": 6_778_455},
    "HAVCR2": {"chrom": "chr5",  "start": 157_085_832, "end": 157_109_044},
    "TIGIT":  {"chrom": "chr3",  "start": 114_294_028, "end": 114_310_288},
    # Craniofacial
    "TCOF1":  {"chrom": "chr5",  "start": 150_357_697, "end": 150_400_293},
    "FGFR2":  {"chrom": "chr10", "start": 121_478_330, "end": 121_598_458},
    "TWIST1": {"chrom": "chr7",  "start": 19_113_047,  "end": 19_117_636},
    "IRF6":   {"chrom": "chr1",  "start": 209_785_617, "end": 209_806_142},
}


# ============================================================================
# DISEASE-TO-GENE MAPPING
# ============================================================================

EXPANDED_DISEASE_MAP = {
    # Monogenic Diabetes
    "MODY1":       {"primary_gene": "HNF4A", "category": "monogenic_diabetes"},
    "MODY2":       {"primary_gene": "GCK",   "category": "monogenic_diabetes"},
    "MODY3":       {"primary_gene": "HNF1A", "category": "monogenic_diabetes"},
    "MODY5":       {"primary_gene": "HNF1B", "category": "monogenic_diabetes"},
    "neonatal_diabetes_KCNJ11": {
        "primary_gene": "KCNJ11", "category": "monogenic_diabetes",
    },
    "neonatal_diabetes_ABCC8": {
        "primary_gene": "ABCC8", "category": "monogenic_diabetes",
    },
    # Monogenic Obesity
    "MC4R_obesity":     {"primary_gene": "MC4R",  "category": "monogenic_obesity"},
    "leptin_deficiency": {"primary_gene": "LEP",   "category": "monogenic_obesity"},
    "LEPR_deficiency":  {"primary_gene": "LEPR",  "category": "monogenic_obesity"},
    "POMC_deficiency":  {"primary_gene": "POMC",  "category": "monogenic_obesity"},
    "PCSK1_deficiency": {"primary_gene": "PCSK1", "category": "monogenic_obesity"},
    # Ophthalmology
    "aniridia":       {"primary_gene": "PAX6",  "category": "ophthalmology"},
    "stargardt":      {"primary_gene": "ABCA4", "category": "ophthalmology"},
    "XLRS":           {"primary_gene": "RS1",   "category": "ophthalmology"},
    "best_disease":   {"primary_gene": "BEST1", "category": "ophthalmology"},
    "PCG":            {"primary_gene": "CYP1B1", "category": "ophthalmology"},
    "achromatopsia_CNGA3": {"primary_gene": "CNGA3", "category": "ophthalmology"},
    "achromatopsia_CNGB3": {"primary_gene": "CNGB3", "category": "ophthalmology"},
    # NCL / Batten Disease
    "juvenile_NCL":   {"primary_gene": "CLN3",  "category": "NCL"},
    "CLN2_disease":   {"primary_gene": "TPP1",  "category": "NCL",
                       "approved_therapy": "cerliponase alfa (Brineura)"},
    "infantile_NCL":  {"primary_gene": "PPT1",  "category": "NCL"},
    "CLN7_disease":   {"primary_gene": "MFSD8", "category": "NCL"},
    "CLN5_disease":   {"primary_gene": "CLN5",  "category": "NCL"},
    "CLN6_disease":   {"primary_gene": "CLN6",  "category": "NCL"},
    "CLN8_disease":   {"primary_gene": "CLN8",  "category": "NCL"},
    # Bone Marrow Failure
    "FPD_AML":   {"primary_gene": "RUNX1", "category": "bone_marrow_failure"},
    "GATA2_deficiency": {"primary_gene": "GATA2", "category": "bone_marrow_failure"},
    "CAMT":      {"primary_gene": "MPL",   "category": "bone_marrow_failure"},
    "SDS":       {"primary_gene": "SBDS",  "category": "bone_marrow_failure"},
    # Immune Checkpoints
    "PD1_immunotherapy": {"primary_gene": "PDCD1", "category": "immune_checkpoint"},
    "PDL1_immunotherapy": {"primary_gene": "CD274", "category": "immune_checkpoint"},
    "LAG3_immunotherapy": {"primary_gene": "LAG3", "category": "immune_checkpoint"},
    "TIM3_immunotherapy": {"primary_gene": "HAVCR2", "category": "immune_checkpoint"},
    "TIGIT_immunotherapy": {"primary_gene": "TIGIT", "category": "immune_checkpoint"},
    # Craniofacial
    "treacher_collins": {"primary_gene": "TCOF1", "category": "craniofacial"},
    "crouzon_apert":    {"primary_gene": "FGFR2", "category": "craniofacial"},
    "saethre_chotzen":  {"primary_gene": "TWIST1", "category": "craniofacial"},
    "van_der_woude":    {"primary_gene": "IRF6",  "category": "craniofacial"},
}


# ============================================================================
# APPROVED / PIPELINE THERAPIES SUMMARY
# ============================================================================

EXPANDED_APPROVED_THERAPIES = {
    "Brineura": {
        "generic_name": "cerliponase alfa",
        "target_gene": "TPP1",
        "disease": "CLN2 disease (late infantile NCL)",
        "type": "Recombinant enzyme replacement therapy",
        "route": "intracerebroventricular (ICV)",
        "approval_year": 2017,
        "manufacturer": "BioMarin",
    },
    "IMCIVREE": {
        "generic_name": "setmelanotide",
        "target_pathway": "leptin-melanocortin (MC4R agonist)",
        "diseases": [
            "POMC deficiency obesity (2020)",
            "LEPR deficiency obesity (2020)",
            "PCSK1 deficiency obesity (2022)",
            "Bardet-Biedl syndrome (2022)",
        ],
        "type": "MC4R agonist peptide",
        "route": "subcutaneous daily injection",
        "manufacturer": "Rhythm Pharmaceuticals",
    },
    "MYALEPT": {
        "generic_name": "metreleptin",
        "target_gene": "LEP (used for lipodystrophy; off-label for congenital leptin deficiency)",
        "disease": "Lipodystrophy / congenital leptin deficiency",
        "type": "Recombinant leptin",
        "route": "subcutaneous daily injection",
        "approval_year": 2014,
        "manufacturer": "Amylin / AstraZeneca / Amryt",
    },
    "Opdualag": {
        "generic_name": "relatlimab + nivolumab",
        "target_gene": "LAG3 + PDCD1",
        "disease": "Advanced melanoma",
        "type": "Anti-LAG-3 + anti-PD-1 monoclonal antibody combination",
        "route": "intravenous",
        "approval_year": 2022,
        "manufacturer": "Bristol-Myers Squibb",
    },
}


# ============================================================================
# GENE-TO-CONDITION INDEX (for pipeline lookups)
# ============================================================================

EXPANDED_GENE_CONDITIONS = {
    # Monogenic Diabetes
    "GCK":    ["MODY2", "monogenic_diabetes", "neonatal_diabetes"],
    "HNF1A":  ["MODY3", "monogenic_diabetes"],
    "HNF4A":  ["MODY1", "monogenic_diabetes"],
    "HNF1B":  ["MODY5", "monogenic_diabetes", "RCAD", "17q12_deletion"],
    "KCNJ11": ["neonatal_diabetes", "PNDM", "TNDM", "DEND_syndrome"],
    "ABCC8":  ["neonatal_diabetes", "congenital_hyperinsulinism"],
    # Obesity
    "MC4R":   ["monogenic_obesity", "MC4R_deficiency"],
    "LEP":    ["monogenic_obesity", "leptin_deficiency"],
    "LEPR":   ["monogenic_obesity", "leptin_receptor_deficiency"],
    "POMC":   ["monogenic_obesity", "POMC_deficiency", "adrenal_insufficiency"],
    "PCSK1":  ["monogenic_obesity", "PCSK1_deficiency"],
    # Ophthalmology
    "PAX6":   ["aniridia", "ophthalmology", "corneal_disease"],
    "ABCA4":  ["stargardt_disease", "cone_rod_dystrophy", "ophthalmology"],
    "RS1":    ["X_linked_retinoschisis", "ophthalmology"],
    "BEST1":  ["best_disease", "vitelliform_macular_dystrophy", "ophthalmology"],
    "CYP1B1": ["primary_congenital_glaucoma", "ophthalmology"],
    "CNGA3":  ["achromatopsia", "ophthalmology"],
    "CNGB3":  ["achromatopsia", "ophthalmology"],
    # NCL
    "CLN3":   ["juvenile_NCL", "batten_disease", "NCL"],
    "TPP1":   ["late_infantile_NCL", "CLN2_disease", "batten_disease", "NCL"],
    "PPT1":   ["infantile_NCL", "CLN1_disease", "batten_disease", "NCL"],
    "MFSD8":  ["variant_late_infantile_NCL", "CLN7_disease", "NCL"],
    "CLN5":   ["CLN5_disease", "finnish_variant_NCL", "NCL"],
    "CLN6":   ["CLN6_disease", "kufs_disease", "NCL"],
    "CLN8":   ["CLN8_disease", "northern_epilepsy", "NCL"],
    # Bone Marrow Failure
    "RUNX1":  ["FPD_AML", "bone_marrow_failure", "thrombocytopenia"],
    "GATA2":  ["GATA2_deficiency", "MonoMAC", "bone_marrow_failure"],
    "MPL":    ["CAMT", "bone_marrow_failure", "thrombocytopenia"],
    "SBDS":   ["shwachman_diamond_syndrome", "bone_marrow_failure"],
    # Immune Checkpoints
    "PDCD1":  ["cancer_immunotherapy", "PD1", "immune_checkpoint"],
    "CD274":  ["cancer_immunotherapy", "PDL1", "immune_checkpoint"],
    "LAG3":   ["cancer_immunotherapy", "LAG3", "immune_checkpoint"],
    "HAVCR2": ["cancer_immunotherapy", "TIM3", "immune_checkpoint"],
    "TIGIT":  ["cancer_immunotherapy", "TIGIT", "immune_checkpoint"],
    # Craniofacial
    "TCOF1":  ["treacher_collins_syndrome", "craniofacial"],
    "FGFR2":  ["crouzon_syndrome", "apert_syndrome", "craniosynostosis"],
    "TWIST1": ["saethre_chotzen_syndrome", "craniosynostosis"],
    "IRF6":   ["van_der_woude_syndrome", "cleft_lip_palate", "craniofacial"],
}
