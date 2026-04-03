"""
Metabolic and liver disease gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers urea cycle defects, tyrosinemia, galactosemia, glycogen storage
diseases, maple syrup urine disease, and other inborn errors of hepatic
metabolism amenable to gene replacement or CRISPR correction.

Categories:
  1. Urea Cycle Defects (OTC, ASS1, ASL, CPS1, ARG1)
  2. Tyrosinemia (FAH)
  3. Galactosemia (GALT)
  4. Glycogen Storage Diseases (G6PC, GAA, GBE1)
  5. Maple Syrup Urine Disease (BCKDHA, BCKDHB, DBT)
  6. Wilson Disease (ATP7B)
  7. Hereditary Hemochromatosis (HFE)
  8. Alpha-1 Antitrypsin Deficiency (SERPINA1)

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - ClinicalTrials.gov, PubMed
  - Ultragenyx, Arctus, Intellia, Regeneron gene therapy pipelines
"""


# ============================================================================
# 1. UREA CYCLE DEFECTS
# ============================================================================

UREA_CYCLE_TARGETS = {

    "OTC": {
        "gene_id": 5009,
        "chrom": "chrX",
        "start": 38_211_826,
        "end": 38_280_703,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xp11.4",
        "exon_count": 10,
        "role": (
            "Ornithine transcarbamylase -- mitochondrial enzyme catalysing "
            "the second step of the urea cycle (ornithine + carbamoyl "
            "phosphate -> citrulline).  Most common urea cycle defect.  "
            "Hemizygous males present with neonatal hyperammonemic crisis; "
            "heterozygous females have variable expressivity."
        ),
        "disease": "OTC deficiency (X-linked urea cycle defect)",
        "omim_disease": 311250,
        "omim_gene": 300461,
        "inheritance": "XL",
        "key_variants": [
            {
                "name": "p.Arg109Ter (R109X)",
                "rsid": "rs72554345",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Severe neonatal-onset; no residual enzyme activity.",
            },
            {
                "name": "p.Arg277Trp (R277W)",
                "rsid": "rs72554365",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Partial deficiency; late-onset phenotype.",
            },
        ],
        "strategy": (
            "(1) AAV8 liver-directed gene replacement (cDNA ~1.1 kb, ideal "
            "for AAV) -- Ultragenyx DTX401 approach; "
            "(2) LNP-delivered mRNA (Arctus ARC-AAT approach); "
            "(3) CRISPR insertion of OTC cDNA at albumin safe harbor locus "
            "for durable hepatic expression (Intellia approach); "
            "(4) Base editing for specific missense mutations.  "
            "Delivery: AAV8 IV (liver tropism) or LNP IV."
        ),
        "clinical_programs": (
            "Ultragenyx DTX401 (AAV8-OTC) Phase 1/2 (NCT02991144) -- dose-"
            "dependent ammonia reduction.  "
            "Ultragenyx DTX401 Phase 3 (HARMONY, NCT05345171).  "
            "Arctus ARC-OTC mRNA Phase 1.  "
            "No CRISPR trials for OTC as of 2026."
        ),
        "conditions": ["OTC_deficiency", "urea_cycle_defect", "metabolic_liver",
                        "hyperammonemia", "neonatal_metabolic_crisis"],
    },

    "ASS1": {
        "gene_id": 445,
        "chrom": "chr9",
        "start": 130_444_867,
        "end": 130_501_274,
        "strand": "+",
        "refseq": "NC_000009.12",
        "cytoband": "9q34.11",
        "exon_count": 14,
        "role": (
            "Argininosuccinate synthase 1 -- cytosolic enzyme in urea cycle, "
            "third step (citrulline + aspartate -> argininosuccinate).  "
            "Deficiency causes citrullinemia type I."
        ),
        "disease": "Citrullinemia type I",
        "omim_disease": 215700,
        "omim_gene": 603470,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gly390Arg (G390R)",
                "rsid": "rs121908641",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common worldwide; reduces enzyme activity >90%.",
            },
        ],
        "strategy": (
            "(1) AAV8 liver-directed gene replacement (cDNA ~1.2 kb); "
            "(2) CRISPR-mediated insertion at albumin locus; "
            "(3) Base editing for G390R."
        ),
        "clinical_programs": "Preclinical AAV programs. No CRISPR trials.",
        "conditions": ["citrullinemia", "urea_cycle_defect", "metabolic_liver",
                        "hyperammonemia"],
    },

    "ASL": {
        "gene_id": 435,
        "chrom": "chr7",
        "start": 66_072_631,
        "end": 66_090_360,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q11.21",
        "exon_count": 16,
        "role": (
            "Argininosuccinate lyase -- fourth step of urea cycle "
            "(argininosuccinate -> arginine + fumarate).  Deficiency causes "
            "argininosuccinic aciduria (ASA), the second most common UCD."
        ),
        "disease": "Argininosuccinic aciduria",
        "omim_disease": 207900,
        "omim_gene": 608310,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg379Cys (R379C)",
                "rsid": "rs121964990",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Severe neonatal form; disrupts tetramer interface.",
            },
        ],
        "strategy": (
            "(1) AAV liver gene replacement (cDNA ~1.4 kb, fits AAV); "
            "(2) CRISPR insertion at albumin safe harbor.  "
            "Delivery: AAV8 IV."
        ),
        "clinical_programs": "Preclinical. No gene therapy or CRISPR trials.",
        "conditions": ["argininosuccinic_aciduria", "urea_cycle_defect",
                        "metabolic_liver", "hyperammonemia"],
    },

    "CPS1": {
        "gene_id": 1373,
        "chrom": "chr2",
        "start": 210_421_479,
        "end": 210_543_314,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q34",
        "exon_count": 38,
        "role": (
            "Carbamoyl phosphate synthetase I -- first and rate-limiting step "
            "of the urea cycle (ammonia + CO2 + ATP -> carbamoyl phosphate).  "
            "Deficiency is the rarest but most severe UCD."
        ),
        "disease": "CPS1 deficiency",
        "omim_disease": 237300,
        "omim_gene": 608307,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Thr544Met (T544M)",
                "rsid": "rs121964985",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Catalytic domain; severe neonatal onset.",
            },
        ],
        "strategy": (
            "(1) AAV liver gene replacement -- cDNA ~4.5 kb, tight AAV fit; "
            "(2) Dual-AAV split approach for larger constructs.  "
            "Delivery: AAV8 IV."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["CPS1_deficiency", "urea_cycle_defect", "metabolic_liver",
                        "hyperammonemia"],
    },

    "ARG1": {
        "gene_id": 383,
        "chrom": "chr6",
        "start": 131_549_940,
        "end": 131_561_283,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6q23.2",
        "exon_count": 8,
        "role": (
            "Arginase 1 -- final step of urea cycle (arginine -> urea + "
            "ornithine).  Deficiency causes hyperargininemia with progressive "
            "spastic diplegia rather than acute hyperammonemia."
        ),
        "disease": "Arginase deficiency (hyperargininemia)",
        "omim_disease": 207800,
        "omim_gene": 608313,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gly235Arg (G235R)",
                "rsid": "rs121964991",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Active site adjacent; abolishes arginase activity.",
            },
        ],
        "strategy": (
            "(1) AAV liver gene replacement (cDNA ~1.0 kb, ideal for AAV); "
            "(2) Ex vivo CRISPR correction in hepatocytes."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["arginase_deficiency", "hyperargininemia",
                        "urea_cycle_defect", "metabolic_liver"],
    },
}


# ============================================================================
# 2. TYROSINEMIA TYPE I
# ============================================================================

TYROSINEMIA_TARGETS = {

    "FAH": {
        "gene_id": 2184,
        "chrom": "chr15",
        "start": 80_152_893,
        "end": 80_187_253,
        "strand": "+",
        "refseq": "NC_000015.10",
        "cytoband": "15q25.1",
        "exon_count": 14,
        "role": (
            "Fumarylacetoacetate hydrolase -- terminal enzyme of tyrosine "
            "catabolism.  Deficiency causes hepatorenal tyrosinemia type I "
            "(HT1): accumulation of fumarylacetoacetate/succinylacetone -> "
            "liver failure, hepatocellular carcinoma, renal Fanconi syndrome, "
            "and porphyria-like crises."
        ),
        "disease": "Hereditary tyrosinemia type I (HT1)",
        "omim_disease": 276700,
        "omim_gene": 613871,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "IVS12+5G>A (Quebec founder)",
                "rsid": "rs121964988",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Founder mutation in Quebec/Saguenay-Lac-St-Jean "
                    "(carrier rate 1/20).  Causes exon 12 skipping."
                ),
            },
            {
                "name": "p.Pro261Leu (P261L)",
                "rsid": "rs121964986",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Northern European; associated with pseudodeficiency.",
            },
        ],
        "strategy": (
            "(1) CRISPR-mediated metabolic pathway rerouting: disruption of "
            "HPD (4-hydroxyphenylpyruvate dioxygenase) converts HT1 to benign "
            "tyrosinemia type III -- proven in mouse models (Grompe lab); "
            "(2) AAV liver gene replacement (cDNA ~1.3 kb, ideal for AAV); "
            "(3) Base editing for specific mutations (IVS12+5G>A); "
            "(4) In vivo hepatocyte CRISPR correction provides growth advantage "
            "to corrected cells (self-selecting gene therapy).  "
            "Delivery: AAV8 IV to liver or LNP."
        ),
        "clinical_programs": (
            "Nitisinone (NTBC) -- approved drug blocking HPD upstream, "
            "standard of care.  "
            "Genethon liver-directed AAV-FAH preclinical.  "
            "No CRISPR trials for FAH as of 2026, but CRISPR HPD-disruption "
            "validated in Fah-null mice (Yin et al., Nature Biotech 2016)."
        ),
        "conditions": ["tyrosinemia_type_I", "HT1", "metabolic_liver",
                        "hepatorenal", "liver_failure"],
    },
}


# ============================================================================
# 3. GALACTOSEMIA
# ============================================================================

GALACTOSEMIA_TARGETS = {

    "GALT": {
        "gene_id": 2592,
        "chrom": "chr9",
        "start": 34_636_561,
        "end": 34_640_812,
        "strand": "-",
        "refseq": "NC_000009.12",
        "cytoband": "9p13.3",
        "exon_count": 11,
        "role": (
            "Galactose-1-phosphate uridylyltransferase -- second enzyme in "
            "the Leloir pathway of galactose metabolism.  Deficiency causes "
            "classic galactosemia: neonatal liver failure, E. coli sepsis, "
            "cataracts, and long-term cognitive/ovarian impairment despite "
            "dietary galactose restriction."
        ),
        "disease": "Classic galactosemia",
        "omim_disease": 230400,
        "omim_gene": 606999,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gln188Arg (Q188R)",
                "rsid": "rs75391579",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mutation in Europeans (~70% of alleles).  "
                    "No residual enzyme activity.  CpG hotspot."
                ),
            },
            {
                "name": "p.Lys285Asn (K285N, Duarte-2)",
                "rsid": "rs2070074",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Second most common; partial activity.",
            },
        ],
        "strategy": (
            "(1) AAV liver gene replacement (cDNA ~1.1 kb, ideal for AAV); "
            "(2) Base editing to correct Q188R (single nucleotide change); "
            "(3) mRNA/LNP hepatic delivery for enzyme replacement.  "
            "Delivery: AAV8 IV liver-directed."
        ),
        "clinical_programs": (
            "No gene therapy trials as of 2026.  "
            "Galactose-free diet is sole treatment but does not prevent "
            "long-term complications (speech apraxia, POI, tremor)."
        ),
        "conditions": ["galactosemia", "metabolic_liver", "neonatal_metabolic",
                        "galactose_metabolism"],
    },
}


# ============================================================================
# 4. GLYCOGEN STORAGE DISEASES
# ============================================================================

GSD_TARGETS = {

    "G6PC": {
        "gene_id": 2538,
        "chrom": "chr17",
        "start": 42_900_793,
        "end": 42_914_438,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.31",
        "exon_count": 5,
        "role": (
            "Glucose-6-phosphatase catalytic subunit -- ER membrane enzyme, "
            "final step of gluconeogenesis and glycogenolysis.  Deficiency "
            "causes GSD type Ia (von Gierke disease): severe fasting "
            "hypoglycemia, hepatomegaly, lactic acidosis, hyperuricemia, "
            "and hepatic adenomas."
        ),
        "disease": "Glycogen storage disease type Ia (von Gierke)",
        "omim_disease": 232200,
        "omim_gene": 613742,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg83Cys (R83C)",
                "rsid": "rs80356484",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common in European populations.",
            },
            {
                "name": "c.648G>T (p.Leu216=, splice mutation)",
                "rsid": "rs80356482",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": "Common in Chinese/Asian populations (~55% of alleles).",
            },
        ],
        "strategy": (
            "(1) AAV8 liver gene replacement (cDNA ~1.1 kb, ideal for AAV) -- "
            "Ultragenyx DTX401 and Genethon GSD-Ia programs; "
            "(2) CRISPR insertion at albumin locus; "
            "(3) mRNA/LNP for transient expression.  "
            "Delivery: AAV8 IV to liver."
        ),
        "clinical_programs": (
            "Ultragenyx DTX401 (AAV8-G6PC) Phase 1/2 (NCT03517085) -- "
            "showed improved fasting tolerance.  Phase 3 planned.  "
            "No CRISPR trials as of 2026."
        ),
        "conditions": ["GSD_Ia", "von_gierke_disease", "glycogen_storage_disease",
                        "metabolic_liver", "hypoglycemia"],
    },

    "GAA": {
        "gene_id": 2548,
        "chrom": "chr17",
        "start": 80_101_526,
        "end": 80_119_881,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q25.3",
        "exon_count": 20,
        "role": (
            "Acid alpha-glucosidase (acid maltase) -- lysosomal enzyme "
            "degrading glycogen.  Deficiency causes Pompe disease (GSD II): "
            "infantile-onset with hypertrophic cardiomyopathy and skeletal "
            "myopathy, or late-onset with progressive limb-girdle/respiratory "
            "muscle weakness."
        ),
        "disease": "Pompe disease (GSD II)",
        "omim_disease": 232300,
        "omim_gene": 606800,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.-32-13T>G (IVS1)",
                "rsid": "rs386834236",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common late-onset Pompe mutation in Europeans.  "
                    "Causes leaky splicing with ~10-20% residual activity."
                ),
            },
            {
                "name": "p.Asp645Glu (D645E)",
                "rsid": "rs28940868",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Common in infantile-onset; catalytic domain.",
            },
        ],
        "strategy": (
            "(1) AAV liver-directed gene replacement for cross-correction "
            "(cDNA ~2.9 kb fits AAV) -- Spark SPK-3006 approach; "
            "(2) AAV9 muscle-directed (Sarepta/Genethon); "
            "(3) CRISPR correction of IVS1 splice variant; "
            "(4) Ex vivo HSC gene therapy.  "
            "Delivery: AAV8 IV (liver) or AAV9 IV (muscle+CNS)."
        ),
        "clinical_programs": (
            "Spark/Roche SPK-3006 (AAV-GAA liver) Phase 1/2 (NCT04093349).  "
            "Sarepta SRP-9003 (AAV-GAA muscle) Phase 1/2 (NCT04174105).  "
            "Asklepios AskBio007 Phase 1.  "
            "ERT: alglucosidase alfa (Myozyme) and avalglucosidase alfa "
            "(Nexviazyme) approved."
        ),
        "conditions": ["pompe_disease", "GSD_II", "glycogen_storage_disease",
                        "metabolic_liver", "myopathy", "lysosomal_storage"],
    },

    "GBE1": {
        "gene_id": 2632,
        "chrom": "chr3",
        "start": 81_538_925,
        "end": 81_809_158,
        "strand": "+",
        "refseq": "NC_000003.12",
        "cytoband": "3p12.2",
        "exon_count": 16,
        "role": (
            "Glycogen branching enzyme 1 -- creates alpha-1,6 branch points "
            "in glycogen.  Deficiency causes GSD IV (Andersen disease): "
            "accumulation of abnormal amylopectin-like glycogen (polyglucosan) "
            "in liver, muscle, and CNS."
        ),
        "disease": "GSD IV (Andersen disease) / APBD",
        "omim_disease": 232500,
        "omim_gene": 607839,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Tyr329Ser (Y329S)",
                "rsid": "rs121918368",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Adult polyglucosan body disease (APBD) variant.",
            },
        ],
        "strategy": (
            "(1) AAV liver gene replacement (cDNA ~2.1 kb, fits AAV); "
            "(2) AAV9 for neuromuscular forms.  "
            "Delivery: AAV8 (liver) or AAV9 (systemic)."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["GSD_IV", "andersen_disease", "glycogen_storage_disease",
                        "metabolic_liver", "APBD"],
    },
}


# ============================================================================
# 5. MAPLE SYRUP URINE DISEASE (MSUD)
# ============================================================================

MSUD_TARGETS = {

    "BCKDHA": {
        "gene_id": 593,
        "chrom": "chr19",
        "start": 41_399_798,
        "end": 41_427_142,
        "strand": "-",
        "refseq": "NC_000019.10",
        "cytoband": "19q13.2",
        "exon_count": 9,
        "role": (
            "Branched-chain keto acid dehydrogenase E1-alpha subunit -- "
            "part of the BCKD complex that catabolizes branched-chain amino "
            "acids (leucine, isoleucine, valine).  Deficiency causes MSUD "
            "type Ia."
        ),
        "disease": "Maple syrup urine disease type Ia",
        "omim_disease": 248600,
        "omim_gene": 608348,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Tyr393Asn (Y393N, Mennonite founder)",
                "rsid": "rs121964968",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Old Order Mennonite founder; carrier rate ~1/10.",
            },
        ],
        "strategy": (
            "(1) AAV8 liver gene replacement (cDNA ~1.3 kb); "
            "(2) CRISPR insertion at albumin locus; "
            "(3) mRNA/LNP hepatic delivery.  "
            "Delivery: AAV8 IV."
        ),
        "clinical_programs": "Preclinical AAV programs.",
        "conditions": ["MSUD", "maple_syrup_urine_disease", "metabolic_liver",
                        "branched_chain_amino_acid"],
    },

    "BCKDHB": {
        "gene_id": 594,
        "chrom": "chr6",
        "start": 80_839_404,
        "end": 80_903_768,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6q14.1",
        "exon_count": 10,
        "role": (
            "BCKD E1-beta subunit -- partner of E1-alpha in the BCKD complex.  "
            "Deficiency causes MSUD type Ib."
        ),
        "disease": "Maple syrup urine disease type Ib",
        "omim_disease": 248600,
        "omim_gene": 248611,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gly278Ser (G278S)",
                "rsid": "rs121964969",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Thiamine-responsive MSUD variant.",
            },
        ],
        "strategy": (
            "(1) AAV8 liver gene replacement (cDNA ~1.3 kb); "
            "(2) Thiamine trial for responsive genotypes (G278S)."
        ),
        "clinical_programs": "Preclinical.",
        "conditions": ["MSUD", "maple_syrup_urine_disease", "metabolic_liver"],
    },

    "DBT": {
        "gene_id": 1629,
        "chrom": "chr1",
        "start": 100_187_953,
        "end": 100_240_093,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p21.2",
        "exon_count": 11,
        "role": (
            "Dihydrolipoamide branched-chain transacylase (E2 subunit of "
            "BCKD complex).  Deficiency causes MSUD type II."
        ),
        "disease": "Maple syrup urine disease type II",
        "omim_disease": 248600,
        "omim_gene": 248610,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Asn222Ser",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts E2 core assembly.",
            },
        ],
        "strategy": (
            "(1) AAV8 liver gene replacement (cDNA ~1.4 kb); "
            "(2) Base editing for specific variants."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["MSUD", "maple_syrup_urine_disease", "metabolic_liver"],
    },
}


# ============================================================================
# 6. WILSON DISEASE
# ============================================================================

WILSON_TARGETS = {

    "ATP7B": {
        "gene_id": 540,
        "chrom": "chr13",
        "start": 51_932_797,
        "end": 52_012_374,
        "strand": "+",
        "refseq": "NC_000013.11",
        "cytoband": "13q14.3",
        "exon_count": 21,
        "role": (
            "Copper-transporting P-type ATPase -- hepatic copper transporter "
            "responsible for biliary copper excretion and incorporation of "
            "copper into ceruloplasmin.  Deficiency causes Wilson disease: "
            "hepatic copper overload -> liver disease (cirrhosis, fulminant "
            "failure), neuropsychiatric symptoms (dystonia, tremor), and "
            "Kayser-Fleischer corneal rings."
        ),
        "disease": "Wilson disease",
        "omim_disease": 277900,
        "omim_gene": 606882,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.His1069Gln (H1069Q)",
                "rsid": "rs76151636",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mutation in Europeans (~35-75% of alleles).  "
                    "ATP-binding domain; misfolded protein retained in ER."
                ),
            },
            {
                "name": "p.Arg778Leu (R778L)",
                "rsid": "rs28942074",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common in East Asian populations (~30-50%).",
            },
        ],
        "strategy": (
            "(1) AAV8 liver gene replacement -- cDNA ~4.4 kb is tight for "
            "AAV but feasible with compact promoter (Ultragenyx/Vivet); "
            "(2) Dual-AAV or mini-gene approach for larger regulatory elements; "
            "(3) Base editing for H1069Q (European) or R778L (Asian); "
            "(4) Truncated ATP7B mini-gene retaining copper transport function.  "
            "Delivery: AAV8 IV (liver-directed)."
        ),
        "clinical_programs": (
            "Vivet Therapeutics VTX-801 (AAV-ATP7B liver) Phase 1/2 "
            "(NCT04537377).  "
            "Ultragenyx/GeneTx preclinical.  "
            "Standard of care: chelation (D-penicillamine, trientine) and "
            "zinc acetate."
        ),
        "conditions": ["wilson_disease", "copper_metabolism", "metabolic_liver",
                        "hepatolenticular_degeneration", "liver_disease"],
    },
}


# ============================================================================
# 7. HEREDITARY HEMOCHROMATOSIS
# ============================================================================

HEMOCHROMATOSIS_TARGETS = {

    "HFE": {
        "gene_id": 3077,
        "chrom": "chr6",
        "start": 26_087_281,
        "end": 26_098_343,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6p22.2",
        "exon_count": 7,
        "role": (
            "Hereditary hemochromatosis protein -- regulates hepcidin "
            "expression, the master iron regulator.  C282Y homozygosity "
            "causes type 1 hemochromatosis: reduced hepcidin -> excessive "
            "intestinal iron absorption -> iron overload in liver, heart, "
            "pancreas, and joints."
        ),
        "disease": "Hereditary hemochromatosis type 1",
        "omim_disease": 235200,
        "omim_gene": 613609,
        "inheritance": "AR (variable penetrance)",
        "key_variants": [
            {
                "name": "p.Cys282Tyr (C282Y)",
                "rsid": "rs1800562",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Major hemochromatosis allele; ~85% of clinical cases "
                    "are C282Y homozygous.  Allele frequency ~6% in Northern "
                    "Europeans.  Penetrance ~30-50% for clinical disease."
                ),
            },
            {
                "name": "p.His63Asp (H63D)",
                "rsid": "rs1799945",
                "consequence": "missense",
                "clinical_significance": "risk_factor",
                "notes": (
                    "Mild risk when compound heterozygous with C282Y.  "
                    "H63D homozygosity rarely causes clinical disease."
                ),
            },
        ],
        "strategy": (
            "(1) Base editing to correct C282Y (single nucleotide change) -- "
            "high impact given population frequency; "
            "(2) CRISPRa to upregulate HAMP (hepcidin gene) as downstream "
            "correction; "
            "(3) CRISPR-mediated HFE restoration in hepatocytes.  "
            "Delivery: LNP to liver or AAV8 IV.  "
            "Note: phlebotomy is highly effective, so gene therapy reserved "
            "for patients with organ damage or poor compliance."
        ),
        "clinical_programs": (
            "No gene therapy or CRISPR trials.  "
            "Phlebotomy standard of care.  "
            "Rusfertide (hepcidin mimetic) Phase 2 for polycythemia vera."
        ),
        "conditions": ["hemochromatosis", "iron_overload", "metabolic_liver",
                        "liver_disease"],
    },
}


# ============================================================================
# 8. ALPHA-1 ANTITRYPSIN DEFICIENCY
# ============================================================================

A1AT_TARGETS = {

    "SERPINA1": {
        "gene_id": 5265,
        "chrom": "chr14",
        "start": 94_376_747,
        "end": 94_390_654,
        "strand": "-",
        "refseq": "NC_000014.9",
        "cytoband": "14q32.13",
        "exon_count": 5,
        "role": (
            "Alpha-1 antitrypsin (AAT) -- serine protease inhibitor, primarily "
            "made in liver, protects lungs from neutrophil elastase.  Z allele "
            "(E342K) causes AAT polymerization in hepatocyte ER -> liver disease "
            "(cirrhosis, HCC) AND lung disease (emphysema) from AAT deficiency.  "
            "1 in 2500 Europeans are ZZ homozygous."
        ),
        "disease": "Alpha-1 antitrypsin deficiency",
        "omim_disease": 613490,
        "omim_gene": 107400,
        "inheritance": "AR (codominant expression)",
        "key_variants": [
            {
                "name": "p.Glu342Lys (E342K, Pi*Z allele)",
                "rsid": "rs28929474",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Pi*Z: ~95% of severe AATD.  Causes AAT misfolding and "
                    "ER retention.  Liver disease from toxic gain-of-function "
                    "(polymer accumulation), lung disease from loss-of-function."
                ),
            },
            {
                "name": "p.Glu264Val (E264V, Pi*S allele)",
                "rsid": "rs17580",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Milder deficiency (~60% AAT levels). Common in Southern Europe.",
            },
        ],
        "strategy": (
            "(1) AAV liver gene replacement with M-allele AAT cDNA (~1.2 kb, "
            "ideal for AAV) -- Arctus and others; "
            "(2) CRISPR correction of Z mutation (E342K) in hepatocytes -- "
            "addresses both liver and lung disease; "
            "(3) Base editing to correct Z->M (single nucleotide); "
            "(4) CRISPRi to reduce Z-AAT production (reduces liver toxicity "
            "but worsens lung disease unless combined with augmentation).  "
            "Delivery: AAV8 IV to liver."
        ),
        "clinical_programs": (
            "Intellia NTLA-2002 (CRISPR liver, reduce Z-AAT + augment) "
            "preclinical with IND planned.  "
            "Arctus ARC-AAT (mRNA/LNP AAT augmentation) Phase 1.  "
            "Arrowhead ARO-AAT (siRNA to reduce Z-AAT liver toxicity) "
            "Phase 2/3 (NCT03945292).  "
            "AAT augmentation therapy (Prolastin-C, Zemaira) approved."
        ),
        "conditions": ["alpha1_antitrypsin_deficiency", "AATD", "metabolic_liver",
                        "emphysema", "liver_disease", "COPD"],
    },
}


# ============================================================================
# Combined export
# ============================================================================

ALL_METABOLIC_LIVER_TARGETS = {}
ALL_METABOLIC_LIVER_TARGETS.update(UREA_CYCLE_TARGETS)
ALL_METABOLIC_LIVER_TARGETS.update(TYROSINEMIA_TARGETS)
ALL_METABOLIC_LIVER_TARGETS.update(GALACTOSEMIA_TARGETS)
ALL_METABOLIC_LIVER_TARGETS.update(GSD_TARGETS)
ALL_METABOLIC_LIVER_TARGETS.update(MSUD_TARGETS)
ALL_METABOLIC_LIVER_TARGETS.update(WILSON_TARGETS)
ALL_METABOLIC_LIVER_TARGETS.update(HEMOCHROMATOSIS_TARGETS)
ALL_METABOLIC_LIVER_TARGETS.update(A1AT_TARGETS)
