"""
Curated gene-therapy / CRISPR targets for lysosomal storage diseases, kidney
diseases, endocrine disorders, hearing loss, and connective tissue disorders
-- GRCh38 (hg38) coordinates.

Sources:  NCBI Gene, dbSNP, ClinVar, Ensembl (GRCh38.p14), OMIM,
          ClinicalTrials.gov, PubMed, CRISPR Medicine News,
          Innovative Genomics Institute clinical-trial tracker.

Coordinates pulled from NCBI Gene (RefSeq annotation, GRCh38.p14) on 2026-04-02.
Variant positions pulled from NCBI dbSNP / ClinVar on the same date.

IMPORTANT -- verify every coordinate against current NCBI / Ensembl releases
before production use.  Numbering can shift between patch levels.

All interventions require informed patient consent and IRB / ethics approval.

Categories:
  1. Lysosomal Storage Diseases (GBA1, GLA, GAA, IDUA, IDS, NPC1, HEXA, GALC, ARSA)
  2. Kidney / Renal (PKD1, PKD2, COL4A3, COL4A4, COL4A5)
  3. Endocrine (CYP21A2, INS/HLA, GH1)
  4. Hearing Loss (GJB2, COCH, SLC26A4, OTOF)
  5. Connective Tissue (COL1A1, COL1A2)
"""


# ============================================================================
# 1. LYSOSOMAL STORAGE DISEASES
# ============================================================================

LYSOSOMAL_STORAGE_TARGETS = {

    # -------------------------------------------------------------------
    # 1a.  Gaucher Disease -- GBA1
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
            "Glucocerebrosidase (acid beta-glucosidase) -- lysosomal enzyme "
            "that cleaves the beta-glucosidic linkage of glucosylceramide "
            "(glucocerebroside).  Deficiency causes accumulation of "
            "glucosylceramide in macrophages ('Gaucher cells') in spleen, "
            "liver, bone marrow, and CNS (in neuronopathic forms)."
        ),
        "disease": "Gaucher Disease (types 1, 2, 3)",
        "omim_disease": 230800,  # type 1
        "omim_gene": 606463,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1226A>G (p.Asn409Ser / N370S -- legacy nomenclature)",
                "rsid": "rs76763715",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mutation in Ashkenazi Jewish population. "
                    "Accounts for ~70% of alleles in AJ patients, ~25% in "
                    "non-Jewish populations.  Homozygotes always have type 1 "
                    "(non-neuronopathic).  Presence of at least one N370S "
                    "allele is protective against neuronopathic forms."
                ),
            },
            {
                "name": "c.1448T>C (p.Leu483Pro / L444P -- legacy nomenclature)",
                "rsid": "rs421016",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Second most common globally; most common in non-Jewish "
                    "and Asian populations.  Homozygous L444P strongly "
                    "associated with neuronopathic disease (types 2 and 3).  "
                    "Also a major Parkinson disease risk allele."
                ),
            },
            {
                "name": "c.1504C>T (p.Arg502Cys / R463C)",
                "rsid": "rs364897",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Common in type 1 Gaucher disease.",
            },
            {
                "name": "c.84dupG (84GG insertion -- legacy)",
                "rsid": "rs387906315",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Common in Ashkenazi Jewish population (~10% of alleles). "
                    "Null allele; associated with severe disease in compound "
                    "heterozygotes."
                ),
            },
            {
                "name": "c.1297G>T (p.Val433Leu / V394L)",
                "rsid": "rs75548401",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Associated with mild type 1 phenotype.",
            },
        ],
        "mutational_landscape": (
            ">500 pathogenic variants catalogued.  The 5 variants above "
            "account for ~95% of alleles in Ashkenazi Jewish and ~75% in "
            "non-Jewish populations.  CAUTION: GBA1 has a nearby pseudogene "
            "(GBAP1) causing recombination events (RecNciI, RecTL, complex "
            "alleles) that complicate sequencing; long-read or targeted "
            "approaches recommended."
        ),
        "crispr_strategy": (
            "Ex vivo approaches: (1) CRISPR/Cas9 or base editing of "
            "autologous CD34+ HSPCs to correct point mutations (especially "
            "N370S, L444P), followed by transplant -- corrected macrophages "
            "engraft and produce functional GCase.  (2) CRISPR knock-in of "
            "GBA1 cDNA at safe-harbor locus (AAVS1) in HSPCs.  "
            "In vivo: (3) AAV-mediated liver-directed gene addition -- "
            "hepatocytes secrete GCase for cross-correction.  (4) LNP-CRISPR "
            "for liver-targeted gene insertion.  Key challenge: pseudogene "
            "GBAP1 homology requires careful guide design to avoid off-target "
            "editing at the pseudogene."
        ),
        "clinical_trials": [
            {
                "name": "PR001 / LY3884961 (PROCEED trial)",
                "sponsor": "Prevail Therapeutics / Eli Lilly",
                "nct": "NCT04411654",
                "phase": "I/II",
                "approach": (
                    "AAV9-GBA1 gene therapy, intravenous delivery.  Delivers "
                    "functional GBA1 cDNA to produce glucocerebrosidase."
                ),
                "status_2026": (
                    "Phase I/II PROCEED trial for type 1 Gaucher disease.  "
                    "AAV9 vector, IV delivery.  Dose escalation ongoing."
                ),
            },
            {
                "name": "PR001 (PROVIDE trial -- Type 2 Gaucher)",
                "sponsor": "Prevail Therapeutics / Eli Lilly",
                "nct": "NCT04411654",
                "phase": "I/II",
                "approach": (
                    "AAV9-GBA1 intracisternal injection for neuronopathic "
                    "type 2 Gaucher disease in infants.  Pre-treated with "
                    "methylprednisolone and sirolimus."
                ),
                "status_2026": (
                    "Phase I/II enrolling type 2 infants.  CNS-directed AAV9."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1b.  Fabry Disease -- GLA
    # -------------------------------------------------------------------
    "GLA": {
        "gene_id": 2717,
        "chrom": "chrX",
        "start": 101_397_803,
        "end": 101_407_925,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq22.1",
        "exon_count": 7,
        "role": (
            "Alpha-galactosidase A -- lysosomal hydrolase that cleaves "
            "terminal alpha-galactosyl moieties from glycosphingolipids "
            "(primarily globotriaosylceramide, Gb3).  Deficiency causes "
            "progressive Gb3 accumulation in vascular endothelium, kidneys, "
            "heart, and nervous system."
        ),
        "disease": "Fabry Disease",
        "omim_disease": 301500,
        "omim_gene": 300644,
        "inheritance": "XLR (but females can be symptomatic due to X-inactivation)",
        "key_variants": [
            {
                "name": "c.644A>G (p.Asn215Ser / N215S)",
                "rsid": "rs28935490",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common variant in UK Biobank screening.  "
                    "Associated with cardiac-predominant (late-onset) "
                    "phenotype.  Often identified through newborn screening."
                ),
            },
            {
                "name": "c.982G>A (p.Gly328Arg)",
                "rsid": "rs104894834",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Classic Fabry disease phenotype.",
            },
            {
                "name": "c.141G>C (p.Trp47Cys / W47C)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Classic severe phenotype; disrupts enzyme folding.",
            },
            {
                "name": "c.936+1G>T (IVS6+1G>T)",
                "rsid": None,
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "notes": "Splice-site variant causing exon 6 skipping; null allele.",
            },
            {
                "name": "c.335G>A (p.Arg112His)",
                "rsid": "rs104894835",
                "consequence": "missense",
                "clinical_significance": "disputed -- reclassified as benign/VUS by many labs",
                "notes": (
                    "Population-frequent variant; does NOT cause classic "
                    "Fabry but sometimes listed as pathogenic in ClinVar.  "
                    "Exercise caution in interpretation."
                ),
            },
        ],
        "mutational_landscape": (
            ">1000 pathogenic/likely-pathogenic variants reported.  "
            "Majority are 'private' (family-specific): 306+ point mutations, "
            "115+ small indels, 8+ gross rearrangements.  No single common "
            "founder mutation -- most families carry unique variants.  "
            "X-linked, so hemizygous males fully affected; heterozygous "
            "females variably affected depending on X-inactivation patterns."
        ),
        "crispr_strategy": (
            "(1) Ex vivo CRISPR/Cas9 correction in patient iPSCs or HSPCs -- "
            "Karl-Scholler et al. (2025) demonstrated CRISPR/Cas9 restoration "
            "of wild-type GLA in Fabry iPSCs with clearance of Gb3 deposits "
            "and preserved pluripotency.  (2) In vivo liver-directed AAV or "
            "LNP delivery of GLA cDNA for secreted enzyme cross-correction.  "
            "(3) mRNA/LNP transient expression approaches.  "
            "Key consideration: X-linked gene, so female carriers with "
            "skewed X-inactivation also benefit from gene therapy."
        ),
        "clinical_trials": [
            {
                "name": "Elfabrio (pegunigalsidase alfa-iwxj)",
                "sponsor": "Chiesi / Protalix",
                "phase": "APPROVED (FDA May 2023)",
                "approach": (
                    "PEGylated enzyme replacement therapy (not gene therapy).  "
                    "Plant-cell-derived alpha-galactosidase A with extended "
                    "half-life.  IV infusion every 2 weeks."
                ),
                "status_2026": "FDA approved for adult Fabry disease.",
            },
            {
                "name": "4D-310",
                "sponsor": "4D Molecular Therapeutics",
                "nct": "NCT04519749",
                "phase": "I/II",
                "approach": (
                    "AAV-based gene therapy (IV) using a cardiac-tropic "
                    "AAV capsid (4D-C102) delivering GLA transgene.  "
                    "Targets cardiac manifestations specifically."
                ),
                "status_2026": (
                    "Phase I/II; cardiac-targeted AAV gene therapy.  "
                    "Data readouts ongoing."
                ),
            },
            {
                "name": "ST-920 / isaralgagene civaparvovec",
                "sponsor": "Sangamo / Sanofi",
                "nct": "NCT04046224",
                "phase": "I/II",
                "approach": "AAV6-based liver-directed GLA gene therapy, IV.",
                "status_2026": (
                    "Phase I/II STAAR trial; liver-targeted AAV gene therapy "
                    "showing sustained alpha-Gal A enzyme levels in blood."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1c.  Pompe Disease -- GAA
    # -------------------------------------------------------------------
    "GAA": {
        "gene_id": 2548,
        "chrom": "chr17",
        "start": 80_101_581,
        "end": 80_119_881,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q25.3",
        "exon_count": 20,
        "role": (
            "Acid alpha-glucosidase (acid maltase) -- lysosomal enzyme that "
            "hydrolyzes glycogen to glucose.  Deficiency causes glycogen "
            "accumulation in lysosomes of cardiac and skeletal muscle, "
            "leading to progressive myopathy and cardiomyopathy."
        ),
        "disease": "Pompe Disease (Glycogen Storage Disease type II)",
        "omim_disease": 232300,
        "omim_gene": 606800,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.-32-13T>G (IVS1-13T>G)",
                "rsid": "rs386834236",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common variant in late-onset Pompe disease.  "
                    "Allows ~10-20% normal splicing, producing some residual "
                    "enzyme activity.  Present in most Caucasian late-onset "
                    "patients (>50% of alleles)."
                ),
            },
            {
                "name": "c.525delT (p.Glu176Argfs*45)",
                "rsid": "rs386834235",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Common in Dutch and Northern European populations.  "
                    "Null allele; associated with classic infantile-onset "
                    "when in trans with another null allele."
                ),
            },
            {
                "name": "c.2560C>T (p.Arg854Ter)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Null allele found across multiple populations.",
            },
            {
                "name": "c.2481+102_2646+31del (exon 18 deletion)",
                "rsid": None,
                "consequence": "gross_deletion",
                "clinical_significance": "pathogenic",
                "notes": "Common in Chinese and Taiwanese patients.",
            },
        ],
        "mutational_landscape": (
            ">600 pathogenic variants reported.  IVS1 c.-32-13T>G is the "
            "most common in late-onset; 2 null alleles in trans cause "
            "classic infantile-onset.  Genotype-phenotype correlation: "
            "residual enzyme activity <1% = infantile onset; 1-30% = "
            "late-onset with variable severity."
        ),
        "crispr_strategy": (
            "(1) Ex vivo CRISPR correction of patient HSPCs or muscle "
            "progenitors.  (2) In vivo AAV-mediated GAA gene addition -- "
            "AAV vectors with ubiquitous or muscle-specific promoters.  "
            "(3) Base editing to correct the common IVS1 splice variant "
            "(single T>G change amenable to ABE correction).  "
            "Key challenge: GAA gene is large (3.6 kb CDS) but fits in AAV."
        ),
        "clinical_trials": [
            {
                "name": "Nexviazyme (avalglucosidase alfa-ngpt)",
                "sponsor": "Sanofi Genzyme",
                "phase": "APPROVED (FDA Aug 2021)",
                "approach": (
                    "Enzyme replacement therapy (not gene therapy).  "
                    "Enhanced uptake via bis-M6P glycoengineering.  "
                    "IV infusion every 2 weeks."
                ),
                "status_2026": (
                    "FDA approved for late-onset Pompe disease (age >= 1 yr).  "
                    "COMET Phase III showed superiority over alglucosidase alfa."
                ),
            },
            {
                "name": "GC301 (AAV9 gene therapy)",
                "sponsor": "Genethon / Sarepta",
                "nct": "NCT05793307",
                "phase": "I/II",
                "approach": (
                    "AAV9 vector with ubiquitous promoter delivering GAA "
                    "cDNA.  IV administration.  Targets both cardiac and "
                    "skeletal muscle."
                ),
                "status_2026": (
                    "Phase I/II for infantile-onset Pompe disease.  "
                    "Enrolling; 41 active clinical trials for Pompe as of "
                    "February 2025."
                ),
            },
            {
                "name": "SPK-3006 (AAV liver-directed)",
                "sponsor": "Spark Therapeutics",
                "nct": "NCT04093349",
                "phase": "I/II",
                "approach": (
                    "AAV-based liver-directed gene therapy.  Liver produces "
                    "and secretes GAA for systemic cross-correction of muscle."
                ),
                "status_2026": "Phase I/II; liver-directed secreted GAA approach.",
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1d.  Mucopolysaccharidosis I (MPS I / Hurler Syndrome) -- IDUA
    # -------------------------------------------------------------------
    "IDUA": {
        "gene_id": 3425,
        "chrom": "chr4",
        "start": 986_997,
        "end": 1_004_564,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4p16.3",
        "exon_count": 14,
        "role": (
            "Alpha-L-iduronidase -- lysosomal enzyme required for "
            "degradation of heparan sulfate and dermatan sulfate "
            "glycosaminoglycans (GAGs).  Deficiency causes GAG accumulation "
            "in all tissues leading to progressive multi-organ disease."
        ),
        "disease": (
            "Mucopolysaccharidosis type I (MPS I): Hurler (severe), "
            "Hurler-Scheie (intermediate), Scheie (attenuated)"
        ),
        "omim_disease": 607014,  # Hurler
        "omim_gene": 252800,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1205G>A (p.Trp402Ter / W402X)",
                "rsid": "rs121965034",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mutation in Northern European populations "
                    "(~40% of alleles).  Null allele; homozygotes have "
                    "severe Hurler phenotype."
                ),
            },
            {
                "name": "c.208C>T (p.Gln70Ter / Q70X)",
                "rsid": "rs121965033",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Second most common in European populations (~15% of "
                    "alleles).  Null allele; severe phenotype."
                ),
            },
            {
                "name": "c.1469T>C (p.Leu490Pro)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Associated with attenuated (Scheie) phenotype.",
            },
        ],
        "mutational_landscape": (
            ">200 pathogenic variants.  W402X and Q70X together account "
            "for >50% of alleles in European patients.  Genotype-phenotype: "
            "two null alleles = severe Hurler; at least one missense allele "
            "with residual activity = attenuated phenotypes."
        ),
        "crispr_strategy": (
            "(1) Ex vivo lentiviral gene addition to HSPCs (current lead "
            "approach -- OTL-203).  (2) In vivo AAV-mediated gene addition "
            "to liver (RGX-111) or CNS.  (3) CRISPR base editing to correct "
            "W402X (C>T, amenable to ABE) or Q70X.  (4) In vivo CRISPR at "
            "AAVS1 safe-harbor for IDUA knock-in.  Key challenge: CNS "
            "penetration -- blood-brain barrier limits systemic ERT efficacy "
            "for Hurler neurological disease."
        ),
        "clinical_trials": [
            {
                "name": "OTL-203 (HURCULES trial)",
                "sponsor": "Orchard Therapeutics",
                "phase": "Registrational (Phase III equivalent)",
                "approach": (
                    "Ex vivo lentiviral gene therapy in autologous HSPCs.  "
                    "Transduced cells engraft and produce supraphysiological "
                    "IDUA levels; microglial engraftment provides CNS enzyme."
                ),
                "status_2026": (
                    "Last patient treated July 2025 in HURCULES registrational "
                    "trial comparing OTL-203 to allo-HSCT.  Updated data at "
                    "WORLDSymposium Feb 2025 showed improvements in neuro, "
                    "skeletal, and clinical outcomes.  Filing anticipated."
                ),
            },
            {
                "name": "RGX-111",
                "sponsor": "REGENXBIO / Nippon Shinyaku",
                "nct": "NCT05514249",
                "phase": "I/II",
                "approach": (
                    "AAV9-IDUA intracisternal gene therapy for CNS.  "
                    "Orphan Drug, Rare Pediatric Disease, and Fast Track "
                    "designations from FDA."
                ),
                "status_2026": (
                    "Phase I/II interim results show encouraging biological "
                    "activity and favorable safety.  REGENXBIO partnered "
                    "with Nippon Shinyaku (March 2025) for US/Asia."
                ),
            },
            {
                "name": "ISP-001",
                "sponsor": "Immusoft",
                "phase": "I",
                "approach": (
                    "Autologous B cell therapy -- patient B cells engineered "
                    "ex vivo to express and secrete IDUA, then reinfused."
                ),
                "status_2026": (
                    "Results from first human trial presented at "
                    "WORLDSymposium Jan 2025."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1e.  MPS II (Hunter Syndrome) -- IDS
    # -------------------------------------------------------------------
    "IDS": {
        "gene_id": 3423,
        "chrom": "chrX",
        "start": 149_476_988,
        "end": 149_505_306,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq28",
        "exon_count": 9,
        "role": (
            "Iduronate-2-sulfatase -- lysosomal enzyme that hydrolyzes the "
            "2-sulfate groups of dermatan sulfate and heparan sulfate GAGs.  "
            "Deficiency causes progressive multi-system GAG accumulation.  "
            "Approximately two-thirds of patients have neuronopathic form."
        ),
        "disease": "Mucopolysaccharidosis type II (Hunter Syndrome)",
        "omim_disease": 309900,
        "omim_gene": 300823,
        "inheritance": "XLR",
        "key_variants": [
            {
                "name": "IDS-IDS2 recombination / inversion",
                "rsid": None,
                "consequence": "gross_rearrangement",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Inversion involving IDS and its pseudogene IDS2, caused "
                    "by intrachromosomal recombination.  Found in ~13% of "
                    "patients.  Disrupts gene structure; null allele."
                ),
            },
            {
                "name": "Complete/partial gene deletions",
                "rsid": None,
                "consequence": "gross_deletion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~10-15% of patients; always associated with severe "
                    "neuronopathic phenotype.  Includes contiguous gene "
                    "deletion syndromes."
                ),
            },
            {
                "name": "c.1122C>T (p.Arg374Cys)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Recurrent missense in attenuated form.",
            },
        ],
        "mutational_landscape": (
            ">650 pathogenic variants.  No single common mutation -- highly "
            "allelic heterogeneous.  ~20% are gross rearrangements "
            "(deletions, IDS-IDS2 recombination).  Missense, nonsense, "
            "splice, and indel variants distributed across all 9 exons.  "
            "X-linked: males hemizygous, females rarely symptomatic."
        ),
        "crispr_strategy": (
            "(1) Ex vivo lentiviral transduction of HSPCs with IDS cDNA -- "
            "autologous HSCT approach under clinical investigation.  "
            "(2) In vivo AAV-mediated CNS-directed gene therapy (RGX-121) "
            "via intracisternal injection to address neuronopathic disease.  "
            "(3) CRISPR-mediated IDS cDNA knock-in at safe-harbor loci in "
            "HSPCs.  Key challenge: like MPS I, CNS penetration is critical "
            "for neuronopathic patients."
        ),
        "clinical_trials": [
            {
                "name": "ELAPRASE (idursulfase)",
                "sponsor": "Takeda",
                "phase": "APPROVED (FDA 2006)",
                "approach": (
                    "Enzyme replacement therapy, IV weekly.  Does NOT cross "
                    "blood-brain barrier; only addresses somatic disease."
                ),
                "status_2026": "Standard of care for somatic MPS II.",
            },
            {
                "name": "RGX-121 (CAMPSIITE trial)",
                "sponsor": "REGENXBIO",
                "nct": "NCT03566043",
                "phase": "I/II/III",
                "approach": (
                    "AAV9-IDS intracisternal gene therapy.  One-time injection "
                    "targeting CNS for neuronopathic MPS II.  Enrolled boys "
                    "4-59 months with neuronopathic form."
                ),
                "status_2026": (
                    "Phase I/II/III (CAMPSIITE).  Also Phase II/III trial "
                    "NCT07236606 for children >= 5 years.  Dose escalation "
                    "showing CNS GAG reduction."
                ),
            },
            {
                "name": "HSC Gene Therapy (Manchester/UCL)",
                "sponsor": "University of Manchester / UCL",
                "phase": "I/II",
                "approach": (
                    "Ex vivo lentiviral HSPC gene therapy.  Five children "
                    "under 1 year treated.  Autologous HSCs transduced with "
                    "IDS-expressing lentiviral vector."
                ),
                "status_2026": (
                    "Phase I/II; groundbreaking trial cleared to begin in UK "
                    "(2024-2025).  Treating pre-symptomatic infants."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1f.  Niemann-Pick Type C -- NPC1
    # -------------------------------------------------------------------
    "NPC1": {
        "gene_id": 4864,
        "chrom": "chr18",
        "start": 23_506_184,
        "end": 23_586_506,
        "strand": "-",
        "refseq": "NC_000018.10",
        "cytoband": "18q11.2",
        "exon_count": 25,
        "role": (
            "NPC intracellular cholesterol transporter 1 -- large "
            "transmembrane protein in late endosomes/lysosomes that "
            "mediates intracellular cholesterol trafficking.  Deficiency "
            "causes cholesterol and sphingolipid accumulation in lysosomes, "
            "leading to progressive neurodegeneration and visceral disease."
        ),
        "disease": "Niemann-Pick Disease Type C1",
        "omim_disease": 257220,
        "omim_gene": 607623,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.3182T>C (p.Ile1061Thr / I1061T)",
                "rsid": "rs80358259",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common NPC1 mutation worldwide, especially in "
                    "Western European and Hispanic populations (~20% of "
                    "alleles).  Protein misfolding leads to ER retention "
                    "and degradation."
                ),
            },
            {
                "name": "c.3019C>G (p.Pro1007Ala / P1007A)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Blocks cholesterol egress from lysosomes but protein "
                    "reaches lysosome; used for phenotypic classification "
                    "(variant biochemical phenotype)."
                ),
            },
        ],
        "mutational_landscape": (
            ">400 pathogenic NPC1 variants.  I1061T is the most common; "
            "most other variants are private.  95% of NPC cases caused by "
            "NPC1 mutations (vs. 5% NPC2).  NPC1 is large (1278 aa, 25 "
            "exons), complicating gene therapy approaches due to AAV "
            "packaging limits."
        ),
        "crispr_strategy": (
            "(1) AAV9-mediated gene therapy with truncated NPC1 promoter -- "
            "novel small promoter (2023) permits packaging of the large NPC1 "
            "cDNA (3.8 kb) in AAV9.  CSF delivery for CNS.  (2) Base editing "
            "to correct I1061T (T>C, amenable to CBE).  (3) ASO approaches "
            "to enhance NPC1 protein folding.  Key challenge: large gene "
            "size near AAV packaging limit (~4.7 kb); requires compact "
            "promoter/regulatory elements."
        ),
        "clinical_trials": [
            {
                "name": "BGT-NPC (Bloomsbury Genetic Therapies / UCL)",
                "sponsor": "Bloomsbury Genetic Therapies",
                "phase": "Preclinical (IND-enabling)",
                "approach": (
                    "AAV9-NPC1 with novel small endogenous NPC1 promoter, "
                    "intrathecal delivery.  Enhanced therapeutic efficacy "
                    "shown in Npc1nih and Npc1nmf164 mouse models."
                ),
                "status_2026": (
                    "Preclinical / IND-enabling studies.  Received Rare "
                    "Pediatric Disease designation.  Spinout from UCL "
                    "(October 2022).  First-in-human anticipated."
                ),
            },
            {
                "name": "Levacetylleucine",
                "sponsor": "IntraBio",
                "phase": "APPROVED (FDA)",
                "approach": (
                    "Small molecule (not gene therapy).  Modulates lipid "
                    "trafficking; standalone symptomatic treatment."
                ),
                "status_2026": "FDA approved; disease-modifying but not curative.",
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1g.  Tay-Sachs Disease -- HEXA
    # -------------------------------------------------------------------
    "HEXA": {
        "gene_id": 3073,
        "chrom": "chr15",
        "start": 72_340_924,
        "end": 72_376_014,
        "strand": "-",
        "refseq": "NC_000015.10",
        "cytoband": "15q23",
        "exon_count": 14,
        "role": (
            "Hexosaminidase subunit alpha -- alpha subunit of the lysosomal "
            "enzyme beta-hexosaminidase A (HexA = alpha+beta heterodimer).  "
            "HexA degrades GM2 ganglioside in the brain.  Deficiency causes "
            "catastrophic GM2 ganglioside accumulation in neurons."
        ),
        "disease": "Tay-Sachs Disease (GM2 gangliosidosis, variant B)",
        "omim_disease": 272800,
        "omim_gene": 606869,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1274_1277dupTATC (p.Tyr427IlefsTer5 -- '1278insTATC')",
                "rsid": "rs387906309",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common in Ashkenazi Jewish population (~80% of "
                    "AJ carrier alleles).  Exon 11 TATC insertion.  Null "
                    "allele; infantile onset."
                ),
            },
            {
                "name": "c.1421+1G>C (IVS12+1G>C)",
                "rsid": "rs387906310",
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Second most common in Ashkenazi Jewish.  Intron 12 "
                    "splice donor site.  Null allele."
                ),
            },
            {
                "name": "c.805G>A (p.Gly269Ser / G269S)",
                "rsid": "rs121907966",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Late-onset (chronic/adult) Tay-Sachs.  Retains some "
                    "residual enzyme activity (~2-4%).  Common in non-Jewish "
                    "populations."
                ),
            },
            {
                "name": "c.78G>A (p.Trp26Ter -- pseudodeficiency allele R247W/R249W)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic/pseudodeficiency",
                "notes": (
                    "Pseudodeficiency alleles show reduced enzyme activity "
                    "in vitro but are NOT pathogenic.  Important for "
                    "genetic counseling to avoid false positives."
                ),
            },
        ],
        "mutational_landscape": (
            ">180 pathogenic variants.  3 mutations account for 98% of "
            "carrier alleles in Ashkenazi Jewish: 1278insTATC, IVS12+1G>C, "
            "and G269S.  Carrier frequency ~1/30 in AJ.  Infantile onset "
            "from two null alleles; late-onset from at least one allele "
            "with residual activity."
        ),
        "crispr_strategy": (
            "(1) In vivo dual-AAV delivery of HEXA + HEXB to CNS -- "
            "intrathalamic and intrathecal injection to provide both "
            "subunits for functional HexA enzyme.  (2) Base editing "
            "to correct the 1278insTATC frameshift (complex -- may "
            "require prime editing).  (3) Ex vivo HSPC gene therapy.  "
            "Key challenge: dual-subunit enzyme (alpha from HEXA, beta "
            "from HEXB) -- gene therapy for Tay-Sachs must supply HEXA "
            "and leverage endogenous HEXB.  Sandhoff disease requires HEXB."
        ),
        "clinical_trials": [
            {
                "name": "Dual-vector rAAVrh8 gene therapy (UMass)",
                "sponsor": "UMass Medical School / Auburn University",
                "nct": "NCT04669535",
                "phase": "I/II",
                "approach": (
                    "Dual intrathalamic and intrathecal AAVrh8 vectors "
                    "delivering HEXA and HEXB cDNAs.  Single-dose, "
                    "dose-escalation study in infantile and juvenile "
                    "GM2 gangliosidosis (Tay-Sachs and Sandhoff)."
                ),
                "status_2026": (
                    "Phase I/II published in Nature Medicine (2025).  "
                    "9 participants (6 infantile, 3 juvenile) received "
                    "single dose.  All showed >2x normal HexA enzyme "
                    "production.  Infantile patients showed clinical "
                    "stabilization, prolonged oral feeding to 3-3.5 yrs.  "
                    "Juvenile dystonia worsening noted (excluded from "
                    "further enrollment).  First-in-human gene therapy "
                    "for Tay-Sachs."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1h.  Krabbe Disease (Globoid Cell Leukodystrophy) -- GALC
    # -------------------------------------------------------------------
    "GALC": {
        "gene_id": 2581,
        "chrom": "chr14",
        "start": 87_933_014,
        "end": 87_993_667,
        "strand": "-",
        "refseq": "NC_000014.9",
        "cytoband": "14q31.3",
        "exon_count": 17,
        "role": (
            "Galactosylceramidase (galactocerebroside beta-galactosidase) "
            "-- lysosomal enzyme that degrades galactosylceramide and "
            "psychosine.  Deficiency causes toxic psychosine accumulation "
            "leading to oligodendrocyte death and severe demyelination."
        ),
        "disease": "Krabbe Disease (Globoid Cell Leukodystrophy)",
        "omim_disease": 245200,
        "omim_gene": 606890,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.857G>A (p.Gly286Asp) -- 30 kb deletion",
                "rsid": "rs398123095",
                "consequence": "gross_deletion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~30 kb deletion spanning exons 11-17; most common "
                    "mutation in European populations (~45% of alleles in "
                    "Scandinavians).  Null allele; infantile onset."
                ),
            },
            {
                "name": "c.1700A>C (p.Tyr567Ser / Y567S)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Associated with late-onset Krabbe disease.",
            },
            {
                "name": "c.1901T>C (p.Leu634Ser / L634S)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Found in multiple late-onset patients.",
            },
        ],
        "mutational_landscape": (
            ">200 pathogenic variants.  The ~30 kb deletion is the "
            "predominant European mutation.  Most other variants are "
            "private.  Disease incidence ~1:100,000; severe infantile "
            "form is most common (onset 3-6 months, fatal by 2 years "
            "without treatment)."
        ),
        "crispr_strategy": (
            "(1) In vivo AAV-mediated gene therapy -- AAV9 or AAVrh10 "
            "with GALC cDNA; IV delivery for peripheral + CNS expression.  "
            "(2) Ex vivo HSPC gene therapy combined with pre-transplant "
            "conditioning.  (3) CRISPR correction not practical for ~30 kb "
            "deletion; gene addition is preferred.  Key challenge: timing "
            "-- treatment must begin before irreversible demyelination, "
            "ideally neonatal if detected by newborn screening."
        ),
        "clinical_trials": [
            {
                "name": "FBX-101 (RESKUE trial)",
                "sponsor": "Forge Biologics",
                "nct": "NCT04693598",
                "phase": "I/II",
                "approach": (
                    "AAV-GALC gene therapy, IV infusion.  Targets both "
                    "CNS and PNS.  Shown to correct neuropathy, improve "
                    "myelination, and extend survival in animal models."
                ),
                "status_2026": (
                    "Phase I/II RESKUE trial -- active, not recruiting.  "
                    "Early data promising.  FBX-101 received UK Innovation "
                    "Passport and EMA PRIME designation."
                ),
            },
            {
                "name": "PBKR03",
                "sponsor": "Passage Bio",
                "nct": "NCT04771416",
                "phase": "I/II",
                "approach": (
                    "AAV-based GALC gene therapy, intracisternal or "
                    "intravenous delivery."
                ),
                "status_2026": (
                    "SUSPENDED -- trial paused due to changes in company "
                    "strategy.  Previously active in US and Canada."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1i.  Metachromatic Leukodystrophy (MLD) -- ARSA
    # -------------------------------------------------------------------
    "ARSA": {
        "gene_id": 410,
        "chrom": "chr22",
        "start": 50_622_754,
        "end": 50_628_152,
        "strand": "-",
        "refseq": "NC_000022.11",
        "cytoband": "22q13.33",
        "exon_count": 8,
        "role": (
            "Arylsulfatase A -- lysosomal enzyme that desulfates "
            "3-O-sulfogalactosylceramide (sulfatide).  Deficiency causes "
            "sulfatide accumulation in myelin-producing oligodendrocytes "
            "and Schwann cells, leading to progressive central and "
            "peripheral demyelination."
        ),
        "disease": "Metachromatic Leukodystrophy (MLD)",
        "omim_disease": 250100,
        "omim_gene": 607574,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.465+1G>A (IVS2+1G>A)",
                "rsid": "rs121907983",
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mutation in European populations.  "
                    "Null allele (complete loss of correctly spliced mRNA).  "
                    "Homozygotes have late infantile onset."
                ),
            },
            {
                "name": "c.1283C>T (p.Pro426Leu / P426L)",
                "rsid": "rs121907985",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Second most common in Europeans.  Retains ~2-5% "
                    "residual enzyme activity.  Associated with juvenile "
                    "and adult-onset MLD."
                ),
            },
            {
                "name": "c.1210+1G>A (IVS7+1G>A)",
                "rsid": None,
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "notes": "Null allele; less common splice variant.",
            },
            {
                "name": "ARSA pseudodeficiency alleles (N350S, *96A>G)",
                "rsid": "rs2071421 / rs1065757",
                "consequence": "regulatory/missense",
                "clinical_significance": "benign (pseudodeficiency)",
                "notes": (
                    "Common in general population (~1-2% of Caucasians).  "
                    "Reduce in-vitro enzyme activity to 5-20% but do NOT "
                    "cause disease.  Critical for genetic counseling to "
                    "distinguish from true pathogenic alleles."
                ),
            },
        ],
        "mutational_landscape": (
            ">250 pathogenic variants.  IVS2+1G>A and P426L account for "
            "~50% of European alleles.  Genotype-phenotype: two null "
            "alleles = late infantile onset; one null + one residual "
            "activity allele = juvenile/adult onset."
        ),
        "crispr_strategy": (
            "(1) Ex vivo lentiviral gene addition to HSPCs -- APPROVED "
            "approach (Lenmeldy/atidarsagene autotemcel).  Corrected HSPCs "
            "differentiate into microglia that cross BBB and deliver ARSA.  "
            "(2) In vivo AAV-mediated gene therapy (intrathecal or IV).  "
            "(3) Base editing correction of IVS2+1G>A splice variant "
            "(G>A, amenable to ABE).  Gene therapy is most advanced of all "
            "LSDs for this disease."
        ),
        "clinical_trials": [
            {
                "name": "Lenmeldy (atidarsagene autotemcel / arsa-cel)",
                "sponsor": "Orchard Therapeutics",
                "phase": "APPROVED (FDA March 18, 2024)",
                "approach": (
                    "Ex vivo lentiviral ARSA gene therapy in autologous "
                    "CD34+ HSPCs.  One-time treatment.  Busulfan "
                    "myeloablative conditioning, then modified cell infusion.  "
                    "Engrafted cells produce supraphysiological ARSA levels; "
                    "microglia deliver enzyme to CNS."
                ),
                "status_2026": (
                    "FDA APPROVED March 2024 for pre-symptomatic late "
                    "infantile (PSLI), pre-symptomatic early juvenile "
                    "(PSEJ), and early symptomatic early juvenile (ESEJ) "
                    "MLD.  Cost: $4.25M (most expensive gene therapy).  "
                    "12+ years follow-up data: preserved motor and cognitive "
                    "function; survival benefit vs natural history.  "
                    "Also approved in EU as Libmeldy (Dec 2020).  "
                    "FIRST APPROVED GENE THERAPY FOR A LYSOSOMAL STORAGE "
                    "DISEASE."
                ),
            },
        ],
    },
}


# ============================================================================
# 2. KIDNEY / RENAL DISEASES
# ============================================================================

KIDNEY_RENAL_TARGETS = {

    # -------------------------------------------------------------------
    # 2a.  Autosomal Dominant Polycystic Kidney Disease -- PKD1
    # -------------------------------------------------------------------
    "PKD1": {
        "gene_id": 5310,
        "chrom": "chr16",
        "start": 2_088_708,
        "end": 2_135_898,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16p13.3",
        "exon_count": 46,
        "role": (
            "Polycystin-1 -- large transmembrane glycoprotein (~4303 aa) "
            "involved in cell-cell and cell-matrix interactions, calcium "
            "signalling, and mechanosensation via primary cilia.  Forms "
            "complex with polycystin-2 (PKD2) at the cilium.  "
            "Haploinsufficiency or loss of function disrupts tubular cell "
            "proliferation/differentiation leading to cyst formation."
        ),
        "disease": "Autosomal Dominant Polycystic Kidney Disease (ADPKD)",
        "omim_disease": 173900,
        "omim_gene": 601313,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Truncating variants (nonsense/frameshift) distributed across gene",
                "rsid": None,
                "consequence": "protein_truncating",
                "clinical_significance": "pathogenic",
                "notes": (
                    "PKD1 truncating variants (nonsense, frameshift, "
                    "canonical splice) cause more severe disease and earlier "
                    "ESRD (mean age ~55 yr) than non-truncating.  No single "
                    "common founder mutation -- most are private."
                ),
            },
            {
                "name": "Non-truncating missense variants",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "variable -- many VUS",
                "notes": (
                    "Non-truncating variants associated with later ESRD "
                    "(mean ~67 yr).  Hypomorphic alleles may require "
                    "'second hit' for cystogenesis.  Large number of VUS "
                    "complicates interpretation."
                ),
            },
        ],
        "mutational_landscape": (
            ">2500 pathogenic/likely-pathogenic variants.  PKD1 accounts "
            "for ~78% of ADPKD.  CAUTION: PKD1 exons 1-33 have 6 "
            "pseudogenes (PKD1P1-P6) with ~97% sequence identity, making "
            "sequencing extremely challenging.  Long-read sequencing or "
            "LR-PCR required for accurate genotyping.  No common founder "
            "mutation; each family typically has a unique variant."
        ),
        "crispr_strategy": (
            "(1) CRISPR base editing to correct individual missense "
            "mutations -- Mayo Clinic (Feb 2026) demonstrated single-dose "
            "CRISPR base editing in preclinical ADPKD model, correcting "
            "single-letter PKD1 mutation, slowing cyst growth and extending "
            "survival.  (2) AAV-delivered CRISPR for somatic correction of "
            "cystic epithelium.  (3) miR-17 ASO (RGLS8429) to downregulate "
            "cyst-promoting miR-17 pathway.  Key challenge: kidney delivery "
            "-- renal tubular cells are difficult to target; PKD1 mRNA is "
            "~14 kb (too large for single AAV); pseudogene homology."
        ),
        "clinical_trials": [
            {
                "name": "RGLS8429 (miR-17 ASO)",
                "sponsor": "Regulus Therapeutics",
                "phase": "I",
                "approach": (
                    "Anti-miR-17 oligonucleotide.  Targets cyst growth "
                    "pathway rather than gene correction.  Subcutaneous."
                ),
                "status_2026": (
                    "Phase I FDA-cleared; currently enrolling.  First "
                    "ASO/oligonucleotide specifically for ADPKD."
                ),
            },
            {
                "name": "Tolvaptan (Jynarque)",
                "sponsor": "Otsuka",
                "phase": "APPROVED (FDA 2018)",
                "approach": (
                    "V2 vasopressin receptor antagonist.  Slows cyst growth "
                    "by reducing cAMP.  Not gene therapy."
                ),
                "status_2026": "Standard of care for rapidly progressing ADPKD.",
            },
        ],
    },

    # -------------------------------------------------------------------
    # 2b.  ADPKD -- PKD2
    # -------------------------------------------------------------------
    "PKD2": {
        "gene_id": 5311,
        "chrom": "chr4",
        "start": 88_007_635,
        "end": 88_077_777,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4q22.1",
        "exon_count": 15,
        "role": (
            "Polycystin-2 (TRPP2) -- cation channel of the TRP family, "
            "located at primary cilia, ER, and plasma membrane.  Forms "
            "heteromeric complex with polycystin-1 for calcium signalling.  "
            "Loss of function leads to the same cystic phenotype as PKD1 "
            "mutations but with milder/later presentation."
        ),
        "disease": "Autosomal Dominant Polycystic Kidney Disease (ADPKD)",
        "omim_disease": 613095,
        "omim_gene": 173910,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Truncating and missense variants throughout gene",
                "rsid": None,
                "consequence": "protein_truncating/missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "PKD2 accounts for ~15% of ADPKD.  Milder than PKD1: "
                    "ESRD at mean age ~79 yr for PKD2 vs ~55-67 for PKD1.  "
                    "Like PKD1, variants are private to families."
                ),
            },
        ],
        "mutational_landscape": (
            ">300 pathogenic variants.  No pseudogenes (unlike PKD1), "
            "so sequencing is straightforward.  PKD2 represents ~15% of "
            "ADPKD cases.  Milder phenotype with later ESRD."
        ),
        "crispr_strategy": (
            "Same strategies as PKD1 but simpler due to smaller gene size "
            "(2.9 kb CDS, fits in AAV) and no pseudogene interference.  "
            "(1) AAV-mediated gene supplementation to kidney.  "
            "(2) CRISPR base/prime editing for individual variants.  "
            "Kidney delivery remains the key bottleneck."
        ),
        "clinical_trials": [
            {
                "name": "No PKD2-specific gene therapy trials as of 2026",
                "status_2026": (
                    "PKD2-specific gene therapy is in preclinical stages.  "
                    "RGLS8429 miR-17 ASO (Phase I) applies to both PKD1 "
                    "and PKD2 ADPKD as it targets a common downstream "
                    "pathway."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 2c.  Alport Syndrome -- COL4A3
    # -------------------------------------------------------------------
    "COL4A3": {
        "gene_id": 1285,
        "chrom": "chr2",
        "start": 227_164_624,
        "end": 227_314_792,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2q36.3",
        "exon_count": 52,
        "role": (
            "Collagen type IV alpha-3 chain -- component of the alpha-3/4/5 "
            "collagen IV network in glomerular basement membrane (GBM), "
            "cochlea, and lens capsule.  Mutations disrupt the collagen "
            "network leading to progressive nephropathy, hearing loss, and "
            "ocular abnormalities."
        ),
        "disease": "Alport Syndrome (autosomal recessive or digenic forms)",
        "omim_disease": 203780,
        "omim_gene": 120070,
        "inheritance": "AR or AD (thin basement membrane nephropathy heterozygotes)",
        "key_variants": [
            {
                "name": "Glycine substitutions in Gly-X-Y collagenous domain",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Glycine residues in the collagenous domain are critical "
                    "for triple helix formation.  Gly>non-Gly substitutions "
                    "disrupt helix assembly.  Most are private."
                ),
            },
        ],
        "mutational_landscape": (
            ">100 pathogenic variants.  Autosomal Alport: biallelic COL4A3 "
            "or COL4A4 variants.  Heterozygous carriers may have thin "
            "basement membrane nephropathy (TBMN)."
        ),
        "crispr_strategy": (
            "Gene editing in podocyte-lineage cells to correct COL4A3/A5 "
            "variants (2019 proof of concept by Daga et al).  Challenges: "
            "large gene, tissue-specific expression, collagen assembly "
            "requires stoichiometric balance of alpha-3/4/5 chains."
        ),
        "clinical_trials": [
            {
                "name": "No approved gene therapies as of 2026",
                "status_2026": (
                    "Preclinical.  Current treatment limited to off-label "
                    "ACEi/ARB to delay ESRD.  Gene therapy research ongoing."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 2d.  Alport Syndrome -- COL4A4
    # -------------------------------------------------------------------
    "COL4A4": {
        "gene_id": 1286,
        "chrom": "chr2",
        "start": 226_967_360,
        "end": 227_164_488,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q36.3",
        "exon_count": 48,
        "role": (
            "Collagen type IV alpha-4 chain -- essential component of the "
            "alpha-3/4/5 collagen IV heterotrimer in GBM.  Adjacent to "
            "COL4A3 in a head-to-head arrangement sharing a bidirectional "
            "promoter on chromosome 2."
        ),
        "disease": "Alport Syndrome (autosomal forms)",
        "omim_disease": 203780,
        "omim_gene": 120131,
        "inheritance": "AR or AD",
        "key_variants": [
            {
                "name": "Glycine substitutions and truncating variants",
                "rsid": None,
                "consequence": "protein_truncating/missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Similar spectrum to COL4A3.  Biallelic mutations in "
                    "COL4A4 cause autosomal recessive Alport syndrome."
                ),
            },
        ],
        "mutational_landscape": (
            ">100 pathogenic variants.  COL4A3 and COL4A4 share a "
            "bidirectional promoter; both contribute to autosomal Alport."
        ),
        "crispr_strategy": (
            "Same approaches as COL4A3.  Podocyte-directed gene correction "
            "or supplementation.  Requires kidney-targeted delivery."
        ),
        "clinical_trials": [
            {
                "name": "No approved gene therapies as of 2026",
                "status_2026": "Preclinical.  ACEi/ARB remain standard of care.",
            },
        ],
    },

    # -------------------------------------------------------------------
    # 2e.  Alport Syndrome (X-linked) -- COL4A5
    # -------------------------------------------------------------------
    "COL4A5": {
        "gene_id": 1287,
        "chrom": "chrX",
        "start": 108_439_838,
        "end": 108_697_545,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xq22.3",
        "exon_count": 51,
        "role": (
            "Collagen type IV alpha-5 chain -- critical component of the "
            "alpha-3/4/5 collagen IV network in GBM, cochlea, and lens.  "
            "Mutations cause the most common form of Alport syndrome (~80% "
            "of cases are X-linked due to COL4A5)."
        ),
        "disease": "Alport Syndrome (X-linked -- most common form)",
        "omim_disease": 301050,
        "omim_gene": 303630,
        "inheritance": "XLD",
        "key_variants": [
            {
                "name": "Glycine substitutions in collagenous domain",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~40% of X-linked Alport mutations are Gly>X "
                    "substitutions.  Position of glycine substitution "
                    "correlates with age at ESRD."
                ),
            },
            {
                "name": "Large deletions (multi-exon)",
                "rsid": None,
                "consequence": "gross_deletion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~10% of mutations.  Associated with severe juvenile-onset "
                    "ESRD and often anti-GBM disease after transplant.  "
                    "Includes contiguous COL4A5-COL4A6 deletions (with "
                    "leiomyomatosis)."
                ),
            },
        ],
        "mutational_landscape": (
            ">900 pathogenic COL4A5 variants.  Most common form of Alport "
            "(~80% of cases).  Large gene (~258 kb, 51 exons).  No common "
            "founder mutation -- private variants predominate.  Hemizygous "
            "males have severe ESRD by age 20-40; heterozygous females "
            "variably affected."
        ),
        "crispr_strategy": (
            "(1) SonoThera ultrasound-mediated delivery (UMD) of full-length "
            "Col4a5 to kidney -- presented at ASN 2025 demonstrating full "
            "gene expression in murine and NHP models.  Novel non-viral "
            "delivery technology.  (2) CRISPR/Cas9 editing of podocytes "
            "ex vivo (proof of concept 2019).  (3) Antisense oligonucleotide "
            "splice correction for deep-intronic variants (JCI Insight 2025).  "
            "Key challenge: very large gene (258 kb genomic, 6.4 kb CDS); "
            "exceeds AAV packaging capacity.  Non-viral delivery critical."
        ),
        "clinical_trials": [
            {
                "name": "SonoThera UMD-COL4A5",
                "sponsor": "SonoThera",
                "phase": "Preclinical",
                "approach": (
                    "Ultrasound-mediated delivery of full-length COL4A5 "
                    "genetic medicine to kidney.  Demonstrated expression "
                    "in murine and NHP models."
                ),
                "status_2026": (
                    "Preclinical.  Data presented at ASN Kidney Week "
                    "November 2025.  First technology to demonstrate "
                    "full-length COL4A5 expression in kidney in vivo."
                ),
            },
        ],
    },
}


# ============================================================================
# 3. ENDOCRINE DISORDERS
# ============================================================================

ENDOCRINE_TARGETS = {

    # -------------------------------------------------------------------
    # 3a.  Congenital Adrenal Hyperplasia -- CYP21A2
    # -------------------------------------------------------------------
    "CYP21A2": {
        "gene_id": 1589,
        "chrom": "chr6",
        "start": 32_038_415,
        "end": 32_041_644,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6p21.33",
        "exon_count": 10,
        "role": (
            "Steroid 21-hydroxylase -- cytochrome P450 enzyme in adrenal "
            "cortex that converts 17-hydroxyprogesterone to "
            "11-deoxycortisol (cortisol pathway) and progesterone to "
            "deoxycorticosterone (aldosterone pathway).  Deficiency causes "
            "cortisol deficiency, aldosterone deficiency (salt-wasting), "
            "and androgen excess."
        ),
        "disease": "Congenital Adrenal Hyperplasia (CAH) due to 21-hydroxylase deficiency",
        "omim_disease": 201910,
        "omim_gene": 613815,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "IVS2 splice variant (c.293-13A/C>G or c.293-13C>G)",
                "rsid": "rs6467",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common CYP21A2 mutation (~25-35% of alleles).  "
                    "Intron 2 splice variant; allows some residual enzyme "
                    "activity.  Associated with simple virilizing form."
                ),
            },
            {
                "name": "c.955C>T (p.Gln319Ter / Q318X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Null allele.  Found in salt-wasting form.  Relatively "
                    "common (~10% of alleles)."
                ),
            },
            {
                "name": "c.844G>T (p.Val282Leu / V281L)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common variant in non-classic (late-onset) CAH.  "
                    "~50% residual enzyme activity.  Very common in "
                    "Ashkenazi Jewish population."
                ),
            },
            {
                "name": "Large gene deletion / gene conversion with CYP21A1P",
                "rsid": None,
                "consequence": "gross_deletion/gene_conversion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~20-25% of alleles.  CYP21A2 lies in the HLA class III "
                    "region near its pseudogene CYP21A1P; non-allelic "
                    "homologous recombination causes deletions and gene "
                    "conversions transferring pseudogene mutations to the "
                    "active gene."
                ),
            },
        ],
        "mutational_landscape": (
            "~95% of CYP21A2 mutations arise from gene conversions with "
            "the adjacent pseudogene CYP21A1P or large deletions due to "
            "NAHR in the RCCX module of the MHC.  Only ~5% are de novo.  "
            "About 10 mutations account for >90% of disease alleles.  "
            "CAUTION: CYP21A2 lies in the HLA class III region at 6p21.33, "
            "making sequencing and CRISPR design challenging due to "
            "pseudogene CYP21A1P (~98% sequence identity)."
        ),
        "crispr_strategy": (
            "(1) In vivo AAV5-CYP21A2 gene therapy (BBP-631) -- deliver "
            "functional CYP21A2 transgene to adrenal glands.  "
            "(2) Genomic editing of adrenal cortical cells -- challenging "
            "due to pseudogene homology.  (3) Ex vivo iPSC-derived adrenal "
            "cell replacement.  (4) CRISPR base editing for point mutations "
            "like V281L (G>T, could use CBE on antisense strand).  "
            "Key challenge: pseudogene CYP21A1P homology complicates guide "
            "design; adrenal-specific delivery is underdeveloped."
        ),
        "clinical_trials": [
            {
                "name": "BBP-631 (Adrenas Therapeutics)",
                "sponsor": "Adrenas Therapeutics",
                "nct": "NCT04783181",
                "phase": "I/II (HALTED)",
                "approach": (
                    "AAV5 vector containing human CYP21A2 transgene, "
                    "intravenous delivery.  Single dose designed to provide "
                    "endogenous cortisol production."
                ),
                "status_2026": (
                    "Phase I/II enrolled 7 adults.  BBP-631 was well "
                    "tolerated with increases in 11-deoxycortisol and some "
                    "cortisol production at higher doses.  However, study "
                    "closed to enrollment and further development of "
                    "BBP-631 has been HALTED.  Higher doses will not be "
                    "tested in this program."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 3b.  Type 1 Diabetes -- INS, HLA region
    # -------------------------------------------------------------------
    "INS": {
        "gene_id": 3630,
        "chrom": "chr11",
        "start": 2_159_779,
        "end": 2_161_209,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11p15.5",
        "exon_count": 3,
        "role": (
            "Insulin -- peptide hormone secreted by pancreatic beta cells.  "
            "Essential for glucose homeostasis.  In type 1 diabetes, "
            "autoimmune destruction of beta cells leads to absolute insulin "
            "deficiency.  NOTE: T1D is predominantly autoimmune/polygenic, "
            "not monogenic; INS is listed for completeness."
        ),
        "disease": "Type 1 Diabetes (autoimmune, polygenic)",
        "omim_disease": 222100,
        "omim_gene": 176730,
        "inheritance": "Complex/polygenic (HLA + INS + >60 other loci)",
        "key_variants": [
            {
                "name": "INS-VNTR (variable number tandem repeat, 5' of INS)",
                "rsid": None,
                "consequence": "regulatory",
                "clinical_significance": "risk_allele",
                "notes": (
                    "Class I VNTR alleles (26-63 repeats) associated with "
                    "increased T1D susceptibility; Class III (140-210 "
                    "repeats) protective.  Affects INS expression in thymus "
                    "and immune tolerance."
                ),
            },
            {
                "name": "HLA-DRB1*03:01, HLA-DRB1*04:01 (HLA class II)",
                "rsid": None,
                "consequence": "HLA_haplotype",
                "clinical_significance": "risk_allele",
                "notes": (
                    "HLA region accounts for ~50% of T1D genetic risk.  "
                    "DR3/DR4 heterozygosity confers highest risk (OR ~15).  "
                    "NOT a single-gene target for CRISPR -- complex region."
                ),
            },
        ],
        "mutational_landscape": (
            "T1D is NOT a monogenic disorder.  >60 susceptibility loci "
            "identified by GWAS.  HLA accounts for ~50% of heritability.  "
            "INS-VNTR is the second strongest genetic locus.  Monogenic "
            "neonatal diabetes (MODY) involves different genes (KCNJ11, "
            "ABCC8, INS mutations).  CRISPR approaches for T1D focus on "
            "cell replacement rather than gene correction."
        ),
        "crispr_strategy": (
            "(1) CRISPR-edited stem cell-derived beta cells -- VCTX210 "
            "(CRISPR Therapeutics/ViaCyte): iPSC-derived beta cells with "
            "CRISPR-mediated immune evasion (B2M knockout, CIITA knockout, "
            "CD47 overexpression) to prevent rejection without "
            "immunosuppression.  (2) VCTX211 with additional A20/TNFAIP3 "
            "and MANF gene inserts for cytokine protection.  "
            "(3) Encapsulated cell therapy.  "
            "Key consideration: T1D gene therapy is about creating "
            "immune-evasive insulin-producing cells, not correcting a "
            "single gene mutation."
        ),
        "clinical_trials": [
            {
                "name": "VCTX210",
                "sponsor": "CRISPR Therapeutics / ViaCyte (now Vertex)",
                "phase": "I",
                "approach": (
                    "Gene-edited, allogeneic, stem cell-derived pancreatic "
                    "islet cells.  CRISPR edits make cells hypoimmunogenic "
                    "(evade immune destruction without immunosuppression)."
                ),
                "status_2026": (
                    "Phase I; first patient dosed January 2026.  "
                    "Breakthrough: 2025 first-in-human data showed "
                    "gene-edited cells engrafted in one T1D patient without "
                    "immune rejection and without immunosuppressants "
                    "(published in Nature 2025).  Paradigm-shifting result."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 3c.  Growth Hormone Deficiency -- GH1
    # -------------------------------------------------------------------
    "GH1": {
        "gene_id": 2688,
        "chrom": "chr17",
        "start": 63_917_203,
        "end": 63_918_839,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q23.3",
        "exon_count": 5,
        "role": (
            "Growth hormone 1 (somatotropin) -- pituitary peptide hormone "
            "that stimulates growth, cell reproduction, and IGF-1 "
            "production.  Part of the GH gene cluster (GH1, GH2, CSH1, "
            "CSH2, CSHL1) on chromosome 17."
        ),
        "disease": "Isolated Growth Hormone Deficiency (IGHD)",
        "omim_disease": 262400,  # type IA
        "omim_gene": 139250,
        "inheritance": "AR (type IA/IB), AD (type II), XLR (type III)",
        "key_variants": [
            {
                "name": "GH1 gene deletion (6.7 kb -- type IA)",
                "rsid": None,
                "consequence": "gross_deletion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Homozygous deletion of GH1 due to NAHR between "
                    "flanking repeat sequences in the GH gene cluster.  "
                    "Causes type IA IGHD (most severe; absent GH, anti-GH "
                    "antibody formation with exogenous GH treatment)."
                ),
            },
            {
                "name": "c.291+1G>A (IVS3+1G>A -- type II AD)",
                "rsid": None,
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common cause of autosomal dominant IGHD type II.  "
                    "Causes exon 3 skipping; 17.5 kDa isoform acts as "
                    "dominant-negative, disrupting GH secretory granules."
                ),
            },
            {
                "name": "c.172G>A (p.Asp58Asn)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Bioinactive GH; normal immunoreactive GH but no bioactivity.",
            },
        ],
        "mutational_landscape": (
            "IGHD is rare (~1:4000-10,000 births).  Most cases are "
            "idiopathic/sporadic.  Monogenic forms: type IA (AR, GH1 "
            "deletion -- severest), type IB (AR, point mutations with "
            "some residual GH), type II (AD, splice mutations), "
            "type III (XLR, BTK region -- not GH1).  GH gene cluster "
            "homology makes NAHR-mediated deletions common."
        ),
        "crispr_strategy": (
            "(1) AAV-mediated GH1 gene replacement to pituitary or liver "
            "(for secreted GH).  Small gene (0.8 kb CDS) easily fits AAV.  "
            "(2) Base editing to correct splice mutations like IVS3+1G>A "
            "(G>A, amenable to ABE).  (3) Ex vivo iPSC-derived "
            "somatotrophs.  Key challenge: pituitary delivery is difficult; "
            "liver-directed secreted GH approach may be more practical "
            "but lacks pulsatile secretion physiology."
        ),
        "clinical_trials": [
            {
                "name": "No gene therapy clinical trials as of 2026",
                "status_2026": (
                    "GH1 gene therapy is in preclinical stages.  "
                    "Recombinant human GH (somatropin -- multiple brands) "
                    "and long-acting GH analogs (somapacitan/Sogroya, "
                    "lonapegsomatropin/Skytrofa) are standard of care.  "
                    "Gene therapy less urgent given effective replacement "
                    "therapies available."
                ),
            },
        ],
    },
}


# ============================================================================
# 4. HEARING LOSS
# ============================================================================

HEARING_LOSS_TARGETS = {

    # -------------------------------------------------------------------
    # 4a.  DFNB1 (Most common genetic deafness) -- GJB2
    # -------------------------------------------------------------------
    "GJB2": {
        "gene_id": 2706,
        "chrom": "chr13",
        "start": 20_187_470,
        "end": 20_192_938,
        "strand": "-",
        "refseq": "NC_000013.11",
        "cytoband": "13q12.11",
        "exon_count": 2,
        "role": (
            "Connexin 26 -- gap junction protein forming intercellular "
            "channels in the cochlea between supporting cells.  Critical "
            "for potassium recycling from hair cells back to endolymph.  "
            "Mutations disrupt cochlear potassium homeostasis causing "
            "sensorineural hearing loss."
        ),
        "disease": "DFNB1 Non-syndromic Hearing Loss (most common genetic deafness)",
        "omim_disease": 220290,
        "omim_gene": 121011,
        "inheritance": "AR (DFNB1A) or AD (DFNA3A -- rare)",
        "key_variants": [
            {
                "name": "c.35delG (p.Gly12ValfsTer2)",
                "rsid": "rs80338939",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common deafness-causing mutation worldwide.  "
                    "Predominant in European, North African, and Middle "
                    "Eastern populations.  Carrier frequency ~2-4% in "
                    "Southern Europe, ~1-2% in Northern Europe.  "
                    "Null allele; homozygotes typically have severe-profound "
                    "hearing loss."
                ),
            },
            {
                "name": "c.235delC",
                "rsid": "rs80338943",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common GJB2 mutation in East Asian populations "
                    "(~69% of Chinese GJB2 alleles).  Null allele; "
                    "homozygotes have severe-profound hearing loss."
                ),
            },
            {
                "name": "c.167delT",
                "rsid": "rs80338942",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Common in Ashkenazi Jewish population (carrier "
                    "frequency ~4%).  Null allele."
                ),
            },
            {
                "name": "c.109G>A (p.Val37Ile / V37I)",
                "rsid": "rs72474224",
                "consequence": "missense",
                "clinical_significance": "pathogenic (mild)",
                "notes": (
                    "Most common in East/Southeast Asian populations.  "
                    "Associated with mild-to-moderate hearing loss.  "
                    "Some residual connexin function retained."
                ),
            },
        ],
        "mutational_landscape": (
            "GJB2 mutations are the #1 cause of genetic sensorineural "
            "hearing loss worldwide, accounting for ~50% of all prelingual "
            "genetic deafness.  c.35delG dominates in Europeans; c.235delC "
            "in East Asians; c.167delT in Ashkenazi Jews.  Truncating "
            "variants cause severe-profound loss; non-truncating (V37I) "
            "cause mild-moderate.  GJB2 has only 1 coding exon (exon 2), "
            "making it compact for gene therapy."
        ),
        "crispr_strategy": (
            "(1) AAV-mediated gene addition with cell-specific promoters -- "
            "2025: co-administration of AAV1 and AAV-ie with SCpro (supporting "
            "cell-specific promoter) restored hearing in Gjb2-deficient mice.  "
            "(2) Adenine base editing (ABE) delivered via AAV-Sia6e (Japan, "
            "April 2025) -- single-base correction small enough for AAV.  "
            "(3) CRISPR-Cas9 allele-specific disruption of dominant-negative "
            "alleles (for DFNA3).  Key challenge: GJB2 expression is required "
            "in supporting cells (not hair cells) -- cell-type-specific "
            "targeting needed.  Timing critical: must treat before cochlear "
            "development is complete."
        ),
        "clinical_trials": [
            {
                "name": "AAV-GJB2 gene therapy (preclinical/early clinical)",
                "sponsor": "Multiple academic groups (China, Japan, UC Irvine)",
                "phase": "Preclinical / Phase I (emerging)",
                "approach": (
                    "AAV-mediated delivery of functional GJB2 to cochlear "
                    "supporting cells.  Various serotypes and promoters "
                    "under investigation."
                ),
                "status_2026": (
                    "Preclinical in 2025-2026.  First clinical trial for "
                    "hereditary deafness gene therapy (10 participants ages "
                    "1.5-23.9 yrs) demonstrated safety and efficacy, though "
                    "this was for OTOF/DFNB9, not GJB2 specifically.  "
                    "GJB2 trials expected to follow given strong preclinical "
                    "data.  ABE approach (Japan 2025) designed specifically "
                    "for AAV-compatible size."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 4b.  DFNA9 -- COCH
    # -------------------------------------------------------------------
    "COCH": {
        "gene_id": 1690,
        "chrom": "chr14",
        "start": 30_874_559,
        "end": 30_895_615,
        "strand": "+",
        "refseq": "NC_000014.9",
        "cytoband": "14q12",
        "exon_count": 12,
        "role": (
            "Cochlin -- secreted extracellular matrix protein; most "
            "abundant protein in the inner ear.  Contains LCCL domain and "
            "two vWFA domains.  Expressed in spiral ligament, spiral limbus, "
            "and vestibular organs.  Dominant-negative mutations cause "
            "accumulation of misfolded cochlin leading to progressive "
            "hearing loss and vestibular dysfunction."
        ),
        "disease": "DFNA9 (Autosomal Dominant Non-syndromic Hearing Loss + vestibular dysfunction)",
        "omim_disease": 601369,
        "omim_gene": 603196,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "c.151C>T (p.Pro51Ser / P51S)",
                "rsid": "rs121908026",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "LCCL domain mutation; one of the most commonly "
                    "reported DFNA9 variants.  Belgian founder mutation.  "
                    "Onset 3rd-4th decade."
                ),
            },
            {
                "name": "c.317T>C (p.Ile109Asn / I109N)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "LCCL domain; progressive high-frequency loss + vestibular.",
            },
            {
                "name": "c.259G>A (p.Gly87Val / G87V)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "LCCL domain; identified in Chinese family.",
            },
        ],
        "mutational_landscape": (
            ">30 pathogenic COCH variants, mostly missense in the LCCL "
            "domain (exons 4-5).  Dominant gain-of-function/dominant-negative "
            "mechanism: misfolded cochlin aggregates in the inner ear.  "
            "Progressive hearing loss starting in adulthood (20s-50s), "
            "often with vestibular symptoms."
        ),
        "crispr_strategy": (
            "(1) Allele-specific knockdown using ASOs targeting the mutant "
            "COCH mRNA while preserving wild-type -- lead therapeutic "
            "approach (RNID-funded research 2024-2025).  (2) CRISPR "
            "allele-specific silencing using SNP-guided Cas9 or CRISPRi.  "
            "(3) CRISPR-mediated disruption of the mutant allele only.  "
            "Key consideration: dominant gain-of-function mechanism means "
            "gene addition alone won't work; must silence/destroy the "
            "mutant allele."
        ),
        "clinical_trials": [
            {
                "name": "No clinical trials as of 2026",
                "status_2026": (
                    "Preclinical.  RNID-funded research developing ASOs "
                    "for DFNA9.  First preclinical study published.  "
                    "Results expected to support future clinical trial "
                    "applications.  No approved treatments."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 4c.  Pendred Syndrome -- SLC26A4
    # -------------------------------------------------------------------
    "SLC26A4": {
        "gene_id": 5172,
        "chrom": "chr7",
        "start": 107_660_828,
        "end": 107_717_809,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q22.3",
        "exon_count": 21,
        "role": (
            "Pendrin -- transmembrane anion transporter (Cl-/HCO3-/I- "
            "exchanger) expressed in inner ear (endolymphatic sac, cochlea, "
            "vestibule), thyroid, and kidney.  In the ear, maintains "
            "endolymphatic pH and ionic homeostasis.  In thyroid, mediates "
            "iodide efflux into follicular lumen."
        ),
        "disease": (
            "Pendred Syndrome (SNHL + goiter) / DFNB4 (non-syndromic "
            "hearing loss with enlarged vestibular aqueduct)"
        ),
        "omim_disease": 274600,  # Pendred
        "omim_gene": 605646,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.919-2A>G (IVS7-2A>G)",
                "rsid": "rs111033313",
                "consequence": "splice_acceptor",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common SLC26A4 mutation worldwide.  Particularly "
                    "prevalent in East Asian populations.  Null allele."
                ),
            },
            {
                "name": "c.2168A>G (p.His723Arg / H723R)",
                "rsid": "rs121908362",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Very common in Japanese/Korean populations.  ~50% of "
                    "SLC26A4 alleles in Japanese Pendred/EVA patients."
                ),
            },
            {
                "name": "c.1246A>C (p.Thr416Pro / T416P)",
                "rsid": "rs111033206",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Common in European and Middle Eastern populations.",
            },
            {
                "name": "c.1001+1G>A (IVS8+1G>A)",
                "rsid": None,
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "notes": "Common in Northern European populations.",
            },
        ],
        "mutational_landscape": (
            ">500 pathogenic/likely-pathogenic SLC26A4 variants.  "
            "Second most common cause of genetic deafness (after GJB2) "
            "in many Asian populations.  Hallmark radiological finding: "
            "enlarged vestibular aqueduct (EVA) on CT/MRI.  Goiter in "
            "Pendred (not in DFNB4).  Hearing loss may be fluctuating "
            "or progressive."
        ),
        "crispr_strategy": (
            "(1) AAV-mediated gene therapy -- 2019 mouse study demonstrated "
            "feasibility of SLC26A4 gene delivery to inner ear, preventing "
            "endolymphatic hydrops and restoring hearing (but not vestibular "
            "function).  (2) Base editing for common splice variants "
            "(IVS7-2A>G amenable to ABE).  (3) AAV delivery to "
            "endolymphatic sac.  Key challenge: spatially and temporally "
            "limited expression from current vectors; pendrin must be "
            "expressed in specific inner ear cell types.  Timing critical: "
            "treatment before inner ear structural damage."
        ),
        "clinical_trials": [
            {
                "name": "No clinical trials as of 2026",
                "status_2026": (
                    "Preclinical.  Mouse proof-of-concept (2019) showed "
                    "hearing restoration but not vestibular rescue.  "
                    "Multiple preclinical programs advancing.  "
                    "Cochlear implants remain standard of care."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 4d.  OTOF-Related Deafness (DFNB9) -- OTOF
    # -------------------------------------------------------------------
    "OTOF": {
        "gene_id": 9381,
        "chrom": "chr2",
        "start": 26_457_203,
        "end": 26_558_756,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2p23.3",
        "exon_count": 48,
        "role": (
            "Otoferlin -- multi-C2 domain calcium sensor protein critical "
            "for synaptic vesicle exocytosis at inner hair cell (IHC) "
            "ribbon synapses.  Acts as the Ca2+-sensor for synaptic "
            "vesicle fusion at the IHC-spiral ganglion neuron synapse.  "
            "Loss of function abolishes neurotransmitter release, causing "
            "'auditory synaptopathy' -- OAEs preserved but ABR absent."
        ),
        "disease": "DFNB9 (Autosomal Recessive Non-syndromic Deafness -- auditory neuropathy/synaptopathy)",
        "omim_disease": 601071,
        "omim_gene": 603681,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.5098G>C (p.Glu1700Gln -- founder in Spanish population)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Identified in Spanish families.  One of the more "
                    "common recurrent variants."
                ),
            },
            {
                "name": "c.2485C>T (p.Gln829Ter / Q829X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Truncating; null allele.  Found in multiple populations.",
            },
            {
                "name": "c.5473C>G (p.Pro1825Ala)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "C2F domain; disrupts calcium-sensing function.",
            },
            {
                "name": "c.5567G>A (p.Arg1856Gln)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "C2F domain variant.",
            },
        ],
        "mutational_landscape": (
            "~220 pathogenic/likely-pathogenic OTOF variants.  Distribution: "
            "84 missense, 44 frameshift, 43 nonsense, 36 splice, 7 in-frame "
            "indels, 3 CNVs.  Variants span all 48 exons.  Large gene "
            "(~6 kb CDS) exceeds single AAV packaging capacity (~4.7 kb), "
            "requiring dual-AAV strategies.  DFNB9 accounts for ~2-8% of "
            "prelingual deafness; up to 91% of auditory neuropathy spectrum "
            "disorder (ANSD) in Korean cohorts."
        ),
        "crispr_strategy": (
            "(1) Dual-AAV gene replacement (DB-OTO) -- split OTOF cDNA "
            "across two AAV1 vectors that recombine in the target cell via "
            "overlap or trans-splicing.  Intracochlear delivery.  THIS IS "
            "THE LEADING GENE THERAPY FOR GENETIC HEARING LOSS.  "
            "(2) Base/prime editing for individual missense variants.  "
            "(3) Single-AAV with truncated/mini OTOF constructs.  "
            "Key advantage: only inner hair cells need correction (small "
            "target population); intracochlear injection is feasible and "
            "established surgically."
        ),
        "clinical_trials": [
            {
                "name": "DB-OTO (CHORD trial)",
                "sponsor": "Regeneron Pharmaceuticals (acquired Decibel Therapeutics 2023 for $109M)",
                "nct": "NCT05788536",
                "phase": "I/II (pivotal)",
                "approach": (
                    "Dual-AAV1 gene therapy delivering full-length OTOF cDNA "
                    "via intracochlear infusion.  Single administration.  "
                    "Two AAV1 vectors that recombine intracellularly to "
                    "produce functional otoferlin protein in inner hair cells."
                ),
                "status_2026": (
                    "PHASE I/II CHORD TRIAL -- BREAKTHROUGH RESULTS.  "
                    "Published in NEJM (October 2025).  12 pediatric "
                    "participants (10 months to 16 years) with biallelic "
                    "OTOF variants.  KEY RESULTS: 11/12 had clinically "
                    "meaningful hearing improvement; 3 achieved NORMAL "
                    "hearing levels; 8 with longer follow-up showed "
                    "stability or continued improvement.  Met primary "
                    "endpoint: 9/12 achieved <=70 dBHL at week 24.  "
                    "3 bilateral recipients also showed improvement.  "
                    "Among 3 with speech assessments, all showed significant "
                    "improvement.  Regeneron plans FDA submission by end "
                    "of 2025/early 2026.  POTENTIALLY FIRST APPROVED "
                    "GENE THERAPY FOR HEARING LOSS."
                ),
            },
            {
                "name": "Chinese dual-AAV OTOF trials (multiple sponsors)",
                "sponsor": "Fudan University / multiple Chinese companies",
                "phase": "I/II",
                "approach": (
                    "Dual-AAV OTOF gene therapy, intracochlear delivery.  "
                    "Multiple independent trials in China."
                ),
                "status_2026": (
                    "Multiple Phase I/II trials in China reporting positive "
                    "results.  Bilateral gene therapy in children with DFNB9 "
                    "published in Nature Medicine (2024).  China is leading "
                    "in parallel with Regeneron's CHORD trial."
                ),
            },
        ],
    },
}


# ============================================================================
# 5. CONNECTIVE TISSUE DISORDERS
# ============================================================================

CONNECTIVE_TISSUE_TARGETS = {

    # -------------------------------------------------------------------
    # 5a.  Osteogenesis Imperfecta -- COL1A1
    # -------------------------------------------------------------------
    "COL1A1": {
        "gene_id": 1277,
        "chrom": "chr17",
        "start": 50_184_101,
        "end": 50_201_631,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.33",
        "exon_count": 51,
        "role": (
            "Collagen type I alpha-1 chain -- major structural protein of "
            "bone, skin, tendons, and ligaments.  Forms heterotrimers "
            "(2x alpha-1 + 1x alpha-2) that assemble into fibrils.  "
            "Mutations disrupt collagen triple helix formation, leading "
            "to bone fragility."
        ),
        "disease": "Osteogenesis Imperfecta (types I-IV, primarily)",
        "omim_disease": 166200,  # OI type I
        "omim_gene": 120150,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Glycine substitutions in Gly-X-Y triple helix domain",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common OI mutation class.  Gly>Ser, Gly>Cys, "
                    "Gly>Arg, Gly>Asp, etc.  Glycine substitutions disrupt "
                    "triple helix formation (dominant-negative).  Position "
                    "and substitution type determine severity: Gly>Asp more "
                    "severe than Gly>Ser; C-terminal substitutions generally "
                    "more severe in COL1A1."
                ),
            },
            {
                "name": "Null alleles (haploinsufficiency -- OI type I)",
                "rsid": None,
                "consequence": "protein_truncating",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Nonsense, frameshift, or splice mutations causing "
                    "nonsense-mediated mRNA decay of one allele.  "
                    "Haploinsufficiency produces quantitatively reduced "
                    "but structurally normal collagen.  Mildest OI (type I): "
                    "blue sclerae, moderate fractures, normal stature."
                ),
            },
        ],
        "mutational_landscape": (
            ">1500 pathogenic COL1A1 variants.  ~80-90% of OI caused by "
            "COL1A1 or COL1A2.  Two main mechanisms: (1) Dominant-negative "
            "glycine substitutions (structural defect, types II-IV); "
            "(2) Haploinsufficiency (quantitative defect, type I).  No "
            "common founder mutation; most are private."
        ),
        "crispr_strategy": (
            "(1) CRISPR/Cas9 allele-specific disruption of dominant-negative "
            "mutant allele -- convert OI type II-IV to milder type I "
            "(haploinsufficiency).  Multiple proof-of-concept studies in "
            "iPSCs (2021-2023).  (2) Base editing for specific Gly>X "
            "substitutions.  (3) AAV-delivered CRISPR gene editing in "
            "osteoblasts.  Key challenge: must target osteoblasts in bone; "
            "systemic delivery to bone is difficult.  Dominant-negative "
            "mechanism means gene addition alone insufficient -- must "
            "silence or correct the mutant allele."
        ),
        "clinical_trials": [
            {
                "name": "No gene therapy clinical trials as of 2026",
                "status_2026": (
                    "Preclinical.  CRISPR correction of COL1A1 in patient "
                    "iPSCs demonstrated in multiple studies (2021-2025).  "
                    "Bisphosphonates and anti-RANKL antibodies are current "
                    "standard of care.  Clinical translation limited by "
                    "bone-targeted delivery challenges."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 5b.  Osteogenesis Imperfecta -- COL1A2
    # -------------------------------------------------------------------
    "COL1A2": {
        "gene_id": 1278,
        "chrom": "chr7",
        "start": 94_394_895,
        "end": 94_431_227,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q21.3",
        "exon_count": 52,
        "role": (
            "Collagen type I alpha-2 chain -- one copy per collagen I "
            "heterotrimer (2x alpha-1 + 1x alpha-2).  Mutations cause "
            "similar OI spectrum as COL1A1 but with some genotype-phenotype "
            "differences."
        ),
        "disease": "Osteogenesis Imperfecta",
        "omim_disease": 166210,  # OI type II
        "omim_gene": 120160,
        "inheritance": "AD (dominant-negative) or AR (rare biallelic null = OI-like EDS)",
        "key_variants": [
            {
                "name": "Glycine substitutions in Gly-X-Y triple helix domain",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Same mechanism as COL1A1.  In COL1A2, Gly>Asp "
                    "substitutions are particularly severe.  N-terminal "
                    "substitutions tend to be more severe in COL1A2 "
                    "(opposite pattern from COL1A1)."
                ),
            },
            {
                "name": "Homozygous null COL1A2 alleles",
                "rsid": None,
                "consequence": "protein_truncating",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Rare AR form.  Complete absence of alpha-2 chains "
                    "causes collagen I to form alpha-1 homotrimers.  "
                    "Phenotype: EDS/OI overlap (joint hypermobility, "
                    "skin fragility, moderate bone fragility)."
                ),
            },
        ],
        "mutational_landscape": (
            ">900 pathogenic COL1A2 variants.  Similar spectrum to COL1A1 "
            "but fewer overall variants.  COL1A2 glycine substitutions "
            "generally cause milder OI than equivalent COL1A1 positions "
            "because there is only one alpha-2 chain per trimer vs. two "
            "alpha-1 chains."
        ),
        "crispr_strategy": (
            "Same strategies as COL1A1: (1) allele-specific CRISPR "
            "disruption of dominant-negative allele; (2) base editing for "
            "glycine substitutions; (3) gene correction in osteoblast "
            "progenitors.  AAV-delivered rAAV gene editing (Molecular "
            "Therapy Nucleic Acids, 2023) showed correction of Col1a2 "
            "mutation in osteoblast-lineage cells in vivo in mouse models."
        ),
        "clinical_trials": [
            {
                "name": "No gene therapy clinical trials as of 2026",
                "status_2026": (
                    "Preclinical.  In vivo rAAV-CRISPR correction of "
                    "Col1a2 mutations demonstrated in mouse models (2023).  "
                    "Improved bone microarchitecture, mechanical strength, "
                    "and collagen quality documented.  Clinical translation "
                    "pending bone-targeted delivery development."
                ),
            },
        ],
    },
}


# ============================================================================
# COMBINED TARGETS DICTIONARY
# ============================================================================

ALL_TARGETS = {}
ALL_TARGETS.update(LYSOSOMAL_STORAGE_TARGETS)
ALL_TARGETS.update(KIDNEY_RENAL_TARGETS)
ALL_TARGETS.update(ENDOCRINE_TARGETS)
ALL_TARGETS.update(HEARING_LOSS_TARGETS)
ALL_TARGETS.update(CONNECTIVE_TISSUE_TARGETS)


# ============================================================================
# QUICK-REFERENCE: BED-format coordinates for pipeline use
# ============================================================================

def get_bed_regions(padding: int = 0) -> list[tuple[str, int, int, str]]:
    """Return (chrom, start, end, gene_name) tuples in BED format.

    Parameters
    ----------
    padding : int
        Number of bases to add upstream and downstream of each gene.

    Returns
    -------
    list of (chrom, start, end, gene_name) tuples sorted by chrom+start.
    """
    rows = []
    for gene, info in ALL_TARGETS.items():
        rows.append((
            info["chrom"],
            max(0, info["start"] - padding),
            info["end"] + padding,
            gene,
        ))
    return sorted(rows, key=lambda r: (r[0], r[1]))


# Add conditions arrays for pipeline condition-based filtering
_CONDITIONS = {
    "GBA1":    ["gaucher_disease", "lysosomal_storage", "metabolic"],
    "GLA":     ["fabry_disease", "lysosomal_storage", "metabolic"],
    "GAA":     ["pompe_disease", "lysosomal_storage", "metabolic"],
    "IDUA":    ["MPS_I", "hurler_syndrome", "lysosomal_storage", "metabolic"],
    "IDS":     ["MPS_II", "hunter_syndrome", "lysosomal_storage", "metabolic"],
    "NPC1":    ["niemann_pick_C", "lysosomal_storage", "metabolic"],
    "HEXA":    ["tay_sachs", "lysosomal_storage", "metabolic"],
    "GALC":    ["krabbe_disease", "lysosomal_storage", "leukodystrophy"],
    "ARSA":    ["MLD", "metachromatic_leukodystrophy", "lysosomal_storage",
                "leukodystrophy"],
    "PKD1":    ["polycystic_kidney", "ADPKD", "kidney_disease"],
    "PKD2":    ["polycystic_kidney", "ADPKD", "kidney_disease"],
    "COL4A3":  ["alport_syndrome", "kidney_disease"],
    "COL4A4":  ["alport_syndrome", "kidney_disease"],
    "COL4A5":  ["alport_syndrome", "kidney_disease"],
    "CYP21A2": ["congenital_adrenal_hyperplasia", "CAH", "endocrine"],
    "INS":     ["type_1_diabetes", "T1D", "endocrine"],
    "GH1":     ["growth_hormone_deficiency", "endocrine"],
    "GJB2":    ["DFNB1", "hearing_loss", "deafness"],
    "COCH":    ["DFNA9", "hearing_loss", "deafness"],
    "SLC26A4": ["pendred_syndrome", "hearing_loss", "deafness"],
    "OTOF":    ["DFNB9", "hearing_loss", "deafness"],
    "COL1A1":  ["osteogenesis_imperfecta", "OI", "connective_tissue"],
    "COL1A2":  ["osteogenesis_imperfecta", "OI", "connective_tissue"],
}
for _gene, _conds in _CONDITIONS.items():
    if _gene in ALL_TARGETS:
        ALL_TARGETS[_gene]["conditions"] = _conds


def get_gene_info(gene_symbol: str) -> dict | None:
    """Look up a gene by symbol (case-insensitive)."""
    for key, val in ALL_TARGETS.items():
        if key.upper() == gene_symbol.upper():
            return val
    return None


# ============================================================================
# Summary statistics (for quick reference)
# ============================================================================

if __name__ == "__main__":
    print(f"Total genes in database: {len(ALL_TARGETS)}")
    print(f"\nLysosomal Storage Diseases: {len(LYSOSOMAL_STORAGE_TARGETS)} genes")
    for g in LYSOSOMAL_STORAGE_TARGETS:
        info = LYSOSOMAL_STORAGE_TARGETS[g]
        print(f"  {g}: {info['chrom']}:{info['start']}-{info['end']} ({info.get('disease', 'N/A')})")
    print(f"\nKidney/Renal: {len(KIDNEY_RENAL_TARGETS)} genes")
    for g in KIDNEY_RENAL_TARGETS:
        info = KIDNEY_RENAL_TARGETS[g]
        print(f"  {g}: {info['chrom']}:{info['start']}-{info['end']} ({info.get('disease', 'N/A')})")
    print(f"\nEndocrine: {len(ENDOCRINE_TARGETS)} genes")
    for g in ENDOCRINE_TARGETS:
        info = ENDOCRINE_TARGETS[g]
        print(f"  {g}: {info['chrom']}:{info['start']}-{info['end']} ({info.get('disease', 'N/A')})")
    print(f"\nHearing Loss: {len(HEARING_LOSS_TARGETS)} genes")
    for g in HEARING_LOSS_TARGETS:
        info = HEARING_LOSS_TARGETS[g]
        print(f"  {g}: {info['chrom']}:{info['start']}-{info['end']} ({info.get('disease', 'N/A')})")
    print(f"\nConnective Tissue: {len(CONNECTIVE_TISSUE_TARGETS)} genes")
    for g in CONNECTIVE_TISSUE_TARGETS:
        info = CONNECTIVE_TISSUE_TARGETS[g]
        print(f"  {g}: {info['chrom']}:{info['start']}-{info['end']} ({info.get('disease', 'N/A')})")
