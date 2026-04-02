"""
Curated gene-therapy targets for primary immunodeficiencies, cystic fibrosis,
and metabolic / liver diseases -- GRCh38 (hg38) coordinates.

Sources:  NCBI Gene, dbSNP, ClinVar, Ensembl (GRCh38.p14), OMIM, CFF,
          Innovative Genomics Institute clinical-trial tracker, PubMed.

Coordinates pulled from NCBI Gene (RefSeq annotation RS_2025_08) on 2026-04-02.
Variant positions pulled from NCBI dbSNP / ClinVar on the same date.

IMPORTANT -- verify every coordinate against current NCBI / Ensembl releases
before production use.  Numbering can shift between patch levels.

All interventions require informed patient consent and IRB / ethics approval.
"""

# ============================================================================
# 1. PRIMARY IMMUNODEFICIENCIES
# ============================================================================

PRIMARY_IMMUNODEFICIENCY_TARGETS = {

    # -------------------------------------------------------------------
    # 1a.  X-linked Severe Combined Immunodeficiency (X-SCID)
    # -------------------------------------------------------------------
    "IL2RG": {
        "gene_id": 3561,
        "chrom": "chrX",
        "start": 71_107_404,
        "end": 71_111_577,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq13.1",
        "exon_count": 8,
        "role": (
            "Interleukin-2 receptor subunit gamma (common gamma chain). "
            "Shared signalling subunit for IL-2, IL-4, IL-7, IL-9, IL-15, "
            "and IL-21 receptors.  Loss-of-function causes X-SCID (T-B+NK-)."
        ),
        "disease": "X-linked Severe Combined Immunodeficiency (X-SCID)",
        "omim": "300400",
        "key_variants": [
            {
                "name": "c.676C>T (p.Arg226Cys)",
                "hgvs_genomic": "NC_000023.11:g.71109765G>A",  # complement
                "rsid": "rs104895492",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Extracellular domain; variable expressivity from "
                    "classic SCID to atypical/leaky presentations."
                ),
            },
            {
                "name": "c.677G>A (p.Arg226His)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Adjacent to R226C; also extracellular domain.",
            },
            {
                "name": "c.670C>T (p.Arg224Trp)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Extracellular domain hotspot region.",
            },
        ],
        "mutational_landscape": (
            ">140 pathogenic/likely-pathogenic variants in ClinVar. "
            "~50% cluster in exon 5 (WSXWS motif, extracellular domain) "
            "and exon 3.  Mutation types: missense, nonsense, splice-site, "
            "small insertions/deletions.  Null alleles cause classic SCID; "
            "hypomorphic missense alleles cause atypical/leaky SCID."
        ),
        "crispr_strategy": (
            "Ex vivo adenine base editing (ABE) of autologous CD34+ HSPCs "
            "to correct point mutations, followed by re-infusion after "
            "reduced-intensity conditioning.  Correction of the common "
            "gamma chain restores signalling for all six cytokine receptors."
        ),
        "clinical_trials": [
            {
                "nct": "NCT06851767",
                "sponsor": "NIAID / NIH",
                "phase": "I/II",
                "approach": "Adenine base editing of autologous HSPCs (BE-HSPC IL2RG)",
                "status_2026": (
                    "Phase I/II open-label trial; 18 planned participants "
                    "(age >= 3). First patient dosed June 2025.  As of late "
                    "2025, 3 X-SCID participants + 1 CD40L-deficient "
                    "Hyper-IgM participant dosed; no severe adverse events "
                    "reported.  This is the first clinical demonstration of "
                    "CRISPR-mediated correction (not disruption) of a "
                    "disease-causing mutation."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1b.  ADA-SCID
    # -------------------------------------------------------------------
    "ADA": {
        "gene_id": 100,
        "chrom": "chr20",
        "start": 44_619_522,
        "end": 44_651_699,
        "strand": "-",
        "refseq": "NC_000020.11",
        "cytoband": "20q13.12",
        "exon_count": 12,
        "role": (
            "Adenosine deaminase -- catalyses irreversible deamination of "
            "adenosine and deoxyadenosine in the purine salvage pathway.  "
            "Deficiency causes toxic accumulation of deoxyadenosine and "
            "dATP leading to lymphocyte apoptosis (T-B-NK- SCID)."
        ),
        "disease": "ADA-Severe Combined Immunodeficiency (ADA-SCID)",
        "omim": "608958",
        "key_variants": [
            {
                "name": "c.646G>A (p.Gly216Arg / G216R)",
                "rsid": "rs121908723",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Identified in Amish kindred; exon 7; found in ~2/21 "
                    "additional ADA-SCID patients in one study."
                ),
            },
            {
                "name": "c.302T>C (p.Leu101Pro)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Loss of catalytic activity; severe phenotype.",
            },
            {
                "name": "c.467C>T (p.Ala156Val)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Classic ADA-SCID; near catalytic site.",
            },
        ],
        "mutational_landscape": (
            "~100 pathogenic variants in ClinVar.  Mostly missense affecting "
            "catalytic residues.  Also splice-site, nonsense, and small indels.  "
            "Genotype-phenotype correlation: null alleles = early-onset SCID; "
            "partial activity = delayed/partial ADA deficiency."
        ),
        "crispr_strategy": (
            "Ex vivo gene correction or gene addition in autologous CD34+ HSPCs.  "
            "Historically treated with PEG-ADA enzyme replacement as a bridge, "
            "then gamma-retroviral (Strimvelis) or lentiviral (OTL-101) gene "
            "addition.  CRISPR HDR-based correction or targeted cDNA knock-in "
            "at safe-harbour locus being explored preclinically."
        ),
        "gene_therapy_history": [
            {
                "product": "Strimvelis",
                "approach": "Ex vivo gamma-retroviral vector gene addition in CD34+ HSPCs",
                "status": (
                    "EMA-approved 2016; first approved ex vivo gene therapy.  "
                    "Orchard Therapeutics discontinued in 2022; marketing "
                    "authorisation transferred to Fondazione Telethon (Italy).  "
                    "One case of treatment-related leukaemia reported "
                    "(insertional oncogenesis risk with gamma-retroviral vectors)."
                ),
            },
            {
                "product": "OTL-101",
                "approach": "Ex vivo lentiviral vector gene addition (safer self-inactivating vector)",
                "status": (
                    "Orchard Therapeutics; pivoted away from ADA-SCID to "
                    "focus on neurometabolic diseases.  Lentiviral approach "
                    "considered safer than gamma-retroviral (lower insertional "
                    "oncogenesis risk)."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1c.  RAG1-SCID / RAG2-SCID
    # -------------------------------------------------------------------
    "RAG1": {
        "gene_id": 5896,
        "chrom": "chr11",
        "start": 36_510_353,
        "end": 36_579_762,
        "strand": "+",
        "refseq": "NC_000011.10",
        "cytoband": "11p12",
        "exon_count": 7,
        "role": (
            "Recombination activating gene 1 -- essential for V(D)J "
            "recombination of immunoglobulin and T-cell receptor genes.  "
            "Null mutations cause T-B-NK+ SCID; hypomorphic mutations "
            "cause Omenn syndrome or leaky SCID."
        ),
        "disease": "RAG1-SCID / Omenn Syndrome",
        "omim": "179615",
        "key_variants": [
            {
                "name": "c.2333G>A (p.Arg778Gln / R778Q)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Core RAG1 domain; causes classic T-B- SCID.",
            },
            {
                "name": "c.1681C>T (p.Arg561His / R561H)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Hypomorphic; associated with Omenn syndrome. "
                    "Residual V(D)J activity produces oligoclonal T cells."
                ),
            },
            {
                "name": "c.1187G>A (p.Arg396His / R396H)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Core RAG1 catalytic domain.",
            },
        ],
        "mutational_landscape": (
            "345 variants in ClinVar; ~92 pathogenic/likely-pathogenic, "
            "majority are missense in the core catalytic domain (residues "
            "384-1008).  Over 200 distinct RAG1 pathogenic variants reported.  "
            "Entire coding sequence is in a single large exon."
        ),
        "crispr_strategy": (
            "Full coding-sequence replacement via CRISPR-Cas9-mediated HDR "
            "at the endogenous locus (demonstrated in Nature Communications "
            "2023).  This 'universal' approach corrects all RAG1 mutations "
            "with a single construct.  Ex vivo in autologous HSPCs."
        ),
    },

    "RAG2": {
        "gene_id": 5897,
        "chrom": "chr11",
        "start": 36_590_996,
        "end": 36_598_236,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11p12",
        "exon_count": 2,
        "role": (
            "Recombination activating gene 2 -- forms a complex with RAG1 "
            "for V(D)J recombination.  Located ~11 kb from RAG1, convergently "
            "transcribed.  Entire protein encoded in a single exon."
        ),
        "disease": "RAG2-SCID / Omenn Syndrome",
        "omim": "179616",
        "key_variants": [
            {
                "name": "c.686G>A (p.Arg229Gln / R229Q)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Compound heterozygous with other RAG2 variants; "
                    "associated with T-B- SCID."
                ),
            },
            {
                "name": "c.1G>A (p.Met1?)",
                "rsid": None,
                "consequence": "start_loss",
                "clinical_significance": "pathogenic",
                "notes": "Loss of initiation codon; null allele.",
            },
        ],
        "mutational_landscape": (
            "Fewer known variants than RAG1.  Most pathogenic variants "
            "are missense in the Kelch-repeat / PHD domain regions.  "
            "Similar genotype-phenotype spectrum: null = SCID, "
            "hypomorphic = Omenn syndrome or granuloma/autoimmunity."
        ),
        "crispr_strategy": (
            "Complete coding-sequence replacement via CRISPR-Cas9-mediated "
            "HDR (same approach as RAG1).  Single exon simplifies construct "
            "design.  Ex vivo correction in autologous HSPCs."
        ),
    },

    # -------------------------------------------------------------------
    # 1d.  Chronic Granulomatous Disease (CGD)
    # -------------------------------------------------------------------
    "CYBB": {
        "gene_id": 1536,
        "chrom": "chrX",
        "start": 37_780_059,
        "end": 37_813_461,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xp21.1-p11.4",
        "exon_count": 14,
        "role": (
            "Cytochrome b-245 beta chain (gp91-phox) -- membrane-bound "
            "catalytic subunit of the phagocyte NADPH oxidase.  Hemizygous "
            "loss-of-function causes X-linked CGD (~65% of all CGD cases)."
        ),
        "disease": "X-linked Chronic Granulomatous Disease",
        "omim": "300481",
        "key_variants": [
            {
                "name": "c.676C>T (p.Arg226Ter / R226X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Premature stop; WASP-negative; NMD or truncated protein.  "
                    "Identified in multiple CGD families."
                ),
            },
            {
                "name": "c.252G>A (splice / p.?)",
                "rsid": None,
                "consequence": "splice_donor / missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Target of CRISPR D10A-nickase correction with no "
                    "detectable off-target activity (research demonstration)."
                ),
            },
            {
                "name": "c.376C>T (p.Arg126Cys / C126R de novo)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "De novo; recurrent pneumonia and BCGitis.",
            },
        ],
        "mutational_landscape": (
            "Hundreds of pathogenic variants across all 14 exons.  Mutation "
            "types include missense, nonsense, splice-site, small indels, "
            "and large deletions/CNVs.  No single common founder mutation; "
            "highly heterogeneous."
        ),
        "crispr_strategy": (
            "Base editing: NIAID Phase I/II trial of base-edited autologous "
            "HSPCs to repair CYBB missense mutations.  One participant dosed; "
            "a severe adverse event was reported.  Also: targeted cDNA insertion "
            "near CYBB exon 1 as a 'near-universal' approach correcting "
            "all downstream mutations."
        ),
        "clinical_trials": [
            {
                "nct": "NCT06325709",
                "sponsor": "NIAID / NIH",
                "phase": "I/II",
                "approach": "Adenine base editing of autologous HSPCs for CYBB repair",
                "status_2026": (
                    "Open-label; 1 participant dosed; a severe adverse event "
                    "was reported (details pending)."
                ),
            },
        ],
    },

    "CYBA": {
        "gene_id": 1535,
        "chrom": "chr16",
        "start": 88_643_289,
        "end": 88_651_053,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16q24.2",
        "exon_count": 6,
        "role": (
            "Cytochrome b-245 alpha chain (p22-phox) -- small membrane "
            "subunit of phagocyte NADPH oxidase.  Biallelic loss-of-function "
            "causes autosomal recessive CGD (~5% of cases)."
        ),
        "disease": "Autosomal Recessive CGD (p22-phox deficiency)",
        "omim": "608508",
        "crispr_strategy": (
            "Targeted cDNA insertion near CYBA exon 1 demonstrated in "
            "preclinical studies (Nature Communications 2025).  Similar "
            "'near-universal' approach as for CYBB."
        ),
    },

    "NCF1": {
        "gene_id": 653361,
        "chrom": "chr7",
        "start": 74_774_011,
        "end": 74_789_315,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q11.23",
        "exon_count": 11,
        "role": (
            "Neutrophil cytosolic factor 1 (p47-phox) -- cytoplasmic "
            "regulatory subunit of NADPH oxidase.  Biallelic loss-of-function "
            "causes AR-CGD (~20% of cases).  NCF1 locus complicated by "
            "pseudogenes NCF1B and NCF1C."
        ),
        "disease": "Autosomal Recessive CGD (p47-phox deficiency)",
        "omim": "608512",
        "crispr_strategy": (
            "Prime editing: Prime Medicine PM359 trial.  First-ever prime-editing "
            "clinical data reported May 2025 in a CGD patient.  December 2025 "
            "report: 2 participants treated, robust restoration of NADPH "
            "oxidase function, resolution of infections and CGD colitis."
        ),
        "clinical_trials": [
            {
                "nct": "NCT06559176",
                "sponsor": "Prime Medicine",
                "product": "PM359",
                "phase": "I/II",
                "approach": "Prime editing in autologous HSPCs for NCF1 correction",
                "status_2026": (
                    "First-ever prime-editing human efficacy data (May 2025).  "
                    "2 patients treated; robust restoration of neutrophil "
                    "oxidase function; resolution of skin, lung, and GI "
                    "infections and CGD colitis.  Prime Medicine seeking "
                    "FDA approval pathway as of March 2026."
                ),
            },
        ],
    },

    "NCF2": {
        "gene_id": 4688,
        "chrom": "chr1",
        "start": 183_555_562,
        "end": 183_601_849,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1q25.3",
        "exon_count": 20,
        "role": (
            "Neutrophil cytosolic factor 2 (p67-phox) -- cytoplasmic "
            "activating subunit of NADPH oxidase.  Biallelic loss-of-function "
            "causes AR-CGD (~5% of cases)."
        ),
        "disease": "Autosomal Recessive CGD (p67-phox deficiency)",
        "omim": "608515",
    },

    "NCF4": {
        "gene_id": 4689,
        "chrom": "chr22",
        "start": 36_861_006,
        "end": 36_878_015,
        "strand": "+",
        "refseq": "NC_000022.11",
        "cytoband": "22q12.3",
        "exon_count": 12,
        "role": (
            "Neutrophil cytosolic factor 4 (p40-phox) -- regulatory "
            "component of NADPH oxidase.  Biallelic mutations cause "
            "AR-CGD with selective defects in neutrophil oxidase activity "
            "(~1% of cases)."
        ),
        "disease": "Autosomal Recessive CGD (p40-phox deficiency)",
        "omim": "601488",
    },

    # -------------------------------------------------------------------
    # 1e.  Wiskott-Aldrich Syndrome (WAS)
    # -------------------------------------------------------------------
    "WAS": {
        "gene_id": 7454,
        "chrom": "chrX",
        "start": 48_676_636,
        "end": 48_691_427,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xp11.23",
        "exon_count": 13,
        "role": (
            "Wiskott-Aldrich syndrome protein (WASp) -- actin cytoskeleton "
            "regulator in haematopoietic cells.  Loss-of-function causes "
            "WAS triad: thrombocytopenia, eczema, immunodeficiency.  "
            "Milder mutations cause X-linked thrombocytopenia (XLT)."
        ),
        "disease": "Wiskott-Aldrich Syndrome / X-linked Thrombocytopenia",
        "omim": "301000",
        "key_variants": [
            {
                "name": "c.665C>T (p.Arg211Ter / R211X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Hotspot mutation; found in 10 families in large series.  "
                    "WASp-negative; high disease severity score.  Classic WAS."
                ),
            },
            {
                "name": "c.291G>A (p.Arg86His or nearby missense at codon 86)",
                "rsid": None,
                "consequence": "missense / splice",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Hotspot at codon 86 (R86H, R86C, R86S, R86G, R86L).  "
                    "17 unrelated families.  WASp-positive; mild phenotype "
                    "(XLT); low disease severity score."
                ),
            },
            {
                "name": "c.290C>T (p.Arg86Cys / R86C)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common at codon-86 hotspot; mild XLT phenotype.",
            },
        ],
        "mutational_landscape": (
            "Study of 577 patients identified six hotspot mutations.  "
            "Genotype is a biomarker for disease severity and survival.  "
            "Null (WASp-negative) = classic WAS; missense/reduced WASp = XLT."
        ),
        "crispr_strategy": (
            "Ex vivo lentiviral gene addition already approved (see below).  "
            "CRISPR-Cas9 targeted gene correction in HSPCs demonstrated "
            "preclinically (Rai et al., 2020, Nature Communications)."
        ),
        "gene_therapy_history": [
            {
                "product": "Waskyra (etuvetidigene autotemcel)",
                "approach": "Ex vivo lentiviral vector gene addition in CD34+ HSPCs",
                "status": (
                    "FDA-approved December 2025 for WAS patients >= 6 months "
                    "without suitable HLA-matched related donor.  Single "
                    "administration.  Long-term data (median 7.6 yr follow-up): "
                    "severe infections and eczema resolved; autoimmune episodes "
                    "and bleeding significantly reduced.  No treatment-related "
                    "serious adverse events or insertional oncogenesis reported."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1f.  Hyper-IgM Syndrome
    # -------------------------------------------------------------------
    "CD40LG": {
        "gene_id": 959,
        "chrom": "chrX",
        "start": 136_648_158,
        "end": 136_660_390,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xq26.3",
        "exon_count": 5,
        "role": (
            "CD40 ligand (CD154 / TNFSF5) -- expressed on activated T cells; "
            "engages CD40 on B cells to drive immunoglobulin class-switch "
            "recombination.  Hemizygous loss-of-function causes X-linked "
            "Hyper-IgM syndrome (HIGM1): low/absent IgG, IgA, IgE with "
            "normal/elevated IgM."
        ),
        "disease": "X-linked Hyper-IgM Syndrome (HIGM1)",
        "omim": "300386",
        "key_variants": [
            {
                "name": "c.607C>T (p.Arg203Ile / R203I)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Atypical HIGM1 phenotype; TNF-homology domain.",
            },
        ],
        "mutational_landscape": (
            "Most mutations cluster in the TNF-homology domain (exon 5).  "
            "Types: missense, nonsense, splice-site, small indels.  "
            "Highly heterogeneous; few recurrent mutations."
        ),
        "crispr_strategy": (
            "Base editing in the NIAID X-SCID trial has also enrolled a "
            "CD40L-deficient Hyper-IgM patient (same NCT06851767 trial "
            "with expanded eligibility).  Ex vivo ABE correction of "
            "CD40LG point mutations in autologous HSPCs."
        ),
    },

    "AICDA": {
        "gene_id": 57379,
        "chrom": "chr12",
        "start": 8_602_170,
        "end": 8_612_859,
        "strand": "-",
        "refseq": "NC_000012.12",
        "cytoband": "12p13.31",
        "exon_count": 5,
        "role": (
            "Activation-induced cytidine deaminase (AID) -- essential for "
            "somatic hypermutation and class-switch recombination in B cells.  "
            "Biallelic loss-of-function causes autosomal recessive Hyper-IgM "
            "syndrome type 2 (HIGM2)."
        ),
        "disease": "Autosomal Recessive Hyper-IgM Syndrome Type 2 (HIGM2)",
        "omim": "605257",
        "mutational_landscape": (
            "Rare; mutations throughout the gene including catalytic domain.  "
            "Missense, nonsense, and splice-site variants reported."
        ),
    },
}

# ============================================================================
# 2. CYSTIC FIBROSIS
# ============================================================================

CYSTIC_FIBROSIS_TARGETS = {
    "CFTR": {
        "gene_id": 1080,
        "chrom": "chr7",
        "start": 117_480_025,
        "end": 117_668_665,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q31.2",
        "exon_count": 27,
        "gene_size_bp": 188_640,
        "role": (
            "Cystic fibrosis transmembrane conductance regulator -- cAMP-gated "
            "chloride/bicarbonate channel on epithelial cell surfaces.  "
            "Loss-of-function causes cystic fibrosis: progressive lung disease, "
            "pancreatic insufficiency, elevated sweat chloride."
        ),
        "disease": "Cystic Fibrosis",
        "omim": "602421",

        # ---------------------------------------------------------------
        # Key pathogenic variants with exact GRCh38 positions
        # ---------------------------------------------------------------
        "key_variants": [
            # ---- F508del (Class II) -- most common, ~70% of CF alleles ----
            {
                "name": "F508del / p.Phe508del",
                "hgvs_coding": "NM_000492.4:c.1521_1523delCTT",
                "hgvs_genomic": "NC_000007.14:g.117559592_117559594del",
                "grch38_position": "chr7:117,559,591-117,559,594",
                "rsid": "rs113993960",
                "ref_allele": "TCTT",
                "alt_allele": "T",  # deletion of CTT
                "consequence": "inframe_deletion",
                "clinical_significance": "pathogenic",
                "mutation_class": "II",
                "frequency": "~70% of CF alleles worldwide",
                "notes": (
                    "Most common CF mutation globally.  Protein misfolding / "
                    "trafficking defect -- CFTR does not reach cell surface.  "
                    "Targeted by lumacaftor+ivacaftor (Orkambi), "
                    "tezacaftor+ivacaftor (Symdeko), and "
                    "elexacaftor+tezacaftor+ivacaftor (Trikafta)."
                ),
            },
            # ---- G551D (Class III) -- gating mutation, Kalydeco target ----
            {
                "name": "G551D / p.Gly551Asp",
                "hgvs_coding": "NM_000492.4:c.1652G>A",
                "hgvs_genomic": "NC_000007.14:g.117587806G>A",
                "grch38_position": "chr7:117,587,806",
                "rsid": "rs75527207",
                "ref_allele": "G",
                "alt_allele": "A",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "mutation_class": "III",
                "frequency": "~1.5% of CF alleles",
                "notes": (
                    "Gating mutation; CFTR reaches surface but channel does "
                    "not open properly.  Ivacaftor (Kalydeco) potentiates "
                    "channel gating -- FDA-approved January 2012 specifically "
                    "for G551D.  Also drug-response annotation in ClinVar."
                ),
            },
            # ---- G542X (Class I) -- nonsense / premature stop ----
            {
                "name": "G542X / p.Gly542Ter",
                "hgvs_coding": "NM_000492.4:c.1624G>T",
                "hgvs_genomic": "NC_000007.14:g.117587778G>T",
                "grch38_position": "chr7:117,587,778",
                "rsid": "rs113993959",
                "ref_allele": "G",
                "alt_allele": "T",
                "consequence": "stop_gained",
                "clinical_significance": "pathogenic",
                "mutation_class": "I",
                "frequency": "~2.6% of CF alleles",
                "notes": (
                    "Premature stop codon in exon 12.  No functional CFTR "
                    "protein produced.  Candidate for adenine base editing "
                    "(ABE) to correct the premature stop.  Preclinical "
                    "demonstration: ABE rescued G542X in patient-derived "
                    "intestinal organoids (iScience, 2025)."
                ),
            },
            # ---- W1282X (Class I) -- nonsense ----
            {
                "name": "W1282X / p.Trp1282Ter",
                "hgvs_coding": "NM_000492.4:c.3846G>A",
                "hgvs_genomic": "NC_000007.14:g.117642566G>A",
                "grch38_position": "chr7:117,642,566",
                "rsid": "rs77010898",
                "ref_allele": "G",
                "alt_allele": "A",
                "consequence": "stop_gained",
                "clinical_significance": "pathogenic",
                "mutation_class": "I",
                "frequency": "~1.0% of CF alleles (higher in Ashkenazi Jewish)",
                "notes": (
                    "Premature stop in exon 22.  Target of base-editing "
                    "approaches (ABE to revert stop codon).  One of the most "
                    "studied targets for CRISPR-based CF research."
                ),
            },
            # ---- N1303K (Class II) -- missense ----
            {
                "name": "N1303K / p.Asn1303Lys",
                "hgvs_coding": "NM_000492.4:c.3909C>G",
                "hgvs_genomic": "NC_000007.14:g.117652877C>G",
                "grch38_position": "chr7:117,652,877",
                "rsid": "rs80034486",
                "ref_allele": "C",
                "alt_allele": "G",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "mutation_class": "II",
                "frequency": "~1.6% of CF alleles",
                "notes": (
                    "Missense in nucleotide-binding domain 2 (NBD2).  "
                    "Processing/trafficking defect similar to F508del."
                ),
            },
            # ---- R553X (Class I) -- nonsense ----
            {
                "name": "R553X / p.Arg553Ter",
                "hgvs_coding": "NM_000492.4:c.1657C>T",
                "hgvs_genomic": "NC_000007.14:g.117587811C>T",
                "grch38_position": "chr7:117,587,811",
                "rsid": "rs74597325",
                "ref_allele": "C",
                "alt_allele": "T",
                "consequence": "stop_gained",
                "clinical_significance": "pathogenic",
                "mutation_class": "I",
                "frequency": "~0.7% of CF alleles",
                "notes": (
                    "Premature stop in exon 12; also causes aberrant splicing "
                    "(exon 11 skipping) leading to frameshift.  >96% of "
                    "patients are pancreatic-insufficient."
                ),
            },
        ],

        # ---------------------------------------------------------------
        # CFTR mutation classification system
        # ---------------------------------------------------------------
        "mutation_classes": {
            "I": "No functional CFTR protein (nonsense, frameshift, splice-site). "
                 "Examples: G542X, W1282X, R553X, 621+1G>T.",
            "II": "CFTR protein misfolded / trafficking defect; degraded before "
                  "reaching cell surface.  Examples: F508del, N1303K.",
            "III": "CFTR reaches surface but channel gating defective.  "
                   "Examples: G551D, S549N.",
            "IV": "CFTR reaches surface and opens, but reduced chloride "
                  "conductance.  Examples: R117H, R334W.",
            "V": "Reduced quantity of normal CFTR (splicing defects, "
                 "promoter variants).  Examples: 3849+10kbC>T, A455E.",
            "VI": "Normal CFTR produced but accelerated turnover / instability "
                  "at surface.  Examples: Q1412X, 4326delTC.",
        },

        # ---------------------------------------------------------------
        # CRISPR therapeutic strategies for CF
        # ---------------------------------------------------------------
        "crispr_strategies": {
            "base_editing_nonsense": (
                "Adenine base editing (ABE) to correct Class I premature stop "
                "codons (e.g., G542X, W1282X, R553X) by converting the "
                "pathogenic A back to G in the antisense strand.  Demonstrated "
                "in patient-derived organoids and airway epithelia."
            ),
            "prime_editing_F508del": (
                "Prime editing to precisely re-insert the deleted CTT trinucleotide "
                "for F508del.  Broad Institute demonstrated efficient correction "
                "in human lung cells.  Requires pegRNA design targeting the "
                "3-bp deletion.  More versatile than base editing for indels."
            ),
            "crispr_hdr": (
                "Traditional CRISPR-Cas9 + HDR template for precise correction.  "
                "Less efficient in post-mitotic airway cells; more applicable "
                "ex vivo in basal stem cells or iPSC-derived models."
            ),
            "delivery_challenges": (
                "Major barrier: in vivo delivery to airway epithelium.  "
                "LNP delivery to primate airway achieved ~5.3% editing "
                "efficiency in rhesus epithelia.  In mice, edited epithelia "
                "persisted for 12+ months.  Clinical trials for CRISPR-based "
                "CF therapy are NOT yet underway as of early 2026, but "
                "preclinical data are advancing rapidly."
            ),
        },
        "clinical_trial_status_2026": (
            "No CRISPR/gene-editing clinical trials for CF as of early 2026.  "
            "Many challenges remain: efficient in vivo delivery to airways, "
            "editing efficiency in basal stem cells, immune response to "
            "delivery vehicles.  Systematic review (2025) analyzed 27 "
            "preclinical studies; F508del and W1282X are the most studied.  "
            "Current clinical treatment relies on CFTR modulators (Trikafta) "
            "and investigational mRNA/gene transfer approaches."
        ),
    },
}

# ============================================================================
# 3. METABOLIC / LIVER DISEASES
# ============================================================================

METABOLIC_LIVER_TARGETS = {

    # -------------------------------------------------------------------
    # 3a.  Phenylketonuria (PKU)
    # -------------------------------------------------------------------
    "PAH": {
        "gene_id": 5053,
        "chrom": "chr12",
        "start": 102_836_889,
        "end": 102_958_441,
        "strand": "-",
        "refseq": "NC_000012.12",
        "cytoband": "12q23.2",
        "exon_count": 15,
        "role": (
            "Phenylalanine hydroxylase -- biopterin-dependent enzyme catalysing "
            "conversion of phenylalanine to tyrosine.  Deficiency causes PKU: "
            "toxic phenylalanine accumulation leading to intellectual disability "
            "if untreated."
        ),
        "disease": "Phenylketonuria (PKU)",
        "omim": "612349",
        "key_variants": [
            {
                "name": "R408W / p.Arg408Trp",
                "hgvs_coding": "NM_000277.3:c.1222C>T",
                "hgvs_genomic": "NC_000012.12:g.102840493G>A",
                "grch38_position": "chr12:102,840,493",
                "rsid": "rs5030858",
                "ref_allele": "G",
                "alt_allele": "A",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "frequency": "gnomAD 0.2%; most prevalent PKU allele in Europe",
                "notes": (
                    "Exon 12; CGG>TGG transition.  Associated with classic "
                    "severe PKU.  Highest frequency in Eastern Europe "
                    "(Bulgaria, Lithuania, eastern Germany); Balto-Slavic "
                    "origin on haplotype 2.  Homozygous R408W found in "
                    "~4.8% of European PKU patients."
                ),
            },
            {
                "name": "R252W / p.Arg252Trp",
                "hgvs_coding": "NM_000277.3:c.754C>T",
                "rsid": "rs5030849",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Mediterranean kindreds; CpG deamination.  Negligible "
                    "enzyme activity; undetectable PAH protein."
                ),
            },
            {
                "name": "R261Q / p.Arg261Gln",
                "hgvs_coding": "NM_000277.3:c.782G>A",
                "rsid": "rs5030850",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "BH4-responsive variant; mild/moderate PKU.  One of "
                    "the most common BH4-responsive mutations."
                ),
            },
            {
                "name": "IVS12+1G>A (c.1315+1G>A)",
                "hgvs_coding": "NM_000277.3:c.1315+1G>A",
                "rsid": "rs5030855",
                "consequence": "splice_donor",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Splice-donor variant; most frequent allele in Norway "
                    "(19% of PKU chromosomes).  Also common in Denmark/England."
                ),
            },
            {
                "name": "Y414C / p.Tyr414Cys",
                "hgvs_coding": "NM_000277.3:c.1241A>G",
                "rsid": "rs5030856",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Mild PKU; common in Scandinavia.  BH4 acts as "
                    "chemical chaperone preventing protein degradation."
                ),
            },
        ],
        "crispr_strategy": (
            "Preclinical: in vivo liver-directed base editing or prime editing "
            "to correct common PAH missense mutations.  PAH is expressed in "
            "liver hepatocytes, making LNP delivery feasible.  No clinical "
            "trials for gene editing as of early 2026."
        ),
    },

    # -------------------------------------------------------------------
    # 3b.  Wilson Disease
    # -------------------------------------------------------------------
    "ATP7B": {
        "gene_id": 540,
        "chrom": "chr13",
        "start": 51_932_669,
        "end": 52_012_132,
        "strand": "-",
        "refseq": "NC_000013.11",
        "cytoband": "13q14.3",
        "exon_count": 26,
        "role": (
            "ATPase copper transporting beta -- P-type ATPase that transports "
            "copper into the trans-Golgi network for incorporation into "
            "ceruloplasmin, and into bile for excretion.  Biallelic "
            "loss-of-function causes Wilson disease: hepatic and neurologic "
            "copper accumulation."
        ),
        "disease": "Wilson Disease",
        "omim": "606882",
        "key_variants": [
            {
                "name": "H1069Q / p.His1069Gln",
                "hgvs_coding": "NM_000053.4:c.3207C>A",
                "rsid": "rs76151636",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common Wilson disease mutation in European "
                    "populations (~35-45% of alleles).  ATP-binding domain; "
                    "protein misfolding / trafficking defect."
                ),
            },
            {
                "name": "R778L / p.Arg778Leu",
                "hgvs_coding": "NM_000053.4:c.2333G>T",
                "rsid": "rs28942074",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mutation in East Asian populations.  "
                    "Transmembrane domain; loss of copper transport."
                ),
            },
        ],
        "crispr_strategy": (
            "Preclinical: liver-directed in vivo base editing or dual-AAV "
            "gene replacement.  Liver-expressed gene suitable for LNP or "
            "AAV delivery.  No clinical trials for gene editing as of "
            "early 2026."
        ),
    },

    # -------------------------------------------------------------------
    # 3c.  Alpha-1 Antitrypsin Deficiency (AATD)
    # -------------------------------------------------------------------
    "SERPINA1": {
        "gene_id": 5265,
        "chrom": "chr14",
        "start": 94_376_747,
        "end": 94_390_635,
        "strand": "-",
        "refseq": "NC_000014.9",
        "cytoband": "14q32.13",
        "exon_count": 7,
        "role": (
            "Serpin family A member 1 (alpha-1 antitrypsin, AAT) -- serine "
            "protease inhibitor produced by hepatocytes.  Inhibits neutrophil "
            "elastase in lungs.  Deficiency causes emphysema / COPD and "
            "liver disease (toxic polymer accumulation in hepatocytes)."
        ),
        "disease": "Alpha-1 Antitrypsin Deficiency (AATD)",
        "omim": "613490",
        "key_variants": [
            {
                "name": "Z allele / E342K (p.Glu366Lys in current numbering)",
                "hgvs_coding": "NM_000295.5:c.1096G>A",
                "hgvs_genomic": "NC_000014.9:g.94378610C>T",
                "grch38_position": "chr14:94,378,610",
                "rsid": "rs28929474",
                "ref_allele": "C",
                "alt_allele": "T",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most clinically significant AATD allele.  Accounts for "
                    "~95% of clinically recognised AATD.  ZZ homozygotes have "
                    "~15% of normal AAT levels.  Protein misfolds and "
                    "polymerises in hepatocyte ER causing liver disease, while "
                    "lung disease results from uninhibited elastase activity.  "
                    "5x less effective as neutrophil elastase inhibitor than "
                    "normal M allele."
                ),
            },
            {
                "name": "S allele / E264V (p.Glu288Val in current numbering)",
                "hgvs_coding": "NM_000295.5:c.863A>T",
                "hgvs_genomic": "NC_000014.9:g.94380925T>A",
                "grch38_position": "chr14:94,380,925",
                "rsid": "rs17580",
                "ref_allele": "T",
                "alt_allele": "A",
                "consequence": "missense",
                "clinical_significance": "pathogenic / risk_factor",
                "notes": (
                    "Milder deficiency.  SS homozygotes have ~40-50% of "
                    "normal AAT levels.  SZ compound heterozygotes have "
                    "~35% of normal levels and intermediate clinical risk."
                ),
            },
        ],
        "crispr_strategy": (
            "Liver-directed in vivo base editing to correct the Z allele "
            "(G>A at coding level, C>T at genomic level on minus strand).  "
            "ABE could directly revert the Z mutation.  Also: CRISPR knockout "
            "of mutant allele to prevent toxic polymer accumulation in "
            "liver (reduces liver disease risk while AAT augmentation "
            "therapy covers lung protection).  No clinical trials as of "
            "early 2026."
        ),
    },

    # -------------------------------------------------------------------
    # 3d.  Familial Hypercholesterolemia (FH) -- multi-gene
    # -------------------------------------------------------------------
    "LDLR": {
        "gene_id": 3949,
        "chrom": "chr19",
        "start": 11_089_463,
        "end": 11_133_820,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19p13.2",
        "exon_count": 19,
        "role": (
            "Low-density lipoprotein receptor -- mediates hepatic LDL "
            "cholesterol uptake.  Heterozygous loss-of-function causes "
            "heterozygous FH (~1:250 prevalence); homozygous/compound "
            "heterozygous causes severe homozygous FH (~1:300,000)."
        ),
        "disease": "Familial Hypercholesterolemia (FH)",
        "omim": "606945",
        "mutational_landscape": (
            ">2,000 known LDLR pathogenic variants.  Highly heterogeneous: "
            "missense, nonsense, splice-site, large deletions/duplications.  "
            "No single common mutation; population-specific founder effects."
        ),
        "crispr_strategy": (
            "Indirect approach: rather than correcting the vast array of "
            "LDLR mutations, CRISPR knockout of PCSK9 or ANGPTL3 (see below) "
            "to upregulate remaining LDLR activity or reduce LDL via "
            "alternative pathways."
        ),
    },

    "PCSK9": {
        "gene_id": 255738,
        "chrom": "chr1",
        "start": 55_039_548,
        "end": 55_064_852,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1p32.3",
        "exon_count": 15,
        "role": (
            "Proprotein convertase subtilisin/kexin type 9 -- secreted by "
            "hepatocytes; binds LDLR and promotes its lysosomal degradation.  "
            "PCSK9 gain-of-function = hypercholesterolemia; loss-of-function "
            "= low LDL and cardioprotection.  CRISPR knockout target."
        ),
        "disease": "Target for Familial Hypercholesterolemia therapy",
        "omim": "607786",
        "crispr_strategy": (
            "In vivo hepatic PCSK9 knockout via base editing (VERVE-101/102) "
            "or CRISPR-Cas9 nuclease.  Permanently silences PCSK9, mimicking "
            "the cardioprotective natural loss-of-function variants.  "
            "Single-dose, 'one-and-done' LDL reduction."
        ),
        "clinical_trials": [
            {
                "product": "VERVE-101",
                "nct": "Heart-1",
                "sponsor": "Verve Therapeutics",
                "phase": "I",
                "approach": (
                    "LNP-delivered adenine base editing (ABE) targeting "
                    "a single base in PCSK9 to permanently inactivate "
                    "hepatic PCSK9 production."
                ),
                "status_2026": (
                    "Interim data showed dose-dependent reductions: at "
                    "0.45 mg/kg, PCSK9 reduced 59-84%, LDL reduced 39-48%.  "
                    "At 0.6 mg/kg, PCSK9 reduced 47%, LDL reduced 55% "
                    "sustained >180 days.  Trial paused April 2024 after "
                    "one participant developed elevated liver enzymes and "
                    "low platelets."
                ),
            },
            {
                "product": "VERVE-102",
                "nct": "Heart-2",
                "sponsor": "Verve Therapeutics",
                "phase": "I (transitioning to Phase II)",
                "approach": (
                    "Same ABE payload as VERVE-101 but with improved LNP "
                    "delivery system for better tolerability."
                ),
                "status_2026": (
                    "April 2025 interim data: well-tolerated, no treatment-"
                    "related SAEs.  14 patients with HeFH dosed across 3 "
                    "cohorts (0.3, 0.45, 0.6 mg/kg).  At 0.6 mg/kg: mean "
                    "LDL reduction 53%, maximum 69%.  FDA IND cleared; "
                    "Phase II enrollment initiating in 2026.  Final "
                    "dose-escalation data expected H2 2025."
                ),
            },
        ],
    },

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
            "Angiopoietin-like 3 -- secreted by liver; inhibits lipoprotein "
            "lipase and endothelial lipase.  Loss-of-function causes "
            "familial combined hypolipidaemia (low LDL, TG, and HDL) with "
            "no apparent adverse health effects.  CRISPR knockout target."
        ),
        "disease": "Target for dyslipidaemia / atherosclerosis therapy",
        "omim": "604774",
        "crispr_strategy": (
            "In vivo hepatic ANGPTL3 knockout via CRISPR-Cas9 nuclease "
            "(CTX310, CRISPR Therapeutics).  LNP-delivered Cas9 mRNA + "
            "sgRNA targeting ANGPTL3.  Single IV dose."
        ),
        "clinical_trials": [
            {
                "product": "CTX310",
                "sponsor": "CRISPR Therapeutics",
                "phase": "I",
                "approach": (
                    "LNP-encapsulated CRISPR-Cas9 mRNA + guide RNA "
                    "targeting ANGPTL3 in hepatocytes."
                ),
                "status_2026": (
                    "Phase I ascending-dose trial (0.1-0.8 mg/kg) published "
                    "NEJM November 2025.  At highest dose (0.8 mg/kg): mean "
                    "LDL reduction 48.9%, triglyceride reduction 55.2% "
                    "sustained through 60+ days.  Two serious adverse events "
                    "in 15 participants: one spinal disc herniation (unrelated), "
                    "one sudden death 179 days post-treatment at lowest dose "
                    "(0.1 mg/kg, relationship to treatment unclear)."
                ),
            },
        ],
    },
}


# ============================================================================
# HELPER: Flat gene lookup across all categories
# ============================================================================

ALL_GENE_THERAPY_TARGETS = {}
ALL_GENE_THERAPY_TARGETS.update(PRIMARY_IMMUNODEFICIENCY_TARGETS)
ALL_GENE_THERAPY_TARGETS.update(CYSTIC_FIBROSIS_TARGETS)
ALL_GENE_THERAPY_TARGETS.update(METABOLIC_LIVER_TARGETS)


def get_target(gene_symbol: str) -> dict | None:
    """Return the target dict for *gene_symbol*, or None."""
    return ALL_GENE_THERAPY_TARGETS.get(gene_symbol.upper())


def get_bed_regions() -> list[tuple[str, int, int, str]]:
    """Return a BED-format list of (chrom, start, end, name) for all targets."""
    regions = []
    for name, info in ALL_GENE_THERAPY_TARGETS.items():
        regions.append((info["chrom"], info["start"], info["end"], name))
    return sorted(regions, key=lambda r: (r[0], r[1]))


# Add conditions arrays for pipeline condition filtering
_CONDITIONS_MAP = {
    "IL2RG": ["SCID", "X-SCID", "immunodeficiency"],
    "ADA": ["SCID", "ADA-SCID", "immunodeficiency"],
    "RAG1": ["SCID", "RAG-SCID", "immunodeficiency"],
    "RAG2": ["SCID", "RAG-SCID", "immunodeficiency"],
    "CYBB": ["CGD", "chronic_granulomatous_disease", "immunodeficiency"],
    "CYBA": ["CGD", "chronic_granulomatous_disease", "immunodeficiency"],
    "NCF1": ["CGD", "chronic_granulomatous_disease", "immunodeficiency"],
    "NCF2": ["CGD", "chronic_granulomatous_disease", "immunodeficiency"],
    "NCF4": ["CGD", "chronic_granulomatous_disease", "immunodeficiency"],
    "WAS": ["wiskott_aldrich", "immunodeficiency"],
    "CD40LG": ["hyper_igm", "immunodeficiency"],
    "AICDA": ["hyper_igm", "immunodeficiency"],
    "CFTR": ["cystic_fibrosis", "pulmonary"],
    "PAH": ["phenylketonuria", "metabolic"],
    "ATP7B": ["wilson_disease", "metabolic", "liver"],
    "SERPINA1": ["alpha1_antitrypsin", "metabolic", "pulmonary", "liver"],
    "LDLR": ["familial_hypercholesterolemia", "cardiovascular", "metabolic"],
    "PCSK9": ["familial_hypercholesterolemia", "cardiovascular", "metabolic"],
    "ANGPTL3": ["dyslipidemia", "cardiovascular", "metabolic"],
}
for _gene, _conds in _CONDITIONS_MAP.items():
    if _gene in ALL_GENE_THERAPY_TARGETS:
        ALL_GENE_THERAPY_TARGETS[_gene]["conditions"] = _conds


def get_variant_vcf_lines() -> list[str]:
    """
    Return VCF-style lines (CHROM POS ID REF ALT) for all variants that
    have an rsID and a grch38_position defined.
    """
    lines = []
    for gene_name, info in ALL_GENE_THERAPY_TARGETS.items():
        for var in info.get("key_variants", []):
            rsid = var.get("rsid")
            pos_str = var.get("grch38_position")
            ref = var.get("ref_allele")
            alt = var.get("alt_allele")
            if rsid and pos_str and ref and alt:
                # Parse "chr7:117,559,591-117,559,594" or "chr7:117,587,806"
                chrom_part, coord_part = pos_str.split(":")
                coord_clean = coord_part.replace(",", "")
                if "-" in coord_clean:
                    pos = coord_clean.split("-")[0]
                else:
                    pos = coord_clean
                lines.append(f"{chrom_part}\t{pos}\t{rsid}\t{ref}\t{alt}")
    return lines
