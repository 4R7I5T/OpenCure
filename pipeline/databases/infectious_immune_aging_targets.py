"""
Curated gene-therapy / CRISPR targets for infectious diseases, immune
dysregulation, telomere/aging disorders, pulmonary, hepatic, skeletal,
and dermatologic diseases -- GRCh38 (hg38) coordinates.

Sources:  NCBI Gene (RefSeq annotation RS_2025_08, GRCh38.p14), dbSNP,
          ClinVar, Ensembl, OMIM, ClinicalTrials.gov, PubMed,
          CRISPR Medicine News, Innovative Genomics Institute clinical-trial
          tracker.

Coordinates pulled from NCBI Gene (GRCh38.p14) on 2026-04-02.
Variant positions pulled from NCBI dbSNP / ClinVar on the same date.

IMPORTANT -- verify every coordinate against current NCBI / Ensembl releases
before production use.  Numbering can shift between patch levels.

All interventions require informed patient consent and IRB / ethics approval.

Categories:
  1. Infectious Disease Targets (CCR5, CXCR4, HBV cccDNA, HPV E6/E7)
  2. Immune Dysregulation (FOXP3, PRF1, UNC13D, STX11, STXBP2, FAS/TNFRSF6,
     FASLG, CTLA4, LRBA, STAT3, STAT1)
  3. Telomere / Aging Disorders (TERT, TERC, DKC1, WRN, LMNA)
  4. Pulmonary (BMPR2, SFTPB, SFTPC)
  5. Hepatic (UGT1A1, ATP8B1, ABCB11, ABCB4)
  6. Skeletal (FGFR3, ALPL)
  7. Dermatologic (TGM1, FLG)
"""


# ============================================================================
# 1. INFECTIOUS DISEASE TARGETS
# ============================================================================

INFECTIOUS_DISEASE_TARGETS = {

    # -------------------------------------------------------------------
    # 1a.  HIV Cure -- CCR5 knockout (Berlin patient delta32 strategy)
    # -------------------------------------------------------------------
    "CCR5": {
        "gene_id": 1234,
        "chrom": "chr3",
        "start": 46_370_142,
        "end": 46_376_206,
        "strand": "+",
        "refseq": "NC_000003.12",
        "cytoband": "3p21.31",
        "exon_count": 4,
        "role": (
            "C-C chemokine receptor type 5 -- primary co-receptor for "
            "macrophage-tropic (R5) HIV-1 entry.  Expressed on CD4+ T cells, "
            "macrophages, dendritic cells.  The natural delta32 deletion "
            "(rs333) renders individuals homozygously resistant to R5-tropic "
            "HIV-1 infection."
        ),
        "disease": "HIV-1 infection (functional cure target)",
        "omim_gene": 601373,
        "inheritance": "N/A (somatic / ex vivo editing)",
        "key_variants": [
            {
                "name": "CCR5-delta32 (c.554_585del32 / rs333)",
                "hgvs_genomic": "NC_000003.12:g.46373456_46373487del",
                "rsid": "rs333",
                "consequence": "frameshift / loss-of-function",
                "clinical_significance": "protective (HIV resistance)",
                "allele_frequency": "~10% in European populations, rare in other ancestries",
                "notes": (
                    "32-bp deletion in coding region creates premature stop "
                    "codon.  Homozygotes (~1% of Europeans) are highly resistant "
                    "to R5-tropic HIV-1.  Timothy Ray Brown ('Berlin patient') "
                    "received CCR5-delta32 homozygous stem cell transplant and "
                    "achieved HIV cure.  Adam Castillejo ('London patient') "
                    "replicated this approach."
                ),
            },
        ],
        "mutational_landscape": (
            "The delta32 deletion (rs333) is the primary variant of interest.  "
            "Other rare CCR5 variants (m303, G106R, C101X, FS299) also reduce "
            "surface expression and confer partial HIV resistance.  CCR5-null "
            "is generally well tolerated; homozygous delta32 individuals have "
            "mildly increased susceptibility to West Nile virus and possibly "
            "worse influenza outcomes."
        ),
        "crispr_strategy": (
            "Ex vivo CRISPR-Cas9 knockout of CCR5 in autologous CD4+ T cells "
            "or CD34+ HSPCs to mimic delta32.  Approaches: (1) NHEJ-mediated "
            "disruption of CCR5 exon 3 coding region; (2) precise delta32 "
            "deletion via dual-guide strategy; (3) base editing to introduce "
            "premature stop codon.  Goal is to create HIV-resistant immune "
            "reconstitution without allogeneic transplant."
        ),
        "clinical_trials": [
            {
                "nct": "NCT03164135",
                "sponsor": "Affiliated Hospital to Academy of Military Medical Sciences (China)",
                "phase": "I",
                "approach": "CRISPR-Cas9 CCR5-ablated CD34+ HSPCs",
                "status_2026": (
                    "Published in NEJM 2019: patient with HIV and AML received "
                    "CCR5-edited HSPCs.  Engraftment achieved with 5-8% editing "
                    "efficiency; insufficient for HIV cure but demonstrated "
                    "safety and long-term persistence of edited cells (>19 months)."
                ),
            },
            {
                "nct": "NCT02500849",
                "sponsor": "University of Pennsylvania / Carl June",
                "phase": "I",
                "approach": "Zinc finger nuclease (ZFN) CCR5 disruption in CD4+ T cells (SB-728-T)",
                "status_2026": (
                    "Completed.  Demonstrated safety; transient reduction in "
                    "HIV viral load during ART interruption; CCR5-modified cells "
                    "showed selective survival advantage during viremia."
                ),
            },
            {
                "nct": "NCT05397184",
                "sponsor": "Excision BioTherapeutics",
                "phase": "I/II",
                "approach": "EBT-101: in vivo AAV9-delivered CRISPR targeting HIV proviral DNA",
                "status_2026": (
                    "First in-human in vivo CRISPR gene therapy for HIV.  "
                    "EBT-101 uses dual gRNAs targeting HIV-1 LTR and Gag to "
                    "excise integrated provirus.  Phase I/II enrolling; "
                    "preliminary safety data presented.  Combines with CCR5 "
                    "disruption strategies."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1b.  HIV Cure -- CXCR4 (secondary co-receptor)
    # -------------------------------------------------------------------
    "CXCR4": {
        "gene_id": 7852,
        "chrom": "chr2",
        "start": 136_114_349,
        "end": 136_118_149,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q22.1",
        "exon_count": 2,
        "role": (
            "C-X-C chemokine receptor type 4 -- co-receptor for X4-tropic "
            "and dual-tropic HIV-1 strains.  Also essential receptor for "
            "CXCL12 (SDF-1) signaling involved in hematopoiesis, immune "
            "cell trafficking, and organ development.  CXCR4-tropic virus "
            "typically emerges in late-stage HIV disease."
        ),
        "disease": "HIV-1 infection (X4-tropic / dual-tropic strains)",
        "omim_gene": 162643,
        "inheritance": "N/A (therapeutic target)",
        "key_variants": [
            {
                "name": "WHIM syndrome mutations (e.g., p.R334X truncation)",
                "rsid": None,
                "consequence": "gain-of-function (impaired receptor internalization)",
                "clinical_significance": "pathogenic (WHIM syndrome); note resistance to X4-HIV",
                "notes": (
                    "WHIM syndrome C-terminal truncation mutations cause "
                    "constitutive CXCR4 signaling and neutropenia.  Paradoxically, "
                    "WHIM patients may have altered susceptibility to X4-tropic HIV.  "
                    "Complete CXCR4 knockout is embryonic lethal in mice -- "
                    "cannot be fully ablated like CCR5."
                ),
            },
        ],
        "mutational_landscape": (
            "Unlike CCR5, CXCR4 is essential for hematopoiesis and development.  "
            "No natural protective null alleles exist.  WHIM syndrome gain-of-function "
            "mutations cluster in the C-terminal tail (exon 2).  Full knockout is "
            "not a viable HIV strategy."
        ),
        "crispr_strategy": (
            "NOT a knockout target (essential gene).  Strategies: (1) Transient "
            "CRISPRi downregulation during acute X4-tropic infection; "
            "(2) Combined CCR5 knockout + CXCR4 partial disruption to block "
            "both co-receptors; (3) Engineering CXCR4 to prevent HIV gp120 "
            "binding while preserving CXCL12 signaling (structure-guided "
            "base editing of extracellular loops); (4) Dual CCR5-KO + "
            "anti-CXCR4 peptide expression from safe harbour locus."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No current CRISPR clinical trials directly targeting CXCR4.  "
                    "Most HIV cure strategies focus on CCR5 disruption +/- "
                    "antiretroviral intensification.  AMD3100 (plerixafor) is "
                    "an approved CXCR4 antagonist used for HSC mobilization."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 1c.  Hepatitis B Cure -- HBV cccDNA targeting
    # -------------------------------------------------------------------
    "HBV_cccDNA": {
        "gene_id": None,  # viral, not human gene
        "chrom": "N/A (viral episomal DNA in hepatocyte nucleus)",
        "start": None,
        "end": None,
        "strand": "N/A",
        "refseq": "NC_003977.2 (HBV reference genome, ~3.2 kb)",
        "cytoband": "N/A",
        "role": (
            "Hepatitis B virus covalently closed circular DNA (cccDNA) -- "
            "the viral persistence reservoir.  cccDNA resides as a stable "
            "minichromosome in hepatocyte nuclei and serves as the "
            "transcriptional template for all HBV mRNAs including pregenomic "
            "RNA.  Current antivirals (nucleos(t)ide analogues) suppress "
            "viral replication but cannot eliminate cccDNA, explaining why "
            "HBV is not cured by antiviral therapy alone."
        ),
        "disease": "Chronic Hepatitis B (functional cure)",
        "omim_gene": None,
        "inheritance": "N/A (acquired viral infection)",
        "hbv_genome_targets": [
            {
                "region": "Surface antigen (HBsAg) ORF",
                "coordinates": "nt 155-835 / 3205-3221 (overlapping ORFs)",
                "crispr_target": (
                    "Disruption of HBsAg production; targets S gene to "
                    "eliminate HBsAg secretion (functional cure endpoint)"
                ),
            },
            {
                "region": "Core/precore promoter",
                "coordinates": "nt 1742-1849",
                "crispr_target": (
                    "Disruption of HBeAg and core protein production"
                ),
            },
            {
                "region": "Polymerase ORF (overlaps envelope)",
                "coordinates": "nt 2307-1623 (spans origin)",
                "crispr_target": (
                    "Disrupting polymerase eliminates viral replication "
                    "but must account for overlapping reading frames"
                ),
            },
            {
                "region": "X gene / HBx",
                "coordinates": "nt 1374-1838",
                "crispr_target": (
                    "HBx is essential for cccDNA transcription; disruption "
                    "silences the entire cccDNA minichromosome"
                ),
            },
        ],
        "mutational_landscape": (
            "HBV cccDNA is a compact ~3.2 kb genome with 4 overlapping ORFs.  "
            "The overlapping nature means a single CRISPR cut can disrupt "
            "multiple genes simultaneously.  Key challenge: 8 HBV genotypes "
            "(A-H) with ~8% intergenotype divergence require guide RNA design "
            "against conserved regions.  HBV integrations in host genome "
            "(random, fragmented) also produce HBsAg but are NOT cccDNA-derived "
            "and require separate targeting strategy."
        ),
        "crispr_strategy": (
            "In vivo LNP-delivered CRISPR-Cas9 or Cas12 targeting multiple "
            "conserved sites on HBV cccDNA simultaneously.  Multiplexed "
            "approach: (1) Dual guides flanking HBsAg ORF to excise surface "
            "antigen gene; (2) HBx disruption to silence cccDNA transcription; "
            "(3) Must account for all 8 HBV genotypes -- target conserved "
            "regions.  Alternative: base editing to introduce stop codons in "
            "essential ORFs without creating DSBs (avoids chromosomal "
            "rearrangement risk from integrated HBV DNA)."
        ),
        "clinical_trials": [
            {
                "program": "Excision BioTherapeutics EBT-107",
                "phase": "Preclinical (IND-enabling 2025-2026)",
                "approach": (
                    "AAV-delivered CRISPR targeting HBV cccDNA; dual-guide "
                    "strategy against conserved regions.  Preclinical data "
                    "in HBV transgenic mice showed >95% reduction in HBsAg "
                    "and cccDNA.  Planning first-in-human 2026-2027."
                ),
            },
            {
                "program": "Intellia / Regeneron",
                "phase": "Preclinical",
                "approach": (
                    "LNP-delivered CRISPR-Cas9 targeting HBV cccDNA in "
                    "hepatocytes.  LNP delivery avoids AAV immunogenicity "
                    "and enables repeat dosing.  In development."
                ),
            },
        ],
        "human_host_genes_for_hbv_therapy": {
            "NTCP_SLC10A1": {
                "gene_id": 6554,
                "chrom": "chr14",
                "start": 69_806_443,
                "end": 69_829_296,
                "role": (
                    "Sodium taurocholate co-transporting polypeptide -- the "
                    "HBV/HDV cell entry receptor.  NTCP knockout or editing "
                    "could prevent new hepatocyte infection.  Bulevirtide "
                    "(Hepcludex), an NTCP entry inhibitor, is EMA-approved "
                    "for HDV."
                ),
            },
        },
    },

    # -------------------------------------------------------------------
    # 1d.  HPV -- E6/E7 oncogene targeting
    # -------------------------------------------------------------------
    "HPV_E6E7": {
        "gene_id": None,  # viral genes
        "chrom": "N/A (viral DNA integrated in host genome)",
        "start": None,
        "end": None,
        "strand": "N/A",
        "refseq": "NC_001526.4 (HPV16 reference), NC_001355.1 (HPV18 reference)",
        "role": (
            "HPV E6 and E7 oncoproteins -- the primary drivers of HPV-induced "
            "carcinogenesis.  E6 binds and degrades p53 (TP53) via E6AP; "
            "E7 inactivates Rb (RB1) family proteins.  Together they drive "
            "cell cycle dysregulation, genomic instability, and immortalization.  "
            "HPV16/18 cause ~70% of cervical cancers and majority of "
            "HPV-positive oropharyngeal cancers."
        ),
        "disease": "HPV-driven cancers (cervical, oropharyngeal, anal, penile)",
        "omim_gene": None,
        "inheritance": "N/A (acquired viral oncogene)",
        "hpv_genome_targets": [
            {
                "gene": "E6",
                "hpv16_coordinates": "nt 83-559",
                "effect": (
                    "E6 disruption restores p53 function, triggering apoptosis "
                    "and senescence in HPV-transformed cells"
                ),
            },
            {
                "gene": "E7",
                "hpv16_coordinates": "nt 562-858",
                "effect": (
                    "E7 disruption restores Rb function, causing G1 arrest "
                    "and senescence in HPV-transformed cells"
                ),
            },
        ],
        "crispr_strategy": (
            "CRISPR-Cas9 targeting of HPV E6 and/or E7 in infected/transformed "
            "cells.  Disruption of E6/E7 reactivates p53 and Rb tumour "
            "suppressors, inducing apoptosis and senescence of HPV+ cancer "
            "cells.  Approaches: (1) Dual guides to E6+E7 simultaneously; "
            "(2) Targeting the shared E6/E7 promoter (P97 in HPV16) to "
            "silence both; (3) Cas13 RNA targeting of E6/E7 mRNAs for "
            "transient knockdown.  Delivery: local intratumoral injection of "
            "AAV or lipid nanoparticle formulations.  Must be HPV-type specific "
            "due to sequence divergence between HPV16 and HPV18."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR clinical trials for HPV E6/E7 registered as of "
                    "early 2026.  Extensive preclinical data in cell lines and "
                    "xenograft models showing HPV+ cancer cell death upon E6/E7 "
                    "disruption.  Therapeutic HPV vaccines (VGX-3100 DNA vaccine, "
                    "various mRNA vaccines) are in clinical trials as an "
                    "alternative approach to clear HPV-transformed cells."
                ),
            },
        ],
    },
}


# ============================================================================
# 2. IMMUNE DYSREGULATION TARGETS
# ============================================================================

IMMUNE_DYSREGULATION_TARGETS = {

    # -------------------------------------------------------------------
    # 2a.  IPEX Syndrome -- FOXP3
    # -------------------------------------------------------------------
    "FOXP3": {
        "gene_id": 50943,
        "chrom": "chrX",
        "start": 49_250_436,
        "end": 49_264_710,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xp11.23",
        "exon_count": 12,
        "role": (
            "Forkhead box P3 -- master transcription factor for regulatory "
            "T cell (Treg) development and function.  Loss-of-function "
            "causes IPEX (Immune dysregulation, Polyendocrinopathy, "
            "Enteropathy, X-linked) syndrome -- severe multi-organ "
            "autoimmunity due to Treg deficiency."
        ),
        "disease": "IPEX Syndrome",
        "omim_disease": 304790,
        "omim_gene": 300292,
        "inheritance": "XLR",
        "key_variants": [
            {
                "name": "c.1150G>A (p.Ala384Thr / A384T)",
                "rsid": "rs104894391",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Forkhead domain; disrupts DNA binding.  One of the "
                    "most commonly reported IPEX mutations."
                ),
            },
            {
                "name": "c.1190G>A (p.Arg397Gln / R397Q)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Forkhead domain; impairs dimerization and DNA binding.",
            },
            {
                "name": "c.748T>C (p.Phe250Leu / F250L)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Leucine zipper domain; disrupts FOXP3 oligomerization.",
            },
            {
                "name": "c.227delT (p.Leu76fs)",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Null allele; severe early-onset IPEX.",
            },
        ],
        "mutational_landscape": (
            "Over 70 pathogenic variants reported in ClinVar/HGMD.  Mutation "
            "types: missense (especially in forkhead domain residues 337-421), "
            "nonsense, splice-site, small indels, and rare large deletions.  "
            "Genotype-phenotype correlation is limited; identical mutations can "
            "have variable severity.  Hemizygous in males (X-linked)."
        ),
        "crispr_strategy": (
            "Ex vivo CRISPR-mediated gene correction in autologous CD4+ T cells "
            "or HSPCs.  Approaches: (1) HDR-based correction of specific "
            "mutations in Treg progenitors; (2) Targeted FOXP3 cDNA knock-in "
            "at the endogenous locus to preserve regulatory elements; "
            "(3) Lentiviral FOXP3 gene addition in autologous HSPCs (analogous "
            "to conventional gene therapy).  Challenge: need to preserve "
            "precise FOXP3 expression levels -- overexpression can impair "
            "conventional T cell function."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR clinical trials for IPEX as of early 2026.  "
                    "Standard treatment is allogeneic HSCT.  Preclinical "
                    "lentiviral FOXP3 gene therapy programs being developed "
                    "(UCSF, Great Ormond Street Hospital).  Treg cell therapy "
                    "approaches also in development."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 2b.  Familial HLH -- PRF1 (perforin)
    # -------------------------------------------------------------------
    "PRF1": {
        "gene_id": 5551,
        "chrom": "chr10",
        "start": 70_597_348,
        "end": 70_602_741,
        "strand": "-",
        "refseq": "NC_000010.11",
        "cytoband": "10q22.1",
        "exon_count": 3,
        "role": (
            "Perforin 1 -- pore-forming cytolytic protein stored in "
            "cytotoxic granules of NK cells and cytotoxic T lymphocytes.  "
            "Essential for granule-mediated cytotoxicity.  Biallelic "
            "loss-of-function causes Familial HLH type 2 (FHL2), the "
            "most common genetic form of HLH."
        ),
        "disease": "Familial Hemophagocytic Lymphohistiocytosis type 2 (FHL2)",
        "omim_disease": 603553,
        "omim_gene": 170280,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1122G>A (p.Trp374Ter / W374X)",
                "rsid": "rs104894401",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Null allele; severe early-onset HLH.  Common in Turkish/Kurdish populations.",
            },
            {
                "name": "c.272C>T (p.Ala91Val / A91V)",
                "rsid": "rs35947132",
                "consequence": "missense",
                "clinical_significance": "risk factor / VUS",
                "notes": (
                    "Hypomorphic variant; ~3-5% carrier frequency in Europeans.  "
                    "Reduces perforin activity ~50%.  Predisposes to late-onset "
                    "HLH, especially in combination with other triggers "
                    "(infections, malignancy).  Debated pathogenicity as "
                    "isolated finding."
                ),
            },
            {
                "name": "c.1349C>T (p.Thr450Met / T450M)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "C2 domain; impairs membrane binding and pore formation.",
            },
        ],
        "mutational_landscape": (
            ">110 pathogenic variants in ClinVar.  Missense, nonsense, "
            "splice-site, and small indels throughout the gene.  Most "
            "mutations cluster in exon 2 (C2 domain) and exon 3 (MACPF domain).  "
            "Null alleles cause severe early-onset FHL; hypomorphic alleles "
            "cause late-onset or attenuated disease."
        ),
        "crispr_strategy": (
            "Ex vivo gene correction or gene addition in autologous HSPCs.  "
            "Approaches: (1) Base editing for recurrent missense variants; "
            "(2) Full PRF1 cDNA (~1.7 kb) knock-in at endogenous locus or "
            "safe harbour; (3) Lentiviral gene addition.  Small gene size "
            "is favorable for AAV/lentiviral packaging.  Must restore "
            "perforin expression specifically in cytotoxic lymphocytes."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR clinical trials for FHL2 as of early 2026.  "
                    "Lentiviral PRF1 gene therapy preclinical programs underway "
                    "(UCL / GOSH).  Standard curative treatment is allogeneic "
                    "HSCT; emapalumab (anti-IFN-gamma) approved as bridge."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 2c.  Familial HLH -- UNC13D (Munc13-4)
    # -------------------------------------------------------------------
    "UNC13D": {
        "gene_id": 201294,
        "chrom": "chr17",
        "start": 75_788_918,
        "end": 75_825_099,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q25.1",
        "exon_count": 32,
        "role": (
            "Protein unc-13 homolog D (Munc13-4) -- essential for "
            "priming of cytolytic granules for exocytosis in NK cells "
            "and CTLs.  Biallelic loss-of-function causes FHL3, the "
            "second most common genetic FHL subtype."
        ),
        "disease": "Familial Hemophagocytic Lymphohistiocytosis type 3 (FHL3)",
        "omim_disease": 608898,
        "omim_gene": 608897,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.118-308C>T (deep intronic, creates cryptic exon)",
                "rsid": "rs753767388",
                "consequence": "splice alteration (deep intronic)",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common FHL3 mutation worldwide.  Activates a "
                    "cryptic splice site creating a pseudo-exon that disrupts "
                    "the reading frame.  Often missed by exon-only sequencing."
                ),
            },
            {
                "name": "c.753+1G>T",
                "rsid": None,
                "consequence": "splice_donor_loss",
                "clinical_significance": "pathogenic",
                "notes": "Canonical splice donor disruption; null allele.",
            },
        ],
        "mutational_landscape": (
            "~80 pathogenic variants reported.  Large gene (32 exons); mutations "
            "distributed throughout.  Deep intronic variant c.118-308C>T is the "
            "most common and requires specific sequencing coverage.  Mix of "
            "splice-site, nonsense, missense, and large deletions."
        ),
        "crispr_strategy": (
            "Ex vivo correction in HSPCs.  Approaches: (1) Base editing to "
            "correct the common deep intronic c.118-308C>T variant (ideal ABE "
            "target); (2) CRISPR-mediated deletion of the cryptic exon; "
            "(3) Full cDNA knock-in (~3.2 kb, tight for AAV but feasible with "
            "dual-AAV or lentiviral delivery).  The prevalence of the single "
            "deep intronic founder mutation makes this a compelling base-editing "
            "target."
        ),
    },

    # -------------------------------------------------------------------
    # 2d.  Familial HLH -- STX11 (Syntaxin-11)
    # -------------------------------------------------------------------
    "STX11": {
        "gene_id": 8676,
        "chrom": "chr6",
        "start": 132_109_537,
        "end": 132_112_190,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6q24.2",
        "exon_count": 2,
        "role": (
            "Syntaxin-11 -- SNARE protein mediating cytolytic granule "
            "fusion with the plasma membrane.  Biallelic loss-of-function "
            "causes FHL4.  Small intronless coding region (single coding exon)."
        ),
        "disease": "Familial Hemophagocytic Lymphohistiocytosis type 4 (FHL4)",
        "omim_disease": 603552,
        "omim_gene": 605014,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.173T>G (p.Leu58Arg / L58R)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts SNARE complex formation.",
            },
            {
                "name": "c.369_370insA (p.Glu124fs)",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Founder mutation in Turkish/Kurdish populations.",
            },
        ],
        "crispr_strategy": (
            "Single-exon coding region (~840 bp) makes STX11 an ideal target "
            "for full gene replacement via HDR or prime editing.  Small size "
            "easily fits any delivery vector.  Ex vivo correction in HSPCs."
        ),
    },

    # -------------------------------------------------------------------
    # 2e.  Familial HLH -- STXBP2 (Munc18-2)
    # -------------------------------------------------------------------
    "STXBP2": {
        "gene_id": 6813,
        "chrom": "chr19",
        "start": 7_573_534,
        "end": 7_592_730,
        "strand": "+",
        "refseq": "NC_000019.10",
        "cytoband": "19p13.2",
        "exon_count": 19,
        "role": (
            "Syntaxin-binding protein 2 (Munc18-2) -- chaperone for "
            "Syntaxin-11; essential for cytolytic granule exocytosis.  "
            "Biallelic loss-of-function causes FHL5."
        ),
        "disease": "Familial Hemophagocytic Lymphohistiocytosis type 5 (FHL5)",
        "omim_disease": 613101,
        "omim_gene": 601717,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1247-1G>C",
                "rsid": "rs387906832",
                "consequence": "splice_acceptor_loss",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Canonical splice site variant; common in Saudi Arabian "
                    "families.  Leads to exon skipping and null allele."
                ),
            },
            {
                "name": "c.1430C>T (p.Pro477Leu / P477L)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts Syntaxin-11 binding interface.",
            },
        ],
        "crispr_strategy": (
            "Ex vivo gene correction in HSPCs.  Approaches: (1) Base editing "
            "for recurrent splice-site and missense variants; (2) Full cDNA "
            "knock-in (~1.8 kb, fits AAV and lentiviral vectors easily); "
            "(3) STXBP2 gene addition at safe harbour locus."
        ),
    },

    # -------------------------------------------------------------------
    # 2f.  ALPS -- FAS / TNFRSF6
    # -------------------------------------------------------------------
    "FAS": {
        "gene_id": 355,
        "chrom": "chr10",
        "start": 88_953_029,
        "end": 88_979_255,
        "strand": "+",
        "refseq": "NC_000010.11",
        "cytoband": "10q23.31",
        "exon_count": 9,
        "aliases": ["TNFRSF6", "CD95", "APO-1"],
        "role": (
            "Fas cell surface death receptor (CD95/TNFRSF6) -- death "
            "receptor mediating extrinsic apoptosis pathway.  Engagement "
            "by FasL triggers caspase-8-dependent apoptosis.  Essential for "
            "lymphocyte homeostasis.  Heterozygous dominant-negative "
            "mutations cause ALPS type Ia (autoimmune lymphoproliferative "
            "syndrome), the most common ALPS subtype."
        ),
        "disease": "Autoimmune Lymphoproliferative Syndrome type Ia (ALPS-FAS)",
        "omim_disease": 601859,
        "omim_gene": 134637,
        "inheritance": "AD (dominant-negative)",
        "key_variants": [
            {
                "name": "c.763C>T (p.Arg255Ter / R255X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Intracellular death domain; loss of apoptotic signaling.",
            },
            {
                "name": "c.718A>G (p.Thr240Ala / T240A)",
                "rsid": "rs371816832",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Death domain; dominant-negative.  Mutant FAS incorporates "
                    "into trimeric receptor complex and poisons signaling."
                ),
            },
            {
                "name": "Somatic FAS mutations (ALPS-sFAS)",
                "rsid": None,
                "consequence": "various (mosaic somatic)",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Somatic FAS mutations in double-negative T cells (DNT) "
                    "cause ALPS-sFAS.  Detected by deep sequencing of sorted "
                    "DNT cells.  Present in ~20-30% of ALPS patients."
                ),
            },
        ],
        "mutational_landscape": (
            ">100 pathogenic variants reported.  Most cluster in exon 9 "
            "(intracellular death domain, residues 230-314) with dominant-negative "
            "effect.  Extracellular domain mutations tend to cause "
            "haploinsufficiency.  Many ALPS-FAS cases are somatic mosaics.  "
            "Penetrance is incomplete (~60-70% of heterozygous carriers develop "
            "clinical ALPS)."
        ),
        "crispr_strategy": (
            "Challenging due to dominant-negative mechanism.  Approaches: "
            "(1) Allele-specific silencing of the mutant FAS allele using "
            "SNP-guided CRISPR or ASO; (2) CRISPR-mediated disruption of "
            "the mutant allele only (requires distinguishing WT from mutant); "
            "(3) For somatic ALPS-sFAS, editing is less relevant as the "
            "mutation is confined to DNT cells.  Alternative: sirolimus "
            "(mTOR inhibitor) is standard medical treatment."
        ),
    },

    # -------------------------------------------------------------------
    # 2g.  ALPS -- FASLG (Fas Ligand)
    # -------------------------------------------------------------------
    "FASLG": {
        "gene_id": 356,
        "chrom": "chr1",
        "start": 172_659_106,
        "end": 172_667_102,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q24.3",
        "exon_count": 4,
        "aliases": ["TNFSF6", "CD178", "CD95L"],
        "role": (
            "Fas ligand (FasL/CD178) -- type II transmembrane protein "
            "that binds FAS to trigger apoptosis.  Biallelic loss-of-function "
            "causes ALPS type Ib (very rare)."
        ),
        "disease": "Autoimmune Lymphoproliferative Syndrome type Ib (ALPS-FASLG)",
        "omim_disease": 601859,
        "omim_gene": 134638,
        "inheritance": "AR (homozygous) or AD (rare dominant-negative)",
        "key_variants": [
            {
                "name": "c.438G>A (p.Trp146Ter / W146X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Null allele; reported in rare homozygous ALPS-Ib cases.",
            },
        ],
        "mutational_landscape": (
            "Very rare cause of ALPS (<5 families reported).  Most ALPS is "
            "caused by FAS mutations.  FASLG mutations are mostly in the "
            "extracellular TNF homology domain."
        ),
        "crispr_strategy": (
            "Ex vivo gene correction or gene addition in HSPCs.  Small gene "
            "(~845 bp coding sequence) easily fits delivery vectors.  "
            "Base editing or HDR for specific variants."
        ),
    },

    # -------------------------------------------------------------------
    # 2h.  CTLA-4 Haploinsufficiency
    # -------------------------------------------------------------------
    "CTLA4": {
        "gene_id": 1493,
        "chrom": "chr2",
        "start": 203_867_788,
        "end": 203_873_960,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q33.2",
        "exon_count": 4,
        "role": (
            "Cytotoxic T-lymphocyte associated protein 4 -- critical immune "
            "checkpoint receptor.  Competes with CD28 for B7 ligand binding, "
            "delivers inhibitory signals, and mediates transendocytosis of "
            "B7 from APCs.  Heterozygous loss-of-function causes CTLA-4 "
            "haploinsufficiency (CHAI disease) with multi-organ autoimmunity, "
            "lymphoproliferation, and hypogammaglobulinemia."
        ),
        "disease": "CTLA-4 Haploinsufficiency (CHAI disease)",
        "omim_disease": 616100,
        "omim_gene": 123890,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "c.151C>T (p.Arg51Ter / R51X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Null allele; haploinsufficiency.  Severe phenotype.",
            },
            {
                "name": "c.208T>C (p.Tyr70His / Y70H)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Extracellular IgV domain; disrupts B7 ligand binding "
                    "and transendocytosis."
                ),
            },
            {
                "name": "c.410G>A (p.Arg137Gln / R137Q)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Intracellular domain; impairs signaling.",
            },
        ],
        "mutational_landscape": (
            "~60 pathogenic variants reported.  Penetrance ~60-70% with "
            "variable expressivity.  Mutations distributed across all exons: "
            "extracellular domain (B7 binding), transmembrane, and intracellular "
            "(signaling).  Both null alleles and dominant-negative missense.  "
            "Haploinsufficiency is the primary mechanism for most mutations."
        ),
        "crispr_strategy": (
            "Haploinsufficiency is amenable to CRISPRa (transcriptional "
            "activation of the remaining wild-type allele) or gene supplementation.  "
            "Approaches: (1) CRISPRa to upregulate WT CTLA4 from the intact "
            "allele; (2) Targeted knock-in of functional CTLA4 cDNA (~672 bp) "
            "at the endogenous locus; (3) For dominant-negative missense, "
            "allele-specific disruption of mutant allele.  Medical alternative: "
            "abatacept (CTLA4-Ig fusion) replaces CTLA-4 function systemically."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR trials for CTLA-4 haploinsufficiency as of 2026.  "
                    "Abatacept (CTLA4-Ig) is used off-label with significant "
                    "clinical benefit (restored immune regulation in >85% of "
                    "patients).  Sirolimus also used.  Gene therapy in preclinical."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 2i.  LRBA Deficiency
    # -------------------------------------------------------------------
    "LRBA": {
        "gene_id": 987,
        "chrom": "chr4",
        "start": 150_256_406,
        "end": 150_666_920,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4q31.3",
        "exon_count": 57,
        "role": (
            "Lipopolysaccharide-responsive beige-like anchor protein -- "
            "regulates CTLA-4 trafficking and recycling.  LRBA deficiency "
            "phenocopies CTLA-4 haploinsufficiency because LRBA prevents "
            "lysosomal degradation of CTLA-4.  Loss leads to reduced CTLA-4 "
            "surface expression and immune dysregulation."
        ),
        "disease": "LRBA Deficiency (immune dysregulation with autoimmunity)",
        "omim_disease": 614700,
        "omim_gene": 606453,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various truncating mutations (nonsense, frameshift, splice)",
                "rsid": None,
                "consequence": "loss-of-function",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most reported mutations are private (family-specific) "
                    "truncating variants.  Large gene (57 exons, ~2861 aa "
                    "protein).  No single common mutation."
                ),
            },
        ],
        "mutational_landscape": (
            "~50 pathogenic variants reported.  Primarily truncating (nonsense, "
            "frameshift, splice-site) and large deletions.  Very large gene "
            "(~410 kb genomic, 57 exons, ~8.6 kb CDS).  Most mutations are "
            "private.  Consanguinity is a major risk factor."
        ),
        "crispr_strategy": (
            "Challenging due to very large gene size (~8.6 kb CDS exceeds "
            "AAV packaging limit).  Approaches: (1) Dual-AAV split-intein "
            "strategy for gene replacement; (2) Correction of specific "
            "variants by base editing or prime editing (patient-specific); "
            "(3) CRISPRa to upregulate residual LRBA expression from "
            "hypomorphic alleles; (4) Bypass strategy: CRISPRa to upregulate "
            "CTLA4 directly, compensating for LRBA loss.  Medical "
            "alternative: abatacept (same as CTLA-4 haploinsufficiency)."
        ),
    },

    # -------------------------------------------------------------------
    # 2j.  STAT3 Gain-of-Function
    # -------------------------------------------------------------------
    "STAT3": {
        "gene_id": 6774,
        "chrom": "chr17",
        "start": 42_313_324,
        "end": 42_388_569,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.2",
        "exon_count": 24,
        "role": (
            "Signal transducer and activator of transcription 3 -- key "
            "transcription factor downstream of multiple cytokine receptors "
            "(IL-6, IL-10, IL-21, IL-23, IFN).  STAT3 gain-of-function (GOF) "
            "causes early-onset multi-organ autoimmunity (lymphoproliferation, "
            "cytopenias, enteropathy, diabetes, interstitial lung disease).  "
            "Note: STAT3 loss-of-function causes Hyper-IgE syndrome (AD-HIES), "
            "a distinct immunodeficiency."
        ),
        "disease": "STAT3 Gain-of-Function Disease (early-onset autoimmunity)",
        "omim_disease": 615952,
        "omim_gene": 102582,
        "inheritance": "AD (gain-of-function, dominant active)",
        "key_variants": [
            {
                "name": "c.1144C>T (p.Arg382Trp / R382W)",
                "rsid": None,
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": (
                    "SH2 domain; enhanced STAT3 dimerization and activation.  "
                    "One of the most common STAT3-GOF mutations."
                ),
            },
            {
                "name": "c.1145G>A (p.Arg382Gln / R382Q)",
                "rsid": None,
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": "SH2 domain hotspot; same codon as R382W.",
            },
            {
                "name": "c.833G>A (p.Arg278His / R278H)",
                "rsid": None,
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": "DNA-binding domain; enhanced transcriptional activity.",
            },
            {
                "name": "c.2144A>G (p.Tyr715Cys / Y715C)",
                "rsid": None,
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": "Transactivation domain; enhanced STAT3 signaling.",
            },
        ],
        "mutational_landscape": (
            "~40 distinct GOF mutations reported.  Hotspot in SH2 domain "
            "(residues 380-465) but also in DNA-binding and transactivation "
            "domains.  All are heterozygous dominant-active missense.  "
            "Penetrance is essentially 100% but expressivity is highly variable.  "
            "Distinct from STAT3-LOF (Hyper-IgE) mutations which cluster in "
            "SH2 and DNA-binding domains but cause loss of function."
        ),
        "crispr_strategy": (
            "Allele-specific silencing of the GOF mutant allele.  Approaches: "
            "(1) CRISPR-mediated disruption of mutant STAT3 allele using "
            "variant-specific guide RNA; (2) CRISPRi to reduce total STAT3 "
            "expression (must preserve some WT STAT3 function); (3) Base "
            "editing to revert the GOF mutation to wild-type (ideal for "
            "recurrent hotspot mutations like R382W/Q).  Challenge: STAT3 "
            "is essential -- complete knockout is lethal.  Must selectively "
            "target mutant allele or carefully titrate total expression."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR trials for STAT3-GOF as of 2026.  JAK inhibitors "
                    "(ruxolitinib, tofacitinib, baricitinib) are the primary "
                    "medical treatment and block upstream STAT3 activation.  "
                    "Tocilizumab (anti-IL-6R) also used.  Some patients undergo "
                    "HSCT for severe disease."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 2k.  STAT1 Gain-of-Function
    # -------------------------------------------------------------------
    "STAT1": {
        "gene_id": 6772,
        "chrom": "chr2",
        "start": 191_833_674,
        "end": 191_878_870,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q32.2",
        "exon_count": 25,
        "role": (
            "Signal transducer and activator of transcription 1 -- mediates "
            "interferon (IFN-alpha/beta/gamma) and IL-27 signaling.  STAT1 "
            "GOF causes chronic mucocutaneous candidiasis (CMC) due to "
            "impaired Th17 differentiation (STAT1 GOF suppresses STAT3 "
            "signaling required for IL-17 production).  Also causes "
            "autoimmunity, invasive infections, and cerebral aneurysms.  "
            "Note: STAT1 LOF causes susceptibility to mycobacteria and "
            "viruses (distinct phenotype)."
        ),
        "disease": "STAT1 Gain-of-Function (CMC, autoimmunity, immunodeficiency)",
        "omim_disease": 614162,
        "omim_gene": 600555,
        "inheritance": "AD (gain-of-function)",
        "key_variants": [
            {
                "name": "c.821G>A (p.Arg274Gln / R274Q)",
                "rsid": None,
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Coiled-coil domain; impaired STAT1 dephosphorylation "
                    "leading to prolonged IFN signaling.  One of the most "
                    "common STAT1-GOF mutations."
                ),
            },
            {
                "name": "c.800C>T (p.Thr267Ile / T267I)",
                "rsid": None,
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": "Coiled-coil domain; enhanced STAT1 phosphorylation.",
            },
            {
                "name": "c.1057G>A (p.Asp353Asn / D353N)",
                "rsid": None,
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": "DNA-binding domain; enhanced transcriptional activity.",
            },
            {
                "name": "c.971A>C (p.His324Pro / H324P)",
                "rsid": None,
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": "DNA-binding domain.",
            },
        ],
        "mutational_landscape": (
            "~80 distinct GOF variants reported in >400 patients.  Hotspots in "
            "coiled-coil domain (residues 235-316, impair dephosphorylation) "
            "and DNA-binding domain (residues 317-488, enhance DNA binding).  "
            "All heterozygous dominant-active.  Variable penetrance for "
            "different phenotypic features.  Mechanism: impaired STAT1 "
            "dephosphorylation by SHP2/PTPN11 leads to prolonged IFN signaling "
            "and suppression of STAT3-dependent Th17 responses."
        ),
        "crispr_strategy": (
            "Allele-specific silencing of the GOF mutant allele.  Approaches: "
            "(1) CRISPR-mediated disruption of mutant STAT1 allele using "
            "variant-specific guide RNA (must distinguish single nucleotide "
            "change); (2) Base editing to revert GOF mutation to wild-type; "
            "(3) CRISPRi to reduce total STAT1 expression to ~50% (mimicking "
            "haploinsufficiency, which is tolerated).  Challenge: complete "
            "STAT1 loss causes severe immunodeficiency -- must preserve "
            "some IFN signaling."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR trials for STAT1-GOF as of 2026.  JAK inhibitors "
                    "(ruxolitinib) are increasingly used and show benefit for "
                    "CMC and autoimmune features.  HSCT is considered for severe "
                    "cases but has high mortality (~40% in some series).  "
                    "Anti-IFN therapies in development."
                ),
            },
        ],
    },
}


# ============================================================================
# 3. TELOMERE / AGING DISORDER TARGETS
# ============================================================================

TELOMERE_AGING_TARGETS = {

    # -------------------------------------------------------------------
    # 3a.  Dyskeratosis Congenita -- TERT
    # -------------------------------------------------------------------
    "TERT": {
        "gene_id": 7015,
        "chrom": "chr5",
        "start": 1_253_167,
        "end": 1_295_068,
        "strand": "-",
        "refseq": "NC_000005.10",
        "cytoband": "5p15.33",
        "exon_count": 16,
        "role": (
            "Telomerase reverse transcriptase -- catalytic protein subunit "
            "of telomerase.  Adds TTAGGG telomeric repeats to chromosome "
            "ends using TERC RNA as a template.  Essential for telomere "
            "maintenance in stem cells, germ cells, and activated immune "
            "cells.  Heterozygous loss-of-function causes telomere biology "
            "disorders (TBDs) including dyskeratosis congenita, aplastic "
            "anemia, idiopathic pulmonary fibrosis, and hepatic cirrhosis."
        ),
        "disease": "Dyskeratosis Congenita / Telomere Biology Disorders",
        "omim_disease": 613989,
        "omim_gene": 187270,
        "inheritance": "AD (haploinsufficiency) or AR (biallelic for severe forms)",
        "key_variants": [
            {
                "name": "c.2701C>T (p.Arg901Trp / R901W)",
                "rsid": "rs199422291",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "C-terminal extension (CTE) domain; reduces telomerase "
                    "activity.  Common in pulmonary fibrosis families."
                ),
            },
            {
                "name": "c.2383-2A>G (splice variant)",
                "rsid": None,
                "consequence": "splice_acceptor_loss",
                "clinical_significance": "pathogenic",
                "notes": "Canonical splice site; null allele.",
            },
            {
                "name": "c.1234C>T (p.His412Tyr / H412Y)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Reverse transcriptase domain; reduces catalytic activity.",
            },
        ],
        "mutational_landscape": (
            ">80 pathogenic variants in ClinVar.  Missense in the reverse "
            "transcriptase (RT) and C-terminal extension (CTE) domains "
            "predominate.  Haploinsufficiency mechanism for AD forms: 50% "
            "telomerase activity is insufficient for adequate telomere "
            "maintenance, leading to progressive telomere shortening across "
            "generations (genetic anticipation)."
        ),
        "crispr_strategy": (
            "Haploinsufficiency is amenable to CRISPRa or gene supplementation.  "
            "Approaches: (1) CRISPRa to upregulate endogenous TERT from the "
            "intact allele in HSPCs (must be tightly controlled -- TERT "
            "overexpression risks oncogenesis); (2) Base editing to correct "
            "specific missense variants; (3) Targeted TERT cDNA knock-in "
            "(~3.4 kb, fits AAV with minimal promoter); (4) Transient mRNA "
            "delivery of TERT to extend telomeres without permanent genetic "
            "modification (demonstrated by Ramunas et al., FASEB J 2015).  "
            "CRITICAL SAFETY CONCERN: TERT activation is a hallmark of >90% "
            "of cancers -- any TERT-enhancing therapy must be tightly regulated."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR trials for DC/TBDs as of 2026.  Danazol (androgen) "
                    "shown to increase telomere length in TBD patients (NEJM 2016).  "
                    "Standard treatment is HSCT for bone marrow failure.  "
                    "Gene therapy preclinical programs exploring transient TERT "
                    "activation in HSPCs."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 3b.  Dyskeratosis Congenita -- TERC
    # -------------------------------------------------------------------
    "TERC": {
        "gene_id": 7012,
        "chrom": "chr3",
        "start": 169_764_516,
        "end": 169_765_067,
        "strand": "-",
        "refseq": "NC_000003.12",
        "cytoband": "3q26.2",
        "exon_count": 1,
        "role": (
            "Telomerase RNA component -- the RNA template subunit of "
            "telomerase.  Contains the 5'-CAAUCCCAAUC-3' template that "
            "directs TTAGGG repeat synthesis.  Also contains structural "
            "domains essential for telomerase RNP assembly, Cajal body "
            "localization, and catalytic activity.  Heterozygous deletions "
            "or mutations cause AD dyskeratosis congenita."
        ),
        "disease": "Dyskeratosis Congenita / Telomere Biology Disorders",
        "omim_disease": 127550,
        "omim_gene": 602322,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "n.58G>A (in template region)",
                "rsid": None,
                "consequence": "RNA functional variant",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Template region mutation; alters telomeric repeat sequence "
                    "synthesis.  Severe DC phenotype."
                ),
            },
            {
                "name": "n.110_113delGACT",
                "rsid": None,
                "consequence": "deletion in pseudoknot domain",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts essential pseudoknot structure.",
            },
            {
                "name": "n.37_43del (CR4/CR5 domain)",
                "rsid": None,
                "consequence": "deletion",
                "clinical_significance": "pathogenic",
                "notes": "Disrupts TERT binding domain.",
            },
        ],
        "mutational_landscape": (
            "~50 pathogenic variants reported.  TERC is a small non-coding "
            "RNA gene (451 nt).  Mutations throughout: template region, "
            "pseudoknot, CR4/CR5 (TERT-binding), H/ACA box (stability/processing).  "
            "Large deletions also reported.  Haploinsufficiency mechanism."
        ),
        "crispr_strategy": (
            "Small gene size (~451 bp) is ideal for gene supplementation.  "
            "Approaches: (1) CRISPRa to upregulate WT TERC allele; "
            "(2) TERC gene addition at safe harbour locus (very small insert); "
            "(3) For point mutations, base editing is feasible if within "
            "ABE/CBE editing windows.  Same oncogenesis concerns as TERT "
            "-- enhanced telomerase risks cancer."
        ),
    },

    # -------------------------------------------------------------------
    # 3c.  Dyskeratosis Congenita -- DKC1 (dyskerin)
    # -------------------------------------------------------------------
    "DKC1": {
        "gene_id": 1736,
        "chrom": "chrX",
        "start": 154_762_741,
        "end": 154_777_618,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq28",
        "exon_count": 15,
        "role": (
            "Dyskerin pseudouridine synthase -- component of the H/ACA "
            "ribonucleoprotein complex that stabilizes TERC and catalyzes "
            "pseudouridylation of rRNA and snRNA.  X-linked loss-of-function "
            "causes the classic form of dyskeratosis congenita (X-linked DC), "
            "the most severe and common genetic form of DC."
        ),
        "disease": "X-linked Dyskeratosis Congenita",
        "omim_disease": 305000,
        "omim_gene": 300126,
        "inheritance": "XLR",
        "key_variants": [
            {
                "name": "c.1058C>T (p.Ala353Val / A353V)",
                "rsid": "rs104894700",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common DKC1 mutation; hotspot.  PUA domain; "
                    "impairs TERC binding and stability."
                ),
            },
            {
                "name": "c.106T>G (p.Phe36Val / F36V)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "N-terminal domain; severe phenotype.",
            },
            {
                "name": "c.1213A>G (p.Thr405Ala / T405A)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "PUA domain; Hoyeraal-Hreidarsson syndrome (severe DC).",
            },
        ],
        "mutational_landscape": (
            ">60 pathogenic variants.  Almost all are missense (DKC1 null "
            "is likely embryonic lethal).  Hotspots in N-terminal domain "
            "and PUA (pseudouridine synthase and archaeosine transglycosylase) "
            "domain.  A353V is the most common single variant.  X-linked: "
            "hemizygous males affected, carrier females usually unaffected "
            "due to skewed X-inactivation."
        ),
        "crispr_strategy": (
            "Ex vivo correction in autologous HSPCs.  Approaches: "
            "(1) Base editing for recurrent missense variants (A353V is an "
            "ideal adenine base editor target: C>T on coding strand = G>A); "
            "(2) DKC1 cDNA (~1.5 kb) knock-in at endogenous locus; "
            "(3) Lentiviral gene addition.  Small CDS is favorable for "
            "all delivery platforms."
        ),
    },

    # -------------------------------------------------------------------
    # 3d.  Werner Syndrome -- WRN
    # -------------------------------------------------------------------
    "WRN": {
        "gene_id": 7486,
        "chrom": "chr8",
        "start": 30_887_318,
        "end": 31_028_637,
        "strand": "+",
        "refseq": "NC_000008.11",
        "cytoband": "8p12",
        "exon_count": 35,
        "role": (
            "Werner syndrome ATP-dependent helicase -- RecQ family DNA "
            "helicase with 3'-5' exonuclease activity.  Essential for DNA "
            "repair, replication, recombination, and telomere maintenance.  "
            "Biallelic loss-of-function causes Werner syndrome (adult "
            "progeroid syndrome): premature aging with onset in teens/20s "
            "including bilateral cataracts, scleroderma-like skin, short "
            "stature, premature graying/hair loss, type 2 diabetes, "
            "atherosclerosis, osteoporosis, and high cancer risk."
        ),
        "disease": "Werner Syndrome (adult progeria)",
        "omim_disease": 277700,
        "omim_gene": 604611,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1105C>T (p.Arg369Ter / R369X)",
                "rsid": "rs17847577",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common WRN mutation worldwide (founder in Japanese "
                    "population where Werner syndrome is most prevalent).  "
                    "Truncates protein before helicase domain."
                ),
            },
            {
                "name": "c.3139-1G>C (splice variant)",
                "rsid": None,
                "consequence": "splice_acceptor_loss",
                "clinical_significance": "pathogenic",
                "notes": "Founder mutation in Sardinian population.",
            },
            {
                "name": "c.2089-3024A>G (deep intronic)",
                "rsid": None,
                "consequence": "cryptic splice site activation",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Deep intronic mutation creating pseudo-exon; may be "
                    "missed by exon-only sequencing."
                ),
            },
        ],
        "mutational_landscape": (
            "~80 pathogenic variants reported.  Predominantly truncating "
            "(nonsense, frameshift, splice-site) -- nearly all cause NMD and "
            "absent protein.  Missense variants are rare and generally not "
            "confirmed pathogenic.  c.1105C>T (R369X) accounts for ~50% of "
            "Japanese WRN alleles.  WRN protein is nuclear; all known pathogenic "
            "variants either truncate the protein or disrupt nuclear localization."
        ),
        "crispr_strategy": (
            "Large gene (~4.1 kb CDS, 35 exons) is challenging for full "
            "gene replacement.  Approaches: (1) Base editing for the common "
            "R369X nonsense variant (C>T reversion); (2) Prime editing for "
            "specific truncating variants; (3) Dual-AAV split-intein for "
            "full CDS delivery; (4) CRISPR-mediated exon skipping to restore "
            "reading frame for specific frameshift/splice variants; "
            "(5) For cells with NMD, suppression of NMD + readthrough of "
            "premature stop codons.  Tissue challenge: WRN is ubiquitously "
            "expressed, requiring systemic delivery."
        ),
    },

    # -------------------------------------------------------------------
    # 3e.  Hutchinson-Gilford Progeria -- LMNA
    # -------------------------------------------------------------------
    "LMNA": {
        "gene_id": 4000,
        "chrom": "chr1",
        "start": 156_082_573,
        "end": 156_140_081,
        "strand": "+",
        "refseq": "NC_000001.11",
        "cytoband": "1q22",
        "exon_count": 12,
        "role": (
            "Lamin A/C -- type V intermediate filament protein forming the "
            "nuclear lamina.  LMNA encodes both lamin A and lamin C via "
            "alternative splicing.  Hutchinson-Gilford progeria syndrome "
            "(HGPS) is caused by a specific de novo heterozygous variant "
            "c.1824C>T (p.Gly608Gly) that activates a cryptic splice donor "
            "in exon 11, producing a truncated lamin A ('progerin') with "
            "50 amino acids deleted from the C-terminus.  Progerin is "
            "permanently farnesylated and disrupts nuclear architecture, "
            "DNA repair, and chromatin organization."
        ),
        "disease": "Hutchinson-Gilford Progeria Syndrome (HGPS)",
        "omim_disease": 176670,
        "omim_gene": 150330,
        "inheritance": "AD (de novo, dominant-negative progerin)",
        "note_other_laminopathies": (
            "LMNA mutations also cause Emery-Dreifuss muscular dystrophy, "
            "limb-girdle muscular dystrophy 1B, dilated cardiomyopathy "
            "(DCM1A), familial partial lipodystrophy (FPLD2), Charcot-Marie-Tooth "
            "type 2B1, and restrictive dermopathy.  The cardiac entry for LMNA "
            "exists in cardiac_targets.py."
        ),
        "key_variants": [
            {
                "name": "c.1824C>T (p.Gly608Gly / G608G -- synonymous but pathogenic)",
                "hgvs_genomic": "NC_000001.11:g.156109880G>A",
                "rsid": "rs387906571",
                "consequence": "synonymous BUT activates cryptic splice donor",
                "clinical_significance": "pathogenic",
                "notes": (
                    "THE classic HGPS mutation.  Does not change amino acid "
                    "(G608G) but creates a cryptic splice donor in exon 11 "
                    "that causes deletion of 150 nt from lamin A mRNA, "
                    "producing progerin (50 aa truncation retaining the CAAX "
                    "farnesylation motif).  Present in >90% of HGPS cases.  "
                    "De novo in virtually all cases."
                ),
            },
            {
                "name": "c.1822G>A (p.Gly608Ser / G608S)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Adjacent to c.1824C>T; also activates cryptic splicing "
                    "but less efficiently.  Causes atypical progeria."
                ),
            },
        ],
        "mutational_landscape": (
            "c.1824C>T accounts for >90% of classic HGPS cases.  All are de "
            "novo.  A few other LMNA mutations cause progeroid phenotypes "
            "(c.1822G>A, c.1821G>A).  HGPS is extremely rare (~1:4-8 million "
            "births).  Total LMNA pathogenic variants across all laminopathies "
            "number >400, but HGPS is caused almost exclusively by c.1824C>T."
        ),
        "crispr_strategy": (
            "c.1824C>T is an ideal base-editing target.  Approaches: "
            "(1) Adenine base editing (ABE) to revert the T>C on the coding "
            "strand (A>G on template strand) -- this eliminates the cryptic "
            "splice site and restores normal lamin A splicing.  Demonstrated "
            "by Koblan et al. (Nature 2021) in HGPS mouse model: single IV "
            "injection of AAV9-ABE corrected ~20-60% of cells in aorta and "
            "heart, dramatically extending lifespan.  (2) Antisense "
            "oligonucleotide to block the cryptic splice site.  (3) CRISPR "
            "disruption of the mutant allele.  The uniformity of the c.1824C>T "
            "mutation across nearly all HGPS patients makes this one of the "
            "most compelling single-variant gene therapy targets."
        ),
        "clinical_trials": [
            {
                "program": "Lonafarnib (Zokinvy)",
                "status": (
                    "FDA-approved 2020 -- first treatment for HGPS.  "
                    "Farnesyltransferase inhibitor that reduces progerin "
                    "toxicity.  Extends life ~2.5 years on average but does "
                    "not cure the disease."
                ),
            },
            {
                "notes": (
                    "No CRISPR clinical trials for HGPS as of 2026, but the "
                    "Koblan et al. 2021 mouse study (ABE8e delivered by AAV9 "
                    "to correct c.1824C>T) is considered one of the most "
                    "compelling preclinical demonstrations of base editing for "
                    "genetic disease.  Translation to human clinical trials "
                    "expected to follow.  Progeria Research Foundation actively "
                    "supporting gene therapy development."
                ),
            },
        ],
    },
}


# ============================================================================
# 4. PULMONARY TARGETS
# ============================================================================

PULMONARY_TARGETS = {

    # -------------------------------------------------------------------
    # 4a.  Pulmonary Arterial Hypertension -- BMPR2
    # -------------------------------------------------------------------
    "BMPR2": {
        "gene_id": 659,
        "chrom": "chr2",
        "start": 202_371_949,
        "end": 202_561_032,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2q33.1-q33.2",
        "exon_count": 13,
        "role": (
            "Bone morphogenetic protein receptor type 2 -- serine/threonine "
            "kinase receptor for BMPs.  Heterozygous loss-of-function is the "
            "most common genetic cause of heritable pulmonary arterial "
            "hypertension (HPAH, ~70% of familial PAH and ~20% of idiopathic "
            "PAH).  BMPR2 deficiency leads to uncontrolled pulmonary vascular "
            "smooth muscle proliferation and vasoconstriction."
        ),
        "disease": "Heritable Pulmonary Arterial Hypertension (HPAH)",
        "omim_disease": 178600,
        "omim_gene": 600799,
        "inheritance": "AD (haploinsufficiency, ~20% penetrance)",
        "key_variants": [
            {
                "name": "c.994C>T (p.Arg332Ter / R332X)",
                "rsid": "rs137852679",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Kinase domain; NMD of mutant transcript.  Recurrent.",
            },
            {
                "name": "c.1471C>T (p.Arg491Trp / R491W)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Kinase domain; impairs BMP signaling.",
            },
            {
                "name": "c.2695C>T (p.Arg899Ter / R899X)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Cytoplasmic tail; common truncating variant.",
            },
            {
                "name": "Exonic/whole-gene deletions",
                "rsid": None,
                "consequence": "large deletion / CNV",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Large deletions account for ~6% of HPAH.  Require MLPA "
                    "or CNV analysis for detection."
                ),
            },
        ],
        "mutational_landscape": (
            ">400 pathogenic variants reported.  ~70% are truncating (nonsense, "
            "frameshift, splice), ~25% missense (mostly in kinase domain), "
            "~5% large deletions.  Low penetrance (~20% of carriers develop PAH) "
            "suggests modifying factors.  NMD of most truncating variants leads "
            "to haploinsufficiency."
        ),
        "crispr_strategy": (
            "Haploinsufficiency mechanism supports gene supplementation.  "
            "Approaches: (1) CRISPRa to upregulate WT BMPR2 from the intact "
            "allele; (2) AAV-delivered BMPR2 cDNA (~3.1 kb, fits AAV) to "
            "pulmonary endothelium via IV injection (lung first-pass); "
            "(3) Base editing for recurrent nonsense variants (e.g., R332X); "
            "(4) mRNA/LNP delivery of BMPR2 to lung endothelium.  Target "
            "tissue is pulmonary artery endothelial cells, accessible via "
            "IV delivery."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR trials for HPAH as of 2026.  Sotatercept "
                    "(activin trap, FDA approved 2024) is the first "
                    "disease-modifying therapy for PAH, restoring BMP/TGF-beta "
                    "balance.  Gene therapy approaches in preclinical stage."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 4b.  Surfactant Deficiency -- SFTPB
    # -------------------------------------------------------------------
    "SFTPB": {
        "gene_id": 6439,
        "chrom": "chr2",
        "start": 85_657_214,
        "end": 85_668_493,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2p11.2",
        "exon_count": 11,
        "role": (
            "Surfactant protein B -- essential hydrophobic surfactant protein "
            "that reduces alveolar surface tension and prevents atelectasis.  "
            "Biallelic loss-of-function causes lethal neonatal surfactant "
            "deficiency (fatal within weeks-months without lung transplant).  "
            "SP-B is processed from a 42 kDa proprotein to the mature 8 kDa "
            "active peptide."
        ),
        "disease": "Hereditary Surfactant Protein B Deficiency",
        "omim_disease": 265120,
        "omim_gene": 178640,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.397delC (p.Pro133fs / 121ins2, legacy nomenclature)",
                "rsid": "rs397508256",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "THE most common SFTPB mutation (~60-70% of mutant alleles).  "
                    "Founder mutation.  Also called '121ins2' (legacy GAA insertion "
                    "nomenclature at codon 121).  Causes complete SP-B deficiency "
                    "and lethal neonatal respiratory failure."
                ),
            },
        ],
        "mutational_landscape": (
            "~40 pathogenic variants.  c.397delC (121ins2) dominates (~65% of "
            "alleles).  Remaining mutations are private truncating variants.  "
            "Complete SP-B deficiency is uniformly fatal.  Partial deficiency "
            "(rare hypomorphic alleles) causes chronic interstitial lung disease."
        ),
        "crispr_strategy": (
            "Gene replacement is the primary approach.  SFTPB cDNA (~1.2 kb) "
            "easily fits AAV.  Approaches: (1) AAV-delivered SFTPB to alveolar "
            "type II cells (target cell for surfactant production); "
            "(2) Base editing or prime editing for the common c.397delC; "
            "(3) Intratracheal delivery of LNP-mRNA as a bridge therapy.  "
            "Challenge: must achieve sufficient expression in AT2 cells, "
            "which are difficult to transduce.  Time-critical in neonatal "
            "presentation."
        ),
    },

    # -------------------------------------------------------------------
    # 4c.  Surfactant Deficiency -- SFTPC
    # -------------------------------------------------------------------
    "SFTPC": {
        "gene_id": 6440,
        "chrom": "chr8",
        "start": 22_112_735,
        "end": 22_116_088,
        "strand": "+",
        "refseq": "NC_000008.11",
        "cytoband": "8p21.3",
        "exon_count": 6,
        "role": (
            "Surfactant protein C -- small hydrophobic surfactant protein "
            "essential for surfactant function.  Heterozygous mutations cause "
            "dominant interstitial lung disease via toxic gain-of-function "
            "(misfolded SP-C protein triggers ER stress and AT2 cell injury).  "
            "Variable age of onset from neonatal to adult."
        ),
        "disease": "Surfactant Dysfunction / Interstitial Lung Disease",
        "omim_disease": 610913,
        "omim_gene": 178620,
        "inheritance": "AD (dominant-negative / toxic gain-of-function)",
        "key_variants": [
            {
                "name": "c.218T>C (p.Ile73Thr / I73T)",
                "rsid": "rs121918544",
                "consequence": "missense (toxic gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common SFTPC mutation.  BRICHOS domain; causes "
                    "SP-C misfolding, ER stress, and AT2 cell toxicity.  "
                    "Highly variable penetrance and severity."
                ),
            },
            {
                "name": "c.460+1G>A (exon 4 splice donor)",
                "rsid": None,
                "consequence": "splice_donor_loss",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Leads to exon 4 skipping and production of toxic "
                    "truncated protein (Delta-exon4)."
                ),
            },
        ],
        "mutational_landscape": (
            "~50 pathogenic variants.  Most are missense in the BRICHOS domain "
            "(exon 4-5) causing misfolded protein and ER stress.  Splice variants "
            "causing exon 4 skipping are also common.  Dominant-negative / toxic "
            "gain-of-function mechanism (not haploinsufficiency)."
        ),
        "crispr_strategy": (
            "Must silence the mutant allele (toxic protein mechanism).  "
            "Approaches: (1) Allele-specific CRISPR disruption of the mutant "
            "allele; (2) CRISPRi to selectively silence the mutant allele "
            "using variant-specific guides; (3) For I73T, base editing to "
            "revert C>T (ABE target); (4) Combined mutant silencing + WT "
            "gene supplementation (suppress-and-replace).  Challenge: lung "
            "delivery to AT2 cells."
        ),
    },
}


# ============================================================================
# 5. HEPATIC TARGETS
# ============================================================================

HEPATIC_TARGETS = {

    # -------------------------------------------------------------------
    # 5a.  Crigler-Najjar Syndrome -- UGT1A1
    # -------------------------------------------------------------------
    "UGT1A1": {
        "gene_id": 54658,
        "chrom": "chr2",
        "start": 233_757_013,
        "end": 233_773_299,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q37.1",
        "exon_count": 5,
        "note_gene_structure": (
            "UGT1A1 has a unique gene structure: exon 1 is specific to UGT1A1, "
            "while exons 2-5 are shared with other UGT1A family members "
            "(UGT1A3-UGT1A10).  The UGT1A locus spans ~200 kb with multiple "
            "unique exon 1 sequences spliced to common exons 2-5."
        ),
        "role": (
            "UDP-glucuronosyltransferase 1A1 -- hepatic enzyme that conjugates "
            "bilirubin with glucuronic acid for biliary excretion.  Complete "
            "absence causes Crigler-Najjar syndrome type I (CN-I, severe "
            "unconjugated hyperbilirubinemia, risk of kernicterus).  Reduced "
            "activity causes CN type II (moderate) or Gilbert syndrome (benign)."
        ),
        "disease": "Crigler-Najjar Syndrome (types I and II)",
        "omim_disease": 218800,
        "omim_gene": 191740,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1070A>G (p.Gln357Arg / Q357R)",
                "rsid": "rs200879313",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Abolishes enzyme activity; CN type I.",
            },
            {
                "name": "c.923G>A (p.Gly308Glu / G308E)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Common in CN-I; shared exon 4 (affects all UGT1A isoforms).",
            },
            {
                "name": "TA repeat polymorphism in promoter (UGT1A1*28)",
                "rsid": "rs8175347",
                "consequence": "promoter variant (TA)7 instead of (TA)6",
                "clinical_significance": "benign (Gilbert syndrome)",
                "notes": (
                    "Homozygous (TA)7 causes Gilbert syndrome (~30% reduced "
                    "activity).  Not pathogenic for CN but relevant for "
                    "irinotecan pharmacogenomics."
                ),
            },
        ],
        "mutational_landscape": (
            "~130 pathogenic variants in ClinVar.  CN-I: null variants (nonsense, "
            "splice, severe missense) in exons 1-5; CN-II: missense with residual "
            "activity.  Must distinguish exon 1 variants (UGT1A1-specific) from "
            "exons 2-5 variants (affect multiple UGT1A family members)."
        ),
        "crispr_strategy": (
            "Liver-directed gene therapy is well-suited.  Approaches: "
            "(1) AAV-delivered UGT1A1 cDNA (~1.6 kb, fits AAV easily) to "
            "hepatocytes; this is the most advanced approach.  (2) LNP-CRISPR "
            "for targeted knock-in at a hepatocyte safe harbour (e.g., albumin "
            "intron 1).  (3) Base editing for recurrent missense variants.  "
            "Liver is the ideal target organ; ~5% hepatocyte transduction is "
            "sufficient to normalize bilirubin."
        ),
        "clinical_trials": [
            {
                "nct": "NCT03466463",
                "sponsor": "Genethon / UFGT-CURE",
                "phase": "I/II",
                "approach": "AAV8-UGT1A1 gene therapy (IV infusion)",
                "status_2026": (
                    "Phase I/II for Crigler-Najjar type I.  AAV8 vector "
                    "expressing UGT1A1 under liver-specific promoter.  "
                    "Preliminary data showed bilirubin reduction in treated "
                    "patients.  Ongoing enrollment."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 5b.  PFIC Type 1 -- ATP8B1 (FIC1)
    # -------------------------------------------------------------------
    "ATP8B1": {
        "gene_id": 5205,
        "chrom": "chr18",
        "start": 55_298_969,
        "end": 55_440_346,
        "strand": "+",
        "refseq": "NC_000018.10",
        "cytoband": "18q21.31",
        "exon_count": 28,
        "role": (
            "ATPase phospholipid transporting 8B1 (FIC1) -- P4-type ATPase "
            "aminophospholipid flippase maintaining phospholipid asymmetry "
            "in the canalicular membrane of hepatocytes.  Biallelic "
            "loss-of-function causes PFIC type 1 (Byler disease) with "
            "severe cholestasis, low-GGT, pruritus, and progressive "
            "liver failure.  Also causes benign recurrent intrahepatic "
            "cholestasis type 1 (BRIC1, milder alleles)."
        ),
        "disease": "Progressive Familial Intrahepatic Cholestasis Type 1 (PFIC1)",
        "omim_disease": 211600,
        "omim_gene": 602397,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.2966T>C (p.Ile989Thr / I989T)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Severe PFIC1; disrupts ATPase activity.",
            },
            {
                "name": "c.1982T>C (p.Ile661Thr / I661T)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Founder mutation in Amish/Mennonite Byler kindred.",
            },
        ],
        "crispr_strategy": (
            "Large gene (~3.8 kb CDS, 28 exons) at the edge of AAV packaging.  "
            "Approaches: (1) AAV-delivered ATP8B1 cDNA with minimal promoter "
            "(tight fit but feasible); (2) Dual-AAV for larger constructs; "
            "(3) LNP-delivered mRNA to hepatocytes as bridge/chronic therapy; "
            "(4) Base editing for specific recurrent variants."
        ),
    },

    # -------------------------------------------------------------------
    # 5c.  PFIC Type 2 -- ABCB11 (BSEP)
    # -------------------------------------------------------------------
    "ABCB11": {
        "gene_id": 8647,
        "chrom": "chr2",
        "start": 168_916_580,
        "end": 169_031_185,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2q31.1",
        "exon_count": 28,
        "role": (
            "ATP-binding cassette subfamily B member 11 -- the bile salt "
            "export pump (BSEP).  Primary transporter of conjugated bile "
            "acids across the hepatocyte canalicular membrane.  Biallelic "
            "loss-of-function causes PFIC type 2 (most common form of PFIC), "
            "with severe cholestasis, low GGT, hepatocellular damage, and "
            "high risk of hepatocellular carcinoma and cholangiocarcinoma."
        ),
        "disease": "Progressive Familial Intrahepatic Cholestasis Type 2 (PFIC2)",
        "omim_disease": 601847,
        "omim_gene": 603201,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.890A>G (p.Glu297Gly / E297G)",
                "rsid": "rs72549402",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Common European founder mutation; nucleotide-binding "
                    "domain.  Partially responsive to UDCA."
                ),
            },
            {
                "name": "c.1416T>A (p.Asp472Glu / D482G -- legacy)",
                "rsid": "rs2287622",
                "consequence": "missense",
                "clinical_significance": "pathogenic / pharmacogenomic",
                "notes": (
                    "Also known as V444A polymorphism (rs2287622).  Associated "
                    "with drug-induced cholestasis, intrahepatic cholestasis "
                    "of pregnancy, and risk modifier for PFIC2."
                ),
            },
        ],
        "mutational_landscape": (
            ">200 pathogenic variants.  Missense, nonsense, and splice-site "
            "variants throughout.  E297G is the most common single mutation "
            "in European PFIC2 patients.  Genotype correlates with BSEP "
            "protein expression and treatment response."
        ),
        "crispr_strategy": (
            "Large gene (~3.9 kb CDS) at AAV packaging limit.  Approaches: "
            "(1) AAV-delivered ABCB11 cDNA to hepatocytes (requires compact "
            "promoter); (2) Dual-AAV strategy; (3) LNP-mRNA for chronic "
            "replacement therapy; (4) Base editing for recurrent variants "
            "(E297G).  Maralixibat (Livmarli, FDA approved 2024 for PFIC) "
            "is a medical alternative (IBAT inhibitor)."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR trials for PFIC2 as of 2026.  Maralixibat "
                    "(FDA approved 2024) and odevixibat (FDA approved 2021) "
                    "are IBAT inhibitors for PFIC.  AAV gene therapy programs "
                    "in preclinical development (Vivet Therapeutics VTX-803, "
                    "ultragenyx)."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 5d.  PFIC Type 3 -- ABCB4 (MDR3)
    # -------------------------------------------------------------------
    "ABCB4": {
        "gene_id": 5244,
        "chrom": "chr7",
        "start": 87_401_813,
        "end": 87_477_539,
        "strand": "-",
        "refseq": "NC_000007.14",
        "cytoband": "7q21.12",
        "exon_count": 28,
        "role": (
            "ATP-binding cassette subfamily B member 4 (MDR3) -- phosphatidylcholine "
            "floppase transporting phospholipids across the canalicular membrane.  "
            "Biliary phospholipids protect the bile duct epithelium from bile salt "
            "toxicity.  Biallelic loss causes PFIC3 (high-GGT cholestasis, "
            "cholesterol gallstones, biliary cirrhosis).  Heterozygous variants "
            "cause intrahepatic cholestasis of pregnancy (ICP), low phospholipid-"
            "associated cholelithiasis (LPAC), and drug-induced cholestasis."
        ),
        "disease": "Progressive Familial Intrahepatic Cholestasis Type 3 (PFIC3)",
        "omim_disease": 602347,
        "omim_gene": 171060,
        "inheritance": "AR (PFIC3) / AD with reduced penetrance (ICP, LPAC)",
        "key_variants": [
            {
                "name": "c.1769G>A (p.Gly590Glu / G590E)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "ABC transporter domain; abolishes phospholipid floppase activity.",
            },
            {
                "name": "c.1007G>T (p.Gly336Val / G336V)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Transmembrane domain; disrupts substrate translocation.",
            },
        ],
        "mutational_landscape": (
            ">100 pathogenic variants.  Mix of missense, truncating, and "
            "splice-site.  Unlike PFIC1/2, PFIC3 has high GGT.  Some missense "
            "variants are responsive to UDCA (residual MDR3 expression)."
        ),
        "crispr_strategy": (
            "Similar to ABCB11 -- large CDS (~3.8 kb) at AAV limit.  "
            "Approaches: (1) AAV gene therapy to hepatocytes; (2) Dual-AAV; "
            "(3) LNP-mRNA; (4) Base editing for recurrent variants.  "
            "UDCA is partially effective for missense variants with "
            "residual MDR3 expression."
        ),
    },
}


# ============================================================================
# 6. SKELETAL TARGETS
# ============================================================================

SKELETAL_TARGETS = {

    # -------------------------------------------------------------------
    # 6a.  Achondroplasia -- FGFR3
    # -------------------------------------------------------------------
    "FGFR3": {
        "gene_id": 2261,
        "chrom": "chr4",
        "start": 1_793_293,
        "end": 1_808_872,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4p16.3",
        "exon_count": 19,
        "role": (
            "Fibroblast growth factor receptor 3 -- transmembrane receptor "
            "tyrosine kinase that negatively regulates bone growth.  Gain-of-"
            "function mutations cause constitutive receptor activation, "
            "suppressing chondrocyte proliferation and differentiation in "
            "the growth plate.  The G380R mutation in the transmembrane "
            "domain causes achondroplasia, the most common form of skeletal "
            "dysplasia / disproportionate short stature."
        ),
        "disease": "Achondroplasia",
        "omim_disease": 100800,
        "omim_gene": 134934,
        "inheritance": "AD (gain-of-function, ~97% de novo)",
        "key_variants": [
            {
                "name": "c.1138G>A (p.Gly380Arg / G380R)",
                "hgvs_genomic": "NC_000004.12:g.1801837G>A",
                "rsid": "rs28931614",
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": (
                    "THE achondroplasia mutation; present in >97% of cases.  "
                    "Transmembrane domain; causes constitutive FGFR3 activation.  "
                    "Almost all cases are de novo.  Homozygous G380R is lethal "
                    "(thanatophoric dysplasia-like)."
                ),
            },
            {
                "name": "c.1138G>C (p.Gly380Arg / G380R -- same amino acid change)",
                "rsid": "rs28931614",
                "consequence": "missense (gain-of-function)",
                "clinical_significance": "pathogenic",
                "notes": "Alternative nucleotide change; same amino acid substitution; very rare.",
            },
        ],
        "mutational_landscape": (
            "Achondroplasia is caused almost exclusively by c.1138G>A (p.G380R) "
            "-- one of the most uniform genotype-phenotype relationships in "
            "human genetics.  Other FGFR3 GOF mutations cause different skeletal "
            "dysplasias: hypochondroplasia (N540K), thanatophoric dysplasia "
            "(R248C, K650E), and SADDAN (K650M).  FGFR3 LOF causes camptodactyly, "
            "tall stature, and hearing loss (CATSHL syndrome)."
        ),
        "crispr_strategy": (
            "The uniform G380R mutation makes achondroplasia an ideal "
            "single-variant target.  Approaches: (1) Base editing (ABE) to "
            "revert the A>G at c.1138 (coding strand G>A mutation = template "
            "strand C>T, targetable by CBE); (2) Allele-specific CRISPR "
            "disruption of the mutant allele; (3) CRISPRi to selectively "
            "silence the mutant allele.  Challenge: must be delivered to "
            "growth plate chondrocytes in utero or early childhood (before "
            "growth plate closure).  Timing is critical -- benefit diminishes "
            "after puberty.  AAV delivery to cartilage is being explored."
        ),
        "clinical_trials": [
            {
                "program": "Vosoritide (Voxzogo)",
                "status": (
                    "FDA approved 2021 (EMA 2021).  C-type natriuretic peptide "
                    "analogue that antagonizes FGFR3 signaling downstream.  "
                    "First targeted therapy for achondroplasia.  Requires daily "
                    "subcutaneous injections throughout growth period."
                ),
            },
            {
                "program": "Infigratinib (QED Therapeutics / BridgeBio)",
                "status": (
                    "Phase II PROPEL trial for children with achondroplasia.  "
                    "Oral FGFR1-3 selective tyrosine kinase inhibitor.  "
                    "Preliminary data showed accelerated growth velocity."
                ),
            },
            {
                "notes": (
                    "No CRISPR trials for achondroplasia as of 2026.  In vivo "
                    "base editing to correct G380R demonstrated in mouse models "
                    "of achondroplasia (AAV9-ABE to growth plate chondrocytes).  "
                    "Translation to human trials faces cartilage delivery challenges."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 6b.  Hypophosphatasia -- ALPL
    # -------------------------------------------------------------------
    "ALPL": {
        "gene_id": 249,
        "chrom": "chr1",
        "start": 21_553_005,
        "end": 21_621_849,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1p36.12",
        "exon_count": 12,
        "role": (
            "Alkaline phosphatase, biomineralization associated (tissue-nonspecific "
            "alkaline phosphatase, TNSALP) -- ectoenzyme that hydrolyzes "
            "inorganic pyrophosphate (PPi) to allow hydroxyapatite crystal "
            "deposition in bone.  Loss-of-function causes hypophosphatasia (HPP): "
            "defective bone and tooth mineralization.  Severity ranges from "
            "perinatal lethal to mild adult-onset (odontohypophosphatasia)."
        ),
        "disease": "Hypophosphatasia (HPP)",
        "omim_disease": 241500,
        "omim_gene": 171760,
        "inheritance": "AR (severe forms) / AD (mild forms, dominant-negative)",
        "key_variants": [
            {
                "name": "c.571G>A (p.Glu191Lys / E191K)",
                "rsid": "rs121918007",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Recurrent; common in Japanese HPP.  Active site region; "
                    "abolishes enzyme activity."
                ),
            },
            {
                "name": "c.1559delT (p.Phe520fs)",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Null allele; severe perinatal/infantile HPP.",
            },
            {
                "name": "c.979T>C (p.Phe327Leu / F327L)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Common in European HPP populations.",
            },
        ],
        "mutational_landscape": (
            ">400 pathogenic variants in ClinVar / ALPL mutation database.  "
            "Predominantly missense (~75%), with some nonsense, frameshift, "
            "and splice variants.  Mutations distributed throughout the gene.  "
            "Genotype-phenotype: biallelic severe alleles = perinatal/infantile; "
            "compound heterozygous with one mild allele = childhood/adult; "
            "heterozygous dominant-negative = adult odontohypophosphatasia."
        ),
        "crispr_strategy": (
            "Gene replacement or correction.  ALPL cDNA (~1.6 kb) fits AAV.  "
            "Approaches: (1) AAV-delivered ALPL to osteoblasts (bone-directed "
            "AAV serotypes); (2) Base editing for recurrent missense variants; "
            "(3) Liver-directed AAV to produce secreted TNSALP (enzyme "
            "replacement from liver, similar to ERT approach).  Medical "
            "alternative: asfotase alfa (Strensiq, FDA approved 2015), "
            "enzyme replacement therapy."
        ),
        "clinical_trials": [
            {
                "program": "Asfotase alfa (Strensiq)",
                "status": (
                    "FDA approved 2015 for perinatal/infantile/juvenile HPP.  "
                    "Bone-targeted ERT (TNSALP-Fc-deca-aspartate).  Requires "
                    "frequent SC injections (up to 6x/week).  Gene therapy "
                    "could provide a one-time treatment alternative."
                ),
            },
        ],
    },
}


# ============================================================================
# 7. DERMATOLOGIC TARGETS
# ============================================================================

DERMATOLOGIC_TARGETS = {

    # -------------------------------------------------------------------
    # 7a.  Lamellar Ichthyosis -- TGM1
    # -------------------------------------------------------------------
    "TGM1": {
        "gene_id": 7051,
        "chrom": "chr14",
        "start": 24_265_570,
        "end": 24_280_381,
        "strand": "+",
        "refseq": "NC_000014.9",
        "cytoband": "14q12",
        "exon_count": 15,
        "role": (
            "Transglutaminase 1 (keratinocyte transglutaminase) -- membrane-bound "
            "enzyme that crosslinks proteins (involucrin, loricrin, small proline-"
            "rich proteins) to form the cornified cell envelope, the outermost "
            "barrier layer of the epidermis.  Biallelic loss-of-function causes "
            "autosomal recessive congenital ichthyosis (ARCI), specifically "
            "lamellar ichthyosis type 1 (LI1) -- generalized plate-like scaling "
            "with collodion baby presentation at birth."
        ),
        "disease": "Lamellar Ichthyosis Type 1 (LI1 / ARCI)",
        "omim_disease": 242300,
        "omim_gene": 190195,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.877-2A>G (splice variant)",
                "rsid": None,
                "consequence": "splice_acceptor_loss",
                "clinical_significance": "pathogenic",
                "notes": "Common splice-site variant; exon skipping and null allele.",
            },
            {
                "name": "c.428G>A (p.Arg143His / R143H)",
                "rsid": "rs121909555",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Active site region; common in multiple populations.  "
                    "Abolishes transglutaminase activity."
                ),
            },
            {
                "name": "c.1187G>A (p.Arg396His / R396H)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Catalytic core domain.",
            },
        ],
        "mutational_landscape": (
            ">140 pathogenic variants.  TGM1 is the most common gene for ARCI "
            "(~30-40% of cases).  Missense, nonsense, splice-site, and small "
            "indels.  No single common founder mutation -- ethnically diverse "
            "spectrum."
        ),
        "crispr_strategy": (
            "Skin is accessible for topical/local delivery.  Approaches: "
            "(1) Ex vivo CRISPR correction in autologous keratinocytes / "
            "keratinocyte stem cells followed by engraftment as skin graft "
            "(demonstrated for other skin diseases like epidermolysis bullosa); "
            "(2) Base editing for recurrent missense variants; (3) Topical "
            "delivery of CRISPR components to basal keratinocytes (challenging "
            "due to skin barrier); (4) TGM1 cDNA (~2.5 kb) gene addition via "
            "ex vivo transduction.  Krystal Biotech's beremagene geperpavec "
            "(Vyjuvek, for EB) demonstrates the topical gene therapy paradigm."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No CRISPR or gene therapy trials for lamellar ichthyosis "
                    "as of 2026.  Topical treatments (emollients, retinoids) "
                    "are standard of care.  The success of Vyjuvek for EB "
                    "provides a roadmap for topical gene therapy for ichthyosis."
                ),
            },
        ],
    },

    # -------------------------------------------------------------------
    # 7b.  Ichthyosis Vulgaris -- FLG (filaggrin)
    # -------------------------------------------------------------------
    "FLG": {
        "gene_id": 2312,
        "chrom": "chr1",
        "start": 152_302_166,
        "end": 152_325_246,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1q21.3",
        "exon_count": 3,
        "role": (
            "Filaggrin -- large profilaggrin precursor (~400 kDa) processed "
            "into 10-12 filaggrin repeat units that aggregate keratin filaments "
            "in the stratum corneum.  Filaggrin breakdown products (urocanic "
            "acid, pyrrolidone carboxylic acid) form natural moisturizing factor "
            "(NMF) and maintain skin barrier pH.  Heterozygous loss-of-function "
            "causes ichthyosis vulgaris (most common ichthyosis, ~1:250) and "
            "is the strongest genetic risk factor for atopic dermatitis and "
            "allergic sensitization."
        ),
        "disease": "Ichthyosis Vulgaris / Atopic Dermatitis (major risk factor)",
        "omim_disease": 146700,
        "omim_gene": 135940,
        "inheritance": "semidominant (heterozygous = mild IV; homozygous = moderate IV)",
        "key_variants": [
            {
                "name": "c.2282del4 (p.Ser761fs / 2282del4 -- R501X legacy)",
                "rsid": "rs558269137",
                "consequence": "frameshift / premature stop",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common FLG null allele in European populations "
                    "(~4% carrier frequency).  Also called R501X in old "
                    "nomenclature.  Located in FLG repeat 1."
                ),
            },
            {
                "name": "c.6867delAG (p.Ser2290fs / 2282del4 #2 -- legacy S3247X)",
                "rsid": "rs200995596",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Second most common European null allele.  Located in "
                    "FLG repeat 4.  Also called 2282del4."
                ),
            },
            {
                "name": "c.3321delA (p.Ser1108fs -- legacy 3321delA)",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Common in East Asian populations (Japanese, Chinese, Korean).",
            },
        ],
        "mutational_landscape": (
            ">40 loss-of-function variants identified across populations.  "
            "Almost all are truncating (nonsense, frameshift) within the "
            "repetitive filaggrin repeat units in exon 3.  The repetitive "
            "structure (~12 near-identical ~972 bp repeats) makes FLG "
            "extremely difficult to sequence by short-read NGS and nearly "
            "impossible to analyze by standard bioinformatics pipelines.  "
            "Long-read sequencing (PacBio/ONT) recommended.  Copy number "
            "variation in repeat number (10-12 repeats) is common and "
            "modifies disease severity."
        ),
        "crispr_strategy": (
            "EXTREMELY CHALLENGING target for CRISPR.  Issues: (1) Very large "
            "protein (~400 kDa profilaggrin, >12 kb CDS) -- far exceeds any "
            "viral vector packaging limit; (2) Highly repetitive sequence "
            "makes specific guide RNA design nearly impossible; (3) Short-read "
            "sequencing cannot map edits accurately.  Alternative approaches: "
            "(1) CRISPRa to upregulate the remaining WT allele in heterozygotes; "
            "(2) Protein/peptide replacement therapy (recombinant filaggrin "
            "monomers in topical formulation); (3) Small molecule NMF "
            "supplementation; (4) mRNA delivery of individual filaggrin "
            "repeats (partial restoration).  FLG is not a practical CRISPR "
            "gene correction target with current technology."
        ),
        "clinical_trials": [
            {
                "notes": (
                    "No gene therapy trials for FLG deficiency as of 2026.  "
                    "Dupilumab (anti-IL-4Ra) is approved for atopic dermatitis "
                    "and addresses downstream inflammation but does not correct "
                    "the barrier defect.  Topical therapies restoring skin "
                    "barrier function are the main approach."
                ),
            },
        ],
    },
}


# ============================================================================
# COORDINATE SUMMARY TABLE (for quick pipeline reference)
# ============================================================================

GENE_COORDINATES_GRCH38 = {
    # Infectious Disease
    "CCR5":      {"chrom": "chr3",  "start": 46_370_142,  "end": 46_376_206},
    "CXCR4":     {"chrom": "chr2",  "start": 136_114_349, "end": 136_118_149},
    # HBV/HPV are viral -- no human coordinates
    # Immune Dysregulation
    "FOXP3":     {"chrom": "chrX",  "start": 49_250_436,  "end": 49_264_710},
    "PRF1":      {"chrom": "chr10", "start": 70_597_348,  "end": 70_602_741},
    "UNC13D":    {"chrom": "chr17", "start": 75_788_918,  "end": 75_825_099},
    "STX11":     {"chrom": "chr6",  "start": 132_109_537, "end": 132_112_190},
    "STXBP2":    {"chrom": "chr19", "start": 7_573_534,   "end": 7_592_730},
    "FAS":       {"chrom": "chr10", "start": 88_953_029,  "end": 88_979_255},
    "FASLG":     {"chrom": "chr1",  "start": 172_659_106, "end": 172_667_102},
    "CTLA4":     {"chrom": "chr2",  "start": 203_867_788, "end": 203_873_960},
    "LRBA":      {"chrom": "chr4",  "start": 150_256_406, "end": 150_666_920},
    "STAT3":     {"chrom": "chr17", "start": 42_313_324,  "end": 42_388_569},
    "STAT1":     {"chrom": "chr2",  "start": 191_833_674, "end": 191_878_870},
    # Telomere / Aging
    "TERT":      {"chrom": "chr5",  "start": 1_253_167,   "end": 1_295_068},
    "TERC":      {"chrom": "chr3",  "start": 169_764_516, "end": 169_765_067},
    "DKC1":      {"chrom": "chrX",  "start": 154_762_741, "end": 154_777_618},
    "WRN":       {"chrom": "chr8",  "start": 30_887_318,  "end": 31_028_637},
    "LMNA":      {"chrom": "chr1",  "start": 156_082_573, "end": 156_140_081},
    # Pulmonary
    "BMPR2":     {"chrom": "chr2",  "start": 202_371_949, "end": 202_561_032},
    "SFTPB":     {"chrom": "chr2",  "start": 85_657_214,  "end": 85_668_493},
    "SFTPC":     {"chrom": "chr8",  "start": 22_112_735,  "end": 22_116_088},
    # Hepatic
    "UGT1A1":    {"chrom": "chr2",  "start": 233_757_013, "end": 233_773_299},
    "ATP8B1":    {"chrom": "chr18", "start": 55_298_969,  "end": 55_440_346},
    "ABCB11":    {"chrom": "chr2",  "start": 168_916_580, "end": 169_031_185},
    "ABCB4":     {"chrom": "chr7",  "start": 87_401_813,  "end": 87_477_539},
    # Skeletal
    "FGFR3":     {"chrom": "chr4",  "start": 1_793_293,   "end": 1_808_872},
    "ALPL":      {"chrom": "chr1",  "start": 21_553_005,  "end": 21_621_849},
    # Dermatologic
    "TGM1":      {"chrom": "chr14", "start": 24_265_570,  "end": 24_280_381},
    "FLG":       {"chrom": "chr1",  "start": 152_302_166, "end": 152_325_246},
}


# ============================================================================
# DISEASE-TO-GENE MAPPING (for pipeline routing)
# ============================================================================

INFECTIOUS_IMMUNE_AGING_DISEASES = {
    "HIV": {
        "primary_genes": ["CCR5"],
        "secondary_genes": ["CXCR4"],
        "viral_targets": ["HBV_cccDNA"],
        "category": "infectious",
    },
    "hepatitis_B": {
        "primary_genes": [],
        "viral_targets": ["HBV_cccDNA"],
        "host_genes": ["SLC10A1"],
        "category": "infectious",
    },
    "HPV_cancer": {
        "viral_targets": ["HPV_E6E7"],
        "category": "infectious",
    },
    "IPEX": {
        "primary_gene": "FOXP3",
        "omim": 304790,
        "inheritance": "XLR",
        "category": "immune_dysregulation",
    },
    "familial_HLH": {
        "primary_genes": ["PRF1", "UNC13D", "STX11", "STXBP2"],
        "subtypes": {
            "FHL2": "PRF1",
            "FHL3": "UNC13D",
            "FHL4": "STX11",
            "FHL5": "STXBP2",
        },
        "category": "immune_dysregulation",
    },
    "ALPS": {
        "primary_genes": ["FAS", "FASLG"],
        "subtypes": {
            "ALPS_Ia": "FAS",
            "ALPS_Ib": "FASLG",
        },
        "category": "immune_dysregulation",
    },
    "CTLA4_haploinsufficiency": {
        "primary_gene": "CTLA4",
        "omim": 616100,
        "category": "immune_dysregulation",
    },
    "LRBA_deficiency": {
        "primary_gene": "LRBA",
        "omim": 614700,
        "category": "immune_dysregulation",
    },
    "STAT3_GOF": {
        "primary_gene": "STAT3",
        "omim": 615952,
        "category": "immune_dysregulation",
    },
    "STAT1_GOF": {
        "primary_gene": "STAT1",
        "omim": 614162,
        "category": "immune_dysregulation",
    },
    "dyskeratosis_congenita": {
        "primary_genes": ["TERT", "TERC", "DKC1"],
        "category": "telomere_aging",
    },
    "werner_syndrome": {
        "primary_gene": "WRN",
        "omim": 277700,
        "category": "telomere_aging",
    },
    "HGPS": {
        "primary_gene": "LMNA",
        "omim": 176670,
        "key_variant": "c.1824C>T",
        "category": "telomere_aging",
    },
    "heritable_PAH": {
        "primary_gene": "BMPR2",
        "omim": 178600,
        "category": "pulmonary",
    },
    "surfactant_deficiency": {
        "primary_genes": ["SFTPB", "SFTPC"],
        "category": "pulmonary",
    },
    "crigler_najjar": {
        "primary_gene": "UGT1A1",
        "omim": 218800,
        "category": "hepatic",
    },
    "PFIC": {
        "primary_genes": ["ATP8B1", "ABCB11", "ABCB4"],
        "subtypes": {
            "PFIC1": "ATP8B1",
            "PFIC2": "ABCB11",
            "PFIC3": "ABCB4",
        },
        "category": "hepatic",
    },
    "achondroplasia": {
        "primary_gene": "FGFR3",
        "omim": 100800,
        "key_variant": "c.1138G>A (G380R)",
        "category": "skeletal",
    },
    "hypophosphatasia": {
        "primary_gene": "ALPL",
        "omim": 241500,
        "category": "skeletal",
    },
    "lamellar_ichthyosis": {
        "primary_gene": "TGM1",
        "omim": 242300,
        "category": "dermatologic",
    },
    "ichthyosis_vulgaris": {
        "primary_gene": "FLG",
        "omim": 146700,
        "category": "dermatologic",
    },
}


# ---------------------------------------------------------------------------
# Add conditions arrays for pipeline condition-based filtering
# ---------------------------------------------------------------------------
_CONDITIONS_MAP = {
    "CCR5": ["HIV", "infectious"],
    "CXCR4": ["HIV", "infectious"],
    "HBV_cccDNA": ["hepatitis_B", "HBV", "infectious"],
    "HPV_E6E7": ["HPV", "cervical_cancer", "infectious"],
    "FOXP3": ["IPEX", "immune_dysregulation", "Treg_deficiency"],
    "PRF1": ["HLH", "familial_HLH", "immune_dysregulation"],
    "UNC13D": ["HLH", "familial_HLH", "immune_dysregulation"],
    "STX11": ["HLH", "familial_HLH", "immune_dysregulation"],
    "STXBP2": ["HLH", "familial_HLH", "immune_dysregulation"],
    "FAS": ["ALPS", "immune_dysregulation", "lymphoproliferative"],
    "FASLG": ["ALPS", "immune_dysregulation"],
    "CTLA4": ["CTLA4_haploinsufficiency", "immune_dysregulation"],
    "LRBA": ["LRBA_deficiency", "immune_dysregulation"],
    "STAT3": ["STAT3_GOF", "immune_dysregulation", "autoimmune"],
    "STAT1": ["STAT1_GOF", "CMC", "immune_dysregulation"],
    "TERT": ["dyskeratosis_congenita", "telomere_disorder", "aging"],
    "TERC": ["dyskeratosis_congenita", "telomere_disorder", "aging"],
    "DKC1": ["dyskeratosis_congenita", "telomere_disorder", "aging"],
    "WRN": ["werner_syndrome", "progeroid", "aging"],
    "LMNA": ["progeria", "HGPS", "progeroid", "aging"],
    "BMPR2": ["pulmonary_arterial_hypertension", "PAH", "pulmonary"],
    "SFTPB": ["surfactant_deficiency", "pulmonary"],
    "SFTPC": ["surfactant_deficiency", "pulmonary", "ILD"],
    "UGT1A1": ["crigler_najjar", "hepatic", "hyperbilirubinemia"],
    "ATP8B1": ["PFIC", "PFIC1", "hepatic", "cholestasis"],
    "ABCB11": ["PFIC", "PFIC2", "hepatic", "cholestasis"],
    "ABCB4": ["PFIC", "PFIC3", "hepatic", "cholestasis"],
    "FGFR3": ["achondroplasia", "skeletal_dysplasia", "skeletal"],
    "ALPL": ["hypophosphatasia", "skeletal"],
    "TGM1": ["lamellar_ichthyosis", "ichthyosis", "dermatologic"],
    "FLG": ["ichthyosis_vulgaris", "atopic_dermatitis", "dermatologic"],
}

_ALL_DICTS = [
    INFECTIOUS_DISEASE_TARGETS, IMMUNE_DYSREGULATION_TARGETS,
    TELOMERE_AGING_TARGETS, PULMONARY_TARGETS, HEPATIC_TARGETS,
    SKELETAL_TARGETS, DERMATOLOGIC_TARGETS,
]
for _d in _ALL_DICTS:
    for _gene, _conds in _CONDITIONS_MAP.items():
        if _gene in _d:
            _d[_gene]["conditions"] = _conds
