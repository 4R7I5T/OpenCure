"""
Neurodevelopmental disorder gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers monogenic causes of intellectual disability, autism spectrum disorder,
and related neurodevelopmental conditions with active gene therapy or ASO
research programs.

Categories:
  1. Rett Syndrome (MECP2) -- multiple active gene therapy trials
  2. Fragile X Syndrome (FMR1) -- epigenetic reactivation approaches
  3. Angelman Syndrome (UBE3A) -- ASO to unsilence paternal allele
  4. Phelan-McDermid / SHANK3 (monogenic ASD)
  5. SYNGAP1 Encephalopathy
  6. Kabuki Syndrome (KMT2D)
  7. Rubinstein-Taybi Syndrome (CREBBP)

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - ClinicalTrials.gov
  - Taysha Gene Therapies, Neurogene, GeneTx/Ultragenyx, Ionis
"""


NEURODEVELOPMENTAL_TARGETS = {

    # ===================================================================
    # 1. MECP2 -- Rett Syndrome
    # ===================================================================
    "MECP2": {
        "gene_id": 4204,
        "chrom": "chrX",
        "start": 154_021_573,
        "end": 154_137_103,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xq28",
        "exon_count": 4,
        "role": (
            "Methyl-CpG binding protein 2 -- transcriptional regulator that "
            "binds methylated DNA and modulates chromatin structure.  Essential "
            "for neuronal maturation and synapse function.  Loss-of-function "
            "causes Rett syndrome (almost exclusively in females due to X-linked "
            "lethality in males): normal development for 6-18 months, then "
            "regression with loss of purposeful hand movements, speech loss, "
            "stereotypies, breathing irregularities, seizures.  MECP2 "
            "duplication syndrome (males) causes severe ID + recurrent infections."
        ),
        "disease": "Rett syndrome / MECP2 duplication syndrome",
        "omim_disease": 312750,
        "omim_gene": 300005,
        "inheritance": "XL dominant (Rett); XL (duplication syndrome)",
        "key_variants": [
            {
                "name": "p.Arg168Ter (R168X)",
                "rsid": "rs28934907",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Most common Rett mutation (~12%); severe phenotype.",
            },
            {
                "name": "p.Thr158Met (T158M)",
                "rsid": "rs28934906",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "~10% of Rett; MBD domain; moderate-severe.",
            },
            {
                "name": "p.Arg255Ter (R255X)",
                "rsid": "rs61749721",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Common (~8%); TRD domain; severe.",
            },
            {
                "name": "p.Arg306Cys (R306C)",
                "rsid": "rs28934908",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "~5%; NLS region; milder, preserved speech variant.",
            },
        ],
        "strategy": (
            "(1) AAV9 gene replacement (cDNA ~1.5 kb, fits AAV) with "
            "self-regulatory cassette to prevent overexpression (MECP2 is "
            "dosage-sensitive: too much is as harmful as too little -- MECP2 "
            "duplication syndrome); "
            "(2) Neurogene NGN-401 uses miRNA-responsive element (miR-Tough "
            "Decoy) to auto-regulate MECP2 expression in cells that already "
            "have active X with MECP2; "
            "(3) X-chromosome reactivation -- CRISPR-mediated reactivation "
            "of the silenced wild-type MECP2 on the inactive X; "
            "(4) Readthrough drugs for nonsense mutations (e.g., R168X); "
            "(5) Base editing for specific missense mutations (T158M, R306C).  "
            "Delivery: AAV9 IV or intrathecal.  "
            "CRITICAL: dosage regulation is the key challenge."
        ),
        "clinical_programs": (
            "Neurogene NGN-401 (AAV9-MECP2 with miRNA regulation) Phase 1/2 "
            "(NCT05898620) -- first Rett gene therapy trial with built-in "
            "dose regulation.  "
            "Taysha TSHA-102 (AAV9-MECP2 with miRNA regulation) Phase 1/2 "
            "(NCT05606614).  "
            "Novartis/AveXis OAV201 (AAV9-MECP2) preclinical.  "
            "Acadia trofinetide (IGF-1 analogue) FDA-approved 2023 for Rett "
            "(symptomatic, not gene therapy)."
        ),
        "conditions": ["rett_syndrome", "MECP2_disorder", "neurodevelopmental",
                        "intellectual_disability", "autism_spectrum",
                        "MECP2_duplication"],
    },

    # ===================================================================
    # 2. FMR1 -- Fragile X Syndrome
    # ===================================================================
    "FMR1": {
        "gene_id": 2332,
        "chrom": "chrX",
        "start": 147_911_919,
        "end": 147_951_125,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq27.3",
        "exon_count": 17,
        "role": (
            "Fragile X messenger ribonucleoprotein 1 (FMRP) -- RNA-binding "
            "protein essential for synaptic plasticity, dendritic mRNA "
            "transport, and local translation at synapses.  CGG trinucleotide "
            "repeat expansion in the 5'UTR -> hypermethylation -> FMR1 "
            "silencing -> loss of FMRP.  Most common inherited cause of "
            "intellectual disability and monogenic cause of autism (~2-6% "
            "of ASD)."
        ),
        "disease": "Fragile X syndrome (FXS)",
        "omim_disease": 300624,
        "omim_gene": 309550,
        "inheritance": "XL (trinucleotide repeat expansion)",
        "key_variants": [
            {
                "name": "CGG repeat expansion (>200 = full mutation)",
                "rsid": None,
                "consequence": "repeat_expansion -> epigenetic_silencing",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Normal: 5-44 CGG repeats.  Premutation: 55-200 repeats "
                    "(FXTAS, FXPOI risk).  Full mutation: >200 repeats -> "
                    "CpG methylation -> FMR1 silencing -> no FMRP.  "
                    "Cannot be detected by standard short-read sequencing."
                ),
            },
        ],
        "strategy": (
            "(1) Epigenetic reactivation of silenced FMR1 using CRISPR-dCas9 "
            "fused to TET1 (demethylase) or p300 (acetyltransferase) to remove "
            "CpG methylation at the FMR1 promoter -- shown to restore FMRP "
            "in Fragile X iPSCs (Liu et al., Cell 2018); "
            "(2) CRISPR excision of the expanded CGG repeat tract; "
            "(3) AAV-delivered FMRP gene replacement (cDNA ~1.9 kb fits AAV) "
            "bypassing the silenced locus; "
            "(4) mGluR5 negative allosteric modulators (AFQ056, basimglurant) "
            "as pharmacological rescue (failed in trials).  "
            "Delivery: AAV9 IV or intrathecal."
        ),
        "clinical_programs": (
            "No gene therapy or CRISPR trials as of 2026.  "
            "CRISPR demethylation proof-of-concept in iPSCs published "
            "(Liu et al., Cell 2018).  "
            "mGluR5 antagonists (Novartis, Roche) failed in Phase 2/3.  "
            "BioMarin BMN 111 (not related) in other programs."
        ),
        "conditions": ["fragile_X_syndrome", "FXS", "neurodevelopmental",
                        "intellectual_disability", "autism_spectrum",
                        "trinucleotide_repeat"],
    },

    # ===================================================================
    # 3. UBE3A -- Angelman Syndrome
    # ===================================================================
    "UBE3A": {
        "gene_id": 7337,
        "chrom": "chr15",
        "start": 25_333_728,
        "end": 25_439_056,
        "strand": "-",
        "refseq": "NC_000015.10",
        "cytoband": "15q11.2",
        "exon_count": 16,
        "role": (
            "Ubiquitin-protein ligase E3A -- E3 ubiquitin ligase and "
            "transcriptional co-activator.  Imprinted: only the maternal "
            "allele is expressed in neurons (paternal allele is silenced by "
            "UBE3A-ATS antisense transcript).  Maternal loss causes Angelman "
            "syndrome: severe ID, absent speech, seizures, ataxia, happy "
            "demeanor.  The silenced paternal allele is structurally intact "
            "-> ideal target for unsilencing."
        ),
        "disease": "Angelman syndrome",
        "omim_disease": 105830,
        "omim_gene": 601623,
        "inheritance": "Imprinted (maternal expression in neurons)",
        "key_variants": [
            {
                "name": "Maternal 15q11-q13 deletion (~70%)",
                "rsid": None,
                "consequence": "large_deletion",
                "clinical_significance": "pathogenic",
                "notes": "~5-7 Mb deletion including UBE3A; most severe.",
            },
            {
                "name": "UBE3A point mutations (~10%)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Maternal allele truncating/missense mutations.",
            },
            {
                "name": "Paternal UPD 15 (~3%)",
                "rsid": None,
                "consequence": "uniparental_disomy",
                "clinical_significance": "pathogenic",
                "notes": "Two paternal copies -> no active UBE3A in neurons.",
            },
        ],
        "strategy": (
            "(1) ASO-mediated degradation of UBE3A-ATS antisense transcript "
            "to UNSILENCE the intact paternal UBE3A allele in neurons -- "
            "GeneTx/Ultragenyx GTX-102 approach; "
            "(2) CRISPR-mediated disruption of UBE3A-ATS to permanently "
            "activate paternal UBE3A; "
            "(3) AAV-mediated UBE3A gene replacement for deletion patients "
            "(cDNA ~2.6 kb fits AAV); "
            "(4) Artificial transcription factors (ATFs) to activate "
            "paternal UBE3A.  "
            "Delivery: intrathecal ASO or AAV9 IV."
        ),
        "clinical_programs": (
            "GeneTx/Ultragenyx GTX-102 (anti-UBE3A-ATS ASO) Phase 1/2 "
            "(NCT04259281) -- showed dose-dependent improvement in some "
            "patients; transient lower-limb weakness at high doses (resolved).  "
            "Ionis ION582 (anti-UBE3A-ATS ASO) Phase 1/2.  "
            "Roche RO7248824 (ASO) Phase 1.  "
            "CRISPR-based UBE3A-ATS disruption preclinical (multiple groups)."
        ),
        "conditions": ["angelman_syndrome", "neurodevelopmental",
                        "intellectual_disability", "imprinting_disorder",
                        "UBE3A_deficiency"],
    },

    # ===================================================================
    # 4. SHANK3 -- Phelan-McDermid Syndrome
    # ===================================================================
    "SHANK3": {
        "gene_id": 85358,
        "chrom": "chr22",
        "start": 51_080_103,
        "end": 51_157_843,
        "strand": "-",
        "refseq": "NC_000022.11",
        "cytoband": "22q13.33",
        "exon_count": 24,
        "role": (
            "SH3 and multiple ankyrin repeat domains protein 3 -- master "
            "organizer of the postsynaptic density at glutamatergic synapses.  "
            "Haploinsufficiency causes Phelan-McDermid syndrome: moderate-"
            "severe ID, absent/delayed speech, ASD features, hypotonia, "
            "seizures.  SHANK3 mutations are also found in ~1% of ASD."
        ),
        "disease": "Phelan-McDermid syndrome / SHANK3-related ASD",
        "omim_disease": 606232,
        "omim_gene": 606230,
        "inheritance": "AD (haploinsufficiency, usually de novo)",
        "key_variants": [
            {
                "name": "22q13.3 terminal deletion (most common)",
                "rsid": None,
                "consequence": "large_deletion",
                "clinical_significance": "pathogenic",
                "notes": "Variable deletion size; always includes SHANK3.",
            },
            {
                "name": "SHANK3 truncating mutations",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Point mutations in ~20% without large deletions.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement -- challenging due to multiple isoforms "
            "(largest cDNA ~5.6 kb, tight AAV fit with compact promoter); "
            "(2) CRISPRa to upregulate wild-type SHANK3 allele in neurons; "
            "(3) IGF-1 and insulin-like pathway activation as downstream "
            "rescue (Icahn School trials with IGF-1).  "
            "Delivery: AAV9 IV or intrathecal."
        ),
        "clinical_programs": (
            "No gene therapy trials as of 2026.  "
            "IGF-1 (mecasermin) Phase 2 for Phelan-McDermid showed trends.  "
            "SHANK3 mouse models show reversibility of some phenotypes "
            "with adult gene restoration."
        ),
        "conditions": ["phelan_mcdermid", "SHANK3_disorder", "neurodevelopmental",
                        "intellectual_disability", "autism_spectrum"],
    },

    # ===================================================================
    # 5. SYNGAP1 -- SYNGAP1 Encephalopathy
    # ===================================================================
    "SYNGAP1": {
        "gene_id": 8831,
        "chrom": "chr6",
        "start": 33_386_275,
        "end": 33_421_822,
        "strand": "-",
        "refseq": "NC_000006.12",
        "cytoband": "6p21.32",
        "exon_count": 19,
        "role": (
            "Synaptic Ras GTPase activating protein 1 -- negative regulator "
            "of Ras/ERK signaling at excitatory synapses.  Haploinsufficiency "
            "causes SYNGAP1 encephalopathy: ID, epilepsy (often eyelid myoclonia "
            "with absence), ASD, eating difficulties.  ~1% of all ID."
        ),
        "disease": "SYNGAP1 encephalopathy (MRD5)",
        "omim_disease": 612621,
        "omim_gene": 603384,
        "inheritance": "AD (haploinsufficiency, de novo)",
        "key_variants": [
            {
                "name": "Truncating mutations throughout",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Most are de novo nonsense/frameshift/splice.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~3.3 kb fits AAV) with neuron-"
            "specific promoter; "
            "(2) CRISPRa to upregulate wild-type allele; "
            "(3) ASO readthrough for specific nonsense mutations.  "
            "Delivery: AAV9 IV or intrathecal."
        ),
        "clinical_programs": (
            "Stoke Therapeutics STK-002 (TANGO ASO approach) exploring "
            "upregulation of SYNGAP1 productive splicing.  "
            "No gene therapy trials as of 2026."
        ),
        "conditions": ["SYNGAP1_encephalopathy", "neurodevelopmental",
                        "intellectual_disability", "epileptic_encephalopathy",
                        "autism_spectrum"],
    },

    # ===================================================================
    # 6. KMT2D -- Kabuki Syndrome
    # ===================================================================
    "KMT2D": {
        "gene_id": 8085,
        "chrom": "chr12",
        "start": 49_018_975,
        "end": 49_060_291,
        "strand": "-",
        "refseq": "NC_000012.12",
        "cytoband": "12q13.12",
        "exon_count": 54,
        "role": (
            "Lysine methyltransferase 2D (MLL4) -- histone H3K4 "
            "methyltransferase, key epigenetic activator.  Haploinsufficiency "
            "causes Kabuki syndrome type 1 (~75% of Kabuki): characteristic "
            "facial features, ID, skeletal anomalies, short stature, "
            "immunodeficiency, congenital heart defects."
        ),
        "disease": "Kabuki syndrome type 1",
        "omim_disease": 147920,
        "omim_gene": 602113,
        "inheritance": "AD (mostly de novo)",
        "key_variants": [
            {
                "name": "Truncating mutations throughout (>600 reported)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "No hotspot; gene is very large (19.5 kb cDNA).",
            },
        ],
        "strategy": (
            "(1) Gene too large for AAV (~19.5 kb cDNA); "
            "(2) CRISPRa to upregulate wild-type KMT2D allele; "
            "(3) Histone deacetylase (HDAC) inhibitors to compensate for "
            "reduced H3K4 methylation (AR-42 showed benefit in mouse models); "
            "(4) Epigenetic drugs as pharmacological bypass."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "HDAC inhibitors in preclinical development."
        ),
        "conditions": ["kabuki_syndrome", "neurodevelopmental",
                        "intellectual_disability", "epigenetic_disorder"],
    },

    # ===================================================================
    # 7. CREBBP -- Rubinstein-Taybi Syndrome
    # ===================================================================
    "CREBBP": {
        "gene_id": 1387,
        "chrom": "chr16",
        "start": 3_725_054,
        "end": 3_880_124,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16p13.3",
        "exon_count": 31,
        "role": (
            "CREB-binding protein -- histone acetyltransferase and "
            "transcriptional coactivator.  Haploinsufficiency causes "
            "Rubinstein-Taybi syndrome type 1: broad thumbs/toes, "
            "characteristic facial features, ID, short stature, increased "
            "tumor risk."
        ),
        "disease": "Rubinstein-Taybi syndrome type 1",
        "omim_disease": 180849,
        "omim_gene": 600140,
        "inheritance": "AD (mostly de novo)",
        "key_variants": [
            {
                "name": "Microdeletions and truncating mutations",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "~10% have cytogenetically visible 16p13.3 deletions.",
            },
        ],
        "strategy": (
            "(1) Gene large for AAV (~7.3 kb cDNA) but dual-AAV feasible; "
            "(2) CRISPRa to upregulate wild-type allele; "
            "(3) HDAC inhibitors (SAHA) showed memory improvement in "
            "CREBBP-haploinsufficient mice.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": "No gene therapy or CRISPR trials.",
        "conditions": ["rubinstein_taybi", "RSTS", "neurodevelopmental",
                        "intellectual_disability", "epigenetic_disorder"],
    },
}
