"""
Hereditary hearing loss gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers non-syndromic sensorineural hearing loss genes with active gene therapy
or CRISPR research programs.

Categories:
  1. GJB2 (Connexin 26) -- most common cause of hereditary deafness
  2. OTOF (Otoferlin) -- auditory synaptopathy, active clinical trials
  3. TMC1 -- hair cell mechanotransduction, preclinical CRISPR success
  4. SLC26A4 (Pendrin) -- enlarged vestibular aqueduct / Pendred syndrome
  5. MYO15A -- hair cell stereocilia, autosomal recessive
  6. TMPRSS3 -- cochlear hair cell maturation

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - ClinicalTrials.gov
  - Decibel Therapeutics, Akouos/Lilly, Regeneron hearing gene therapy programs
  - Zheng et al., Nature Medicine 2024 (OTOF gene therapy results)
"""


HEARING_LOSS_TARGETS = {

    # -------------------------------------------------------------------
    # 1.  GJB2 -- Connexin 26
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
            "Connexin 26 -- gap junction protein forming intercellular channels "
            "in the cochlea for potassium recycling.  Mutations cause ~50% of "
            "all non-syndromic autosomal recessive hearing loss (DFNB1).  "
            "The most common cause of congenital deafness worldwide."
        ),
        "disease": "DFNB1 (non-syndromic sensorineural deafness)",
        "omim_disease": 220290,
        "omim_gene": 121011,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.35delG (p.Gly12fs)",
                "rsid": "rs80338939",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common deafness mutation in Europeans (~70% of "
                    "DFNB1 alleles in Mediterranean).  Carrier rate ~2-4% "
                    "in Southern Europe."
                ),
            },
            {
                "name": "c.235delC",
                "rsid": "rs80338943",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Most common GJB2 mutation in East Asians.",
            },
            {
                "name": "c.167delT",
                "rsid": "rs80338942",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Ashkenazi Jewish founder; carrier rate ~4%.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated gene replacement (cDNA ~0.7 kb, ideal for AAV) "
            "via round window cochlear delivery; "
            "(2) CRISPR correction of 35delG in supporting cells; "
            "(3) Dual-function construct with GJB2 + GJB6 for digenic cases.  "
            "Delivery: AAV-ie (inner ear tropic) or AAV1 via round window "
            "membrane.  Critical: early intervention before cochlear maturation."
        ),
        "clinical_programs": (
            "Decibel Therapeutics DB-OTO (AAV1-OTOF) focuses on OTOF.  "
            "Shanghai Eye & ENT Hospital -- AAV-GJB2 in preclinical.  "
            "No approved GJB2 gene therapy as of 2026.  "
            "Cochlear implants remain standard of care."
        ),
        "conditions": ["DFNB1", "hearing_loss", "deafness", "connexin26",
                        "non_syndromic_hearing_loss"],
    },

    # -------------------------------------------------------------------
    # 2.  OTOF -- Otoferlin
    # -------------------------------------------------------------------
    "OTOF": {
        "gene_id": 9381,
        "chrom": "chr2",
        "start": 26_457_734,
        "end": 26_558_419,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2p23.3",
        "exon_count": 48,
        "role": (
            "Otoferlin -- calcium sensor at inner hair cell ribbon synapses, "
            "essential for synaptic vesicle exocytosis and neurotransmitter "
            "release.  Mutations cause auditory neuropathy spectrum disorder / "
            "DFNB9.  Patients have present OAEs but absent ABR."
        ),
        "disease": "DFNB9 (auditory neuropathy / synaptopathy)",
        "omim_disease": 601071,
        "omim_gene": 603681,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gln829Ter (Q829X)",
                "rsid": "rs80356570",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Spanish founder; common in Mediterranean.",
            },
            {
                "name": "c.5098G>C (p.Glu1700Gln)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Temperature-sensitive variant; hearing fluctuates with fever.",
            },
        ],
        "strategy": (
            "(1) Dual-AAV gene replacement (cDNA ~6.0 kb exceeds AAV limit) "
            "using split-intein or overlapping approaches -- Decibel DB-OTO "
            "and Regeneron/Decibel approach; "
            "(2) Chinese groups achieved hearing restoration in OTOF children "
            "(Zheng et al., Nature Medicine 2024) using dual-AAV; "
            "(3) CRISPR correction for specific mutations.  "
            "Delivery: AAV1 or AAVAnc80L65 via round window injection."
        ),
        "clinical_programs": (
            "Regeneron/Decibel DB-OTO Phase 1/2 (NCT05821959) -- dual-AAV OTOF.  "
            "Shanghai/Fudan OTOF gene therapy (AAV1-OTOF) -- first-in-human "
            "results showing hearing restoration in children (2024).  "
            "Akouos/Lilly AK-OTOF Phase 1/2.  "
            "MOST ADVANCED hearing gene therapy program."
        ),
        "conditions": ["DFNB9", "hearing_loss", "auditory_neuropathy",
                        "deafness", "otoferlin_deficiency"],
    },

    # -------------------------------------------------------------------
    # 3.  TMC1 -- Transmembrane Channel-Like 1
    # -------------------------------------------------------------------
    "TMC1": {
        "gene_id": 117531,
        "chrom": "chr9",
        "start": 72_906_547,
        "end": 72_961_168,
        "strand": "+",
        "refseq": "NC_000009.12",
        "cytoband": "9q21.13",
        "exon_count": 24,
        "role": (
            "TMC1 -- pore-forming subunit of the hair cell mechanotransduction "
            "channel.  Recessive mutations cause DFNB7/11 (profound congenital "
            "deafness).  Dominant mutation (M412K, Beethoven mouse) causes "
            "DFNA36 (progressive hearing loss)."
        ),
        "disease": "DFNB7/11 (recessive) / DFNA36 (dominant)",
        "omim_disease": 600974,
        "omim_gene": 606706,
        "inheritance": "AR (DFNB7/11) / AD (DFNA36)",
        "key_variants": [
            {
                "name": "c.100C>T (p.Arg34Ter)",
                "rsid": None,
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Pakistani founder; recessive profound deafness.",
            },
            {
                "name": "p.Met412Lys (M412K, Beethoven)",
                "rsid": "rs80338950",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Dominant progressive hearing loss (DFNA36).  "
                    "Target of Cas9 allele-specific disruption approach "
                    "(Gao et al., Nature 2018)."
                ),
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.3 kb, fits AAV) for recessive "
            "forms; "
            "(2) Allele-specific Cas9 disruption of dominant M412K allele -- "
            "demonstrated hearing rescue in Beethoven mice (Gao et al., "
            "Nature 2018); "
            "(3) Base editing to correct M412K.  "
            "Delivery: AAV-ie or AAV9-PHP.B inner ear injection."
        ),
        "clinical_programs": (
            "No clinical trials as of 2026.  "
            "Preclinical success: Holt lab (Harvard) AAV-TMC1 restored "
            "hearing in Tmc1-null mice.  CRISPR disruption of M412K "
            "allele restored hearing in Beethoven mice."
        ),
        "conditions": ["DFNB7", "DFNB11", "DFNA36", "hearing_loss",
                        "deafness", "mechanotransduction"],
    },

    # -------------------------------------------------------------------
    # 4.  SLC26A4 -- Pendrin
    # -------------------------------------------------------------------
    "SLC26A4": {
        "gene_id": 5172,
        "chrom": "chr7",
        "start": 107_660_828,
        "end": 107_717_802,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q22.3",
        "exon_count": 21,
        "role": (
            "Pendrin -- anion transporter (chloride/iodide/bicarbonate) in "
            "inner ear and thyroid.  Biallelic mutations cause Pendred "
            "syndrome (sensorineural hearing loss + goiter) or non-syndromic "
            "DFNB4 with enlarged vestibular aqueduct (EVA)."
        ),
        "disease": "Pendred syndrome / DFNB4 (EVA)",
        "omim_disease": 274600,
        "omim_gene": 605646,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Leu236Pro (L236P)",
                "rsid": "rs111033313",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Northern European; TM domain; protein misfolding.",
            },
            {
                "name": "p.His723Arg (H723R)",
                "rsid": "rs111033256",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common in East Asians.",
            },
            {
                "name": "c.919-2A>G (IVS7-2A>G)",
                "rsid": "rs111033316",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": "Common in Northern Europeans; exon 8 skipping.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.3 kb, fits AAV); "
            "(2) Base editing for recurrent missense variants; "
            "(3) Early intervention critical before progressive hearing loss.  "
            "Delivery: AAV via round window."
        ),
        "clinical_programs": "No gene therapy trials as of 2026.",
        "conditions": ["pendred_syndrome", "DFNB4", "hearing_loss", "EVA",
                        "enlarged_vestibular_aqueduct", "deafness"],
    },

    # -------------------------------------------------------------------
    # 5.  MYO15A -- Myosin XVA
    # -------------------------------------------------------------------
    "MYO15A": {
        "gene_id": 51168,
        "chrom": "chr17",
        "start": 18_012_176,
        "end": 18_083_581,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17p11.2",
        "exon_count": 66,
        "role": (
            "Unconventional myosin XVA -- motor protein essential for "
            "stereocilia elongation in cochlear hair cells.  Deficiency "
            "causes DFNB3, one of the most common forms of autosomal "
            "recessive deafness worldwide."
        ),
        "disease": "DFNB3 (non-syndromic sensorineural deafness)",
        "omim_disease": 600316,
        "omim_gene": 602666,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Ser1176fs (c.3524delA)",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Pakistani founder mutation.",
            },
        ],
        "strategy": (
            "(1) Dual-AAV gene replacement (cDNA ~10.6 kb, too large for "
            "single AAV) using split-intein or trans-splicing; "
            "(2) Mini-gene approach using truncated MYO15A isoform 2 "
            "(~3.5 kb cDNA fits single AAV).  "
            "Delivery: AAV-ie inner ear injection."
        ),
        "clinical_programs": "Preclinical only as of 2026.",
        "conditions": ["DFNB3", "hearing_loss", "deafness", "stereocilia"],
    },

    # -------------------------------------------------------------------
    # 6.  TMPRSS3 -- Transmembrane Serine Protease 3
    # -------------------------------------------------------------------
    "TMPRSS3": {
        "gene_id": 64699,
        "chrom": "chr21",
        "start": 42_420_654,
        "end": 42_444_753,
        "strand": "+",
        "refseq": "NC_000021.9",
        "cytoband": "21q22.3",
        "exon_count": 13,
        "role": (
            "TMPRSS3 -- type II transmembrane serine protease expressed in "
            "cochlear hair cells.  Mutations cause DFNB8/10 with progressive "
            "or congenital sensorineural hearing loss."
        ),
        "disease": "DFNB8/10 (non-syndromic hearing loss)",
        "omim_disease": 601072,
        "omim_gene": 605511,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Ala306Thr (A306T)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Serine protease domain.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.4 kb, ideal for AAV); "
            "(2) Base editing for specific missense variants.  "
            "Delivery: AAV via round window."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["DFNB8", "DFNB10", "hearing_loss", "deafness"],
    },
}
