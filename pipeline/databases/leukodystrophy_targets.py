"""
Leukodystrophy gene-therapy / CRISPR targets -- GRCh38 (hg38) coordinates.

Leukodystrophies are inherited white-matter disorders caused by defects in
myelin formation or maintenance.  This is one of the most advanced areas
of gene therapy -- two products are already approved (Libmeldy for MLD,
Skysona for cerebral ALD).

Categories:
  1. Metachromatic Leukodystrophy -- MLD (ARSA) -- APPROVED gene therapy
  2. Krabbe Disease (GALC) -- active clinical trials
  3. X-linked Adrenoleukodystrophy (ABCD1) -- APPROVED gene therapy
  4. Canavan Disease (ASPA) -- active clinical trials
  5. Alexander Disease (GFAP) -- ASO approach
  6. Pelizaeus-Merzbacher Disease (PLP1)

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - ClinicalTrials.gov
  - Orchard Therapeutics (Libmeldy), bluebird bio (Skysona/Lenti-D)
  - Passage Bio, Aspa Therapeutics
"""


LEUKODYSTROPHY_TARGETS = {

    # ===================================================================
    # 1. ARSA -- Metachromatic Leukodystrophy (MLD)
    # ===================================================================
    "ARSA": {
        "gene_id": 410,
        "chrom": "chr22",
        "start": 50_622_754,
        "end": 50_628_170,
        "strand": "-",
        "refseq": "NC_000022.11",
        "cytoband": "22q13.33",
        "exon_count": 8,
        "role": (
            "Arylsulfatase A -- lysosomal enzyme hydrolysing cerebroside "
            "3-sulfate (sulfatide).  Deficiency causes metachromatic "
            "leukodystrophy (MLD): sulfatide accumulation in oligodendrocytes "
            "and Schwann cells -> progressive demyelination in CNS and PNS.  "
            "Three clinical forms: late-infantile (most severe, onset 1-2 yr), "
            "juvenile (onset 4-12 yr), and adult (onset >16 yr)."
        ),
        "disease": "Metachromatic leukodystrophy (MLD)",
        "omim_disease": 250100,
        "omim_gene": 607574,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.459+1G>A (IVS2+1G>A, allele I)",
                "rsid": "rs28940893",
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common severe allele in Europeans (~25% of MLD "
                    "alleles).  Causes exon 2 skipping -> no enzyme activity.  "
                    "Homozygosity -> late-infantile MLD."
                ),
            },
            {
                "name": "p.Pro426Leu (P426L, allele A)",
                "rsid": "rs28940894",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mild allele (~20%).  Residual activity "
                    "~3-5%.  Homozygosity -> adult-onset MLD."
                ),
            },
        ],
        "strategy": (
            "(1) Ex vivo lentiviral HSC gene therapy -- APPROVED: Libmeldy "
            "(atidarsagene autotemcel, OTL-200) is the first approved gene "
            "therapy for MLD (EU 2020, US 2024).  Autologous CD34+ cells "
            "transduced with lentiviral ARSA, supraphysiologic enzyme production "
            "enables cross-correction of CNS cells via enzyme secretion; "
            "(2) AAV-mediated CNS-directed gene replacement (intrathecal/IC "
            "AAVrh10-ARSA or AAV9-ARSA) for patients who cannot undergo HSC "
            "transplant; "
            "(3) CRISPR-mediated ARSA insertion at safe harbor in HSCs as "
            "next-gen approach; "
            "(4) Intracerebroventricular enzyme replacement (Takeda TAK-611).  "
            "Delivery: lentiviral HSC (approved) or AAV9/AAVrh10 intrathecal."
        ),
        "clinical_programs": (
            "APPROVED: Libmeldy (OTL-200, Orchard Therapeutics) -- ex vivo "
            "lentiviral HSC gene therapy, EU-approved 2020 (pre-symptomatic "
            "late-infantile and early juvenile MLD), FDA-approved March 2024.  "
            "Passage Bio PBML04 (AAVhu68-ARSA intracisternal) Phase 1/2.  "
            "Takeda TAK-611 (ICV ERT) Phase 2."
        ),
        "conditions": ["metachromatic_leukodystrophy", "MLD", "leukodystrophy",
                        "lysosomal_storage", "demyelinating", "sulfatide"],
    },

    # ===================================================================
    # 2. GALC -- Krabbe Disease (Globoid Cell Leukodystrophy)
    # ===================================================================
    "GALC": {
        "gene_id": 2581,
        "chrom": "chr14",
        "start": 87_933_149,
        "end": 87_992_870,
        "strand": "-",
        "refseq": "NC_000014.9",
        "cytoband": "14q31.3",
        "exon_count": 17,
        "role": (
            "Galactosylceramidase (galactocerebrosidase) -- lysosomal enzyme "
            "degrading galactosylceramide and psychosine.  Deficiency causes "
            "Krabbe disease: psychosine accumulation is toxic to "
            "oligodendrocytes -> rapid demyelination.  Infantile form (90%) "
            "has onset by 6 months with irritability, spasticity, and death "
            "by age 2."
        ),
        "disease": "Krabbe disease (globoid cell leukodystrophy)",
        "omim_disease": 245200,
        "omim_gene": 606890,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.857G>A (p.Gly286Asp, 30-kb deletion)",
                "rsid": "rs398123149",
                "consequence": "large_deletion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~45% of European alleles.  30-kb deletion removing exons "
                    "11-17.  Homozygosity -> severe infantile Krabbe."
                ),
            },
            {
                "name": "p.Tyr551Ser (T551S)",
                "rsid": "rs121908165",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Late-onset variant; residual enzyme activity.",
            },
        ],
        "strategy": (
            "(1) Ex vivo lentiviral HSC gene therapy (analogous to MLD approach) "
            "-- must be delivered pre-symptomatically (newborn screening critical); "
            "(2) AAV9 intrathecal/IV gene replacement (Passage Bio PBKR03); "
            "(3) Combination HSC transplant + AAV for synergistic effect; "
            "(4) Substrate reduction therapy (preclinical).  "
            "Delivery: lentiviral HSC or AAV9 IV/IT."
        ),
        "clinical_programs": (
            "Passage Bio PBKR03 (AAVhu68-GALC) Phase 1/2 (NCT04693598).  "
            "HSCT in presymptomatic infants (identified by newborn screening) "
            "shows benefit but does not fully arrest disease.  "
            "No approved gene therapy as of 2026."
        ),
        "conditions": ["krabbe_disease", "leukodystrophy", "lysosomal_storage",
                        "demyelinating", "globoid_cell"],
    },

    # ===================================================================
    # 3. ABCD1 -- X-linked Adrenoleukodystrophy (X-ALD)
    # ===================================================================
    "ABCD1": {
        "gene_id": 215,
        "chrom": "chrX",
        "start": 153_724_851,
        "end": 153_744_548,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xq28",
        "exon_count": 10,
        "role": (
            "ATP-binding cassette transporter D1 -- peroxisomal membrane "
            "transporter importing very-long-chain fatty acids (VLCFAs) for "
            "beta-oxidation.  Deficiency causes X-linked adrenoleukodystrophy: "
            "VLCFA accumulation -> cerebral inflammatory demyelination "
            "(childhood cerebral ALD, ~35% of males), adrenomyeloneuropathy "
            "(AMN, ~45%), and adrenal insufficiency (~80%).  Most common "
            "peroxisomal disorder."
        ),
        "disease": "X-linked adrenoleukodystrophy (X-ALD)",
        "omim_disease": 300100,
        "omim_gene": 300371,
        "inheritance": "XL",
        "key_variants": [
            {
                "name": ">900 unique mutations (no hotspot)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Extreme allelic heterogeneity.  Same mutation can cause "
                    "different phenotypes within a family -- phenotype NOT "
                    "predicted by genotype (modifier genes + epigenetics)."
                ),
            },
            {
                "name": "c.1415_1416delAG (most common)",
                "rsid": "rs128624224",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "European; complete protein loss.",
            },
        ],
        "strategy": (
            "(1) Ex vivo lentiviral HSC gene therapy -- APPROVED: Skysona "
            "(elivaldogene autotemcel, eli-cel, Lenti-D) by bluebird bio, "
            "FDA-approved 2022 for early cerebral ALD in boys 4-17 with "
            "matched donor unavailable; "
            "(2) AAV9 gene replacement for AMN (targeting spinal cord); "
            "(3) CRISPR-mediated ABCD1 insertion in HSCs (next-gen); "
            "(4) Lorenzo's oil (dietary VLCFA reduction) as adjunct.  "
            "Delivery: lentiviral HSC (approved) or AAV9 intrathecal."
        ),
        "clinical_programs": (
            "APPROVED: Skysona (bluebird bio) -- ex vivo lentiviral HSC gene "
            "therapy for childhood cerebral ALD, FDA-approved Sept 2022.  "
            "However, bluebird bio withdrew from EU and Skysona carries "
            "boxed warning for hematologic malignancy risk (insertional "
            "mutagenesis).  "
            "HSCT from matched donor remains preferred when available.  "
            "Minoryx MIN-102 (PPARgamma agonist) Phase 2/3 for AMN."
        ),
        "conditions": ["adrenoleukodystrophy", "X_ALD", "leukodystrophy",
                        "peroxisomal", "demyelinating", "adrenomyeloneuropathy",
                        "AMN"],
    },

    # ===================================================================
    # 4. ASPA -- Canavan Disease
    # ===================================================================
    "ASPA": {
        "gene_id": 443,
        "chrom": "chr17",
        "start": 3_477_136,
        "end": 3_504_389,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17p13.2",
        "exon_count": 6,
        "role": (
            "Aspartoacylase -- enzyme hydrolysing N-acetylaspartate (NAA) to "
            "aspartate and acetate in oligodendrocytes.  Acetate from NAA is "
            "essential for myelin lipid synthesis.  Deficiency causes Canavan "
            "disease: NAA accumulation (spongiform white matter degeneration), "
            "macrocephaly, hypotonia progressing to spasticity, seizures.  "
            "Especially common in Ashkenazi Jewish population (carrier ~1/40)."
        ),
        "disease": "Canavan disease (spongiform leukodystrophy)",
        "omim_disease": 271900,
        "omim_gene": 608034,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Glu285Ala (E285A)",
                "rsid": "rs28940279",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common mutation worldwide (~85% of Ashkenazi Jewish "
                    "alleles, ~50% of non-Jewish European alleles).  "
                    "Near-complete loss of enzyme activity."
                ),
            },
            {
                "name": "p.Tyr231Ter (Y231X)",
                "rsid": "rs28940280",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Second most common in Ashkenazi Jewish (~15%).",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~0.9 kb, ideal for AAV) -- "
            "one of the earliest CNS gene therapy targets (first AAV trial "
            "by Janson/Leone, 2002); "
            "(2) Aspa Therapeutics (now Myrtelle) rAAV-Olig001-ASPA using "
            "oligodendrocyte-tropic AAV capsid; "
            "(3) Base editing to correct E285A (single nucleotide change); "
            "(4) CRISPR insertion at safe harbor in neural progenitors.  "
            "Delivery: AAV9 IV (crosses BBB in infants) or intracerebroventricular."
        ),
        "clinical_programs": (
            "Myrtelle rAAV-Olig001-ASPA Phase 1/2 (NCT04998396) -- "
            "oligodendrocyte-tropic capsid for improved cell targeting.  "
            "Passage Bio PBCA01 (AAVhu68-ASPA) in development.  "
            "No approved gene therapy as of 2026."
        ),
        "conditions": ["canavan_disease", "leukodystrophy", "spongiform",
                        "NAA_accumulation", "demyelinating"],
    },

    # ===================================================================
    # 5. GFAP -- Alexander Disease
    # ===================================================================
    "GFAP": {
        "gene_id": 2670,
        "chrom": "chr17",
        "start": 44_902_959,
        "end": 44_912_669,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.31",
        "exon_count": 9,
        "role": (
            "Glial fibrillary acidic protein -- type III intermediate filament "
            "protein, principal component of astrocyte cytoskeleton.  "
            "Gain-of-function mutations cause Alexander disease: toxic GFAP "
            "aggregation in astrocytes forms Rosenthal fibers -> astrocyte "
            "dysfunction -> white matter degeneration.  Infantile form is "
            "most severe (megalencephaly, seizures, developmental regression)."
        ),
        "disease": "Alexander disease",
        "omim_disease": 203450,
        "omim_gene": 137780,
        "inheritance": "AD (almost all de novo)",
        "key_variants": [
            {
                "name": "p.Arg79Cys (R79C)",
                "rsid": "rs28929474",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Hotspot in rod domain; ~10% of all Alexander disease.",
            },
            {
                "name": "p.Arg239Cys (R239C)",
                "rsid": "rs267607164",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Hotspot; severe infantile form; ~10% of cases.",
            },
            {
                "name": "p.Arg239His (R239H)",
                "rsid": "rs267607165",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Same codon as R239C; adult form more common.",
            },
        ],
        "strategy": (
            "(1) ASO-mediated GFAP reduction -- Ionis ION373 approach; "
            "reducing total GFAP (mutant is dominant toxic, not LOF) clears "
            "Rosenthal fibers in mouse models; "
            "(2) CRISPRi to silence GFAP expression in astrocytes; "
            "(3) Allele-specific silencing of mutant GFAP allele; "
            "(4) NOT suitable for gene replacement (disease is GOF/toxic, "
            "not loss-of-function).  "
            "Delivery: intrathecal ASO or AAV with astrocyte-specific promoter "
            "(GFAP promoter, ironically)."
        ),
        "clinical_programs": (
            "Ionis ION373 (anti-GFAP ASO) Phase 1/2/3 (NCT04849741) -- "
            "intrathecal ASO, first disease-modifying therapy trial for "
            "Alexander disease.  Preliminary data showed GFAP reduction.  "
            "No CRISPR trials as of 2026."
        ),
        "conditions": ["alexander_disease", "leukodystrophy", "demyelinating",
                        "rosenthal_fibers", "astrocytopathy"],
    },

    # ===================================================================
    # 6. PLP1 -- Pelizaeus-Merzbacher Disease (PMD)
    # ===================================================================
    "PLP1": {
        "gene_id": 5354,
        "chrom": "chrX",
        "start": 103_031_437,
        "end": 103_047_547,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq22.2",
        "exon_count": 7,
        "role": (
            "Proteolipid protein 1 -- major structural protein of CNS myelin "
            "(~50% of myelin protein mass).  Duplications (~60-70% of PMD) "
            "cause overexpression -> ER stress and oligodendrocyte death.  "
            "Point mutations cause misfolded PLP1 -> toxic gain-of-function.  "
            "Null mutations cause milder spastic paraplegia type 2 (SPG2)."
        ),
        "disease": "Pelizaeus-Merzbacher disease (PMD) / SPG2",
        "omim_disease": 312080,
        "omim_gene": 300401,
        "inheritance": "XL",
        "key_variants": [
            {
                "name": "PLP1 duplication (most common)",
                "rsid": None,
                "consequence": "duplication",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~60-70% of PMD.  Gene dosage sensitivity: 1 copy = SPG2, "
                    "2 copies = normal, 3 copies = PMD, 4+ copies = severe PMD."
                ),
            },
            {
                "name": "Various missense mutations",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "~20% of PMD; cause misfolded PLP1 -> ER stress.",
            },
        ],
        "strategy": (
            "(1) ASO-mediated PLP1 reduction for duplication cases "
            "(reduce overexpression to normal levels); "
            "(2) CRISPRi to titrate PLP1 expression to optimal level; "
            "(3) For null mutations (SPG2): AAV gene replacement "
            "(cDNA ~0.8 kb, ideal for AAV); "
            "(4) CRISPR excision of duplicated PLP1 copy.  "
            "Delivery: intrathecal ASO or AAV with oligodendrocyte promoter."
        ),
        "clinical_programs": (
            "Ionis PLP1 ASO preclinical.  "
            "No clinical trials as of 2026.  "
            "Human neural stem cell transplant (StemCells Inc) showed "
            "myelination in Phase 1 but company discontinued."
        ),
        "conditions": ["pelizaeus_merzbacher", "PMD", "SPG2", "leukodystrophy",
                        "hypomyelinating", "demyelinating"],
    },
}
