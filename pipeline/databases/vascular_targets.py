"""
Cerebrovascular and vascular malformation gene-therapy / CRISPR targets -- GRCh38.

Covers CADASIL, hereditary hemorrhagic telangiectasia (HHT), cerebral cavernous
malformations (CCM), and moyamoya disease.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


VASCULAR_TARGETS = {

    # ===================================================================
    # CADASIL -- NOTCH3
    # ===================================================================
    "NOTCH3": {
        "gene_id": 4854,
        "chrom": "chr19",
        "start": 15_159_038,
        "end": 15_200_995,
        "strand": "-",
        "refseq": "NC_000019.10",
        "cytoband": "19p13.12",
        "exon_count": 33,
        "role": (
            "NOTCH3 receptor -- transmembrane receptor in vascular smooth "
            "muscle cells (VSMCs).  Mutations (mostly cysteine-altering in "
            "EGF-like repeats of the extracellular domain) cause CADASIL: "
            "cerebral autosomal dominant arteriopathy with subcortical infarcts "
            "and leukoencephalopathy.  Most common monogenic cause of stroke "
            "and vascular dementia."
        ),
        "disease": "CADASIL",
        "omim_disease": 125310,
        "omim_gene": 600276,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Arg153Cys (R153C)",
                "rsid": "rs28933698",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "EGF-like repeat 4; adds unpaired cysteine -> NOTCH3 aggregation.",
            },
            {
                "name": "p.Arg182Cys (R182C)",
                "rsid": "rs28933699",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common in European cohorts.",
            },
            {
                "name": "p.Arg544Cys (R544C)",
                "rsid": "rs28933706",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "East Asian; relatively mild course.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant NOTCH3 allele using "
            "CRISPRi or Cas13 (mutant protein is toxic, not just LOF); "
            "(2) Base editing to correct cysteine-altering mutations; "
            "(3) Anti-NOTCH3 ECD antibody to clear aggregates (Roche approach).  "
            "Delivery: AAV9 IV with CNS VSMC tropism or intrathecal."
        ),
        "clinical_programs": (
            "No gene therapy or CRISPR trials as of 2026.  "
            "No disease-modifying therapies; antiplatelet agents used.  "
            "Roche anti-NOTCH3 antibody preclinical."
        ),
        "conditions": ["CADASIL", "vascular_dementia", "stroke",
                        "cerebrovascular", "vascular", "leukoencephalopathy"],
    },

    # ===================================================================
    # HHT -- ENG (Endoglin)
    # ===================================================================
    "ENG": {
        "gene_id": 2022,
        "chrom": "chr9",
        "start": 127_815_013,
        "end": 127_854_773,
        "strand": "-",
        "refseq": "NC_000009.12",
        "cytoband": "9q34.11",
        "exon_count": 15,
        "role": (
            "Endoglin -- TGF-beta co-receptor on vascular endothelium.  "
            "Haploinsufficiency causes HHT1 (hereditary hemorrhagic "
            "telangiectasia / Osler-Weber-Rendu): epistaxis, mucocutaneous "
            "telangiectases, pulmonary/hepatic/cerebral arteriovenous "
            "malformations (AVMs)."
        ),
        "disease": "HHT type 1 (Osler-Weber-Rendu)",
        "omim_disease": 187300,
        "omim_gene": 131195,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": ">600 pathogenic variants (no hotspot)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Truncating mutations cause haploinsufficiency.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type ENG allele (restore to "
            "normal dosage); "
            "(2) AAV gene replacement (cDNA ~1.9 kb fits AAV); "
            "(3) mRNA/LNP for transient endoglin augmentation.  "
            "Delivery: endothelial-tropic AAV or IV LNP."
        ),
        "clinical_programs": (
            "Bevacizumab (anti-VEGF) used off-label for HHT bleeding.  "
            "No gene therapy trials as of 2026."
        ),
        "conditions": ["HHT", "HHT1", "hereditary_hemorrhagic_telangiectasia",
                        "vascular", "AVM", "epistaxis"],
    },

    # ===================================================================
    # HHT -- ACVRL1 (ALK1)
    # ===================================================================
    "ACVRL1": {
        "gene_id": 94,
        "chrom": "chr12",
        "start": 51_907_937,
        "end": 51_918_492,
        "strand": "+",
        "refseq": "NC_000012.12",
        "cytoband": "12q13.13",
        "exon_count": 10,
        "role": (
            "ALK1 (activin receptor-like kinase 1) -- TGF-beta/BMP type I "
            "receptor on endothelium.  Haploinsufficiency causes HHT2.  "
            "HHT2 has higher rate of hepatic AVMs than HHT1."
        ),
        "disease": "HHT type 2",
        "omim_disease": 600376,
        "omim_gene": 601284,
        "inheritance": "AD (haploinsufficiency)",
        "key_variants": [
            {
                "name": "Various truncating and missense",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": ">500 variants; kinase domain mutations common.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type allele; "
            "(2) AAV gene replacement (cDNA ~1.5 kb, ideal for AAV).  "
            "Delivery: endothelial-tropic AAV."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["HHT", "HHT2", "vascular", "AVM", "hepatic_AVM"],
    },

    # ===================================================================
    # Cerebral Cavernous Malformation -- KRIT1 (CCM1)
    # ===================================================================
    "KRIT1": {
        "gene_id": 889,
        "chrom": "chr7",
        "start": 91_828_080,
        "end": 91_877_458,
        "strand": "-",
        "refseq": "NC_000007.14",
        "cytoband": "7q21.2",
        "exon_count": 20,
        "role": (
            "KRIT1 (CCM1) -- scaffolding protein at endothelial cell-cell "
            "junctions.  Haploinsufficiency + somatic second hit causes "
            "cerebral cavernous malformations (CCMs): thin-walled vascular "
            "cavities in CNS prone to hemorrhage."
        ),
        "disease": "Cerebral cavernous malformation type 1 (CCM1)",
        "omim_disease": 116860,
        "omim_gene": 604214,
        "inheritance": "AD (two-hit)",
        "key_variants": [
            {
                "name": "p.Gln455Ter (Q455X, Hispanic founder)",
                "rsid": "rs137852759",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Common founder mutation in Hispanic Americans.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type KRIT1 to prevent second hit; "
            "(2) AAV gene replacement (cDNA ~2.1 kb fits AAV); "
            "(3) Rho kinase (ROCK) inhibition as downstream therapy.  "
            "Delivery: brain endothelial-tropic AAV."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "Atorvastatin + vitamin D3 trial for CCM (NCT04345265)."
        ),
        "conditions": ["cerebral_cavernous_malformation", "CCM", "CCM1",
                        "vascular", "cerebrovascular", "brain_hemorrhage"],
    },

    # ===================================================================
    # CCM2
    # ===================================================================
    "CCM2": {
        "gene_id": 83605,
        "chrom": "chr7",
        "start": 45_027_901,
        "end": 45_096_527,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7p13",
        "exon_count": 10,
        "role": (
            "CCM2 (malcavernin/OSM) -- scaffolding protein in KRIT1-CCM2-"
            "PDCD10 complex.  Mutations cause CCM type 2 (~20% of familial CCM)."
        ),
        "disease": "Cerebral cavernous malformation type 2",
        "omim_disease": 603284,
        "omim_gene": 607929,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Truncating mutations",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Various; no hotspot.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.3 kb, ideal for AAV); "
            "(2) CRISPRa of wild-type allele."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["cerebral_cavernous_malformation", "CCM", "CCM2",
                        "vascular", "cerebrovascular"],
    },

    # ===================================================================
    # CCM3 -- PDCD10
    # ===================================================================
    "PDCD10": {
        "gene_id": 11235,
        "chrom": "chr3",
        "start": 167_439_164,
        "end": 167_528_453,
        "strand": "-",
        "refseq": "NC_000003.12",
        "cytoband": "3q26.1",
        "exon_count": 10,
        "role": (
            "PDCD10 (CCM3) -- programmed cell death 10, part of the KRIT1-"
            "CCM2-PDCD10 complex.  CCM3 mutations cause the most aggressive "
            "form of CCM with childhood-onset hemorrhages and meningiomas."
        ),
        "disease": "Cerebral cavernous malformation type 3",
        "omim_disease": 603285,
        "omim_gene": 609118,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Truncating mutations",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Most aggressive CCM form; early pediatric onset.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~0.6 kb, ideal for AAV); "
            "(2) CRISPRa to boost wild-type allele."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["cerebral_cavernous_malformation", "CCM", "CCM3",
                        "vascular", "cerebrovascular", "pediatric"],
    },

    # ===================================================================
    # Moyamoya -- RNF213
    # ===================================================================
    "RNF213": {
        "gene_id": 57674,
        "chrom": "chr17",
        "start": 80_255_479,
        "end": 80_371_516,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q25.3",
        "exon_count": 69,
        "role": (
            "Mysterin (RNF213) -- large AAA+ ATPase/E3 ubiquitin ligase.  "
            "R4810K founder variant is the strongest genetic risk factor for "
            "moyamoya disease (progressive steno-occlusive cerebral vasculopathy) "
            "in East Asian populations.  Carrier frequency ~2% in Japanese."
        ),
        "disease": "Moyamoya disease (susceptibility)",
        "omim_disease": 252350,
        "omim_gene": 613768,
        "inheritance": "complex (low penetrance risk factor)",
        "key_variants": [
            {
                "name": "p.Arg4810Lys (R4810K)",
                "rsid": "rs112735431",
                "consequence": "missense",
                "clinical_significance": "risk_factor",
                "notes": (
                    "East Asian founder; OR ~190 for moyamoya.  ~2% carrier "
                    "frequency in Japan but low penetrance (~1/150 carriers "
                    "develop disease)."
                ),
            },
        ],
        "strategy": (
            "(1) Base editing to correct R4810K in vascular endothelium; "
            "(2) CRISPRi to modulate RNF213 expression; "
            "(3) Gene is too large (~17 kb cDNA) for AAV.  "
            "Delivery: endothelial-tropic LNP or AAV."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "Surgical revascularization (STA-MCA bypass) is standard of care."
        ),
        "conditions": ["moyamoya", "vascular", "cerebrovascular", "stroke",
                        "steno_occlusive_vasculopathy"],
    },
}
