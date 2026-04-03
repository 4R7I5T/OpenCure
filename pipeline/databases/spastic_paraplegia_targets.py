"""
Hereditary spastic paraplegia (HSP) gene-therapy targets -- GRCh38 (hg38).

HSPs are a group of >80 genetic disorders characterised by progressive
spasticity and weakness of the lower limbs due to corticospinal tract
axonal degeneration.  SPG4 (SPAST) is the most common (40% of AD-HSP).

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


SPASTIC_PARAPLEGIA_TARGETS = {

    "SPAST": {
        "gene_id": 6683,
        "chrom": "chr2",
        "start": 32_067_889,
        "end": 32_156_655,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2p22.3",
        "exon_count": 17,
        "role": (
            "Spastin -- AAA ATPase that severs microtubules, essential for "
            "microtubule dynamics in axons.  Haploinsufficiency causes SPG4, "
            "the most common hereditary spastic paraplegia (~40% of AD-HSP, "
            "~1:20,000 prevalence).  Longest axons (corticospinal tract) "
            "are most vulnerable."
        ),
        "disease": "Hereditary spastic paraplegia type 4 (SPG4)",
        "omim_disease": 182601,
        "omim_gene": 604277,
        "inheritance": "AD (haploinsufficiency + some dominant-negative)",
        "key_variants": [
            {
                "name": ">500 unique mutations (extreme heterogeneity)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~40% missense, ~20% splice, ~20% frameshift/nonsense, "
                    "~10% large deletions.  Variable age of onset (1-70 years) "
                    "even within families."
                ),
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type SPAST allele (double "
            "normal spastin dosage is safe in animal models); "
            "(2) AAV gene replacement (cDNA ~1.8 kb, ideal for AAV) to "
            "corticospinal neurons; "
            "(3) Microtubule-stabilizing drugs (vinblastine at low dose) "
            "showed benefit in Drosophila and mouse models.  "
            "Delivery: AAV9 intrathecal or IV."
        ),
        "clinical_programs": (
            "No gene therapy trials as of 2026.  "
            "Noscapine (microtubule-modulating) in preclinical.  "
            "Rehabilitation and antispasticity agents (baclofen, tizanidine) "
            "remain standard of care."
        ),
        "conditions": ["SPG4", "hereditary_spastic_paraplegia", "HSP",
                        "spastic_paraplegia", "corticospinal_degeneration"],
    },

    "ATL1": {
        "gene_id": 51062,
        "chrom": "chr14",
        "start": 50_714_996,
        "end": 50_769_485,
        "strand": "+",
        "refseq": "NC_000014.9",
        "cytoband": "14q22.1",
        "exon_count": 14,
        "role": (
            "Atlastin GTPase 1 -- dynamin-related ER-shaping GTPase, "
            "mediates homotypic ER tubule fusion.  Mutations cause SPG3A, "
            "the most common early-onset AD-HSP (~10% of all AD-HSP).  "
            "Onset typically <10 years, pure spastic paraplegia, slowly "
            "progressive."
        ),
        "disease": "Hereditary spastic paraplegia type 3A (SPG3A)",
        "omim_disease": 182600,
        "omim_gene": 606439,
        "inheritance": "AD (mostly de novo in early-onset)",
        "key_variants": [
            {
                "name": "p.Arg239Cys (R239C)",
                "rsid": "rs121917816",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "GTPase domain; most common SPG3A mutation.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant mutant allele; "
            "(2) AAV gene replacement (cDNA ~1.6 kb, ideal for AAV); "
            "(3) Base editing for R239C."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["SPG3A", "hereditary_spastic_paraplegia", "HSP",
                        "early_onset_spasticity"],
    },

    "SPG7": {
        "gene_id": 6687,
        "chrom": "chr16",
        "start": 89_508_423,
        "end": 89_557_709,
        "strand": "+",
        "refseq": "NC_000016.10",
        "cytoband": "16q24.3",
        "exon_count": 17,
        "role": (
            "Paraplegin -- mitochondrial metalloprotease (m-AAA protease "
            "subunit) in the inner mitochondrial membrane.  Biallelic "
            "mutations cause SPG7: spastic paraplegia + cerebellar ataxia "
            "(spastic ataxia).  Most common AR-HSP in many European cohorts."
        ),
        "disease": "Hereditary spastic paraplegia type 7 (SPG7, spastic ataxia)",
        "omim_disease": 607259,
        "omim_gene": 602783,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Ala510Val (A510V)",
                "rsid": "rs61755348",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common SPG7 mutation in Europeans (~30% of alleles).  "
                    "Hypomorphic; residual protease activity."
                ),
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.4 kb, fits AAV) with "
            "mitochondrial targeting; "
            "(2) Base editing for A510V."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["SPG7", "hereditary_spastic_paraplegia", "HSP",
                        "spastic_ataxia", "mitochondrial"],
    },

    "SPG11": {
        "gene_id": 80208,
        "chrom": "chr15",
        "start": 44_564_766,
        "end": 44_642_828,
        "strand": "+",
        "refseq": "NC_000015.10",
        "cytoband": "15q21.1",
        "exon_count": 40,
        "role": (
            "Spatacsin -- large protein involved in autophagosome-lysosome "
            "reformation and axonal cargo transport.  Mutations cause SPG11 "
            "(most common AR-HSP overall, ~20% of AR-HSP): complicated HSP "
            "with thin corpus callosum, cognitive decline, parkinsonism, "
            "peripheral neuropathy."
        ),
        "disease": "SPG11 (complicated HSP with thin corpus callosum)",
        "omim_disease": 604360,
        "omim_gene": 610844,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.733_734delAT (most common)",
                "rsid": "rs312262729",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Recurrent frameshift; NMD of transcript.",
            },
        ],
        "strategy": (
            "(1) Gene large for AAV (~7.4 kb cDNA); dual-AAV split-intein "
            "approach; "
            "(2) CRISPRa to upregulate residual expression from hypomorphic "
            "alleles."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["SPG11", "hereditary_spastic_paraplegia", "HSP",
                        "thin_corpus_callosum", "complicated_HSP"],
    },

    "KIF5A": {
        "gene_id": 3798,
        "chrom": "chr12",
        "start": 57_481_481,
        "end": 57_525_760,
        "strand": "+",
        "refseq": "NC_000012.12",
        "cytoband": "12q13.3",
        "exon_count": 28,
        "role": (
            "Kinesin family member 5A -- neuronal kinesin-1 heavy chain, "
            "anterograde axonal transport motor.  Motor domain mutations "
            "cause SPG10 (pure or complicated HSP).  C-terminal mutations "
            "cause neonatal intractable myoclonus (NEIMY) and ALS."
        ),
        "disease": "SPG10 / Neonatal intractable myoclonus",
        "omim_disease": 604187,
        "omim_gene": 602821,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "p.Asn256Ser (N256S)",
                "rsid": "rs121918350",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Motor domain; impairs ATPase cycle; SPG10.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant allele; "
            "(2) AAV gene replacement (cDNA ~3.0 kb, fits AAV).  "
            "Delivery: AAV9 intrathecal."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["SPG10", "hereditary_spastic_paraplegia", "HSP",
                        "axonal_transport_defect"],
    },

    "REEP1": {
        "gene_id": 65055,
        "chrom": "chr2",
        "start": 86_186_765,
        "end": 86_249_161,
        "strand": "-",
        "refseq": "NC_000002.12",
        "cytoband": "2p11.2",
        "exon_count": 7,
        "role": (
            "Receptor expression-enhancing protein 1 -- ER-shaping protein "
            "that generates ER tubule curvature (works with ATL1).  "
            "Haploinsufficiency causes SPG31 (~5% of AD-HSP): pure spastic "
            "paraplegia, childhood to adult onset."
        ),
        "disease": "SPG31",
        "omim_disease": 610250,
        "omim_gene": 609139,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Various LOF mutations",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Haploinsufficiency; deletions, nonsense, splice.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type REEP1 allele; "
            "(2) AAV gene replacement (cDNA ~0.6 kb, ideal for AAV)."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["SPG31", "hereditary_spastic_paraplegia", "HSP"],
    },
}
