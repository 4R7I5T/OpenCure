"""
Peroxisomal disorder gene-therapy targets -- GRCh38 (hg38).

Covers Zellweger spectrum disorders (ZSD) and rhizomelic chondrodysplasia
punctata (RCDP).  X-linked adrenoleukodystrophy (ABCD1) is covered
separately in the leukodystrophy pipeline.

Peroxisomes are essential for VLCFA beta-oxidation, plasmalogen synthesis,
bile acid synthesis, and reactive oxygen species metabolism.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


PEROXISOMAL_TARGETS = {

    "PEX1": {
        "gene_id": 5189,
        "chrom": "chr7",
        "start": 92_594_461,
        "end": 92_622_445,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q21.2",
        "exon_count": 24,
        "role": (
            "Peroxin 1 -- AAA ATPase required for peroxisomal matrix protein "
            "import (PEX1-PEX6 complex recycles PEX5 receptor).  Most common "
            "peroxisome biogenesis disorder gene (~68% of ZSD).  Zellweger "
            "spectrum: neonatal Zellweger syndrome (severe, lethal), neonatal "
            "adrenoleukodystrophy (intermediate), and infantile Refsum disease "
            "(mildest ZSD form)."
        ),
        "disease": "Zellweger spectrum disorder (ZSD, most common gene)",
        "omim_disease": 214100,
        "omim_gene": 602136,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Gly843Asp (G843D)",
                "rsid": "rs61750420",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common PEX1 allele worldwide (~30% of PEX1 alleles).  "
                    "Temperature-sensitive: some peroxisomal function at 30C, "
                    "lost at 37C.  Milder phenotype than null alleles."
                ),
            },
            {
                "name": "c.2528G>A (p.Gly843Asp) + null = intermediate",
                "rsid": None,
                "consequence": "compound_heterozygous",
                "clinical_significance": "pathogenic",
                "notes": "G843D/null = intermediate (NALD/IRD phenotype).",
            },
        ],
        "strategy": (
            "(1) AAV liver-directed gene replacement (cDNA ~3.9 kb fits AAV "
            "with compact promoter) -- restoring hepatic peroxisomes corrects "
            "bile acid and VLCFA abnormalities; "
            "(2) Pharmacological chaperones for temperature-sensitive G843D "
            "(betaine/TMAO showed benefit in fibroblasts); "
            "(3) Base editing for G843D (single nucleotide change); "
            "(4) Gene therapy must be delivered early (before irreversible "
            "CNS damage in severe ZSD).  "
            "Delivery: AAV8 IV to liver + AAV9 for CNS."
        ),
        "clinical_programs": (
            "No gene therapy trials as of 2026.  "
            "Cholic acid (Cholbam) FDA-approved for bile acid synthesis "
            "defects in ZSD (supportive, not curative).  "
            "Pharmacological chaperone trials for G843D planned."
        ),
        "conditions": ["zellweger_spectrum", "ZSD", "peroxisomal",
                        "peroxisome_biogenesis_disorder", "VLCFA"],
    },

    "PEX6": {
        "gene_id": 5190,
        "chrom": "chr6",
        "start": 42_928_200,
        "end": 42_948_288,
        "strand": "+",
        "refseq": "NC_000006.12",
        "cytoband": "6p21.1",
        "exon_count": 17,
        "role": (
            "Peroxin 6 -- AAA ATPase partner of PEX1 in the PEX5 receptor "
            "recycling complex.  Second most common ZSD gene (~16% of ZSD).  "
            "Same clinical spectrum as PEX1."
        ),
        "disease": "Zellweger spectrum disorder",
        "omim_disease": 614862,
        "omim_gene": 601498,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various missense and truncating",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "No single common mutation; allelic heterogeneity.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.9 kb fits AAV); "
            "(2) Same approach as PEX1."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["zellweger_spectrum", "ZSD", "peroxisomal"],
    },

    "PEX7": {
        "gene_id": 5191,
        "chrom": "chr6",
        "start": 137_143_991,
        "end": 137_233_052,
        "strand": "-",
        "refseq": "NC_000006.12",
        "cytoband": "6q23.3",
        "exon_count": 10,
        "role": (
            "Peroxin 7 (PTS2 receptor) -- imports PTS2-containing proteins "
            "into peroxisomes.  Mutations cause rhizomelic chondrodysplasia "
            "punctata type 1 (RCDP1): rhizomelic limb shortening, stippled "
            "epiphyses, cataracts, severe ID, failure to thrive.  RCDP1 is "
            "the most common form of RCDP.  Plasmalogen synthesis is "
            "specifically impaired."
        ),
        "disease": "Rhizomelic chondrodysplasia punctata type 1 (RCDP1)",
        "omim_disease": 215100,
        "omim_gene": 601757,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Leu292Ter (L292X)",
                "rsid": "rs121434498",
                "consequence": "nonsense",
                "clinical_significance": "pathogenic",
                "notes": "Common in European populations; null allele.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.0 kb, ideal for AAV); "
            "(2) Plasmalogen replacement therapy (oral alkyl-glycerol "
            "supplementation to bypass the block).  "
            "Delivery: AAV8 IV to liver + systemic."
        ),
        "clinical_programs": (
            "Plasmalogen replacement (PPI-1011) preclinical by Med-Life "
            "Discoveries.  No gene therapy trials."
        ),
        "conditions": ["RCDP", "RCDP1", "peroxisomal",
                        "chondrodysplasia_punctata", "plasmalogen_deficiency"],
    },

    "PEX12": {
        "gene_id": 5193,
        "chrom": "chr17",
        "start": 35_654_371,
        "end": 35_659_505,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q12",
        "exon_count": 3,
        "role": (
            "Peroxin 12 -- E3 ubiquitin ligase component of the peroxisomal "
            "importomer RING finger complex (PEX2-PEX10-PEX12).  Mutations "
            "cause ZSD.  Third most common ZSD gene (~5%)."
        ),
        "disease": "Zellweger spectrum disorder",
        "omim_disease": 614859,
        "omim_gene": 601758,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Ser320Phe (S320F)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "RING domain; impairs PEX5 ubiquitination.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.1 kb, ideal for AAV); "
            "(2) Base editing for specific mutations."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["zellweger_spectrum", "ZSD", "peroxisomal"],
    },

    "ACOX1": {
        "gene_id": 51,
        "chrom": "chr17",
        "start": 75_938_422,
        "end": 75_972_798,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q25.1",
        "exon_count": 14,
        "role": (
            "Acyl-CoA oxidase 1 -- first enzyme of peroxisomal VLCFA "
            "beta-oxidation (straight-chain fatty acids).  Deficiency causes "
            "pseudo-neonatal adrenoleukodystrophy (ACOX1 deficiency): "
            "progressive leukodystrophy, seizures, hearing loss.  Single "
            "peroxisomal enzyme deficiency (peroxisomes are present but "
            "VLCFA accumulates)."
        ),
        "disease": "ACOX1 deficiency (pseudo-NALD)",
        "omim_disease": 264470,
        "omim_gene": 609751,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various missense and splicing",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Rare; <50 cases reported.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.0 kb, fits AAV) to liver "
            "and CNS; "
            "(2) Lorenzo's oil (dietary VLCFA reduction) as adjunct."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["ACOX1_deficiency", "peroxisomal",
                        "leukodystrophy", "VLCFA"],
    },
}
