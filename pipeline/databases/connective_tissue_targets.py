"""
Connective tissue disease gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers osteogenesis imperfecta (OI), Ehlers-Danlos syndrome (EDS) subtypes,
and related heritable connective tissue disorders NOT covered by the cardiac
database (which handles Marfan/FBN1 and vEDS/COL3A1).

Categories:
  1. Osteogenesis Imperfecta (COL1A1, COL1A2)
  2. Classical EDS (COL5A1, COL5A2)
  3. Kyphoscoliotic EDS (PLOD1)
  4. Cutis Laxa (ELN, FBLN5)

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


CONNECTIVE_TISSUE_TARGETS = {

    # ===================================================================
    # OI -- COL1A1
    # ===================================================================
    "COL1A1": {
        "gene_id": 1277,
        "chrom": "chr17",
        "start": 50_183_289,
        "end": 50_201_632,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.33",
        "exon_count": 51,
        "role": (
            "Type I collagen alpha-1 chain -- most abundant protein in bone, "
            "skin, tendon, and dentin.  Forms heterotrimers with COL1A2.  "
            "Mutations cause osteogenesis imperfecta (OI, brittle bone disease).  "
            "Glycine substitutions in Gly-X-Y repeats cause dominant-negative "
            "structural defects (OI types II-IV); haploinsufficiency causes "
            "mild OI type I."
        ),
        "disease": "Osteogenesis imperfecta types I-IV",
        "omim_disease": 166200,
        "omim_gene": 120150,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Glycine substitutions (>800 known)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Gly->Ser most common.  C-terminal substitutions more "
                    "severe (lethal OI type II).  Position-dependent severity."
                ),
            },
            {
                "name": "Null alleles (nonsense/frameshift)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Haploinsufficiency -> mild OI type I (blue sclerae, few fractures).",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant-negative mutant allele "
            "using CRISPRi or siRNA (converts severe OI to mild OI type I); "
            "(2) Suppression-and-replacement: silence both alleles + add "
            "codon-optimized COL1A1 (cDNA ~4.4 kb, tight AAV fit); "
            "(3) Base editing for specific Gly substitutions; "
            "(4) Ex vivo CRISPR correction of patient mesenchymal stem cells "
            "for bone marrow transplant.  "
            "Delivery: AAV9 or AAVrh10; bone-targeted LNP in development."
        ),
        "clinical_programs": (
            "No gene therapy trials for COL1A1 OI as of 2026.  "
            "Ultragenyx UX143 (anti-sclerostin antibody, setrusumab) Phase 2/3.  "
            "Bisphosphonates standard of care."
        ),
        "conditions": ["osteogenesis_imperfecta", "OI", "brittle_bone_disease",
                        "connective_tissue", "bone_fragility"],
    },

    # ===================================================================
    # OI -- COL1A2
    # ===================================================================
    "COL1A2": {
        "gene_id": 1278,
        "chrom": "chr7",
        "start": 94_394_895,
        "end": 94_431_236,
        "strand": "-",
        "refseq": "NC_000007.14",
        "cytoband": "7q21.3",
        "exon_count": 52,
        "role": (
            "Type I collagen alpha-2 chain -- forms heterotrimers with COL1A1.  "
            "Glycine substitutions cause OI types II-IV.  Biallelic null = OI "
            "with hEDS-like joint hypermobility (all-alpha1 homotrimers)."
        ),
        "disease": "Osteogenesis imperfecta types II-IV",
        "omim_disease": 166210,
        "omim_gene": 120160,
        "inheritance": "AD (most); AR (rare)",
        "key_variants": [
            {
                "name": "Glycine substitutions (>400 known)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Same Gly-X-Y disruption mechanism as COL1A1.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant allele; "
            "(2) Base editing for specific Gly substitutions; "
            "(3) Same approaches as COL1A1."
        ),
        "clinical_programs": "Same as COL1A1 programs.",
        "conditions": ["osteogenesis_imperfecta", "OI", "connective_tissue",
                        "bone_fragility"],
    },

    # ===================================================================
    # Classical EDS -- COL5A1
    # ===================================================================
    "COL5A1": {
        "gene_id": 1289,
        "chrom": "chr9",
        "start": 134_641_803,
        "end": 134_845_160,
        "strand": "-",
        "refseq": "NC_000009.12",
        "cytoband": "9q34.3",
        "exon_count": 66,
        "role": (
            "Type V collagen alpha-1 chain -- nucleates type I collagen "
            "fibrillogenesis and regulates fibril diameter.  Haploinsufficiency "
            "causes classical EDS (cEDS): skin hyperextensibility, atrophic "
            "scarring, joint hypermobility, tissue fragility."
        ),
        "disease": "Classical Ehlers-Danlos syndrome (cEDS)",
        "omim_disease": 130000,
        "omim_gene": 120215,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Null alleles (haploinsufficiency)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "~50% of cEDS; NMD of mutant transcript.",
            },
            {
                "name": "Splice variants affecting exon 6 (c.1489G>A)",
                "rsid": None,
                "consequence": "splice_site",
                "clinical_significance": "pathogenic",
                "notes": "Exon skipping creates internally deleted collagen chain.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type allele in haploinsufficiency "
            "cases; "
            "(2) Gene too large for standard AAV (~5.4 kb cDNA); dual-AAV or "
            "lentiviral approach; "
            "(3) Allele-specific silencing for dominant-negative splice mutations."
        ),
        "clinical_programs": "No gene therapy trials as of 2026.",
        "conditions": ["classical_EDS", "cEDS", "ehlers_danlos",
                        "connective_tissue", "joint_hypermobility"],
    },

    # ===================================================================
    # Classical EDS -- COL5A2
    # ===================================================================
    "COL5A2": {
        "gene_id": 1290,
        "chrom": "chr2",
        "start": 189_040_027,
        "end": 189_202_679,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2q32.2",
        "exon_count": 52,
        "role": (
            "Type V collagen alpha-2 chain -- COL5A1/COL5A2 heterotrimer "
            "partner.  Mutations cause classical EDS (minority of cases)."
        ),
        "disease": "Classical EDS",
        "omim_disease": 130000,
        "omim_gene": 120190,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Glycine substitutions in triple helix",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Structural mutations; less common than COL5A1.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant-negative allele; "
            "(2) Base editing for specific Gly substitutions."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["classical_EDS", "cEDS", "connective_tissue"],
    },

    # ===================================================================
    # Kyphoscoliotic EDS -- PLOD1
    # ===================================================================
    "PLOD1": {
        "gene_id": 5351,
        "chrom": "chr1",
        "start": 11_988_097,
        "end": 12_026_657,
        "strand": "-",
        "refseq": "NC_000001.11",
        "cytoband": "1p36.22",
        "exon_count": 19,
        "role": (
            "Lysyl hydroxylase 1 (LH1) -- post-translational modification of "
            "collagen lysine residues for crosslinking.  Deficiency causes "
            "kyphoscoliotic EDS (kEDS): neonatal hypotonia, severe kyphoscoliosis, "
            "joint laxity, skin fragility, ocular fragility."
        ),
        "disease": "Kyphoscoliotic EDS (kEDS)",
        "omim_disease": 225400,
        "omim_gene": 153454,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Alu-Alu recombination duplication (exons 10-16)",
                "rsid": None,
                "consequence": "duplication",
                "clinical_significance": "pathogenic",
                "notes": "Most common in kEDS; founder mutation.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.2 kb fits AAV); "
            "(2) CRISPR excision of duplicated segment.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["kyphoscoliotic_EDS", "kEDS", "connective_tissue",
                        "ehlers_danlos", "kyphoscoliosis"],
    },

    # ===================================================================
    # Cutis Laxa -- ELN (Elastin)
    # ===================================================================
    "ELN": {
        "gene_id": 2006,
        "chrom": "chr7",
        "start": 74_028_172,
        "end": 74_069_905,
        "strand": "-",
        "refseq": "NC_000007.14",
        "cytoband": "7q11.23",
        "exon_count": 34,
        "role": (
            "Elastin -- structural protein of elastic fibers in skin, lungs, "
            "and blood vessels.  Dominant mutations cause cutis laxa type 1 "
            "(loose, inelastic skin) or supravalvular aortic stenosis (SVAS).  "
            "ELN haploinsufficiency causes SVAS; C-terminal frameshift "
            "mutations cause autosomal dominant cutis laxa."
        ),
        "disease": "Cutis laxa AD / Supravalvular aortic stenosis",
        "omim_disease": 123700,
        "omim_gene": 130160,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "C-terminal frameshift mutations",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Dominant-negative; misassembly of elastic fibers.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant-negative allele; "
            "(2) CRISPRa to upregulate wild-type allele for SVAS; "
            "(3) AAV gene replacement (cDNA ~2.3 kb fits AAV) for SVAS.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["cutis_laxa", "SVAS", "connective_tissue",
                        "elastic_fiber_disorder"],
    },

    # ===================================================================
    # Cutis Laxa AR -- FBLN5
    # ===================================================================
    "FBLN5": {
        "gene_id": 10516,
        "chrom": "chr14",
        "start": 91_869_917,
        "end": 91_929_825,
        "strand": "-",
        "refseq": "NC_000014.9",
        "cytoband": "14q32.12",
        "exon_count": 11,
        "role": (
            "Fibulin-5 -- extracellular matrix protein essential for elastic "
            "fiber assembly.  Recessive mutations cause autosomal recessive "
            "cutis laxa type 1A (ARCL1A): severe loose skin, pulmonary "
            "emphysema, vascular tortuosity."
        ),
        "disease": "Autosomal recessive cutis laxa type 1A",
        "omim_disease": 219100,
        "omim_gene": 604580,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Cys217Arg (C217R)",
                "rsid": "rs104893694",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "cbEGF domain; disrupts calcium binding.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.2 kb, ideal for AAV); "
            "(2) Base editing for specific missense mutations."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["cutis_laxa", "ARCL1A", "connective_tissue",
                        "elastic_fiber_disorder"],
    },
}
