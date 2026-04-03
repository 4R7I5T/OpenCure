"""
Bone dysplasia / skeletal disorder gene-therapy targets -- GRCh38 (hg38).

Expands the skeletal pipeline (which covers FGFR3/achondroplasia and
ALPL/hypophosphatasia) with additional bone dysplasias amenable to
gene therapy or targeted molecular therapy.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


BONE_DYSPLASIA_TARGETS = {

    "ACVR1": {
        "gene_id": 90,
        "chrom": "chr2",
        "start": 157_736_092,
        "end": 157_874_082,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2q24.1",
        "exon_count": 11,
        "role": (
            "Activin A receptor type 1 (ALK2) -- BMP type I receptor.  "
            "Gain-of-function mutation R206H causes fibrodysplasia ossificans "
            "progressiva (FOP): progressive heterotopic ossification where "
            "muscle, tendon, and ligament turn to bone after injury/inflammation.  "
            "Most disabling inherited condition of connective tissue.  "
            "~1:1,000,000 prevalence; almost all de novo."
        ),
        "disease": "Fibrodysplasia ossificans progressiva (FOP)",
        "omim_disease": 135100,
        "omim_gene": 102576,
        "inheritance": "AD (almost all de novo)",
        "key_variants": [
            {
                "name": "p.Arg206His (R206H)",
                "rsid": "rs121912678",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "SINGLE recurrent mutation in ~97% of FOP patients.  "
                    "GS domain; constitutive BMP signaling.  Same nucleotide "
                    "change (c.617G>A) in virtually every patient."
                ),
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of R206H mutant ACVR1 allele "
            "(preserving WT essential for normal BMP signaling); "
            "(2) Base editing to correct R206H (single nucleotide, "
            "A>G reversion); "
            "(3) Allele-specific ASO (Regeneron/Ipsen approach); "
            "(4) Palovarotene (RARγ agonist) reduces HO; "
            "(5) Anti-activin A antibody (garetosmab) blocks the trigger.  "
            "Delivery: systemic ASO or AAV."
        ),
        "clinical_programs": (
            "Ipsen palovarotene (Sohonos) FDA-approved 2023 for FOP "
            "(first approved FOP therapy, reduces new HO ~50%).  "
            "Regeneron REGN2477 (anti-activin A, garetosmab) Phase 2 "
            "(NCT03188666).  "
            "Ipsen IPN60130 (anti-ACVR1 ASO) preclinical.  "
            "No CRISPR trials."
        ),
        "conditions": ["fibrodysplasia_ossificans_progressiva", "FOP",
                        "heterotopic_ossification", "bone_dysplasia",
                        "BMP_signaling"],
    },

    "PHEX": {
        "gene_id": 5251,
        "chrom": "chrX",
        "start": 22_032_283,
        "end": 22_251_338,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xp22.11",
        "exon_count": 22,
        "role": (
            "Phosphate-regulating endopeptidase homolog, X-linked -- membrane "
            "metalloprotease expressed in osteocytes/odontoblasts.  Deficiency "
            "causes X-linked hypophosphatemia (XLH): excess FGF23 -> renal "
            "phosphate wasting -> rickets, osteomalacia, short stature, "
            "dental abscesses.  Most common inherited rickets (~1:20,000)."
        ),
        "disease": "X-linked hypophosphatemia (XLH)",
        "omim_disease": 307800,
        "omim_gene": 300550,
        "inheritance": "XL dominant",
        "key_variants": [
            {
                "name": ">400 unique mutations (extreme heterogeneity)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "No hotspot; private mutations predominate.",
            },
        ],
        "strategy": (
            "(1) Burosumab (anti-FGF23 antibody) is FDA-approved and highly "
            "effective, reducing need for gene therapy; "
            "(2) AAV bone-directed gene replacement (cDNA ~2.3 kb fits AAV) "
            "for potential one-time cure vs lifelong burosumab; "
            "(3) CRISPRi to reduce FGF23 as downstream correction.  "
            "Delivery: bone-tropic AAV."
        ),
        "clinical_programs": (
            "Burosumab (Crysvita, Ultragenyx) FDA-approved 2018 for XLH.  "
            "No gene therapy trials (burosumab is effective)."
        ),
        "conditions": ["X_linked_hypophosphatemia", "XLH", "bone_dysplasia",
                        "rickets", "phosphate_wasting"],
    },

    "EXT1": {
        "gene_id": 2131,
        "chrom": "chr8",
        "start": 118_806_745,
        "end": 119_124_076,
        "strand": "-",
        "refseq": "NC_000008.11",
        "cytoband": "8q24.11",
        "exon_count": 11,
        "role": (
            "Exostosin glycosyltransferase 1 -- ER-resident glycosyltransferase "
            "for heparan sulfate (HS) chain elongation.  Haploinsufficiency "
            "causes hereditary multiple exostoses (HME/MO): multiple "
            "osteochondromas (cartilage-capped bony growths) at metaphyses "
            "of long bones, limb length discrepancy, risk of malignant "
            "transformation to chondrosarcoma (~1-5%).  ~1:50,000."
        ),
        "disease": "Hereditary multiple exostoses (HME type 1)",
        "omim_disease": 133700,
        "omim_gene": 608177,
        "inheritance": "AD (haploinsufficiency + LOH in tumors)",
        "key_variants": [
            {
                "name": "Various LOF mutations (>400 in EXT1+EXT2)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "EXT1 mutations ~60-70% of HME; more severe than EXT2.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type EXT1 allele in growth plates; "
            "(2) AAV gene replacement (cDNA ~2.2 kb fits AAV); "
            "(3) Surgical excision of symptomatic osteochondromas remains SOC.  "
            "Delivery: bone-tropic AAV (challenging)."
        ),
        "clinical_programs": (
            "Palovarotene (RARγ agonist) Phase 2 for HME (Clementia/Ipsen) "
            "-- showed reduction in new osteochondromas in preclinical.  "
            "No gene therapy trials."
        ),
        "conditions": ["hereditary_multiple_exostoses", "HME",
                        "osteochondroma", "bone_dysplasia"],
    },

    "EXT2": {
        "gene_id": 2132,
        "chrom": "chr11",
        "start": 44_088_760,
        "end": 44_241_637,
        "strand": "+",
        "refseq": "NC_000011.10",
        "cytoband": "11p11.2",
        "exon_count": 14,
        "role": (
            "Exostosin glycosyltransferase 2 -- EXT1/EXT2 hetero-oligomer "
            "partner for HS biosynthesis.  Mutations cause HME type 2 "
            "(generally milder than EXT1)."
        ),
        "disease": "HME type 2",
        "omim_disease": 133701,
        "omim_gene": 608210,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Various LOF mutations",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "~30-40% of HME; generally milder than EXT1.",
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type allele; "
            "(2) AAV gene replacement (cDNA ~2.1 kb fits AAV)."
        ),
        "clinical_programs": "Same as EXT1.",
        "conditions": ["hereditary_multiple_exostoses", "HME", "bone_dysplasia"],
    },

    "COL2A1": {
        "gene_id": 1280,
        "chrom": "chr12",
        "start": 47_972_964,
        "end": 48_004_507,
        "strand": "-",
        "refseq": "NC_000012.12",
        "cytoband": "12q13.11",
        "exon_count": 54,
        "role": (
            "Type II collagen alpha-1 chain -- major structural protein of "
            "hyaline cartilage, vitreous humor, and intervertebral discs.  "
            "Mutations cause a spectrum of skeletal dysplasias: "
            "achondrogenesis type II (lethal), hypochondrogenesis, "
            "spondyloepiphyseal dysplasia congenita (SED), Kniest dysplasia, "
            "Stickler syndrome type I, and early-onset osteoarthritis."
        ),
        "disease": "Type II collagenopathies (SED, Stickler, Kniest)",
        "omim_disease": 183900,
        "omim_gene": 120140,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Glycine substitutions in Gly-X-Y repeat",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Same mechanism as COL1A1/2 (OI) -- Gly substitutions "
                    "disrupt triple helix.  C-terminal more severe."
                ),
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of dominant-negative mutant allele; "
            "(2) Gene too large for standard AAV (~4.4 kb cDNA) but tight fit; "
            "(3) Base editing for specific Gly substitutions."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["spondyloepiphyseal_dysplasia", "SED", "stickler_syndrome",
                        "kniest_dysplasia", "bone_dysplasia", "collagenopathy"],
    },
}
