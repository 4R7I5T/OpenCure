"""
Mucopolysaccharidosis (MPS) and gangliosidosis gene-therapy targets -- GRCh38.

Dedicated pipeline for lysosomal storage diseases with active gene therapy
clinical programs.  These complement the existing rare_disease pipeline
(which covers IDUA/MPS I, IDS/MPS II, HEXA, HEXB, SMPD1, NPC1) with
MPS/gangliosidosis genes NOT in rare_disease.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - Abeona, Lysogene, Ultragenyx, Passage Bio LSD gene therapy programs
"""


MUCOPOLYSACCHARIDOSIS_TARGETS = {

    "SGSH": {
        "gene_id": 6448,
        "chrom": "chr17",
        "start": 80_208_556,
        "end": 80_219_977,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q25.3",
        "exon_count": 8,
        "role": (
            "N-sulfoglucosamine sulfohydrolase (sulfamidase) -- lysosomal "
            "enzyme cleaving sulfamide bonds in heparan sulfate (HS).  "
            "Deficiency causes MPS IIIA (Sanfilippo A): severe progressive "
            "neurodegeneration with behavioral disturbance, speech loss, "
            "sleep dysfunction, seizures.  Death typically by late teens.  "
            "MPS IIIA is the most common Sanfilippo subtype."
        ),
        "disease": "MPS IIIA (Sanfilippo syndrome type A)",
        "omim_disease": 252900,
        "omim_gene": 605270,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg245His (R245H)",
                "rsid": "rs104894639",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common in Dutch/North European populations.",
            },
            {
                "name": "p.Arg74Cys (R74C)",
                "rsid": "rs104894640",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Polish founder; attenuated phenotype.",
            },
        ],
        "strategy": (
            "(1) AAV-mediated CNS gene replacement (cDNA ~1.5 kb, ideal "
            "for AAV) -- intracerebral or intrathecal delivery; "
            "(2) Lysogene LYS-SAF302 (AAVrh10-SGSH intracerebral) Phase 2/3; "
            "(3) Abeona ABO-102 (AAV9-SGSH IV) Phase 1/2; "
            "(4) Ex vivo lentiviral HSC gene therapy for cross-correction.  "
            "Delivery: AAVrh10 intracerebral or AAV9 IV."
        ),
        "clinical_programs": (
            "Lysogene LYS-SAF302 (AAVrh10-SGSH intracerebral) Phase 2/3 "
            "(NCT03612869).  "
            "Abeona ABO-102 (AAV9-SGSH IV) Phase 1/2 (NCT04088734).  "
            "Orchard Therapeutics (ex vivo lentiviral HSC) Phase 1/2 "
            "(NCT04201405).  MOST ADVANCED MPS III program."
        ),
        "conditions": ["MPS_IIIA", "sanfilippo_A", "mucopolysaccharidosis",
                        "lysosomal_storage", "neurodegeneration"],
    },

    "NAGLU": {
        "gene_id": 4669,
        "chrom": "chr17",
        "start": 42_094_736,
        "end": 42_102_948,
        "strand": "-",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.2",
        "exon_count": 6,
        "role": (
            "N-alpha-acetylglucosaminidase -- lysosomal enzyme cleaving "
            "terminal N-acetylglucosamine from heparan sulfate.  Deficiency "
            "causes MPS IIIB (Sanfilippo B): clinically indistinguishable "
            "from MPS IIIA."
        ),
        "disease": "MPS IIIB (Sanfilippo syndrome type B)",
        "omim_disease": 252920,
        "omim_gene": 609701,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various (>150 mutations, allelic heterogeneity)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "No common mutation; private mutations predominate.",
            },
        ],
        "strategy": (
            "(1) AAV CNS gene replacement (cDNA ~2.2 kb fits AAV); "
            "(2) Abeona ABO-101 (AAVrh10-NAGLU) Phase 1/2; "
            "(3) Intrathecal AAV delivery.  "
            "Delivery: AAVrh10 intracerebral or AAV9 IV."
        ),
        "clinical_programs": (
            "Abeona ABO-101 (AAVrh10-NAGLU) Phase 1/2 (NCT03315182).  "
            "Esteve (formerly Laboratorios del Dr. Esteve) EGT-101 "
            "(AAV9-NAGLU) preclinical."
        ),
        "conditions": ["MPS_IIIB", "sanfilippo_B", "mucopolysaccharidosis",
                        "lysosomal_storage", "neurodegeneration"],
    },

    "GALNS": {
        "gene_id": 2588,
        "chrom": "chr16",
        "start": 88_813_917,
        "end": 88_864_385,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16q24.3",
        "exon_count": 14,
        "role": (
            "N-acetylgalactosamine-6-sulfatase -- lysosomal enzyme cleaving "
            "6-sulfate from keratan sulfate and chondroitin-6-sulfate.  "
            "Deficiency causes MPS IVA (Morquio A): severe skeletal dysplasia, "
            "short stature, odontoid hypoplasia (cervical instability), "
            "corneal clouding, cardiac valve disease.  Normal intelligence."
        ),
        "disease": "MPS IVA (Morquio A syndrome)",
        "omim_disease": 253000,
        "omim_gene": 612222,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg386Cys (R386C)",
                "rsid": "rs121918170",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "British/Irish founder; attenuated phenotype.",
            },
        ],
        "strategy": (
            "(1) AAV liver-directed gene replacement (cDNA ~1.6 kb, ideal "
            "for AAV) for systemic enzyme delivery via secretion; "
            "(2) AAV to bone/cartilage (challenging delivery); "
            "(3) ERT (elosulfase alfa/Vimizim) FDA-approved 2014.  "
            "Delivery: AAV8 IV to liver for cross-correction."
        ),
        "clinical_programs": (
            "ERT: Elosulfase alfa (Vimizim, BioMarin) FDA-approved 2014.  "
            "No gene therapy trials as of 2026."
        ),
        "conditions": ["MPS_IVA", "morquio_A", "mucopolysaccharidosis",
                        "lysosomal_storage", "skeletal_dysplasia"],
    },

    "GLB1": {
        "gene_id": 2720,
        "chrom": "chr3",
        "start": 33_014_066,
        "end": 33_107_385,
        "strand": "-",
        "refseq": "NC_000003.12",
        "cytoband": "3p22.3",
        "exon_count": 16,
        "role": (
            "Beta-galactosidase -- lysosomal enzyme cleaving terminal "
            "galactose from GM1 ganglioside, glycoproteins, and keratan "
            "sulfate.  Deficiency causes GM1 gangliosidosis: progressive "
            "neurodegeneration, hepatosplenomegaly, skeletal dysplasia.  "
            "Infantile form (type I) is most severe (death by age 3-4).  "
            "Also causes MPS IVB (Morquio B, skeletal predominant)."
        ),
        "disease": "GM1 gangliosidosis / MPS IVB (Morquio B)",
        "omim_disease": 230500,
        "omim_gene": 611458,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg208Cys (R208C)",
                "rsid": "rs121907958",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Infantile onset; severe neurodegeneration.",
            },
            {
                "name": "p.Ile51Thr (I51T)",
                "rsid": "rs121907959",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Adult/chronic GM1; milder, some residual activity.",
            },
        ],
        "strategy": (
            "(1) AAV CNS gene replacement (cDNA ~2.0 kb, fits AAV) -- "
            "Passage Bio PBGM01 approach; "
            "(2) IV AAV9 for systemic (CNS + visceral) correction; "
            "(3) Miglustat (substrate reduction) as adjunct.  "
            "Delivery: AAV9 IV or AAVhu68 intracisternal."
        ),
        "clinical_programs": (
            "Passage Bio PBGM01 (AAVhu68-GLB1 intracisternal) Phase 1/2 "
            "(NCT04713475) for infantile GM1.  "
            "Sio Gene Therapies AXO-AAV-GM1 (AAV9-GLB1 IV) Phase 1/2 "
            "(NCT03952637).  "
            "MOST ADVANCED GM1 programs."
        ),
        "conditions": ["GM1_gangliosidosis", "MPS_IVB", "morquio_B",
                        "mucopolysaccharidosis", "lysosomal_storage",
                        "gangliosidosis", "neurodegeneration"],
    },

    "ARSB": {
        "gene_id": 411,
        "chrom": "chr5",
        "start": 78_073_054,
        "end": 78_281_329,
        "strand": "+",
        "refseq": "NC_000005.10",
        "cytoband": "5q14.1",
        "exon_count": 8,
        "role": (
            "Arylsulfatase B (N-acetylgalactosamine-4-sulfatase) -- lysosomal "
            "enzyme degrading dermatan sulfate and chondroitin-4-sulfate.  "
            "Deficiency causes MPS VI (Maroteaux-Lamy): skeletal dysplasia, "
            "corneal clouding, cardiac valve disease, normal intelligence "
            "(unlike MPS I/II)."
        ),
        "disease": "MPS VI (Maroteaux-Lamy syndrome)",
        "omim_disease": 253200,
        "omim_gene": 611542,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Various (>180 mutations, high allelic heterogeneity)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "No common mutation globally; some population-specific founders.",
            },
        ],
        "strategy": (
            "(1) AAV liver-directed gene replacement (cDNA ~1.6 kb, ideal "
            "for AAV) for secretion-mediated cross-correction; "
            "(2) ERT (galsulfase/Naglazyme) FDA-approved 2005; "
            "(3) Ex vivo HSC gene therapy.  "
            "Delivery: AAV8 IV to liver."
        ),
        "clinical_programs": (
            "ERT: Galsulfase (Naglazyme, BioMarin) FDA-approved 2005.  "
            "No gene therapy trials as of 2026."
        ),
        "conditions": ["MPS_VI", "maroteaux_lamy", "mucopolysaccharidosis",
                        "lysosomal_storage", "skeletal_dysplasia"],
    },
}
