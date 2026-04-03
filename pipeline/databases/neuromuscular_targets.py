"""
Neuromuscular disease gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers limb-girdle muscular dystrophies (LGMD), congenital myopathies,
myotubular myopathy, and congenital myasthenic syndromes NOT already
covered by the hereditary disease database (which covers DMD, SMA, DM1).

Categories:
  1. LGMD Type R (recessive): DYSF, SGCA, SGCB, SGCG, CAPN3, FKRP, ANO5
  2. Myotubular Myopathy: MTM1
  3. Congenital Myasthenic Syndromes: RAPSN, DOK7, COLQ

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.

Sources:
  - NCBI Gene (GRCh38.p14), ClinVar, OMIM
  - Sarepta, Solid Biosciences, Atamyo gene therapy programs
  - ClinicalTrials.gov
  - Mendell lab (Nationwide Children's) LGMD gene therapy data
"""


NEUROMUSCULAR_TARGETS = {

    # ===================================================================
    # LGMD R2 (DYSF) -- Dysferlinopathy
    # ===================================================================
    "DYSF": {
        "gene_id": 8291,
        "chrom": "chr2",
        "start": 71_453_570,
        "end": 71_686_762,
        "strand": "+",
        "refseq": "NC_000002.12",
        "cytoband": "2p13.2",
        "exon_count": 55,
        "role": (
            "Dysferlin -- membrane repair protein in skeletal muscle.  "
            "Deficiency causes LGMD R2 (formerly LGMD2B) and Miyoshi myopathy.  "
            "Patients lose ability to reseal sarcolemmal tears."
        ),
        "disease": "LGMD R2 / Miyoshi myopathy",
        "omim_disease": 253601,
        "omim_gene": 603009,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "Truncating variants throughout",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": ">500 pathogenic variants; no hotspot.",
            },
        ],
        "strategy": (
            "(1) Dual-AAV gene replacement (cDNA ~6.2 kb exceeds AAV) using "
            "split-intein or overlapping approach; "
            "(2) Exon skipping for specific truncating variants.  "
            "Delivery: AAVrh74 or AAV9 IV/intramuscular."
        ),
        "clinical_programs": (
            "No clinical trials as of 2026.  "
            "Dual-AAV DYSF replacement showed efficacy in dysferlin-null mice."
        ),
        "conditions": ["LGMD_R2", "dysferlinopathy", "miyoshi_myopathy",
                        "neuromuscular", "muscular_dystrophy"],
    },

    # ===================================================================
    # LGMD R3 (SGCA) -- Alpha-sarcoglycanopathy
    # ===================================================================
    "SGCA": {
        "gene_id": 6442,
        "chrom": "chr17",
        "start": 50_157_516,
        "end": 50_171_427,
        "strand": "+",
        "refseq": "NC_000017.11",
        "cytoband": "17q21.33",
        "exon_count": 10,
        "role": (
            "Alpha-sarcoglycan -- component of the sarcoglycan complex in the "
            "dystrophin-associated glycoprotein complex (DGC).  Deficiency "
            "causes LGMD R3 (formerly LGMD2D), a severe childhood-onset "
            "muscular dystrophy."
        ),
        "disease": "LGMD R3 (alpha-sarcoglycanopathy)",
        "omim_disease": 608099,
        "omim_gene": 600119,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg77Cys (R77C)",
                "rsid": "rs121908948",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common in North American/European patients.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.2 kb, ideal for AAV) -- "
            "Sarepta SRP-9003 approach; "
            "(2) Base editing for R77C.  "
            "Delivery: AAVrh74 IV."
        ),
        "clinical_programs": (
            "Sarepta SRP-9003 (AAVrh74-SGCA) Phase 1/2 (NCT01976091) -- "
            "showed sustained SGCA expression at 2+ years.  "
            "Atamyo AAV-SGCA Phase 1/2 in Europe."
        ),
        "conditions": ["LGMD_R3", "sarcoglycanopathy", "neuromuscular",
                        "muscular_dystrophy"],
    },

    # ===================================================================
    # LGMD R5 (SGCG) -- Gamma-sarcoglycanopathy
    # ===================================================================
    "SGCG": {
        "gene_id": 6445,
        "chrom": "chr13",
        "start": 23_562_259,
        "end": 23_675_894,
        "strand": "+",
        "refseq": "NC_000013.11",
        "cytoband": "13q12.12",
        "exon_count": 8,
        "role": (
            "Gamma-sarcoglycan -- sarcoglycan complex component.  "
            "Deficiency causes LGMD R5 (formerly LGMD2C).  Severe, DMD-like "
            "course, common in North Africa (del521T founder mutation)."
        ),
        "disease": "LGMD R5 (gamma-sarcoglycanopathy)",
        "omim_disease": 253700,
        "omim_gene": 608896,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.525delT (del521T)",
                "rsid": "rs80358234",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "North African / Tunisian founder; >90% of alleles in Maghreb.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~0.9 kb, ideal for AAV); "
            "(2) CRISPR correction of del521T.  "
            "Delivery: AAVrh74 IV."
        ),
        "clinical_programs": (
            "Sarepta SRP-9004 (AAVrh74-SGCG) Phase 1 planned.  "
            "Genethon GNT-0006 (AAV-SGCG) Phase 1/2."
        ),
        "conditions": ["LGMD_R5", "sarcoglycanopathy", "neuromuscular",
                        "muscular_dystrophy"],
    },

    # ===================================================================
    # LGMD R4 (SGCB) -- Beta-sarcoglycanopathy
    # ===================================================================
    "SGCB": {
        "gene_id": 6443,
        "chrom": "chr4",
        "start": 52_855_297,
        "end": 52_869_660,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4q12",
        "exon_count": 6,
        "role": (
            "Beta-sarcoglycan -- sarcoglycan complex component.  "
            "Deficiency causes LGMD R4 (LGMD2E).  Can include cardiomyopathy."
        ),
        "disease": "LGMD R4 (beta-sarcoglycanopathy)",
        "omim_disease": 604286,
        "omim_gene": 600900,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Ser114Phe (S114F)",
                "rsid": "rs121908947",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Amish founder mutation.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.0 kb, ideal for AAV) -- "
            "Sarepta SRP-9005 approach; "
            "(2) Base editing for S114F.  "
            "Delivery: AAVrh74 or AAV9 IV."
        ),
        "clinical_programs": (
            "Sarepta SRP-9005 (AAVrh74-SGCB) Phase 1 (NCT03652259) -- "
            "showed beta-sarcoglycan expression."
        ),
        "conditions": ["LGMD_R4", "sarcoglycanopathy", "neuromuscular",
                        "muscular_dystrophy", "cardiomyopathy"],
    },

    # ===================================================================
    # LGMD R1 (CAPN3) -- Calpainopathy
    # ===================================================================
    "CAPN3": {
        "gene_id": 825,
        "chrom": "chr15",
        "start": 42_358_084,
        "end": 42_411_614,
        "strand": "+",
        "refseq": "NC_000015.10",
        "cytoband": "15q15.1",
        "exon_count": 24,
        "role": (
            "Calpain-3 -- muscle-specific calcium-activated cysteine protease.  "
            "Deficiency causes LGMD R1 (LGMD2A), the most common form of LGMD "
            "worldwide (~30% of all LGMD)."
        ),
        "disease": "LGMD R1 (calpainopathy)",
        "omim_disease": 253600,
        "omim_gene": 114240,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Arg490Trp (R490W)",
                "rsid": "rs121908954",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Active site adjacent; European.",
            },
            {
                "name": "c.550delA",
                "rsid": "rs80338803",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Basque founder mutation.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.5 kb fits AAV); "
            "(2) Base editing for R490W.  "
            "Delivery: AAVrh74 or AAV9 IV."
        ),
        "clinical_programs": (
            "Atamyo/Sarepta AAV-CAPN3 preclinical.  "
            "No clinical trials as of 2026."
        ),
        "conditions": ["LGMD_R1", "calpainopathy", "neuromuscular",
                        "muscular_dystrophy"],
    },

    # ===================================================================
    # LGMD R9 (FKRP) -- Dystroglycanopathy
    # ===================================================================
    "FKRP": {
        "gene_id": 79147,
        "chrom": "chr19",
        "start": 46_750_673,
        "end": 46_765_017,
        "strand": "-",
        "refseq": "NC_000019.10",
        "cytoband": "19q13.32",
        "exon_count": 4,
        "role": (
            "Fukutin-related protein -- glycosyltransferase required for "
            "alpha-dystroglycan O-mannosylation.  Mutations cause a spectrum "
            "from severe Walker-Warburg/MEB to mild LGMD R9 (LGMD2I).  "
            "L276I is the common mild mutation."
        ),
        "disease": "LGMD R9 (FKRP-related dystroglycanopathy)",
        "omim_disease": 607155,
        "omim_gene": 606596,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Leu276Ile (L276I)",
                "rsid": "rs28937900",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Most common FKRP mutation worldwide; ~90% of LGMD R9 "
                    "patients carry at least one L276I allele.  Mild phenotype."
                ),
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.5 kb, ideal for AAV) -- "
            "Atamyo/ML Bio approach; "
            "(2) Base editing for L276I.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": (
            "ML Bio Solutions ML-301 (AAV9-FKRP) Phase 1 planned.  "
            "Atamyo GNT-0004 (AAV-FKRP) preclinical."
        ),
        "conditions": ["LGMD_R9", "dystroglycanopathy", "neuromuscular",
                        "muscular_dystrophy", "walker_warburg"],
    },

    # ===================================================================
    # LGMD R12 (ANO5) -- Anoctaminopathy
    # ===================================================================
    "ANO5": {
        "gene_id": 203859,
        "chrom": "chr11",
        "start": 22_212_914,
        "end": 22_285_345,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11p14.3",
        "exon_count": 22,
        "role": (
            "Anoctamin-5 -- putative calcium-activated chloride channel in "
            "muscle membrane repair.  Mutations cause LGMD R12 (LGMD2L) "
            "and distal Miyoshi-like myopathy.  Common in Northern Europeans."
        ),
        "disease": "LGMD R12 (anoctaminopathy)",
        "omim_disease": 611307,
        "omim_gene": 608662,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.191dupA (p.Asn64Lysfs)",
                "rsid": "rs137854526",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Most common in Northern Europeans; founder mutation.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~2.8 kb fits AAV); "
            "(2) CRISPR correction of c.191dupA.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": "Preclinical only.",
        "conditions": ["LGMD_R12", "anoctaminopathy", "neuromuscular",
                        "muscular_dystrophy"],
    },

    # ===================================================================
    # X-linked Myotubular Myopathy (MTM1)
    # ===================================================================
    "MTM1": {
        "gene_id": 4534,
        "chrom": "chrX",
        "start": 150_555_808,
        "end": 150_659_622,
        "strand": "-",
        "refseq": "NC_000023.11",
        "cytoband": "Xq28",
        "exon_count": 15,
        "role": (
            "Myotubularin -- phosphoinositide phosphatase essential for "
            "endosomal/lysosomal trafficking and muscle T-tubule formation.  "
            "X-linked myotubular myopathy (XLMTM): profound neonatal hypotonia, "
            "respiratory failure, ~50% mortality in first 18 months."
        ),
        "disease": "X-linked myotubular myopathy (XLMTM)",
        "omim_disease": 310400,
        "omim_gene": 300415,
        "inheritance": "XLR",
        "key_variants": [
            {
                "name": "Truncating mutations throughout",
                "rsid": None,
                "consequence": "nonsense/frameshift/splice",
                "clinical_significance": "pathogenic",
                "notes": "~75% are truncating; no hotspot.",
            },
        ],
        "strategy": (
            "(1) AAV8 gene replacement (cDNA ~1.8 kb, ideal for AAV) -- "
            "Astellas/Audentes AT132 approach; "
            "(2) CRISPRa to upregulate MTMR2 (functional paralog) as "
            "alternative rescue strategy.  "
            "Delivery: AAV8 IV."
        ),
        "clinical_programs": (
            "Astellas AT132 (AAV8-MTM1) Phase 1/2 ASPIRO trial "
            "(NCT03199469) -- showed dramatic improvement in survivors but "
            "4 deaths from hepatotoxicity at high dose (clinical hold).  "
            "Dose optimization ongoing.  "
            "Dynacure DYN101 (anti-DNM2 ASO, alternative approach) Phase 1/2."
        ),
        "conditions": ["XLMTM", "myotubular_myopathy", "neuromuscular",
                        "congenital_myopathy", "neonatal_hypotonia"],
    },

    # ===================================================================
    # Congenital Myasthenic Syndrome -- RAPSN
    # ===================================================================
    "RAPSN": {
        "gene_id": 5913,
        "chrom": "chr11",
        "start": 47_465_698,
        "end": 47_479_424,
        "strand": "+",
        "refseq": "NC_000011.10",
        "cytoband": "11p11.2",
        "exon_count": 8,
        "role": (
            "Rapsyn -- postsynaptic scaffolding protein that clusters AChRs "
            "at the neuromuscular junction.  Mutations cause CMS with "
            "fatigable weakness."
        ),
        "disease": "Congenital myasthenic syndrome (CMS)",
        "omim_disease": 616326,
        "omim_gene": 601592,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Asn88Lys (N88K)",
                "rsid": "rs28940574",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common CMS mutation; found on diverse haplotypes.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.2 kb, ideal for AAV); "
            "(2) Base editing for N88K.  "
            "Delivery: AAV9 IV or intramuscular."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "3,4-DAP and salbutamol used off-label for CMS."
        ),
        "conditions": ["congenital_myasthenic_syndrome", "CMS", "neuromuscular",
                        "NMJ_disorder"],
    },

    # ===================================================================
    # Congenital Myasthenic Syndrome -- DOK7
    # ===================================================================
    "DOK7": {
        "gene_id": 285489,
        "chrom": "chr4",
        "start": 3_414_020,
        "end": 3_420_985,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4p16.3",
        "exon_count": 7,
        "role": (
            "Downstream of kinase 7 -- activator of MuSK receptor kinase, "
            "essential for NMJ formation and maintenance.  Mutations cause "
            "DOK7-CMS with limb-girdle pattern weakness."
        ),
        "disease": "DOK7 congenital myasthenic syndrome",
        "omim_disease": 254300,
        "omim_gene": 610285,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1124_1127dupTGCC",
                "rsid": "rs387906688",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Most common DOK7 mutation worldwide.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.5 kb fits AAV); "
            "(2) AAV-DOK7 restored NMJ in DOK7-null mice.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": "Preclinical AAV-DOK7 gene therapy.",
        "conditions": ["congenital_myasthenic_syndrome", "CMS", "neuromuscular",
                        "DOK7_CMS", "NMJ_disorder"],
    },
}
