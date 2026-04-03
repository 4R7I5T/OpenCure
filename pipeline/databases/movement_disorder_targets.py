"""
Movement disorder gene-therapy / CRISPR targets -- GRCh38 (hg38).

Covers monogenic dystonias, dopa-responsive dystonia, paroxysmal dyskinesias,
and myoclonus-dystonia.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


MOVEMENT_DISORDER_TARGETS = {

    "TOR1A": {
        "gene_id": 1861,
        "chrom": "chr9",
        "start": 129_812_442,
        "end": 129_824_424,
        "strand": "+",
        "refseq": "NC_000009.12",
        "cytoband": "9q34.11",
        "exon_count": 5,
        "role": (
            "Torsin-1A -- AAA+ ATPase in the ER/nuclear envelope lumen.  "
            "The 3-bp deletion (deltaE, c.907_909delGAG) causes DYT1 early-"
            "onset generalized dystonia, the most common monogenic dystonia.  "
            "~30% penetrance.  Onset typically in childhood (limb dystonia "
            "spreading to generalized)."
        ),
        "disease": "DYT1 early-onset generalized dystonia",
        "omim_disease": 128100,
        "omim_gene": 605204,
        "inheritance": "AD (30% penetrance)",
        "key_variants": [
            {
                "name": "c.907_909delGAG (deltaE, p.Glu303del)",
                "rsid": "rs80358233",
                "consequence": "in_frame_deletion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "Single mutation accounts for >90% of DYT1; removes "
                    "one glutamic acid from C-terminal.  Ashkenazi Jewish "
                    "founder (carrier ~1/3000)."
                ),
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant TOR1A allele using "
            "CRISPRi or Cas13 (requires distinguishing 3-bp deletion from WT); "
            "(2) AAV-mediated overexpression of wild-type TOR1A to outcompete "
            "dominant-negative mutant protein; "
            "(3) Deep brain stimulation (DBS) remains gold-standard treatment.  "
            "Delivery: AAV9 to basal ganglia (stereotactic)."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "DBS (GPi target) provides ~50-70% improvement in DYT1."
        ),
        "conditions": ["DYT1", "dystonia", "movement_disorder",
                        "generalized_dystonia"],
    },

    "GCH1": {
        "gene_id": 2643,
        "chrom": "chr14",
        "start": 54_841_939,
        "end": 54_874_468,
        "strand": "-",
        "refseq": "NC_000014.9",
        "cytoband": "14q22.2",
        "exon_count": 6,
        "role": (
            "GTP cyclohydrolase 1 -- rate-limiting enzyme in tetrahydrobiopterin "
            "(BH4) biosynthesis.  BH4 is essential cofactor for tyrosine "
            "hydroxylase (dopamine synthesis), phenylalanine hydroxylase, and "
            "nitric oxide synthase.  Haploinsufficiency causes dopa-responsive "
            "dystonia (DRD, Segawa syndrome): childhood-onset dystonia with "
            "dramatic and sustained response to low-dose L-DOPA."
        ),
        "disease": "Dopa-responsive dystonia (DRD, Segawa syndrome, DYT5a)",
        "omim_disease": 128230,
        "omim_gene": 600225,
        "inheritance": "AD (DRD); AR (severe BH4 deficiency)",
        "key_variants": [
            {
                "name": "Various LOF mutations (>150)",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "No hotspot; haploinsufficiency mechanism.",
            },
        ],
        "strategy": (
            "(1) L-DOPA therapy is so effective that gene therapy is rarely "
            "needed (patients respond lifelong to low-dose L-DOPA/carbidopa); "
            "(2) AAV gene replacement (cDNA ~0.8 kb, ideal for AAV) for "
            "severe AR forms (complete BH4 deficiency) resistant to L-DOPA; "
            "(3) CRISPRa to upregulate wild-type allele.  "
            "Delivery: AAV2 to striatum for severe cases."
        ),
        "clinical_programs": (
            "No gene therapy trials (L-DOPA is curative for AD form).  "
            "AAV-AADC gene therapy trials for AADC deficiency (related "
            "BH4 pathway disorder) show proof-of-concept."
        ),
        "conditions": ["dopa_responsive_dystonia", "DRD", "segawa_syndrome",
                        "DYT5", "movement_disorder", "BH4_deficiency"],
    },

    "TH": {
        "gene_id": 7054,
        "chrom": "chr11",
        "start": 2_165_914,
        "end": 2_174_259,
        "strand": "-",
        "refseq": "NC_000011.10",
        "cytoband": "11p15.5",
        "exon_count": 14,
        "role": (
            "Tyrosine hydroxylase -- rate-limiting enzyme in catecholamine "
            "biosynthesis (tyrosine -> L-DOPA).  Biallelic mutations cause "
            "tyrosine hydroxylase deficiency: spectrum from DRD to severe "
            "progressive encephalopathy with autonomic dysfunction."
        ),
        "disease": "Tyrosine hydroxylase deficiency / DYT5b",
        "omim_disease": 605407,
        "omim_gene": 191290,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "p.Leu205Pro (L205P)",
                "rsid": "rs80338909",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Severe infantile form; near-complete enzyme loss.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.5 kb, fits AAV) to substantia "
            "nigra/striatum; "
            "(2) L-DOPA for mild forms; severe forms respond poorly.  "
            "Delivery: AAV2 stereotactic to putamen."
        ),
        "clinical_programs": (
            "PTC Therapeutics PTC-AADC (AAV2-AADC for AADC deficiency) -- "
            "related pathway, FDA-approved 2024 as first gene therapy for "
            "a brain disorder delivered directly to the brain."
        ),
        "conditions": ["TH_deficiency", "DYT5b", "movement_disorder",
                        "dopamine_deficiency", "infantile_parkinsonism"],
    },

    "ATP1A3": {
        "gene_id": 478,
        "chrom": "chr19",
        "start": 42_470_733,
        "end": 42_498_315,
        "strand": "-",
        "refseq": "NC_000019.10",
        "cytoband": "19q13.2",
        "exon_count": 23,
        "role": (
            "Na+/K+ ATPase alpha-3 subunit -- neuronal sodium pump maintaining "
            "electrochemical gradient.  Mutations cause a spectrum: rapid-onset "
            "dystonia-parkinsonism (DYT12, RDP), alternating hemiplegia of "
            "childhood (AHC), and CAPOS syndrome."
        ),
        "disease": "RDP (DYT12) / AHC / CAPOS",
        "omim_disease": 128235,
        "omim_gene": 182350,
        "inheritance": "AD (mostly de novo)",
        "key_variants": [
            {
                "name": "p.Asp801Asn (D801N)",
                "rsid": "rs121918314",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common AHC mutation (~30%); phosphorylation domain.",
            },
            {
                "name": "p.Glu815Lys (E815K)",
                "rsid": "rs121918316",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Severe AHC with epilepsy.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant ATP1A3 (GOF/dominant-negative); "
            "(2) AAV gene replacement (cDNA ~3.3 kb fits AAV); "
            "(3) Base editing for recurrent D801N.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "Flunarizine partial benefit for AHC (not gene therapy)."
        ),
        "conditions": ["rapid_onset_dystonia_parkinsonism", "AHC",
                        "alternating_hemiplegia", "CAPOS", "movement_disorder"],
    },

    "PRRT2": {
        "gene_id": 112476,
        "chrom": "chr16",
        "start": 29_812_193,
        "end": 29_816_020,
        "strand": "-",
        "refseq": "NC_000016.10",
        "cytoband": "16p11.2",
        "exon_count": 4,
        "role": (
            "Proline-rich transmembrane protein 2 -- presynaptic protein "
            "regulating neurotransmitter release and synaptic vesicle fusion.  "
            "Mutations cause paroxysmal kinesigenic dyskinesia (PKD), "
            "benign familial infantile epilepsy (BFIE), and infantile "
            "convulsions with choreoathetosis (ICCA)."
        ),
        "disease": "Paroxysmal kinesigenic dyskinesia (PKD) / BFIE / ICCA",
        "omim_disease": 128200,
        "omim_gene": 614386,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "c.649dupC (p.Arg217Profs*8)",
                "rsid": "rs587784244",
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Most common (~80% of PKD families); hot spot.",
            },
        ],
        "strategy": (
            "(1) Carbamazepine is highly effective for PKD (most patients "
            "are attack-free on low-dose CBZ); gene therapy rarely needed; "
            "(2) AAV gene replacement (cDNA ~1.0 kb, ideal for AAV) for "
            "refractory cases.  "
            "Delivery: AAV9 IV."
        ),
        "clinical_programs": "No gene therapy trials (carbamazepine is curative).",
        "conditions": ["paroxysmal_kinesigenic_dyskinesia", "PKD", "BFIE",
                        "movement_disorder", "paroxysmal_dyskinesia"],
    },

    "SGCE": {
        "gene_id": 8910,
        "chrom": "chr7",
        "start": 94_612_459,
        "end": 94_697_476,
        "strand": "+",
        "refseq": "NC_000007.14",
        "cytoband": "7q21.3",
        "exon_count": 12,
        "role": (
            "Epsilon-sarcoglycan -- transmembrane protein in brain (function "
            "uncertain, possibly synaptic).  Maternally imprinted (paternal "
            "allele expressed).  Paternal mutations cause myoclonus-dystonia "
            "(DYT11): childhood-onset myoclonus + dystonia, often with "
            "psychiatric comorbidities (anxiety, OCD, alcohol dependence).  "
            "Alcohol-responsive."
        ),
        "disease": "Myoclonus-dystonia (DYT11)",
        "omim_disease": 159900,
        "omim_gene": 604149,
        "inheritance": "AD (imprinted, paternal allele)",
        "key_variants": [
            {
                "name": "Various LOF mutations (no hotspot)",
                "rsid": None,
                "consequence": "loss_of_function",
                "clinical_significance": "pathogenic",
                "notes": "Only pathogenic when inherited from father (maternal imprinting).",
            },
        ],
        "strategy": (
            "(1) CRISPRa to reactivate silenced maternal SGCE allele in "
            "brain (imprinting reversal); "
            "(2) AAV gene replacement (cDNA ~1.3 kb, ideal for AAV); "
            "(3) DBS for refractory cases.  "
            "Delivery: AAV9 to basal ganglia."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "DBS (GPi or VIM) provides benefit for refractory myoclonus-dystonia."
        ),
        "conditions": ["myoclonus_dystonia", "DYT11", "movement_disorder",
                        "myoclonus", "imprinting_disorder"],
    },

    "ADCY5": {
        "gene_id": 111,
        "chrom": "chr3",
        "start": 123_281_186,
        "end": 123_462_741,
        "strand": "+",
        "refseq": "NC_000003.12",
        "cytoband": "3q21.1",
        "exon_count": 25,
        "role": (
            "Adenylate cyclase 5 -- striatal-enriched cAMP-producing enzyme.  "
            "Gain-of-function mutations cause ADCY5-related dyskinesia: "
            "childhood-onset choreoathetosis, dystonia, and myoclonus "
            "worsened by drowsiness.  Recently recognized entity."
        ),
        "disease": "ADCY5-related dyskinesia",
        "omim_disease": 606703,
        "omim_gene": 600293,
        "inheritance": "AD (mostly de novo, GOF)",
        "key_variants": [
            {
                "name": "p.Arg418Trp (R418W)",
                "rsid": "rs587777643",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Most common; catalytic domain; constitutive activation.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of GOF mutant allele; "
            "(2) CRISPRi targeted to mutant ADCY5; "
            "(3) Caffeine and tetrabenazine provide partial symptom relief.  "
            "Delivery: AAV to striatum."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["ADCY5_dyskinesia", "choreoathetosis", "movement_disorder",
                        "dyskinesia"],
    },
}
