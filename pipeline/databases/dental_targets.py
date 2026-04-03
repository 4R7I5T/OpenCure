"""
Dental / tooth development disorder gene-therapy targets -- GRCh38 (hg38).

Covers amelogenesis imperfecta (AI) and dentinogenesis imperfecta (DI),
hereditary conditions affecting enamel and dentin formation respectively.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


DENTAL_TARGETS = {

    "AMELX": {
        "gene_id": 265,
        "chrom": "chrX",
        "start": 11_311_471,
        "end": 11_321_048,
        "strand": "+",
        "refseq": "NC_000023.11",
        "cytoband": "Xp22.2",
        "exon_count": 7,
        "role": (
            "Amelogenin X-linked -- major enamel matrix protein (~90% of "
            "enamel organic matrix), secreted by ameloblasts during tooth "
            "development.  Controls enamel crystal growth and organization.  "
            "Mutations cause X-linked amelogenesis imperfecta (AI): "
            "hypoplastic or hypomineralized enamel, discolored fragile teeth, "
            "temperature sensitivity."
        ),
        "disease": "X-linked amelogenesis imperfecta (AI)",
        "omim_disease": 301200,
        "omim_gene": 300391,
        "inheritance": "XL",
        "key_variants": [
            {
                "name": "Various missense and truncating",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Hemizygous males more severely affected.",
            },
        ],
        "strategy": (
            "(1) Recombinant amelogenin biomimetic enamel repair (topical "
            "application during tooth development); "
            "(2) Gene therapy to ameloblasts is challenging as they are "
            "terminally differentiated and lost after tooth eruption; "
            "(3) Most practical approach: stem cell-derived bioengineered "
            "tooth enamel or dental crowns/veneers.  "
            "Note: Gene therapy for dental conditions is in very early stages "
            "due to limited access to ameloblasts in vivo."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "Amelogenin-based enamel repair peptides (Curodont, AMEL) in "
            "early development for caries and AI."
        ),
        "conditions": ["amelogenesis_imperfecta", "AI", "dental",
                        "enamel_defect"],
    },

    "ENAM": {
        "gene_id": 10117,
        "chrom": "chr4",
        "start": 70_586_891,
        "end": 70_611_783,
        "strand": "-",
        "refseq": "NC_000004.12",
        "cytoband": "4q13.3",
        "exon_count": 10,
        "role": (
            "Enamelin -- second most abundant enamel matrix protein.  "
            "Mutations cause autosomal dominant AI (hypoplastic-local) or "
            "autosomal recessive AI (severe generalized)."
        ),
        "disease": "Amelogenesis imperfecta (ENAM-related)",
        "omim_disease": 104500,
        "omim_gene": 606585,
        "inheritance": "AD or AR",
        "key_variants": [
            {
                "name": "g.13185-13186insAG (common AD mutation)",
                "rsid": None,
                "consequence": "frameshift",
                "clinical_significance": "pathogenic",
                "notes": "Heterozygous: local hypoplasia.  Homozygous: severe generalized AI.",
            },
        ],
        "strategy": (
            "(1) Similar challenges as AMELX -- ameloblasts not accessible "
            "post-eruption; "
            "(2) Future: gene-corrected autologous dental stem cells for "
            "bioengineered tooth regeneration; "
            "(3) Restorative dentistry remains primary management."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["amelogenesis_imperfecta", "AI", "dental",
                        "enamel_defect"],
    },

    "FAM83H": {
        "gene_id": 286077,
        "chrom": "chr8",
        "start": 143_799_502,
        "end": 143_819_816,
        "strand": "-",
        "refseq": "NC_000008.11",
        "cytoband": "8q24.3",
        "exon_count": 5,
        "role": (
            "Family with sequence similarity 83 member H -- involved in "
            "keratin cytoskeleton organization in ameloblasts.  Mutations "
            "cause autosomal dominant hypocalcified AI (most common AD AI "
            "type): normal enamel thickness but severely undermineralized, "
            "soft and rapidly worn."
        ),
        "disease": "Hypocalcified amelogenesis imperfecta (AD)",
        "omim_disease": 130900,
        "omim_gene": 611927,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "C-terminal truncating mutations (exon 5)",
                "rsid": None,
                "consequence": "nonsense/frameshift",
                "clinical_significance": "pathogenic",
                "notes": "All known mutations cluster in last exon; GOF/dominant-negative.",
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of mutant allele (since mutations "
            "cluster in exon 5, ASO approach targeting mutant mRNA); "
            "(2) Bioengineered tooth approaches."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["amelogenesis_imperfecta", "AI", "dental",
                        "hypocalcified_enamel"],
    },

    "DSPP": {
        "gene_id": 1834,
        "chrom": "chr4",
        "start": 87_566_536,
        "end": 87_575_151,
        "strand": "+",
        "refseq": "NC_000004.12",
        "cytoband": "4q22.1",
        "exon_count": 5,
        "role": (
            "Dentin sialophosphoprotein -- major non-collagenous dentin "
            "matrix protein, cleaved into dentin sialoprotein (DSP) and "
            "dentin phosphoprotein (DPP).  Mutations cause dentinogenesis "
            "imperfecta types II and III (DGI-II/III) and dentin dysplasia "
            "type II: opalescent discolored teeth, obliterated pulp chambers, "
            "periapical radiolucencies."
        ),
        "disease": "Dentinogenesis imperfecta type II-III / Dentin dysplasia II",
        "omim_disease": 125490,
        "omim_gene": 125485,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Signal peptide mutations (p.Val18Phe, p.Ala15Val)",
                "rsid": None,
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Disrupt ER targeting of DSPP; intracellular retention.",
            },
            {
                "name": "DPP domain repeat expansion/deletion",
                "rsid": None,
                "consequence": "in_frame_deletion",
                "clinical_significance": "pathogenic",
                "notes": "Repetitive serine-serine-aspartate domain; difficult to sequence.",
            },
        ],
        "strategy": (
            "(1) Gene therapy to odontoblasts is technically challenging; "
            "(2) Allele-specific silencing of dominant mutant for DGI; "
            "(3) Most practical: stem cell-based tooth regeneration or "
            "restorative/prosthetic dentistry."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["dentinogenesis_imperfecta", "DGI", "dentin_dysplasia",
                        "dental"],
    },

    "DMP1": {
        "gene_id": 1758,
        "chrom": "chr4",
        "start": 87_518_777,
        "end": 87_530_870,
        "strand": "-",
        "refseq": "NC_000004.12",
        "cytoband": "4q22.1",
        "exon_count": 6,
        "role": (
            "Dentin matrix acidic phosphoprotein 1 -- mineralization protein "
            "expressed in odontoblasts and osteocytes.  Regulates "
            "phosphate/FGF23 homeostasis.  Mutations cause autosomal "
            "recessive hypophosphatemic rickets (ARHR1) with dentin defects."
        ),
        "disease": "Autosomal recessive hypophosphatemic rickets / dentin defect",
        "omim_disease": 241520,
        "omim_gene": 600980,
        "inheritance": "AR",
        "key_variants": [
            {
                "name": "c.1A>G (start codon loss)",
                "rsid": None,
                "consequence": "start_loss",
                "clinical_significance": "pathogenic",
                "notes": "Complete DMP1 loss; rickets + dentin dysplasia.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~1.5 kb fits AAV) to bone "
            "and dentin; "
            "(2) Phosphate + calcitriol supplementation for rickets; "
            "(3) Burosumab (anti-FGF23) may address phosphate wasting."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["hypophosphatemic_rickets", "ARHR", "dental",
                        "dentin_defect", "bone_mineralization"],
    },
}
