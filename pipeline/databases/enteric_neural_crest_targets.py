"""
Enteric / neural crest disorder gene-therapy targets -- GRCh38 (hg38).

Covers Hirschsprung disease (aganglionic megacolon), Waardenburg syndrome,
and congenital central hypoventilation syndrome (CCHS/Ondine's curse) --
all arising from neural crest cell development defects.

Coordinates from NCBI Gene, GRCh38.p14.

All interventions require informed patient consent and IRB / ethics approval.
"""


ENTERIC_NEURAL_CREST_TARGETS = {

    "EDNRB": {
        "gene_id": 1910,
        "chrom": "chr13",
        "start": 77_370_278,
        "end": 77_396_937,
        "strand": "+",
        "refseq": "NC_000013.11",
        "cytoband": "13q22.3",
        "exon_count": 7,
        "role": (
            "Endothelin receptor type B -- GPCR on neural crest cells, "
            "essential for enteric nervous system (ENS) and melanocyte "
            "development.  Mutations cause Hirschsprung disease (HSCR) with "
            "or without Waardenburg-Shah syndrome (WS4): aganglionic megacolon "
            "+ sensorineural deafness + pigmentary anomalies.  Incidence of "
            "HSCR ~1:5000."
        ),
        "disease": "Hirschsprung disease / Waardenburg-Shah syndrome (WS4B)",
        "omim_disease": 613265,
        "omim_gene": 131244,
        "inheritance": "AD/AR (HSCR is oligogenic with variable penetrance)",
        "key_variants": [
            {
                "name": "p.Trp276Cys (W276C, Mennonite founder)",
                "rsid": "rs41310644",
                "consequence": "missense",
                "clinical_significance": "pathogenic",
                "notes": "Homozygous -> WS4 (HSCR + deafness + white hair).",
            },
            {
                "name": "Various LOF mutations",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "HSCR is genetically complex; EDNRB is one of many loci.",
            },
        ],
        "strategy": (
            "(1) Surgical pull-through (Soave/Duhamel/Swenson) is curative "
            "for bowel aganglionosis; "
            "(2) Enteric neural stem cell transplant to aganglionic bowel "
            "(preclinical -- restore ENS); "
            "(3) AAV gene replacement (cDNA ~1.3 kb, fits AAV) for research; "
            "(4) Gene therapy mainly relevant for preventing HSCR in "
            "prenatal/neonatal setting (future).  "
            "Note: Surgery is highly effective, reducing gene therapy urgency."
        ),
        "clinical_programs": (
            "No gene therapy trials.  "
            "ENS stem cell transplant preclinical (Burns/Goldstein labs).  "
            "Pull-through surgery SOC."
        ),
        "conditions": ["hirschsprung_disease", "HSCR", "waardenburg_shah",
                        "WS4", "neural_crest", "aganglionic_megacolon"],
    },

    "EDN3": {
        "gene_id": 1908,
        "chrom": "chr20",
        "start": 59_296_360,
        "end": 59_321_276,
        "strand": "-",
        "refseq": "NC_000020.11",
        "cytoband": "20q13.32",
        "exon_count": 5,
        "role": (
            "Endothelin-3 -- EDNRB ligand expressed in gut mesenchyme, "
            "provides chemoattractant signal for enteric neural crest cell "
            "migration.  Mutations cause Waardenburg-Shah syndrome type 4C "
            "or isolated Hirschsprung disease."
        ),
        "disease": "Waardenburg-Shah syndrome type 4C / Hirschsprung",
        "omim_disease": 613266,
        "omim_gene": 131242,
        "inheritance": "AR (WS4C) / AD modifier (HSCR)",
        "key_variants": [
            {
                "name": "Various LOF mutations",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": "Rare; most HSCR is oligogenic.",
            },
        ],
        "strategy": (
            "(1) AAV gene replacement (cDNA ~0.6 kb, ideal for AAV); "
            "(2) Surgery is curative for bowel aganglionosis."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["hirschsprung_disease", "waardenburg_shah", "WS4C",
                        "neural_crest"],
    },

    "SOX10": {
        "gene_id": 6663,
        "chrom": "chr22",
        "start": 38_368_318,
        "end": 38_381_200,
        "strand": "+",
        "refseq": "NC_000022.11",
        "cytoband": "22q13.1",
        "exon_count": 4,
        "role": (
            "SRY-box transcription factor 10 -- master regulator of neural "
            "crest-derived melanocytes, enteric neurons, and Schwann cells.  "
            "Mutations cause Waardenburg syndrome type 2E (WS2E) or type 4C "
            "(WS4C with HSCR): sensorineural deafness, pigmentary anomalies "
            "(white forelock, heterochromia iridis), and/or Hirschsprung."
        ),
        "disease": "Waardenburg syndrome type 2E/4C / PCWH syndrome",
        "omim_disease": 611584,
        "omim_gene": 602229,
        "inheritance": "AD",
        "key_variants": [
            {
                "name": "Various LOF and dominant-negative mutations",
                "rsid": None,
                "consequence": "various",
                "clinical_significance": "pathogenic",
                "notes": (
                    "HMG domain mutations -> WS2E (deafness + pigment).  "
                    "Truncating near C-terminus -> WS4/PCWH (+ HSCR + "
                    "peripheral neuropathy)."
                ),
            },
        ],
        "strategy": (
            "(1) CRISPRa to upregulate wild-type SOX10 allele; "
            "(2) AAV gene replacement (cDNA ~1.4 kb, fits AAV) -- timing "
            "critical (neural crest migration is embryonic); "
            "(3) Cochlear implant for deafness; surgery for HSCR."
        ),
        "clinical_programs": "No gene therapy trials.",
        "conditions": ["waardenburg_syndrome", "WS2E", "WS4C", "PCWH",
                        "neural_crest", "deafness", "hirschsprung_disease"],
    },

    "PHOX2B": {
        "gene_id": 8929,
        "chrom": "chr4",
        "start": 41_744_083,
        "end": 41_749_159,
        "strand": "-",
        "refseq": "NC_000004.12",
        "cytoband": "4p13",
        "exon_count": 3,
        "role": (
            "Paired-like homeobox 2B -- transcription factor essential for "
            "development of autonomic nervous system (ANS).  Polyalanine "
            "repeat expansions (20+5 to 20+13) or non-polyalanine mutations "
            "cause congenital central hypoventilation syndrome (CCHS / "
            "Ondine's curse): failure of automatic breathing during sleep, "
            "Hirschsprung in ~20%, neuroblastoma risk."
        ),
        "disease": "Congenital central hypoventilation syndrome (CCHS, Ondine's curse)",
        "omim_disease": 209880,
        "omim_gene": 603851,
        "inheritance": "AD (mostly de novo)",
        "key_variants": [
            {
                "name": "Polyalanine repeat expansion (20+5 to 20+13 Ala)",
                "rsid": None,
                "consequence": "repeat_expansion",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~90% of CCHS.  Normal: 20 Ala residues.  20+5 (25 Ala) = "
                    "mildest.  20+13 (33 Ala) = severe + HSCR + neuroblastoma.  "
                    "Genotype predicts need for ventilatory support."
                ),
            },
            {
                "name": "Non-polyalanine repeat mutations (NPARM)",
                "rsid": None,
                "consequence": "frameshift/missense",
                "clinical_significance": "pathogenic",
                "notes": (
                    "~10% of CCHS; more severe: higher HSCR and neuroblastoma "
                    "risk; some require continuous ventilation."
                ),
            },
        ],
        "strategy": (
            "(1) Allele-specific silencing of expanded/mutant PHOX2B allele "
            "(preserving WT is essential for ANS function); "
            "(2) CRISPRi targeting the expanded allele; "
            "(3) Gene therapy to restore respiratory chemosensitivity in "
            "brainstem retrotrapezoid nucleus (RTN); "
            "(4) Phox2b-expressing DREADD for chemogenetic rescue of RTN "
            "CO2 sensitivity (preclinical).  "
            "Delivery: AAV to brainstem RTN (challenging) or intrathecal.  "
            "Standard of care: lifetime mechanical ventilation during sleep "
            "or diaphragm pacing."
        ),
        "clinical_programs": (
            "No gene therapy trials as of 2026.  "
            "Diaphragm pacing (Avery Biomedical) reduces ventilator "
            "dependence.  DREADD-based respiratory rescue preclinical "
            "(Bhatt lab, U. Florida)."
        ),
        "conditions": ["congenital_central_hypoventilation", "CCHS",
                        "ondine_curse", "neural_crest", "autonomic_dysfunction",
                        "hirschsprung_disease"],
    },
}
