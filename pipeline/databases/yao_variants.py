"""
Yao Syndrome (NOD2-associated autoinflammatory disease) variant catalog.

Comprehensive database of NOD2 variants with Yao syndrome associations,
functional effects, population frequencies, and compound heterozygosity
patterns.  Also includes key variants in modifier genes (MEFV, NLRP3,
NLRP12, TNFRSF1A) found in digenic YAOS patients.

Amino acid numbering uses the NM_022162.3 / NP_071445.1 convention
(1040 aa isoform) unless otherwise noted, as this is the standard used
in clinical literature.  MANE Select (NM_001370466.1, 1013 aa) numbering
is provided where applicable.

Data sources:
  - Nomani et al. 2024 (PMC11449693) — 194-patient cohort
  - Yao et al. 2018 (PMC6036904) — pathway/cytokine characterization
  - Maekawa et al. 2016 (PMC4906405) — NOD2 crystal structure
  - Mayo Clinic 2024 (PMC11486699) — 22-patient cohort
  - dbSNP, ClinVar, gnomAD v4
"""

from dataclasses import dataclass, field


@dataclass
class KnownVariant:
    """A cataloged disease-associated variant."""
    gene: str
    rsid: str
    name: str                    # common name (e.g. "R702W", "IVS8+158")
    chrom: str
    pos: int                     # GRCh38 1-based position
    ref: str
    alt: str
    consequence: str             # missense, frameshift, intronic, etc.
    amino_acid: str              # e.g. "p.Arg702Trp" (NM_022162.3)
    amino_acid_mane: str         # MANE Select numbering
    domain: str                  # protein domain
    exon: str                    # exon/intron location
    functional_effect: str       # gain_of_function, loss_of_function, mixed
    gnomad_af: float             # gnomAD global allele frequency
    yaos_frequency: float        # fraction of YAOS patients carrying this
    diseases: list[str] = field(default_factory=list)
    clinvar_significance: str = ""
    notes: str = ""


# ============================================================================
# NOD2 VARIANTS ASSOCIATED WITH YAO SYNDROME
# ============================================================================
# Ordered by YAOS frequency (highest first)

NOD2_YAO_VARIANTS: list[KnownVariant] = [
    KnownVariant(
        gene="NOD2",
        rsid="rs5743289",
        name="IVS8+158",
        chrom="chr16",
        pos=50722863,
        ref="C",
        alt="T",
        consequence="intronic",
        amino_acid="N/A (intronic)",
        amino_acid_mane="N/A (intronic)",
        domain="intronic (near LRR)",
        exon="intron 8",
        functional_effect="gain_of_function",
        gnomad_af=0.1018,
        yaos_frequency=0.60,
        diseases=["yao_syndrome"],
        clinvar_significance="uncertain_significance",
        notes="Primary Yao syndrome variant. Elevates NOD2 transcript levels "
              "(2.47x) and basal p38 MAPK phosphorylation. Increases basal "
              "IL-6 secretion from PBMCs. Does NOT affect intron 8 splicing. "
              "Found in 60% of YAOS patients (116/194). CRITICAL: standard "
              "exome panels may MISS this intronic variant.",
    ),
    KnownVariant(
        gene="NOD2",
        rsid="rs5743291",
        name="V955I",
        chrom="chr16",
        pos=50723365,
        ref="G",
        alt="A",
        consequence="missense",
        amino_acid="p.Val955Ile",
        amino_acid_mane="p.Val928Ile",
        domain="LRR",
        exon="exon 9",
        functional_effect="uncertain",
        gnomad_af=0.0595,
        yaos_frequency=0.26,
        diseases=["yao_syndrome"],
        clinvar_significance="benign/uncertain_significance",
        notes="Second most common YAOS variant (50/194, 26%). Located in "
              "LRR domain. Marginal Crohn's disease association.",
    ),
    KnownVariant(
        gene="NOD2",
        rsid="rs2066844",
        name="R702W",
        chrom="chr16",
        pos=50712015,
        ref="C",
        alt="T",
        consequence="missense",
        amino_acid="p.Arg702Trp",
        amino_acid_mane="p.Arg675Trp",
        domain="NBD/HD1 junction",
        exon="exon 4",
        functional_effect="loss_of_function",
        gnomad_af=0.0399,
        yaos_frequency=0.20,
        diseases=["yao_syndrome", "crohn_disease"],
        clinvar_significance="risk_factor",
        notes="Third most common YAOS variant (38/194, 20%). Usually found "
              "as compound het with IVS8+158. Impairs NF-kB p65 phosphorylation "
              "(loss-of-function for NF-kB). Combined IVS8+158/R702W carriers "
              "show suppressed MDP-induced TNF-alpha. Also a major Crohn's "
              "disease risk allele.",
    ),
    KnownVariant(
        gene="NOD2",
        rsid="rs2066847",
        name="L1007fs",
        chrom="chr16",
        pos=50729870,
        ref="C",
        alt="CC",
        consequence="frameshift",
        amino_acid="p.Leu1007ProfsX2",
        amino_acid_mane="p.Leu980ProfsX2",
        domain="LRR (truncating)",
        exon="exon 11",
        functional_effect="loss_of_function",
        gnomad_af=0.0854,
        yaos_frequency=0.11,
        diseases=["yao_syndrome", "crohn_disease"],
        clinvar_significance="risk_factor",
        notes="Frameshift due to C duplication (3019dupC). Truncates LRR "
              "domain, abolishing MDP sensing. Found in 11% of YAOS patients "
              "(21/194). Strongest Crohn's disease risk allele when homozygous.",
    ),
    KnownVariant(
        gene="NOD2",
        rsid="rs2066845",
        name="G908R",
        chrom="chr16",
        pos=50722629,
        ref="G",
        alt="C",
        consequence="missense",
        amino_acid="p.Gly908Arg",
        amino_acid_mane="p.Gly881Arg",
        domain="LRR",
        exon="exon 8",
        functional_effect="loss_of_function",
        gnomad_af=0.0127,
        yaos_frequency=0.06,
        diseases=["yao_syndrome", "crohn_disease"],
        clinvar_significance="risk_factor",
        notes="Found in 6% of YAOS patients (11/194). Reduces MDP sensing "
              "in LRR domain. Major Crohn's disease risk allele.",
    ),
    KnownVariant(
        gene="NOD2",
        rsid="rs5743277",
        name="R703C",
        chrom="chr16",
        pos=50712018,
        ref="C",
        alt="T",
        consequence="missense",
        amino_acid="p.Arg703Cys",
        amino_acid_mane="p.Arg676Cys",
        domain="NBD",
        exon="exon 4",
        functional_effect="uncertain",
        gnomad_af=0.0030,
        yaos_frequency=0.05,
        diseases=["yao_syndrome"],
        clinvar_significance="uncertain_significance",
        notes="Found in 5% of YAOS patients (10/194). Adjacent to R702W "
              "in the NBD domain. Rare in general population.",
    ),
    KnownVariant(
        gene="NOD2",
        rsid="rs104895467",
        name="N852S",
        chrom="chr16",
        pos=50716899,
        ref="A",
        alt="G",
        consequence="missense",
        amino_acid="p.Asn852Ser",
        amino_acid_mane="p.Asn825Ser",
        domain="LRR",
        exon="exon 6",
        functional_effect="uncertain",
        gnomad_af=0.0010,
        yaos_frequency=0.03,
        diseases=["yao_syndrome"],
        clinvar_significance="uncertain_significance",
        notes="Found in compound het with IVS8+158 in YAOS patients.",
    ),
    KnownVariant(
        gene="NOD2",
        rsid="rs2066842",
        name="P268S",
        chrom="chr16",
        pos=50710713,
        ref="C",
        alt="T",
        consequence="missense",
        amino_acid="p.Pro268Ser",
        amino_acid_mane="p.Pro241Ser",
        domain="CARD2/NBD junction",
        exon="exon 4",
        functional_effect="benign_likely",
        gnomad_af=0.2500,
        yaos_frequency=0.0,
        diseases=["crohn_disease"],
        clinvar_significance="benign",
        notes="Common polymorphism, often in LD with R702W. Not independently "
              "associated with Yao syndrome. Included for completeness and "
              "compound het assessment.",
    ),
]


# ============================================================================
# NOD2 VARIANTS ASSOCIATED WITH BLAU SYNDROME (differential diagnosis)
# Gain-of-function variants in NACHT domain — AD inheritance
# ============================================================================

NOD2_BLAU_VARIANTS: list[KnownVariant] = [
    KnownVariant(
        gene="NOD2",
        rsid="rs104895462",
        name="R334W",
        chrom="chr16",
        pos=50710911,
        ref="C",
        alt="T",
        consequence="missense",
        amino_acid="p.Arg334Trp",
        amino_acid_mane="p.Arg307Trp",
        domain="NBD (NACHT)",
        exon="exon 4",
        functional_effect="gain_of_function",
        gnomad_af=0.0,
        yaos_frequency=0.0,
        diseases=["blau_syndrome", "early_onset_sarcoidosis"],
        clinvar_significance="pathogenic",
        notes="Most common Blau syndrome variant. Constitutive NF-kB "
              "activation. Disrupts autoinhibitory NACHT domain interface.",
    ),
    KnownVariant(
        gene="NOD2",
        rsid="rs104895463",
        name="R334Q",
        chrom="chr16",
        pos=50710912,
        ref="G",
        alt="A",
        consequence="missense",
        amino_acid="p.Arg334Gln",
        amino_acid_mane="p.Arg307Gln",
        domain="NBD (NACHT)",
        exon="exon 4",
        functional_effect="gain_of_function",
        gnomad_af=0.0,
        yaos_frequency=0.0,
        diseases=["blau_syndrome"],
        clinvar_significance="pathogenic",
        notes="Second most common Blau variant. Same residue as R334W.",
    ),
    KnownVariant(
        gene="NOD2",
        rsid="rs587778952",
        name="E383K",
        chrom="chr16",
        pos=50711058,
        ref="G",
        alt="A",
        consequence="missense",
        amino_acid="p.Glu383Lys",
        amino_acid_mane="p.Glu356Lys",
        domain="HD1",
        exon="exon 4",
        functional_effect="gain_of_function",
        gnomad_af=0.0,
        yaos_frequency=0.0,
        diseases=["blau_syndrome"],
        clinvar_significance="pathogenic",
        notes="Blau syndrome. Gain-of-function in helical domain 1.",
    ),
]


# ============================================================================
# MODIFIER GENE VARIANTS (found in digenic YAOS patients)
# ============================================================================

MODIFIER_VARIANTS: list[KnownVariant] = [
    # MEFV
    KnownVariant(
        gene="MEFV",
        rsid="rs3743930",
        name="E148Q",
        chrom="chr16",
        pos=3243800,
        ref="C",
        alt="G",
        consequence="missense",
        amino_acid="p.Glu148Gln",
        amino_acid_mane="p.Glu148Gln",
        domain="exon 2",
        exon="exon 2",
        functional_effect="low_penetrance",
        gnomad_af=0.0800,
        yaos_frequency=0.0,
        diseases=["FMF", "yao_syndrome_modifier"],
        clinvar_significance="uncertain_significance",
        notes="Most common MEFV variant found in digenic YAOS patients "
              "(NOD2+MEFV combination). Low-penetrance FMF variant.",
    ),
    KnownVariant(
        gene="MEFV",
        rsid="rs61752717",
        name="M694V",
        chrom="chr16",
        pos=3243631,
        ref="T",
        alt="C",
        consequence="missense",
        amino_acid="p.Met694Val",
        amino_acid_mane="p.Met694Val",
        domain="exon 10",
        exon="exon 10",
        functional_effect="gain_of_function",
        gnomad_af=0.0050,
        yaos_frequency=0.0,
        diseases=["FMF"],
        clinvar_significance="pathogenic",
        notes="Most severe FMF variant. May contribute to digenic disease "
              "when combined with NOD2 variants.",
    ),

    # NLRP3
    KnownVariant(
        gene="NLRP3",
        rsid="rs35829419",
        name="Q705K",
        chrom="chr1",
        pos=247422814,
        ref="A",
        alt="C",
        consequence="missense",
        amino_acid="p.Gln705Lys",
        amino_acid_mane="p.Gln705Lys",
        domain="LRR",
        exon="exon 6",
        functional_effect="low_penetrance",
        gnomad_af=0.0500,
        yaos_frequency=0.0,
        diseases=["CAPS_low_penetrance", "yao_syndrome_modifier"],
        clinvar_significance="uncertain_significance",
        notes="Low-penetrance CAPS variant. Found in digenic NOD2+NLRP3 "
              "YAOS patients. May enhance inflammasome activation.",
    ),

    # NLRP12
    KnownVariant(
        gene="NLRP12",
        rsid="rs199475867",
        name="F402L",
        chrom="chr19",
        pos=53804950,
        ref="T",
        alt="C",
        consequence="missense",
        amino_acid="p.Phe402Leu",
        amino_acid_mane="p.Phe402Leu",
        domain="NACHT",
        exon="exon 3",
        functional_effect="loss_of_function",
        gnomad_af=0.0010,
        yaos_frequency=0.0,
        diseases=["FCAS2", "yao_syndrome_modifier"],
        clinvar_significance="uncertain_significance",
        notes="Found in digenic NOD2+NLRP12 YAOS patients. NLRP12 is a "
              "negative regulator of NF-kB; loss-of-function may amplify "
              "NOD2-driven inflammation.",
    ),

    # TNFRSF1A
    KnownVariant(
        gene="TNFRSF1A",
        rsid="rs4149584",
        name="R92Q",
        chrom="chr12",
        pos=6335401,
        ref="G",
        alt="A",
        consequence="missense",
        amino_acid="p.Arg92Gln",
        amino_acid_mane="p.Arg92Gln",
        domain="CRD2",
        exon="exon 4",
        functional_effect="low_penetrance",
        gnomad_af=0.0200,
        yaos_frequency=0.0,
        diseases=["TRAPS_low_penetrance", "yao_syndrome_modifier"],
        clinvar_significance="uncertain_significance",
        notes="Low-penetrance TRAPS variant. Common modifier in digenic "
              "autoinflammatory disease.",
    ),
]


# ============================================================================
# COMPOUND HETEROZYGOSITY PATTERNS IN YAO SYNDROME
# ============================================================================
# 49% of YAOS patients carry >= 2 NOD2 variants

COMPOUND_HET_PATTERNS = [
    {
        "pattern": "IVS8+158 / R702W",
        "variants": ["rs5743289", "rs2066844"],
        "frequency": "most common compound het",
        "functional_note": "Mixed: IVS8+158 gain-of-function (p38/IL-6) + "
                           "R702W loss-of-function (NF-kB). Suppressed "
                           "MDP-induced TNF-alpha.",
    },
    {
        "pattern": "IVS8+158 / L1007fs",
        "variants": ["rs5743289", "rs2066847"],
        "frequency": "common",
        "functional_note": "IVS8+158 gain-of-function + L1007fs LRR "
                           "truncation (no MDP sensing). One functional "
                           "allele with elevated expression.",
    },
    {
        "pattern": "IVS8+158 / V955I",
        "variants": ["rs5743289", "rs5743291"],
        "frequency": "common",
        "functional_note": "Both in LRR/intronic region. Compound effect "
                           "on MDP recognition and signaling.",
    },
    {
        "pattern": "IVS8+158 / N852S",
        "variants": ["rs5743289", "rs104895467"],
        "frequency": "occasional",
        "functional_note": "Both affect LRR-proximal region.",
    },
    {
        "pattern": "IVS8+158 / G908R",
        "variants": ["rs5743289", "rs2066845"],
        "frequency": "occasional",
        "functional_note": "IVS8+158 gain-of-function + G908R LRR "
                           "loss-of-function.",
    },
    {
        "pattern": "R702W / G908R",
        "variants": ["rs2066844", "rs2066845"],
        "frequency": "rare in YAOS (more common in Crohn's)",
        "functional_note": "Double loss-of-function. Primarily Crohn's "
                           "disease risk, but some YAOS overlap.",
    },
]


# ============================================================================
# YAO SYNDROME DIAGNOSTIC CRITERIA (Yao 2017)
# For automated assessment against patient phenotype data
# ============================================================================

YAO_DIAGNOSTIC_CRITERIA = {
    "major_criteria": {
        "required_count": 2,
        "criteria": [
            {
                "id": "periodic",
                "description": "Periodic occurrence of symptoms (>= 2 episodes)",
            },
            {
                "id": "fever_or_dermatitis",
                "description": "Recurrent fever OR dermatitis (or both)",
            },
        ],
    },
    "minor_criteria": {
        "required_count": 1,
        "criteria": [
            {
                "id": "arthritis",
                "description": "Oligo/polyarthralgia, inflammatory arthritis, "
                               "or distal extremity swelling",
            },
            {
                "id": "GI_symptoms",
                "description": "Abdominal pain or diarrhea",
            },
            {
                "id": "sicca",
                "description": "Sicca-like symptoms (dry eyes/mouth)",
            },
            {
                "id": "serositis",
                "description": "Pericarditis or pleurisy",
            },
        ],
    },
    "molecular_criterion": {
        "description": "Specific NOD2 genotype (at least one Yao-associated "
                       "variant: IVS8+158, R702W, G908R, L1007fs, V955I, "
                       "R703C, or other cataloged YAOS variant)",
        "primary_variants": ["rs5743289", "rs2066844", "rs2066845",
                             "rs2066847", "rs5743291", "rs5743277"],
    },
    "exclusion_criteria": [
        "Crohn's disease (confirmed by endoscopy/histology)",
        "Blau syndrome (pediatric onset, granulomatous triad)",
        "Primary Sjogren's syndrome",
        "Other established monogenic SAID",
    ],
}


# ============================================================================
# Treatment response by genotype (for pharmacogenomic guidance)
# ============================================================================

TREATMENT_RESPONSE = {
    "glucocorticoid": {
        "overall_response": 0.75,
        "notes": "First-line. 75% respond, but many require chronic low-dose.",
    },
    "IL-1_inhibitor": {
        "overall_response": 0.65,
        "drugs": ["anakinra", "canakinumab", "rilonacept"],
        "best_responders": "IVS8+158 carriers (inflammasome-driven)",
        "notes": "Canakinumab showed 78% improvement. Best biologic option.",
    },
    "IL-6_inhibitor": {
        "overall_response": 0.33,
        "drugs": ["tocilizumab", "sarilumab"],
        "best_responders": "IVS8+158-only genotype (elevated basal IL-6)",
        "notes": "Tocilizumab 85% improvement in IVS8+158 subgroup.",
    },
    "JAK_inhibitor": {
        "overall_response": 0.75,
        "drugs": ["tofacitinib", "upadacitinib"],
        "notes": "Emerging. Upadacitinib + leflunomide reported effective.",
    },
    "hydroxychloroquine": {
        "overall_response": 0.625,
    },
    "anti-TNF": {
        "overall_response": 0.31,
        "drugs": ["adalimumab", "infliximab", "etanercept"],
        "notes": "Limited sustained benefit in YAOS.",
    },
    "colchicine": {
        "overall_response": 0.05,
        "notes": "Not effective for YAOS (unlike FMF).",
    },
}


# ============================================================================
# Helper: build lookup dictionaries
# ============================================================================

def get_all_yao_variants() -> list[KnownVariant]:
    """Return all NOD2 variants associated with Yao syndrome."""
    return list(NOD2_YAO_VARIANTS)


def get_all_known_variants() -> list[KnownVariant]:
    """Return all cataloged variants (YAOS + Blau + modifiers)."""
    return NOD2_YAO_VARIANTS + NOD2_BLAU_VARIANTS + MODIFIER_VARIANTS


def build_rsid_lookup() -> dict[str, KnownVariant]:
    """Build rsID -> KnownVariant lookup for fast matching."""
    return {v.rsid: v for v in get_all_known_variants()}


def build_position_lookup() -> dict[tuple[str, int], list[KnownVariant]]:
    """Build (chrom, pos) -> [KnownVariant] lookup for coordinate matching."""
    lookup: dict[tuple[str, int], list[KnownVariant]] = {}
    for v in get_all_known_variants():
        key = (v.chrom, v.pos)
        lookup.setdefault(key, []).append(v)
    return lookup
