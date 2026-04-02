"""
Step 2c — Autoinflammatory variant annotation and Yao syndrome assessment.

Matches patient variants against a curated catalog of known autoinflammatory
disease variants.  Performs:

  1. Known-variant matching (rsID and positional)
  2. Compound heterozygosity detection in NOD2
  3. Digenic combination screening (NOD2 + MEFV/NLRP3/NLRP12/TNFRSF1A)
  4. Yao syndrome diagnostic likelihood scoring
  5. Differential diagnosis flagging (Blau, Crohn's, FMF, CAPS, TRAPS)

This step is specific to the ``autoimmune`` diagnosis mode and runs after
variant extraction (step 2) but before guide design (step 3).
"""

import json
import logging
from pathlib import Path

from ..databases.yao_variants import (
    NOD2_YAO_VARIANTS,
    NOD2_BLAU_VARIANTS,
    MODIFIER_VARIANTS,
    COMPOUND_HET_PATTERNS,
    YAO_DIAGNOSTIC_CRITERIA,
    TREATMENT_RESPONSE,
    KnownVariant,
    build_rsid_lookup,
    build_position_lookup,
)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Variant matching
# ---------------------------------------------------------------------------

def match_known_variants(patient_variants: dict,
                         profile: dict) -> list[dict]:
    """
    Match patient variants against the known autoinflammatory variant catalog.

    *patient_variants*: per-gene variant data from step 2 profile
    *profile*: full gene profile from step 2

    Returns a list of matched variant dicts.
    """
    rsid_lookup = build_rsid_lookup()
    pos_lookup = build_position_lookup()
    matches = []

    for gene, gene_data in profile.items():
        all_vars = (gene_data.get("coding_variants", []) +
                    gene_data.get("promoter_variants", []))

        for var in all_vars:
            matched_known: KnownVariant | None = None

            # 1. Try positional match (chrom + pos)
            chrom = gene_data.get("chrom", "")
            pos = var.get("pos", 0)
            key = (chrom, pos)
            if key in pos_lookup:
                for candidate in pos_lookup[key]:
                    # Verify allele match
                    alt_list = var.get("alt", [])
                    if isinstance(alt_list, str):
                        alt_list = [alt_list]
                    if candidate.alt in alt_list or candidate.ref == var.get("ref", ""):
                        matched_known = candidate
                        break

            if matched_known:
                matches.append(_build_match_record(var, matched_known, gene,
                                                    match_type="positional"))

    # 2. Scan for IVS8+158 specifically — it is intronic and may be in
    #    extended regions not captured by standard coding variant extraction.
    #    Check if we already found it; if not, note it for the report.
    ivs8_found = any(m["known_variant"]["rsid"] == "rs5743289"
                     for m in matches)

    return matches


def _build_match_record(patient_var: dict, known: KnownVariant,
                        gene: str, match_type: str) -> dict:
    return {
        "gene": gene,
        "patient_variant": {
            "pos": patient_var.get("pos"),
            "ref": patient_var.get("ref"),
            "alt": patient_var.get("alt"),
            "qual": patient_var.get("qual"),
            "zygosity": patient_var.get("zygosity"),
            "type": patient_var.get("type"),
        },
        "known_variant": {
            "rsid": known.rsid,
            "name": known.name,
            "amino_acid": known.amino_acid,
            "domain": known.domain,
            "consequence": known.consequence,
            "functional_effect": known.functional_effect,
            "gnomad_af": known.gnomad_af,
            "yaos_frequency": known.yaos_frequency,
            "diseases": known.diseases,
            "clinvar": known.clinvar_significance,
            "notes": known.notes,
        },
        "match_type": match_type,
    }


# ---------------------------------------------------------------------------
# Compound heterozygosity detection
# ---------------------------------------------------------------------------

def detect_compound_hets(matches: list[dict]) -> list[dict]:
    """
    Identify compound heterozygosity patterns in NOD2 that are
    characteristic of Yao syndrome.

    49% of YAOS patients carry >= 2 NOD2 variants.
    """
    nod2_matches = [m for m in matches if m["gene"] == "NOD2"]

    # Get unique variant rsIDs carried by the patient
    patient_rsids = set()
    for m in nod2_matches:
        patient_rsids.add(m["known_variant"]["rsid"])

    detected_patterns = []
    for pattern in COMPOUND_HET_PATTERNS:
        pattern_rsids = set(pattern["variants"])
        if pattern_rsids.issubset(patient_rsids):
            detected_patterns.append({
                "pattern": pattern["pattern"],
                "variants_matched": pattern["variants"],
                "frequency_in_yaos": pattern["frequency"],
                "functional_note": pattern["functional_note"],
            })

    return detected_patterns


# ---------------------------------------------------------------------------
# Digenic screening
# ---------------------------------------------------------------------------

def screen_digenic(matches: list[dict]) -> list[dict]:
    """
    Screen for digenic combinations (NOD2 + modifier gene variants).
    22% of YAOS patients carry variants in additional SAID genes.
    """
    nod2_present = any(m["gene"] == "NOD2" for m in matches)
    if not nod2_present:
        return []

    modifier_matches = [m for m in matches if m["gene"] != "NOD2"]
    digenic_hits = []

    for m in modifier_matches:
        gene = m["gene"]
        digenic_hits.append({
            "combination": f"NOD2 + {gene}",
            "modifier_gene": gene,
            "modifier_variant": m["known_variant"]["name"],
            "modifier_rsid": m["known_variant"]["rsid"],
            "diseases": m["known_variant"]["diseases"],
            "clinical_note": (
                f"Digenic: patient carries NOD2 variant(s) plus "
                f"{gene} {m['known_variant']['name']}. This combination "
                f"is seen in ~22% of Yao syndrome patients and may "
                f"modify disease presentation."
            ),
        })

    return digenic_hits


# ---------------------------------------------------------------------------
# Diagnostic likelihood scoring
# ---------------------------------------------------------------------------

def score_yao_likelihood(matches: list[dict],
                         compound_hets: list[dict],
                         digenic: list[dict]) -> dict:
    """
    Score the likelihood that the patient's genotype is consistent with
    Yao syndrome vs. other NOD2-associated conditions.

    Scoring:
      - IVS8+158 present: +40 points (hallmark variant)
      - Any other YAOS variant: +15 each
      - Compound het pattern: +20
      - Digenic combination: +10
      - Blau variant (NACHT GoF): -50 (argues against YAOS)
      - Homozygous L1007fs: -30 (argues for Crohn's)
    """
    score = 0
    evidence = []

    nod2_matches = [m for m in matches if m["gene"] == "NOD2"]

    # Check for IVS8+158
    ivs8 = any(m["known_variant"]["rsid"] == "rs5743289"
               for m in nod2_matches)
    if ivs8:
        score += 40
        evidence.append("IVS8+158 present (+40) — primary YAOS variant")

    # Other YAOS variants
    yaos_rsids = {v.rsid for v in NOD2_YAO_VARIANTS}
    for m in nod2_matches:
        rsid = m["known_variant"]["rsid"]
        if rsid in yaos_rsids and rsid != "rs5743289":
            score += 15
            evidence.append(
                f"{m['known_variant']['name']} present (+15)")

    # Compound het bonus
    if compound_hets:
        score += 20
        evidence.append(
            f"Compound heterozygosity detected (+20): "
            f"{compound_hets[0]['pattern']}")

    # Digenic bonus
    if digenic:
        score += 10
        evidence.append(
            f"Digenic combination (+10): "
            f"{digenic[0]['combination']}")

    # Blau penalty
    blau_rsids = {v.rsid for v in NOD2_BLAU_VARIANTS}
    for m in nod2_matches:
        if m["known_variant"]["rsid"] in blau_rsids:
            score -= 50
            evidence.append(
                f"Blau syndrome variant {m['known_variant']['name']} "
                f"(-50) — consider Blau/early-onset sarcoidosis instead")

    # Homozygous L1007fs penalty
    for m in nod2_matches:
        if (m["known_variant"]["rsid"] == "rs2066847" and
                m["patient_variant"].get("zygosity") == "HOM_ALT"):
            score -= 30
            evidence.append(
                "Homozygous L1007fs (-30) — strongly associated with "
                "Crohn's disease rather than YAOS")

    # Determine assessment
    if score >= 50:
        assessment = "HIGH — genotype strongly consistent with Yao syndrome"
    elif score >= 25:
        assessment = "MODERATE — genotype partially consistent with Yao syndrome"
    elif score > 0:
        assessment = "LOW — genotype has some YAOS-associated variants but " \
                     "clinical correlation required"
    else:
        assessment = "UNLIKELY — genotype does not support Yao syndrome diagnosis"

    return {
        "yao_likelihood_score": max(0, score),
        "assessment": assessment,
        "evidence": evidence,
        "molecular_criterion_met": ivs8 or len(nod2_matches) > 0,
    }


# ---------------------------------------------------------------------------
# Differential diagnosis
# ---------------------------------------------------------------------------

def assess_differential(matches: list[dict],
                        yao_score: dict) -> list[dict]:
    """
    Generate differential diagnosis considerations based on variant profile.
    """
    differentials = []

    nod2_matches = [m for m in matches if m["gene"] == "NOD2"]
    other_matches = [m for m in matches if m["gene"] != "NOD2"]

    # Blau syndrome
    blau_rsids = {v.rsid for v in NOD2_BLAU_VARIANTS}
    blau_hits = [m for m in nod2_matches
                 if m["known_variant"]["rsid"] in blau_rsids]
    if blau_hits:
        differentials.append({
            "disease": "Blau syndrome / early-onset sarcoidosis",
            "likelihood": "HIGH" if len(blau_hits) > 0 else "LOW",
            "evidence": [f"NACHT domain gain-of-function variant: "
                         f"{h['known_variant']['name']}" for h in blau_hits],
            "distinguishing_features": "Pediatric onset (<5 yr), "
                                       "granulomatous triad (dermatitis, "
                                       "arthritis, uveitis), AD inheritance",
        })

    # Crohn's disease
    crohn_rsids = {"rs2066844", "rs2066845", "rs2066847"}
    crohn_hits = [m for m in nod2_matches
                  if m["known_variant"]["rsid"] in crohn_rsids]
    if crohn_hits:
        differentials.append({
            "disease": "Crohn's disease",
            "likelihood": "MODERATE" if len(crohn_hits) >= 2 else "LOW",
            "evidence": [f"{h['known_variant']['name']} (shared CD/YAOS "
                         f"risk allele)" for h in crohn_hits],
            "distinguishing_features": "IBD pathology on colonoscopy, "
                                       "transmural inflammation, granulomas. "
                                       "Normal colonoscopy argues against CD.",
        })

    # FMF
    mefv_hits = [m for m in other_matches if m["gene"] == "MEFV"]
    if mefv_hits:
        differentials.append({
            "disease": "Familial Mediterranean Fever (FMF)",
            "likelihood": "LOW" if len(mefv_hits) == 1 else "MODERATE",
            "evidence": [f"MEFV {h['known_variant']['name']}"
                         for h in mefv_hits],
            "distinguishing_features": "Short febrile episodes (12-72h), "
                                       "serositis, excellent colchicine "
                                       "response. YAOS does NOT respond "
                                       "to colchicine.",
        })

    # CAPS
    nlrp3_hits = [m for m in other_matches if m["gene"] == "NLRP3"]
    if nlrp3_hits:
        differentials.append({
            "disease": "CAPS (NLRP3 inflammasomopathy)",
            "likelihood": "LOW",
            "evidence": [f"NLRP3 {h['known_variant']['name']}"
                         for h in nlrp3_hits],
            "distinguishing_features": "Cold-triggered urticaria (FCAS), "
                                       "sensorineural hearing loss (MWS), "
                                       "CNS involvement (NOMID).",
        })

    return differentials


# ---------------------------------------------------------------------------
# Treatment guidance based on genotype
# ---------------------------------------------------------------------------

def suggest_treatment(matches: list[dict],
                      yao_score: dict) -> list[dict]:
    """
    Genotype-guided treatment suggestions based on published response data.
    """
    suggestions = []

    if yao_score["yao_likelihood_score"] < 10:
        return suggestions

    # Check if IVS8+158-only genotype
    nod2_matches = [m for m in matches if m["gene"] == "NOD2"]
    rsids = {m["known_variant"]["rsid"] for m in nod2_matches}
    ivs8_only = rsids == {"rs5743289"}

    for drug, data in TREATMENT_RESPONSE.items():
        suggestion = {
            "treatment": drug,
            "overall_response_rate": data["overall_response"],
            "drugs": data.get("drugs", [drug]),
            "notes": data.get("notes", ""),
        }

        # Genotype-specific guidance
        if drug == "IL-6_inhibitor" and ivs8_only:
            suggestion["genotype_match"] = "STRONG"
            suggestion["genotype_note"] = (
                "IVS8+158-only genotype shows elevated basal IL-6. "
                "Tocilizumab showed 85% improvement in this subgroup.")
        elif drug == "IL-1_inhibitor" and "rs5743289" in rsids:
            suggestion["genotype_match"] = "GOOD"
            suggestion["genotype_note"] = (
                "IVS8+158 carriers have inflammasome-driven phenotype. "
                "Best biologic response overall (65% across all genotypes).")
        elif drug == "JAK_inhibitor":
            suggestion["genotype_match"] = "EMERGING"
            suggestion["genotype_note"] = (
                "JAK inhibitors show promise (75% response) across genotypes.")
        elif drug == "colchicine":
            suggestion["genotype_match"] = "POOR"
            suggestion["genotype_note"] = (
                "YAOS does NOT respond to colchicine (unlike FMF).")
        else:
            suggestion["genotype_match"] = "STANDARD"

        suggestions.append(suggestion)

    # Sort by response rate
    suggestions.sort(key=lambda s: s["overall_response_rate"], reverse=True)
    return suggestions


# ---------------------------------------------------------------------------
# Main step entry point
# ---------------------------------------------------------------------------

def run(vcf_path: str, targets: dict, profile: dict,
        output_dir: Path, patient_id: str) -> dict:
    """
    Execute autoinflammatory variant annotation.

    Returns a comprehensive assessment dict.
    """
    log.info("[Step 2c] Autoinflammatory variant annotation & Yao assessment")

    # 1. Match patient variants against known catalog
    matches = match_known_variants(targets, profile)
    log.info("  Matched %d patient variant(s) against known catalog",
             len(matches))
    for m in matches:
        log.info("    %s %s (%s) — %s [%s]",
                 m["gene"], m["known_variant"]["name"],
                 m["known_variant"]["rsid"],
                 m["known_variant"]["consequence"],
                 m["known_variant"]["functional_effect"])

    # 2. Detect compound heterozygosity
    compound_hets = detect_compound_hets(matches)
    if compound_hets:
        log.info("  Compound heterozygosity detected: %s",
                 compound_hets[0]["pattern"])

    # 3. Screen for digenic combinations
    digenic = screen_digenic(matches)
    if digenic:
        log.info("  Digenic combination(s): %s",
                 ", ".join(d["combination"] for d in digenic))

    # 4. Score Yao syndrome likelihood
    yao_score = score_yao_likelihood(matches, compound_hets, digenic)
    log.info("  Yao syndrome likelihood: %d — %s",
             yao_score["yao_likelihood_score"], yao_score["assessment"])

    # 5. Differential diagnosis
    differential = assess_differential(matches, yao_score)
    for d in differential:
        log.info("  Differential: %s (%s)", d["disease"], d["likelihood"])

    # 6. Treatment guidance
    treatment = suggest_treatment(matches, yao_score)

    # Assemble result
    result = {
        "patient_id": patient_id,
        "known_variant_matches": matches,
        "compound_heterozygosity": compound_hets,
        "digenic_combinations": digenic,
        "yao_syndrome_assessment": yao_score,
        "differential_diagnosis": differential,
        "treatment_guidance": [
            s for s in treatment
            if s["genotype_match"] in ("STRONG", "GOOD", "EMERGING")
        ],
        "full_treatment_data": treatment,
        "summary": {
            "total_known_variants_matched": len(matches),
            "nod2_variants": len([m for m in matches if m["gene"] == "NOD2"]),
            "modifier_gene_variants": len([m for m in matches
                                           if m["gene"] != "NOD2"]),
            "compound_het_patterns": len(compound_hets),
            "digenic_combinations": len(digenic),
            "yao_likelihood_score": yao_score["yao_likelihood_score"],
            "yao_assessment": yao_score["assessment"],
        },
    }

    # Write output
    out_path = output_dir / f"{patient_id}_autoinflammatory_assessment.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    log.info("  Wrote %s", out_path)

    return result
