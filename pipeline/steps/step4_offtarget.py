"""
Step 4 — Off-target analysis.

Integrates with Cas-OFFinder for genome-wide off-target prediction.
Falls back to heuristic analysis if Cas-OFFinder is not installed.
"""

import json
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from ..config import CAS_OFFINDER, MAX_OFFTARGET_MISMATCHES, REFERENCE_GENOME

log = logging.getLogger(__name__)


def _cas_offinder_available() -> bool:
    return shutil.which(CAS_OFFINDER) is not None


def run_cas_offinder(guides: dict[str, str], genome_path: str,
                     max_mm: int = MAX_OFFTARGET_MISMATCHES) -> list[dict]:
    """
    Run Cas-OFFinder and return a list of off-target hit dicts.
    *guides*: {name: 20bp_sequence}
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        inp = os.path.join(tmpdir, "input.txt")
        out = os.path.join(tmpdir, "output.txt")

        with open(inp, "w") as f:
            f.write(f"{genome_path}\n")
            f.write("NNNNNNNNNNNNNNNNNNNNNGG\n")
            for name, seq in guides.items():
                f.write(f"{seq}NNN {max_mm}\n")

        proc = subprocess.run(
            [CAS_OFFINDER, inp, "C", out],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            log.warning("Cas-OFFinder failed: %s", proc.stderr.strip())
            return []

        hits = []
        if os.path.isfile(out):
            with open(out) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        hits.append({
                            "pattern": parts[0],
                            "chromosome": parts[1],
                            "position": int(parts[2]),
                            "actual_sequence": parts[3],
                            "strand": parts[4],
                            "mismatches": int(parts[5]),
                        })
        return hits


def classify_hit(hit: dict) -> str:
    mm = hit["mismatches"]
    if mm == 0:
        return "ON_TARGET"
    if mm <= 1:
        return "HIGH_RISK"
    if mm <= 2:
        return "MODERATE_RISK"
    return "LOW_RISK"


def heuristic_offtarget_check(sequence: str) -> dict:
    """
    Heuristic off-target risk assessment when Cas-OFFinder is unavailable.

    Evaluates the PAM-proximal seed region (positions 1-12 from PAM),
    overall GC content, and sequence complexity.  The seed region is the
    primary determinant of Cas9 binding specificity -- mismatches in the
    seed are far less tolerated than in the PAM-distal region.

    NOT a substitute for genome-wide off-target search, but catches the
    highest-risk guides.
    """
    from ..utils.sequence_utils import gc_content, reverse_complement

    seq = sequence.upper()
    gc = gc_content(seq)
    risk = "LOW"
    notes = []

    # ── 1. Overall GC content ──
    if gc > 0.75:
        risk = "MODERATE"
        notes.append(f"High GC ({gc:.0%}) increases off-target binding stability")
    if gc < 0.30:
        risk = "MODERATE"
        notes.append(f"Low GC ({gc:.0%}) may reduce on-target specificity")

    # ── 2. Seed region analysis (PAM-proximal 12nt = last 12 bases) ──
    seed = seq[-12:]  # PAM-proximal seed
    seed_gc = gc_content(seed)

    # Low-complexity seed is the strongest off-target predictor
    if len(set(seed)) <= 2:
        risk = "HIGH"
        notes.append("Low-complexity seed region (≤2 distinct bases in 12nt)")
    elif len(set(seed)) == 3 and seed_gc > 0.80:
        risk = "MODERATE"
        notes.append("GC-biased seed with limited base diversity")

    # Seed homopolymers increase off-target risk
    for base in "ACGT":
        if base * 4 in seed:
            if risk != "HIGH":
                risk = "MODERATE"
            notes.append(f"Seed contains {base}×4 homopolymer")

    # ── 3. Repetitive motifs (dinucleotide repeats) ──
    for dinuc in ["AT", "TA", "GC", "CG", "GT", "TG", "AC", "CA"]:
        if (dinuc * 4) in seq:  # 4x dinucleotide repeat (8bp)
            if risk != "HIGH":
                risk = "MODERATE"
            notes.append(f"Contains {dinuc}×4 dinucleotide repeat")

    # ── 4. Self-complementarity in seed (potential for off-target at
    #        inverted repeat genomic sites) ──
    seed_rc = reverse_complement(seed)
    if seed_rc[:8] in seq:
        if risk != "HIGH":
            risk = "MODERATE"
        notes.append("Seed region is partially self-complementary")

    return {"risk": risk, "notes": "; ".join(notes) if notes else "Pass"}


def run(guide_designs: dict, output_dir: Path, patient_id: str) -> dict:
    """Run off-target analysis on all designed guides."""
    log.info("[Step 4] Off-target analysis")

    use_cas_offinder = _cas_offinder_available()
    if not use_cas_offinder:
        log.warning("  Cas-OFFinder not found — using heuristic fallback")

    # Collect all guide sequences
    all_guides: dict[str, str] = {}
    for gene, data in guide_designs.items():
        for i, g in enumerate(data.get("top_guides", [])):
            key = f"{gene}_guide{i}"
            all_guides[key] = g["sequence"]

    # Run analysis
    offtarget_results: dict = {}

    if use_cas_offinder and os.path.isfile(REFERENCE_GENOME):
        hits = run_cas_offinder(all_guides, REFERENCE_GENOME)
        for hit in hits:
            hit["risk_level"] = classify_hit(hit)

        # Aggregate per guide
        for key, seq in all_guides.items():
            guide_hits = [h for h in hits if seq[:15] in h.get("pattern", "")]
            high_risk = [h for h in guide_hits if h["risk_level"] == "HIGH_RISK"]
            offtarget_results[key] = {
                "total_hits": len(guide_hits),
                "high_risk_hits": len(high_risk),
                "status": "FAIL" if high_risk else "PASS",
            }
    else:
        # Heuristic fallback
        for key, seq in all_guides.items():
            check = heuristic_offtarget_check(seq)
            offtarget_results[key] = {
                "total_hits": None,
                "high_risk_hits": None,
                "status": "PASS_HEURISTIC" if check["risk"] == "LOW" else "REVIEW",
                "notes": check["notes"],
            }

    # Annotate guide designs with off-target results
    for gene, data in guide_designs.items():
        safe_guides = []
        for i, g in enumerate(data.get("top_guides", [])):
            key = f"{gene}_guide{i}"
            ot = offtarget_results.get(key, {})
            g["offtarget_status"] = ot.get("status", "UNKNOWN")
            g["offtarget_hits"] = ot.get("total_hits")
            if ot.get("status") in ("PASS", "PASS_HEURISTIC"):
                safe_guides.append(g)
        data["safe_guides"] = safe_guides
        data["guides_passing_offtarget"] = len(safe_guides)

    # Write outputs
    ot_dir = output_dir / "offtarget_analysis"
    ot_dir.mkdir(exist_ok=True)
    with open(ot_dir / "offtarget_results.json", "w") as f:
        json.dump(offtarget_results, f, indent=2)
    with open(ot_dir / "filtered_safe_guides.json", "w") as f:
        safe = {gene: data.get("safe_guides", [])
                for gene, data in guide_designs.items()}
        json.dump(safe, f, indent=2)

    log.info("  Total guides analysed: %d", len(all_guides))
    passing = sum(1 for v in offtarget_results.values()
                  if v["status"] in ("PASS", "PASS_HEURISTIC"))
    log.info("  Guides passing off-target filter: %d", passing)

    return guide_designs
