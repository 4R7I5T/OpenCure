"""
Step 6 — Generate pipeline summary report.

Produces a JSON summary and a human-readable text report.
"""

import json
import logging
from datetime import datetime
from pathlib import Path

log = logging.getLogger(__name__)


def run(patient_id: str, diagnosis: str, output_dir: Path,
        validation: dict, strategy: dict, guide_designs: dict,
        constructs: list[dict],
        vep_annotations: list[dict] | None = None,
        autoinflammatory: dict | None = None) -> Path:
    """Generate final pipeline report."""
    log.info("[Step 6] Generating pipeline report")

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    vep_annotations = vep_annotations or []

    # Count statistics
    total_genes = len(guide_designs)
    total_guides = sum(
        len(d.get("top_guides", [])) for d in guide_designs.values()
    )
    safe_guides = sum(
        d.get("guides_passing_offtarget", 0) for d in guide_designs.values()
    )
    genes_with_variants = sum(
        1 for d in guide_designs.values()
        if d.get("patient_variants_in_region", 0) > 0
    )
    oversized = [c for c in constructs if "warning" in c]
    dreadd_constructs = [c for c in constructs
                         if c.get("construct_type") == "DREADD_insertion"]
    crispr_constructs = [c for c in constructs
                         if c.get("construct_type") != "DREADD_insertion"]

    summary = {
        "patient_id": patient_id,
        "diagnosis_type": diagnosis,
        "pipeline_version": "3.0.0",
        "timestamp": timestamp,
        "validation": {
            "vcf_valid": validation.get("vcf", {}).get("valid", False),
            "bam_valid": validation.get("bam", {}).get("valid", False),
            "errors": validation.get("errors", []),
        },
        "analysis": {
            "target_genes_analysed": total_genes,
            "genes_with_patient_variants": genes_with_variants,
            "total_candidate_guides": total_guides,
            "guides_passing_offtarget": safe_guides,
            "constructs_designed": len(constructs),
            "crispr_constructs": len(crispr_constructs),
            "dreadd_constructs": len(dreadd_constructs),
            "constructs_exceeding_aav_limit": len(oversized),
        },
        "per_gene_summary": {},
    }

    # VEP summary
    if vep_annotations:
        vep_genes: dict[str, int] = {}
        for ann in vep_annotations:
            vep_genes[ann["gene"]] = vep_genes.get(ann["gene"], 0) + 1
        summary["vep_annotation"] = {
            "total_high_impact_variants": len(vep_annotations),
            "genes_with_high_impact": len(vep_genes),
            "per_gene_counts": vep_genes,
        }

    for gene, data in guide_designs.items():
        strat = strategy.get(gene, {})
        construct = next((c for c in constructs if c["gene"] == gene), None)
        gene_summary: dict = {
            "role": strat.get("role", data.get("purpose", "")),
            "approach": strat.get("approach", ""),
            "patient_variants": data.get("patient_variants_in_region", 0),
            "top_guide_score": (data["top_guides"][0]["score"]
                                if data.get("top_guides") else None),
            "safe_guides": data.get("guides_passing_offtarget", 0),
            "construct_type": construct["construct_type"] if construct else None,
            "construct_size_bp": construct["total_size_bp"] if construct else None,
        }
        # Add VEP high-impact count for this gene
        gene_vep = [a for a in vep_annotations if a["gene"] == gene]
        if gene_vep:
            gene_summary["vep_high_impact_variants"] = len(gene_vep)
        summary["per_gene_summary"][gene] = gene_summary

    # DREADD summary
    if dreadd_constructs:
        summary["dreadd_summary"] = {
            "constructs": len(dreadd_constructs),
            "activation_ligand": "deschloroclozapine (DCZ)",
            "delivery_method": "stereotactic AAV9 injection",
            "details": [
                {
                    "name": c["name"],
                    "receptor": c.get("receptor", ""),
                    "target_regions": c.get("target_regions", []),
                    "purpose": c.get("purpose", ""),
                    "size_bp": c["total_size_bp"],
                }
                for c in dreadd_constructs
            ],
        }

    # Autoinflammatory / Yao syndrome assessment
    if autoinflammatory:
        summary["autoinflammatory_assessment"] = {
            "yao_likelihood_score": autoinflammatory.get(
                "yao_syndrome_assessment", {}).get("yao_likelihood_score", 0),
            "yao_assessment": autoinflammatory.get(
                "yao_syndrome_assessment", {}).get("assessment", ""),
            "known_variants_matched": autoinflammatory.get(
                "summary", {}).get("total_known_variants_matched", 0),
            "nod2_variants": autoinflammatory.get(
                "summary", {}).get("nod2_variants", 0),
            "compound_het_patterns": autoinflammatory.get(
                "summary", {}).get("compound_het_patterns", 0),
            "digenic_combinations": autoinflammatory.get(
                "summary", {}).get("digenic_combinations", 0),
            "differential_diagnoses": [
                d["disease"] for d in
                autoinflammatory.get("differential_diagnosis", [])
            ],
        }

    # Write JSON summary
    summary_path = output_dir / "pipeline_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # Write human-readable report
    report_path = output_dir / "pipeline_report.txt"
    with open(report_path, "w") as f:
        f.write("=" * 64 + "\n")
        f.write("  OpenCure Pipeline Report\n")
        f.write("=" * 64 + "\n")
        f.write(f"  Patient:    {patient_id}\n")
        f.write(f"  Diagnosis:  {diagnosis}\n")
        f.write(f"  Date:       {timestamp}\n")
        f.write(f"  Pipeline:   v3.0.0\n")
        f.write("=" * 64 + "\n\n")

        f.write("INPUT VALIDATION\n")
        f.write("-" * 40 + "\n")
        f.write(f"  VCF valid:  {validation.get('vcf', {}).get('valid')}\n")
        f.write(f"  BAM valid:  {validation.get('bam', {}).get('valid')}\n")
        if validation.get("errors"):
            for err in validation["errors"]:
                f.write(f"  ERROR: {err}\n")
        f.write("\n")

        # VEP section
        if vep_annotations:
            f.write("VARIANT IMPACT ASSESSMENT (VEP)\n")
            f.write("-" * 40 + "\n")
            f.write(f"  High-impact variants:  {len(vep_annotations)}\n")
            vep_genes_set = set(a["gene"] for a in vep_annotations)
            f.write(f"  Genes affected:        {len(vep_genes_set)}\n")
            for ann in vep_annotations[:10]:
                sift_str = f" SIFT={ann['sift']}" if ann.get("sift") else ""
                pphen_str = (f" PolyPhen={ann['polyphen']}"
                             if ann.get("polyphen") else "")
                f.write(f"    {ann['gene']} {ann['chromosome']}:"
                        f"{ann['position']} {ann['reference']}>"
                        f"{ann['alternative']} [{ann['impact']}] "
                        f"{ann['consequence']}{sift_str}{pphen_str}\n")
            if len(vep_annotations) > 10:
                f.write(f"    ... and {len(vep_annotations) - 10} more\n")
            f.write("\n")

        f.write("ANALYSIS SUMMARY\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Target genes:             {total_genes}\n")
        f.write(f"  Genes w/ patient variants: {genes_with_variants}\n")
        f.write(f"  Total candidate guides:    {total_guides}\n")
        f.write(f"  Guides passing off-target: {safe_guides}\n")
        f.write(f"  CRISPR constructs:         {len(crispr_constructs)}\n")
        if dreadd_constructs:
            f.write(f"  DREADD constructs:         {len(dreadd_constructs)}\n")
        f.write(f"  Total constructs:          {len(constructs)}\n")
        if oversized:
            f.write(f"  Oversized constructs:      {len(oversized)} "
                    f"(need dual-AAV)\n")
        f.write("\n")

        f.write("PER-GENE RESULTS\n")
        f.write("-" * 40 + "\n")
        for gene, gs in summary["per_gene_summary"].items():
            f.write(f"\n  {gene}\n")
            f.write(f"    Role:           {gs['role'][:60]}\n")
            f.write(f"    Approach:       {gs['approach']}\n")
            f.write(f"    Patient vars:   {gs['patient_variants']}\n")
            if gs.get("vep_high_impact_variants"):
                f.write(f"    VEP high-impact: "
                        f"{gs['vep_high_impact_variants']}\n")
            f.write(f"    Top guide score: {gs['top_guide_score']}\n")
            f.write(f"    Safe guides:    {gs['safe_guides']}\n")
            if gs["construct_type"]:
                f.write(f"    Construct:      {gs['construct_type']} "
                        f"({gs['construct_size_bp']} bp)\n")

        # DREADD section
        if dreadd_constructs:
            f.write("\nDREADD CHEMOGENETIC CONSTRUCTS\n")
            f.write("-" * 40 + "\n")
            f.write("  Activation ligand: deschloroclozapine (DCZ)\n")
            f.write("  Delivery: stereotactic AAV9 injection\n")
            f.write("  Estimated onset: 2-4 weeks post-injection\n\n")
            for c in dreadd_constructs:
                f.write(f"  {c['name']}\n")
                f.write(f"    Receptor:   {c.get('receptor', '')}\n")
                f.write(f"    Promoter:   {c.get('promoter', '')}\n")
                f.write(f"    Regions:    "
                        f"{', '.join(c.get('target_regions', []))}\n")
                f.write(f"    Size:       {c['total_size_bp']} bp\n")
                f.write(f"    Purpose:    {c.get('purpose', '')[:70]}\n")
                if c.get("insertion_site"):
                    site = c["insertion_site"]
                    f.write(f"    Safe harbor: {site.get('locus', '')[:50]}\n")
                f.write("\n")

        # Autoinflammatory / Yao syndrome section
        if autoinflammatory:
            yao = autoinflammatory.get("yao_syndrome_assessment", {})
            f.write("\nAUTOINFLAMMATORY / YAO SYNDROME ASSESSMENT\n")
            f.write("-" * 40 + "\n")
            f.write(f"  Yao likelihood score:  "
                    f"{yao.get('yao_likelihood_score', 0)}\n")
            f.write(f"  Assessment:            "
                    f"{yao.get('assessment', 'N/A')}\n")
            f.write(f"  Molecular criterion:   "
                    f"{'MET' if yao.get('molecular_criterion_met') else 'NOT MET'}\n")

            # Evidence
            evidence = yao.get("evidence", [])
            if evidence:
                f.write("\n  Evidence:\n")
                for e in evidence:
                    f.write(f"    - {e}\n")

            # Known variant matches
            matches = autoinflammatory.get("known_variant_matches", [])
            if matches:
                f.write(f"\n  Known variants matched: {len(matches)}\n")
                for m in matches:
                    kv = m["known_variant"]
                    pv = m["patient_variant"]
                    f.write(f"    {m['gene']} {kv['name']} ({kv['rsid']}) "
                            f"[{kv['consequence']}] "
                            f"zyg={pv.get('zygosity', '?')} "
                            f"gnomAD={kv['gnomad_af']:.4f}\n")
                    if kv.get("notes"):
                        note = kv["notes"][:100]
                        f.write(f"      {note}\n")

            # Compound het
            comp = autoinflammatory.get("compound_heterozygosity", [])
            if comp:
                f.write(f"\n  Compound heterozygosity:\n")
                for c in comp:
                    f.write(f"    Pattern: {c['pattern']}\n")
                    f.write(f"    Freq:    {c['frequency_in_yaos']}\n")
                    f.write(f"    Effect:  {c['functional_note'][:80]}\n")

            # Digenic
            dig = autoinflammatory.get("digenic_combinations", [])
            if dig:
                f.write(f"\n  Digenic combinations:\n")
                for d in dig:
                    f.write(f"    {d['combination']}: "
                            f"{d['modifier_variant']} ({d['modifier_rsid']})\n")

            # Differential diagnosis
            diff = autoinflammatory.get("differential_diagnosis", [])
            if diff:
                f.write(f"\n  Differential diagnosis:\n")
                for d in diff:
                    f.write(f"    {d['disease']} — {d['likelihood']}\n")
                    for e in d.get("evidence", []):
                        f.write(f"      {e}\n")
                    f.write(f"      Distinguish: "
                            f"{d.get('distinguishing_features', '')[:70]}\n")

            # Treatment guidance
            tx = autoinflammatory.get("treatment_guidance", [])
            if tx:
                f.write(f"\n  Genotype-guided treatment:\n")
                for t in tx:
                    match = t.get("genotype_match", "")
                    drugs = ", ".join(t.get("drugs", [t["treatment"]]))
                    f.write(f"    [{match}] {drugs} "
                            f"(response: {t['overall_response_rate']:.0%})\n")
                    if t.get("genotype_note"):
                        f.write(f"      {t['genotype_note'][:80]}\n")

            f.write("\n")

        f.write("=" * 64 + "\n")
        f.write("  Output directory: " + str(output_dir) + "\n")
        f.write("=" * 64 + "\n")

    log.info("  Wrote %s", summary_path)
    log.info("  Wrote %s", report_path)

    return report_path
