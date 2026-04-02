"""
Step 2b — Variant impact assessment via Ensembl VEP.

Runs Ensembl VEP (Variant Effect Predictor) on patient variants to annotate
coding consequences, impact severity, and functional predictions (SIFT,
PolyPhen, LoF).  Filters for high-impact missense variants in target genes.

This step slots between variant extraction (step 2) and guide design (step 3)
to inform which variants are functionally significant before designing
therapeutic constructs.

Requires: ensembl-vep (conda install -c bioconda ensembl-vep)
Falls back to a lightweight annotation heuristic if VEP is unavailable.
"""

import json
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from ..config import VEP_EXECUTABLE, VEP_CACHE_DIR

log = logging.getLogger(__name__)


def _vep_available() -> bool:
    return shutil.which(VEP_EXECUTABLE) is not None


def _extract_target_vcf(vcf_path: str, targets: dict,
                         output_vcf: str) -> bool:
    """
    Extract only variants in target gene regions from the full VCF.
    This keeps VEP runtime manageable for whole-genome inputs.
    """
    # Build a BED-like region list for bcftools view
    regions = []
    for gene, coords in targets.items():
        chrom = coords["chrom"]
        start = max(1, coords["start"] - 2000)  # include promoter
        end = coords["end"]
        regions.append(f"{chrom}:{start}-{end}")

    bcftools = shutil.which("bcftools")
    if not bcftools:
        log.warning("bcftools not available — passing full VCF to VEP")
        # Just copy/symlink the original
        if not os.path.exists(output_vcf):
            os.symlink(os.path.abspath(vcf_path), output_vcf)
        return True

    region_str = ",".join(regions)
    cmd = [bcftools, "view", "-r", region_str, "-O", "v",
           "-o", output_vcf, vcf_path]
    log.info("  Extracting target regions from VCF: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log.warning("bcftools view failed: %s", proc.stderr.strip())
        return False
    return True


def run_vep(vcf_path: str, output_vcf: str) -> bool:
    """
    Run Ensembl VEP on a VCF file.
    Returns True on success.
    """
    cmd = [
        VEP_EXECUTABLE,
        "--input_file", vcf_path,
        "--output_file", output_vcf,
        "--format", "vcf",
        "--vcf",
        "--offline",
        "--assembly", "GRCh38",
        "--symbol",
        "--canonical",
        "--protein",
        "--sift", "b",       # SIFT prediction + score
        "--polyphen", "b",   # PolyPhen prediction + score
        "--af",              # allele frequencies
        "--max_af",          # max population AF
        "--force_overwrite",
    ]

    # Use cache dir if configured
    if VEP_CACHE_DIR and os.path.isdir(VEP_CACHE_DIR):
        cmd.extend(["--dir_cache", VEP_CACHE_DIR])

    log.info("  Running VEP: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode != 0:
        log.error("VEP failed (exit %d): %s", proc.returncode,
                  proc.stderr.strip()[:500])
        return False

    log.info("  VEP completed successfully")
    return True


def parse_vep_output(vep_vcf_path: str, target_genes: list[str],
                     impact_filter: set[str] | None = None
                     ) -> list[dict]:
    """
    Parse VEP-annotated VCF and extract high-impact variants in target genes.

    VEP stores consequences in the INFO/CSQ field as pipe-delimited entries.
    """
    if impact_filter is None:
        impact_filter = {"HIGH", "MODERATE"}

    results = []
    csq_fields: list[str] = []

    with open(vep_vcf_path) as fh:
        for line in fh:
            # Parse the CSQ format header to know field positions
            if line.startswith("##INFO=<ID=CSQ"):
                # Format: ...Format: Allele|Consequence|IMPACT|SYMBOL|...">
                fmt_start = line.find("Format: ")
                if fmt_start != -1:
                    fmt_str = line[fmt_start + 8:].rstrip('">\n')
                    csq_fields = fmt_str.split("|")
                continue

            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            info = parts[7]

            # Extract CSQ annotations
            csq_entries = []
            for field in info.split(";"):
                if field.startswith("CSQ="):
                    csq_raw = field[4:]
                    csq_entries = csq_raw.split(",")
                    break

            if not csq_entries or not csq_fields:
                continue

            for entry in csq_entries:
                values = entry.split("|")
                if len(values) < len(csq_fields):
                    values.extend([""] * (len(csq_fields) - len(values)))

                annotation = dict(zip(csq_fields, values))

                symbol = annotation.get("SYMBOL", "")
                impact = annotation.get("IMPACT", "")
                consequence = annotation.get("Consequence", "")

                # Filter: must be in target genes and meet impact threshold
                if symbol not in target_genes:
                    continue
                if impact not in impact_filter:
                    continue

                results.append({
                    "chromosome": chrom,
                    "position": pos,
                    "reference": ref,
                    "alternative": alt,
                    "gene": symbol,
                    "consequence": consequence,
                    "impact": impact,
                    "sift": annotation.get("SIFT", ""),
                    "polyphen": annotation.get("PolyPhen", ""),
                    "protein_position": annotation.get("Protein_position", ""),
                    "amino_acids": annotation.get("Amino_acids", ""),
                    "codons": annotation.get("Codons", ""),
                    "canonical": annotation.get("CANONICAL", ""),
                    "max_af": annotation.get("MAX_AF", ""),
                })

    return results


def heuristic_impact_annotation(profile: dict,
                                target_genes: list[str]) -> list[dict]:
    """
    Lightweight fallback when VEP is unavailable.
    Flags variants based on type and quality alone — no functional prediction.
    """
    results = []

    for gene in target_genes:
        gene_data = profile.get(gene, {})
        for var in gene_data.get("coding_variants", []):
            # Simple heuristic: high-qual SNPs in coding regions
            impact = "UNKNOWN"
            if var.get("type") in ("snp", "mnp"):
                impact = "MODERATE" if var.get("qual", 0) > 30 else "LOW"
            elif var.get("type") in ("ins", "del", "complex"):
                impact = "HIGH" if var.get("qual", 0) > 20 else "MODERATE"

            if impact in ("HIGH", "MODERATE"):
                results.append({
                    "chromosome": gene_data.get("chrom", ""),
                    "position": var["pos"],
                    "reference": var["ref"],
                    "alternative": var["alt"][0] if var["alt"] else "",
                    "gene": gene,
                    "consequence": f"heuristic_{var.get('type', 'unknown')}",
                    "impact": impact,
                    "sift": "",
                    "polyphen": "",
                    "protein_position": "",
                    "amino_acids": "",
                    "codons": "",
                    "canonical": "",
                    "max_af": "",
                })

    return results


def run(vcf_path: str, targets: dict, profile: dict,
        output_dir: Path, patient_id: str) -> list[dict]:
    """
    Execute VEP annotation step.

    Returns a list of high-impact variant annotations.
    """
    log.info("[Step 2b] Variant impact assessment (VEP)")

    target_genes = list(targets.keys())

    if _vep_available():
        with tempfile.TemporaryDirectory() as tmpdir:
            # Extract target-region variants to reduce VEP input size
            subset_vcf = os.path.join(tmpdir, f"{patient_id}_targets.vcf")
            if not _extract_target_vcf(vcf_path, targets, subset_vcf):
                subset_vcf = vcf_path  # fall through to full VCF

            vep_output = str(output_dir / f"{patient_id}_vep_annotated.vcf")
            if run_vep(subset_vcf, vep_output):
                annotations = parse_vep_output(vep_output, target_genes)
            else:
                log.warning("  VEP failed — falling back to heuristic")
                annotations = heuristic_impact_annotation(profile,
                                                          target_genes)
    else:
        log.warning("  VEP not installed — using heuristic impact annotation")
        annotations = heuristic_impact_annotation(profile, target_genes)

    # Write output
    out_path = output_dir / f"{patient_id}_vep_high_impact.json"
    with open(out_path, "w") as f:
        json.dump(annotations, f, indent=2)

    log.info("  Found %d high-impact variants across %d target genes",
             len(annotations), len(set(a["gene"] for a in annotations)))
    log.info("  Wrote %s", out_path)

    # Per-gene summary
    gene_counts: dict[str, int] = {}
    for a in annotations:
        gene_counts[a["gene"]] = gene_counts.get(a["gene"], 0) + 1
    for gene, count in sorted(gene_counts.items()):
        log.info("    %s: %d high-impact variant(s)", gene, count)

    return annotations
