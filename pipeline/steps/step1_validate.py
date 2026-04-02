"""
Step 1 — Validate patient input files (VCF, BAM, BAI).

Checks file integrity, verifies index pairing, and reports coverage
at target gene loci.
"""

import logging
import os
import shutil
import subprocess

from ..config import SAMTOOLS, BCFTOOLS

log = logging.getLogger(__name__)


def _tool_available(name: str) -> bool:
    return shutil.which(name) is not None


def _run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    log.info("Running: %s", " ".join(cmd))
    return subprocess.run(cmd, capture_output=True, text=True, check=check)


def validate_vcf(vcf_path: str) -> dict:
    """Basic VCF validation."""
    result = {"path": vcf_path, "valid": False, "samples": [], "error": None}

    if not os.path.isfile(vcf_path):
        result["error"] = f"VCF file not found: {vcf_path}"
        return result

    # Check we can read the first non-comment line
    import gzip
    opener = gzip.open if vcf_path.endswith(".gz") else open
    mode = "rt" if vcf_path.endswith(".gz") else "r"
    header_found = False
    try:
        with opener(vcf_path, mode) as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header_found = True
                    cols = line.strip().split("\t")
                    if len(cols) > 9:
                        result["samples"] = cols[9:]
                    break
        if not header_found:
            result["error"] = "No #CHROM header line found in VCF"
            return result
    except Exception as exc:
        result["error"] = f"Cannot read VCF: {exc}"
        return result

    # bcftools stats if available
    if _tool_available(BCFTOOLS):
        proc = _run([BCFTOOLS, "stats", vcf_path], check=False)
        if proc.returncode == 0:
            for line in proc.stdout.splitlines():
                if line.startswith("SN") and "number of records" in line:
                    result["record_count"] = line.strip().split("\t")[-1]

    result["valid"] = True
    return result


def validate_bam(bam_path: str, bai_path: str | None = None) -> dict:
    """Validate BAM file and its index."""
    result = {"path": bam_path, "valid": False, "bai": bai_path, "error": None}

    if not os.path.isfile(bam_path):
        result["error"] = f"BAM file not found: {bam_path}"
        return result

    # Look for index
    if bai_path and os.path.isfile(bai_path):
        pass
    elif os.path.isfile(bam_path + ".bai"):
        bai_path = bam_path + ".bai"
    elif os.path.isfile(bam_path.replace(".bam", ".bai")):
        bai_path = bam_path.replace(".bam", ".bai")
    else:
        bai_path = None

    result["bai"] = bai_path

    if not bai_path:
        result["error"] = (
            "BAM index (.bai) not found. Run: samtools index <bam>"
        )
        return result

    # samtools quickcheck if available
    if _tool_available(SAMTOOLS):
        proc = _run([SAMTOOLS, "quickcheck", bam_path], check=False)
        if proc.returncode != 0:
            result["error"] = "BAM file failed samtools quickcheck"
            return result

    result["valid"] = True
    return result


def check_coverage(bam_path: str, chrom: str, start: int, end: int,
                   gene: str = "") -> dict:
    """Return mean coverage over a region using samtools depth."""
    cov = {"gene": gene, "region": f"{chrom}:{start}-{end}", "mean_coverage": None}

    if not _tool_available(SAMTOOLS):
        cov["error"] = "samtools not available"
        return cov

    proc = _run(
        [SAMTOOLS, "depth", "-r", f"{chrom}:{start}-{end}", bam_path],
        check=False,
    )
    if proc.returncode != 0:
        cov["error"] = proc.stderr.strip()
        return cov

    total = 0
    count = 0
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            total += int(parts[2])
            count += 1

    cov["mean_coverage"] = round(total / count, 1) if count else 0.0
    return cov


def run(vcf_path: str, bam_path: str, bai_path: str | None,
        targets: dict) -> dict:
    """
    Full validation step.  Returns a dict with vcf/bam validation results
    and per-gene coverage.
    """
    log.info("[Step 1] Validating input files")

    vcf_result = validate_vcf(vcf_path)
    bam_result = validate_bam(bam_path, bai_path)

    errors = []
    if not vcf_result["valid"]:
        errors.append(f"VCF: {vcf_result['error']}")
    if not bam_result["valid"]:
        errors.append(f"BAM: {bam_result['error']}")

    coverage = {}
    if bam_result["valid"]:
        for gene, coords in targets.items():
            cov = check_coverage(
                bam_path, coords["chrom"], coords["start"], coords["end"],
                gene=gene,
            )
            coverage[gene] = cov

    return {
        "vcf": vcf_result,
        "bam": bam_result,
        "coverage": coverage,
        "errors": errors,
        "passed": len(errors) == 0,
    }
