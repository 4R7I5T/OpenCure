"""
VCF file parsing and variant extraction utilities.
"""

import gzip
import logging
from dataclasses import dataclass, field

log = logging.getLogger(__name__)


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: list[str]
    qual: float
    filter_status: str
    info: dict = field(default_factory=dict)
    genotype: list[int] | None = None
    var_type: str = "snp"
    zygosity: str = "UNKNOWN"
    depth: int = 0           # total read depth (DP) at this site
    alt_depth: int = 0       # alternate allele read depth (AD[1])
    genotype_qual: int = 0   # genotype quality (GQ)

    @property
    def is_pass(self) -> bool:
        return self.filter_status in ("PASS", ".")


def parse_vcf(vcf_path: str, region: tuple[str, int, int] | None = None):
    """
    Parse a VCF/VCF.gz file and yield Variant objects.
    If *region* is provided as (chrom, start, end), only variants in that
    region are yielded.
    """
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    mode = "rt" if str(vcf_path).endswith(".gz") else "r"

    with opener(vcf_path, mode) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(",")
            qual = float(fields[5]) if fields[5] != "." else 0.0
            filt = fields[6]

            # Region filter
            if region:
                r_chrom, r_start, r_end = region
                # Normalize chrom names for comparison
                c1 = chrom.replace("chr", "")
                c2 = r_chrom.replace("chr", "")
                if c1 != c2 or pos < r_start or pos > r_end:
                    continue

            # Parse INFO field
            info = {}
            for entry in fields[7].split(";"):
                if "=" in entry:
                    k, v = entry.split("=", 1)
                    info[k] = v
                else:
                    info[entry] = True

            # Determine variant type
            var_type = info.get("TYPE", "snp")
            if isinstance(var_type, str) and "," in var_type:
                var_type = var_type.split(",")[0]

            # Parse genotype and quality fields from FORMAT/sample columns
            genotype = None
            zygosity = "UNKNOWN"
            depth = 0
            alt_depth = 0
            genotype_qual = 0
            if len(fields) >= 10:
                fmt = fields[8].split(":")
                sample = fields[9].split(":")
                fmt_map = {}
                for fi, fv in zip(fmt, sample):
                    fmt_map[fi] = fv

                # Genotype (GT)
                if "GT" in fmt_map:
                    gt_str = fmt_map["GT"]
                    sep = "/" if "/" in gt_str else "|"
                    alleles = gt_str.split(sep)
                    try:
                        genotype = [int(a) for a in alleles]
                        if genotype[0] == genotype[1]:
                            zygosity = "HOM_REF" if genotype[0] == 0 else "HOM_ALT"
                        else:
                            zygosity = "HET"
                    except (ValueError, IndexError):
                        pass

                # Read depth (DP)
                if "DP" in fmt_map:
                    try:
                        depth = int(fmt_map["DP"])
                    except (ValueError, TypeError):
                        pass

                # Allele depths (AD) -- format: ref_depth,alt_depth[,alt2_depth]
                if "AD" in fmt_map:
                    try:
                        ad_parts = fmt_map["AD"].split(",")
                        if len(ad_parts) >= 2:
                            alt_depth = int(ad_parts[1])
                    except (ValueError, TypeError):
                        pass

                # Genotype quality (GQ)
                if "GQ" in fmt_map:
                    try:
                        genotype_qual = int(fmt_map["GQ"])
                    except (ValueError, TypeError):
                        pass

            yield Variant(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                qual=qual,
                filter_status=filt,
                info=info,
                genotype=genotype,
                var_type=var_type,
                zygosity=zygosity,
                depth=depth,
                alt_depth=alt_depth,
                genotype_qual=genotype_qual,
            )


def extract_variants_in_region(vcf_path, chrom, start, end, pass_only=True):
    """Return a list of Variants within a genomic region."""
    # Try pysam for indexed random-access first
    result = _try_pysam_fetch(vcf_path, chrom, start, end, pass_only)
    if result is not None:
        return result
    # Fallback: linear scan (slow for large files)
    variants = []
    for v in parse_vcf(vcf_path, region=(chrom, start, end)):
        if pass_only and not v.is_pass:
            continue
        variants.append(v)
    return variants


def _try_pysam_fetch(vcf_path, chrom, start, end, pass_only):
    """Use pysam.VariantFile for indexed tabix access if available."""
    try:
        import pysam
    except ImportError:
        return None

    try:
        vf = pysam.VariantFile(str(vcf_path))
    except Exception:
        return None

    # Try both chr-prefixed and bare chrom names
    chroms_to_try = [chrom, chrom.replace("chr", ""), "chr" + chrom.replace("chr", "")]
    variants = []
    fetched = False

    for c in chroms_to_try:
        try:
            for rec in vf.fetch(c, max(0, start - 1), end):
                fetched = True
                filt = "PASS" if not rec.filter or "PASS" in rec.filter else ",".join(rec.filter)
                if pass_only and filt not in ("PASS", "."):
                    continue

                gt = None
                zyg = "UNKNOWN"
                if rec.samples and len(rec.samples) > 0:
                    sample = rec.samples[0]
                    if "GT" in sample:
                        alleles = sample["GT"]
                        if alleles and len(alleles) >= 2:
                            try:
                                gt = [int(a) if a is not None else 0 for a in alleles]
                                if gt[0] == gt[1]:
                                    zyg = "HOM_REF" if gt[0] == 0 else "HOM_ALT"
                                else:
                                    zyg = "HET"
                            except (TypeError, ValueError):
                                pass

                var_type = rec.info.get("TYPE", "snp") if rec.info else "snp"
                if isinstance(var_type, tuple):
                    var_type = var_type[0]

                variants.append(Variant(
                    chrom=rec.chrom,
                    pos=rec.pos,
                    ref=rec.ref,
                    alt=list(rec.alts) if rec.alts else [],
                    qual=rec.qual if rec.qual else 0.0,
                    filter_status=filt,
                    genotype=gt,
                    var_type=str(var_type),
                    zygosity=zyg,
                ))
            if fetched:
                break
        except (ValueError, KeyError):
            continue

    vf.close()
    return variants if fetched else None


def extract_variants_multi_region(vcf_path, regions, pass_only=True):
    """
    Extract variants for multiple regions in a single pass.
    *regions*: {name: (chrom, start, end)}
    Returns: {name: [Variant, ...]}
    """
    # Try indexed access first (fast)
    results = {}
    any_found = False
    for name, (chrom, start, end) in regions.items():
        fetched = _try_pysam_fetch(vcf_path, chrom, start, end, pass_only)
        if fetched is not None:
            results[name] = fetched
            any_found = True
        else:
            results[name] = []

    if any_found:
        return results

    # Fallback: single-pass linear scan matching all regions at once
    log.info("Using single-pass linear scan for %d regions", len(regions))
    results = {name: [] for name in regions}

    # Build lookup: normalized_chrom -> [(name, start, end), ...]
    lookup: dict[str, list] = {}
    for name, (chrom, start, end) in regions.items():
        norm = chrom.replace("chr", "")
        lookup.setdefault(norm, []).append((name, start, end))

    for v in parse_vcf(vcf_path):
        if pass_only and not v.is_pass:
            continue
        norm_chrom = v.chrom.replace("chr", "")
        if norm_chrom not in lookup:
            continue
        for name, start, end in lookup[norm_chrom]:
            if start <= v.pos <= end:
                results[name].append(v)

    return results
