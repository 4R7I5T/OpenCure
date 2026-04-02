import os
import pandas as pd


CROPSR_OUTPUT = "./output.csv"
CAS_OFFINDER_RESULTS = "./cas_offinder_results.tsv"
ANNOTATED_OUTPUT = "./output_annotated.csv"


def _detect_column(df, candidates, source_name):
    for col in candidates:
        if col in df.columns:
            return col
    raise ValueError(
        f"None of the expected columns {candidates} found in {source_name}. "
        f"Available columns: {list(df.columns)}"
    )


def _fallback_copy(cropsr_df: pd.DataFrame) -> None:
    if "on_site_score" not in cropsr_df.columns:
        raise ValueError("Expected 'on_site_score' column in CROPSR output.")

    cropsr_df["on_site_score"] = pd.to_numeric(
        cropsr_df["on_site_score"], errors="coerce"
    )
    cropsr_df["offtarget_risk"] = 0.0
    cropsr_df["combined_score"] = cropsr_df["on_site_score"]
    cropsr_df.to_csv(ANNOTATED_OUTPUT, index=False)
    print(
        f"Cas-OFFinder results missing or unusable; "
        f"wrote fallback annotated file with offtarget_risk=0 to {ANNOTATED_OUTPUT}"
    )


def main():
    if not os.path.exists(CROPSR_OUTPUT):
        raise FileNotFoundError(f"CROPSR output not found: {CROPSR_OUTPUT}")

    cropsr_df = pd.read_csv(CROPSR_OUTPUT)

    if not os.path.exists(CAS_OFFINDER_RESULTS):
        _fallback_copy(cropsr_df)
        return

    try:
        cas_df = pd.read_csv(CAS_OFFINDER_RESULTS, sep="\t")
    except Exception:
        _fallback_copy(cropsr_df)
        return

    try:
        guide_col_cropsr = _detect_column(
            cropsr_df,
            ["guide_sequence", "sequence", "gRNA", "protospacer", "target_sequence"],
            CROPSR_OUTPUT,
        )
        guide_col_cas = _detect_column(
            cas_df,
            ["Id", "guide_sequence", "sequence", "RNA sequence", "RNA_sequence"],
            CAS_OFFINDER_RESULTS,
        )
        mismatch_col = _detect_column(
            cas_df,
            ["mismatches", "mismatch_count", "num_mismatches"],
            CAS_OFFINDER_RESULTS,
        )
    except ValueError:
        _fallback_copy(cropsr_df)
        return

    cas_df[mismatch_col] = pd.to_numeric(cas_df[mismatch_col], errors="coerce")
    cas_df = cas_df.dropna(subset=[mismatch_col])
    if cas_df.empty:
        _fallback_copy(cropsr_df)
        return

    agg = (
        cas_df.groupby(guide_col_cas)[mismatch_col]
        .agg(
            offtarget_hits="size",
            offtarget_min_mismatches="min",
            offtarget_mean_mismatches="mean",
        )
        .reset_index()
    )
    agg["offtarget_risk"] = agg["offtarget_hits"] / (agg["offtarget_min_mismatches"] + 1)

    annotated = cropsr_df.merge(
        agg,
        how="left",
        left_on=guide_col_cropsr,
        right_on=guide_col_cas,
    )

    if guide_col_cas != guide_col_cropsr and guide_col_cas in annotated.columns:
        annotated = annotated.drop(columns=[guide_col_cas])

    if "offtarget_risk" in annotated.columns:
        max_risk = annotated["offtarget_risk"].max()
        if pd.notna(max_risk):
            annotated["offtarget_risk"] = annotated["offtarget_risk"].fillna(max_risk)
        else:
            annotated["offtarget_risk"] = annotated["offtarget_risk"].fillna(0.0)
    else:
        annotated["offtarget_risk"] = 0.0

    if "on_site_score" not in annotated.columns:
        raise ValueError("Expected 'on_site_score' column in CROPSR output.")

    annotated["on_site_score"] = pd.to_numeric(
        annotated["on_site_score"], errors="coerce"
    )
    annotated["combined_score"] = (
        annotated["on_site_score"] - annotated["offtarget_risk"]
    )

    annotated.to_csv(ANNOTATED_OUTPUT, index=False)
    print(
        f"Annotated output with off-target risk written to {ANNOTATED_OUTPUT}"
    )


if __name__ == "__main__":
    main()
