from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd

from paths import PATHS


def _single_row_df(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing input file: {csv_path}")

    df = pd.read_csv(csv_path)
    if df.empty:
        raise ValueError(f"Input file has no rows: {csv_path}")

    return df.iloc[[0]].copy()


def _discover_custom_feature_csvs(features_dir: Path) -> List[Path]:
    return sorted(p for p in features_dir.rglob("*_custom_features.csv") if p.is_file())


def merge_custom_features_for_ids(
    complex_ids: List[str],
    output_csv: str | None = None,
    append: bool = False,
    dedupe_by_id: bool = False,
) -> Path:
    rows = []

    for complex_id in complex_ids:
        in_csv = PATHS.features_dir / f"{complex_id}_custom_features.csv"
        one = _single_row_df(in_csv)

        if "ID" not in one.columns:
            one.insert(0, "ID", complex_id)
        else:
            one["ID"] = complex_id

        rows.append(one)

    merged = pd.concat(rows, axis=0, ignore_index=True)

    if output_csv:
        out_path = Path(output_csv)
    elif append:
        out_path = PATHS.features_dir / "custom_features_merged_all.csv"
    elif len(complex_ids) == 1:
        complex_id = complex_ids[0]
        out_path = PATHS.features_dir / f"{complex_id}_custom_features_merged.csv"
    else:
        out_path = PATHS.features_dir / "custom_features_merged_all.csv"

    out_path.parent.mkdir(parents=True, exist_ok=True)

    if append and out_path.exists():
        existing = pd.read_csv(out_path)
        combined = pd.concat([existing, merged], axis=0, ignore_index=True)
        if dedupe_by_id and "ID" in combined.columns:
            combined = combined.drop_duplicates(subset=["ID"], keep="last")
        merged = combined

    merged.to_csv(out_path, index=False)
    return out_path


def merge_all_custom_features(
    output_csv: str | None = None,
    append: bool = False,
    dedupe_by_id: bool = False,
) -> Path:
    in_csvs = _discover_custom_feature_csvs(PATHS.features_dir)
    if not in_csvs:
        raise FileNotFoundError(f"No *_custom_features.csv files found under: {PATHS.features_dir}")

    rows = []
    for in_csv in in_csvs:
        one = _single_row_df(in_csv)
        if "ID" not in one.columns:
            inferred_id = in_csv.name.replace("_custom_features.csv", "")
            one.insert(0, "ID", inferred_id)
        rows.append(one)

    merged = pd.concat(rows, axis=0, ignore_index=True)

    out_path = Path(output_csv) if output_csv else PATHS.features_dir / "custom_features_merged_all.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if append and out_path.exists():
        existing = pd.read_csv(out_path)
        combined = pd.concat([existing, merged], axis=0, ignore_index=True)
        if dedupe_by_id and "ID" in combined.columns:
            combined = combined.drop_duplicates(subset=["ID"], keep="last")
        merged = combined

    merged.to_csv(out_path, index=False)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Merge all *_custom_features.csv files found under features dir."
        )
    )
    parser.add_argument(
        "--output_csv",
        type=str,
        default=None,
        help="Optional explicit output CSV path",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help=(
            "Append to existing output CSV (preserve all rows). "
            "If --output_csv is not given, appends to <features_dir>/custom_features_merged_all.csv"
        ),
    )
    parser.add_argument(
        "--dedupe_by_id",
        action="store_true",
        help="When used with --append, de-duplicate by ID and keep latest row per ID",
    )
    args = parser.parse_args()

    out = merge_all_custom_features(
        output_csv=args.output_csv,
        append=args.append,
        dedupe_by_id=args.dedupe_by_id,
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
