"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_trfqtl.py

Description:
    Smoke test for the Cancer-tRFQTL validation pipeline.

    This script validates that the real Cancer-tRFQTL supplementary
    workbook can be loaded, validated, summarized, and serialized
    into a QC report without modifying the raw source data.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import json
from pathlib import Path

from pcatrfqtl.validation.runners.trfqtl_workbook import (
    TRFQTLWorkbookValidator,
)


PROJECT_ROOT = Path(__file__).resolve().parent.parent

WORKBOOK = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "trfqtl"
    / "can-25-1282_supplementary_data_suppst1-11.xlsx"
)

OUTPUT = (
    PROJECT_ROOT
    / "data"
    / "interim"
    / "qc"
    / "cancer_trfqtl_validation.json"
)


def main() -> int:
    """Run the Cancer-tRFQTL validation smoke test."""

    print("=" * 72)
    print("pcatRFQTL — Cancer-tRFQTL Smoke Test")
    print("=" * 72)

    if not WORKBOOK.exists():
        print(f"[FAIL] Workbook not found: {WORKBOOK}")
        return 1

    print(f"[PASS] Workbook exists: {WORKBOOK.name}")

    validator = TRFQTLWorkbookValidator(WORKBOOK)

    try:
        report = validator.validate_all()
    except Exception as exc:
        print(f"[FAIL] Validation raised exception: {exc}")
        return 1

    print("[PASS] Workbook validation completed")

    if "sheets" not in report:
        print("[FAIL] Report does not contain 'sheets'")
        return 1

    for required_sheet in ("S2", "S10"):
        if required_sheet not in report["sheets"]:
            print(
                f"[FAIL] Required sheet {required_sheet} "
                "missing from validation report"
            )
            return 1

        print(
            f"[PASS] {required_sheet} present in validation report"
        )

    s2 = report["sheets"]["S2"]
    s10 = report["sheets"]["S10"]

    print()
    print("-" * 72)
    print("S2 — GWAS-associated tRFQTLs")
    print("-" * 72)
    print(f"Status     : {s2['status']}")
    print(f"Rows       : {s2['rows']}")
    print(f"Columns    : {s2['columns']}")
    print(f"Errors     : {s2['validation_error_counts']}")
    print(
        f"Anomalies  : {len(s2.get('anomalies', []))}"
    )

    print()
    print("-" * 72)
    print("S10 — Cancer-risk-associated tRFQTLs")
    print("-" * 72)
    print(f"Status     : {s10['status']}")
    print(f"Rows       : {s10['rows']}")
    print(f"Columns    : {s10['columns']}")
    print(f"Errors     : {s10['validation_error_counts']}")
    print(
        f"Anomalies  : {len(s10.get('anomalies', []))}"
    )

    # --------------------------------------------------------------
    # Known source anomaly
    # --------------------------------------------------------------

    s10_snp_errors = (
        s10["validation_error_counts"]
        .get("snp_ids", 0)
    )

    if s10_snp_errors == 1:
        print(
            "[PASS] Known S10 SNP-ID anomaly detected "
            "(expected count: 1)"
        )
    else:
        print(
            "[WARN] Expected one known S10 SNP-ID anomaly, "
            f"but detected {s10_snp_errors}"
        )

    # --------------------------------------------------------------
    # Anomaly ID sanity check
    # --------------------------------------------------------------

    anomalies = report.get("anomalies", [])

    anomaly_ids = [
        anomaly.get("anomaly_id")
        for anomaly in anomalies
    ]

    if len(anomaly_ids) != len(set(anomaly_ids)):
        print("[FAIL] Duplicate anomaly IDs detected")
        return 1

    print(
        "[PASS] Anomaly IDs are unique within validation run"
    )

    # --------------------------------------------------------------
    # JSON serialization
    # --------------------------------------------------------------

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    try:
        with OUTPUT.open(
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
                default=str,
            )
    except Exception as exc:
        print(
            f"[FAIL] Could not serialize validation report: {exc}"
        )
        return 1

    print(f"[PASS] Validation report written: {OUTPUT}")

    # --------------------------------------------------------------
    # Verify generated JSON
    # --------------------------------------------------------------

    try:
        with OUTPUT.open(
            "r",
            encoding="utf-8",
        ) as handle:
            loaded = json.load(handle)
    except Exception as exc:
        print(
            f"[FAIL] Generated JSON could not be read: {exc}"
        )
        return 1

    if loaded.get("dataset_id") != "cancer_trfqtl":
        print(
            "[FAIL] Unexpected dataset_id in generated report"
        )
        return 1

    print("[PASS] Generated JSON can be read successfully")

    print()
    print("=" * 72)
    print("SMOKE TEST PASSED")
    print("=" * 72)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())