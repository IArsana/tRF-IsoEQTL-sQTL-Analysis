"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_moradi_s3_s4.py

Description:
    Smoke test for Moradi Supplementary Tables S3 and S4.

    The script validates real Moradi trans-QTL and GWAS-linked
    trans-QTL datasets using MoradiWorkbookValidator.

    The smoke test verifies:
        - Source workbooks exist.
        - Expected worksheets are present.
        - Validation completes without unexpected exceptions.
        - Sheet-level QC statistics are generated.
        - Standardized anomalies are produced through AnomalyFactory.
        - The validation report is JSON serializable and readable.

    Raw source files are never modified.

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
from typing import Any

from pcatrfqtl.validation.runners.moradi_workbook import (
    MoradiWorkbookValidator,
)


# ===========================================================================
# Paths
# ===========================================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

MORADI_DIR = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "moradi"
)

S3_WORKBOOK = (
    MORADI_DIR
    / "Supplementary Table S3.xlsx"
)

S4_WORKBOOK = (
    MORADI_DIR
    / "Supplementary Table S4.xlsx"
)

OUTPUT = (
    PROJECT_ROOT
    / "data"
    / "interim"
    / "qc"
    / "moradi_s3_s4_validation.json"
)


# ===========================================================================
# Expected worksheet structure
# ===========================================================================

EXPECTED_S3_SHEETS = {
    "Trans-intron-sQTL",
    "Trans-exon-sQTl",
    "Trans-iso-eQTL",
}

EXPECTED_S4_SHEETS = {
    "Trans-intron-sQTL-LD_0.5",
    "Trans-exon-sQTL-LD_0.5",
    "Trans-iso-eQTL-LD_0.5",
}


# ===========================================================================
# Display helpers
# ===========================================================================


def separator(
    character: str = "=",
    width: int = 80,
) -> None:
    """Print a visual separator."""

    print(character * width)


def print_sheet_result(
    sheet_name: str,
    result: dict[str, Any],
) -> None:
    """Print a concise validation summary for one worksheet."""

    print()
    print(f"  {sheet_name}")

    print(
        f"    Status     : "
        f"{result.get('status', 'UNKNOWN')}"
    )

    print(
        f"    Rows       : "
        f"{result.get('rows', 'N/A')}"
    )

    print(
        f"    Columns    : "
        f"{result.get('columns', 'N/A')}"
    )

    print(
        "    Errors     : "
        f"{result.get('validation_error_counts', {})}"
    )

    print(
        "    Duplicates : "
        f"{result.get('duplicates', {})}"
    )

    notes = result.get(
        "source_schema_notes",
        [],
    )

    if notes:
        print(
            "    Source notes:"
        )

        for note in notes:
            print(
                f"      - {note}"
            )


# ===========================================================================
# Validation helpers
# ===========================================================================


def validate_source_files() -> bool:
    """Verify S3/S4 workbooks exist."""

    success = True

    for workbook in [
        S3_WORKBOOK,
        S4_WORKBOOK,
    ]:
        if workbook.exists():
            print(
                "[PASS] Source workbook exists: "
                f"{workbook.name}"
            )
        else:
            print(
                "[FAIL] Source workbook missing: "
                f"{workbook}"
            )
            success = False

    return success


def validate_expected_sheets(
    report: dict[str, Any],
) -> bool:
    """Verify expected worksheets occur in the report."""

    files = report.get(
        "files",
        {},
    )

    s3 = files.get(
        "supplementary_s3",
        {},
    )

    s4 = files.get(
        "supplementary_s4",
        {},
    )

    actual_s3 = set(
        s3.get(
            "sheets",
            {},
        )
    )

    actual_s4 = set(
        s4.get(
            "sheets",
            {},
        )
    )

    missing_s3 = (
        EXPECTED_S3_SHEETS
        - actual_s3
    )

    missing_s4 = (
        EXPECTED_S4_SHEETS
        - actual_s4
    )

    if missing_s3:
        print(
            "[FAIL] Missing S3 sheets: "
            + ", ".join(
                sorted(missing_s3)
            )
        )
        return False

    if missing_s4:
        print(
            "[FAIL] Missing S4 sheets: "
            + ", ".join(
                sorted(missing_s4)
            )
        )
        return False

    print(
        "[PASS] All expected S3/S4 sheets are present"
    )

    return True


def save_report(
    report: dict[str, Any],
) -> bool:
    """Save validation report to JSON."""

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
            "[FAIL] Could not serialize report: "
            f"{exc}"
        )
        return False

    print(
        "[PASS] Validation report written: "
        f"{OUTPUT}"
    )

    return True


def verify_saved_report() -> bool:
    """Verify generated JSON can be read."""

    try:
        with OUTPUT.open(
            "r",
            encoding="utf-8",
        ) as handle:
            report = json.load(
                handle
            )

    except Exception as exc:
        print(
            "[FAIL] Could not read generated report: "
            f"{exc}"
        )
        return False

    if (
        report.get("dataset_id")
        != "moradi"
    ):
        print(
            "[FAIL] Unexpected dataset_id "
            "in saved report"
        )
        return False

    print(
        "[PASS] Generated JSON can be read successfully"
    )

    return True


# ===========================================================================
# Main
# ===========================================================================


def main() -> int:
    """Run S3/S4 Moradi smoke test."""

    separator()

    print(
        "pcatRFQTL — Moradi S3/S4 Smoke Test"
    )

    separator()

    print()
    print(
        "Checking source datasets..."
    )

    if not validate_source_files():
        return 1

    validator = MoradiWorkbookValidator(
        MORADI_DIR
    )

    # --------------------------------------------------------------
    # Reset anomaly state manually because we are validating only S3/S4.
    # --------------------------------------------------------------

    validator.anomaly_factory = (
        validator.anomaly_factory.__class__()
    )

    validator.anomalies = []

    print()
    print(
        "Running Moradi S3 validation..."
    )

    try:
        s3 = validator.validate_s3()
    except Exception as exc:
        print(
            "[FAIL] S3 validation raised exception: "
            f"{exc}"
        )
        return 1

    print(
        "[PASS] S3 validation completed"
    )

    print()
    print(
        "Running Moradi S4 validation..."
    )

    try:
        s4 = validator.validate_s4()
    except Exception as exc:
        print(
            "[FAIL] S4 validation raised exception: "
            f"{exc}"
        )
        return 1

    print(
        "[PASS] S4 validation completed"
    )

    status = (
        validator._aggregate_status(
            [
                s3["status"],
                s4["status"],
            ]
        )
    )

    report = {
        "dataset_id": "moradi",
        "status": status,
        "files": {
            "supplementary_s3": s3,
            "supplementary_s4": s4,
        },
        "anomalies": list(
            validator.anomalies
        ),
        "validation_summary": (
            validator._summarize_anomalies(
                validator.anomalies
            )
        ),
    }

    if not validate_expected_sheets(
        report
    ):
        return 1

    # --------------------------------------------------------------
    # S3 summary
    # --------------------------------------------------------------

    print()
    separator("-")

    print(
        "Supplementary Table S3 — Trans-QTL associations"
    )

    separator("-")

    print(
        f"Workbook status: "
        f"{s3.get('status')}"
    )

    for (
        sheet_name,
        result,
    ) in s3.get(
        "sheets",
        {},
    ).items():

        print_sheet_result(
            sheet_name,
            result,
        )

    # --------------------------------------------------------------
    # S4 summary
    # --------------------------------------------------------------

    print()
    separator("-")

    print(
        "Supplementary Table S4 — "
        "GWAS-linked trans-QTL associations"
    )

    separator("-")

    print(
        f"Workbook status: "
        f"{s4.get('status')}"
    )

    for (
        sheet_name,
        result,
    ) in s4.get(
        "sheets",
        {},
    ).items():

        print_sheet_result(
            sheet_name,
            result,
        )

    # --------------------------------------------------------------
    # Overall status
    # --------------------------------------------------------------

    print()
    separator("-")

    print(
        "Overall validation status"
    )

    separator("-")

    print(
        f"Dataset status: "
        f"{report.get('status')}"
    )

    print(
        "Total anomalies: "
        f"{len(report['anomalies'])}"
    )

    # --------------------------------------------------------------
    # Serialization
    # --------------------------------------------------------------

    print()

    if not save_report(
        report
    ):
        return 1

    if not verify_saved_report():
        return 1

    print()
    separator()

    print(
        "SMOKE TEST PASSED"
    )

    separator()

    print()
    print(
        "A PASS smoke test means the real source data "
        "was processed successfully. Individual validation "
        "statuses may still be WARNING or FAIL when source "
        "anomalies are detected."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )