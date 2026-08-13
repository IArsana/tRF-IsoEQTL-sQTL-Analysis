"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_moradi.py

Description:
    Smoke test for the Moradi core QTL validation pipeline.

    This script validates the real Moradi Supplementary Tables S1 and S2
    using MoradiWorkbookValidator and verifies that:

        - Source workbooks exist.
        - Core QTL worksheets can be loaded.
        - Validation completes without unexpected exceptions.
        - Expected worksheets are present in the validation report.
        - Validation statistics are generated.
        - The resulting report is JSON serializable.
        - The generated JSON report can be read back successfully.

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

S1_WORKBOOK = (
    MORADI_DIR
    / "Supplementary Table S1.xlsx"
)

S2_WORKBOOK = (
    MORADI_DIR
    / "Supplementary Table S2.xlsx"
)

OUTPUT = (
    PROJECT_ROOT
    / "data"
    / "interim"
    / "qc"
    / "moradi_validation.json"
)


# ===========================================================================
# Expected source structure
# ===========================================================================

EXPECTED_S1_SHEETS = {
    "Cis-intron-retention-sQTL-0.05",
    "Cis-cassette-exon-sQTL-0.05",
    "Cis-iso-eQTL-0.05",
}

EXPECTED_S2_SHEETS = {
    "Cis-intron-sQTL-GWAS-LD_0.5",
    "Cis-exon-sQTL-GWAS-LD_0.5",
    "Cis-iso-eQTl-GWAS-LD_0.5",
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

    errors = result.get(
        "validation_error_counts",
        {},
    )

    print(
        f"    Errors     : {errors}"
    )

    duplicates = result.get(
        "duplicates",
        {},
    )

    print(
        f"    Duplicates : {duplicates}"
    )


# ===========================================================================
# Validation checks
# ===========================================================================


def validate_source_files() -> bool:
    """Verify required Moradi workbooks exist."""

    success = True

    for workbook in [
        S1_WORKBOOK,
        S2_WORKBOOK,
    ]:
        if workbook.exists():
            print(
                f"[PASS] Source workbook exists: "
                f"{workbook.name}"
            )
        else:
            print(
                f"[FAIL] Source workbook missing: "
                f"{workbook}"
            )

            success = False

    return success


def validate_report_structure(
    report: dict[str, Any],
) -> bool:
    """Verify the core validation-report structure."""

    required_keys = {
        "dataset_id",
        "status",
        "files",
    }

    missing = (
        required_keys
        - set(report)
    )

    if missing:
        print(
            "[FAIL] Validation report missing keys: "
            + ", ".join(sorted(missing))
        )

        return False

    if (
        report.get("dataset_id")
        != "moradi"
    ):
        print(
            "[FAIL] Unexpected dataset_id: "
            f"{report.get('dataset_id')!r}"
        )

        return False

    print(
        "[PASS] Validation report structure is valid"
    )

    return True


def validate_expected_sheets(
    report: dict[str, Any],
) -> bool:
    """Verify expected S1/S2 worksheets occur in the report."""

    files = report.get(
        "files",
        {},
    )

    s1 = files.get(
        "supplementary_s1",
        {},
    )

    s2 = files.get(
        "supplementary_s2",
        {},
    )

    s1_sheets = set(
        s1.get(
            "sheets",
            {},
        )
    )

    s2_sheets = set(
        s2.get(
            "sheets",
            {},
        )
    )

    missing_s1 = (
        EXPECTED_S1_SHEETS
        - s1_sheets
    )

    missing_s2 = (
        EXPECTED_S2_SHEETS
        - s2_sheets
    )

    if missing_s1:
        print(
            "[FAIL] Missing S1 validation sheets: "
            + ", ".join(
                sorted(missing_s1)
            )
        )

        return False

    if missing_s2:
        print(
            "[FAIL] Missing S2 validation sheets: "
            + ", ".join(
                sorted(missing_s2)
            )
        )

        return False

    print(
        "[PASS] All expected S1/S2 sheets "
        "are present in the report"
    )

    return True


def save_report(
    report: dict[str, Any],
) -> bool:
    """Serialize the validation report as JSON."""

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
            "[FAIL] Could not serialize "
            f"validation report: {exc}"
        )

        return False

    print(
        f"[PASS] Validation report written: "
        f"{OUTPUT}"
    )

    return True


def verify_saved_report() -> bool:
    """Read the generated JSON report back from disk."""

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
            "[FAIL] Could not read generated "
            f"JSON report: {exc}"
        )

        return False

    if (
        report.get("dataset_id")
        != "moradi"
    ):
        print(
            "[FAIL] Saved report contains "
            "unexpected dataset_id"
        )

        return False

    print(
        "[PASS] Generated JSON can be read successfully"
    )

    return True


# ===========================================================================
# Main smoke test
# ===========================================================================


def main() -> int:
    """Run the Moradi core QTL smoke test."""

    separator()

    print(
        "pcatRFQTL — Moradi Core QTL Smoke Test"
    )

    separator()

    print()
    print("Checking source datasets...")

    if not validate_source_files():
        return 1

    print()
    print("Initializing MoradiWorkbookValidator...")

    validator = MoradiWorkbookValidator(
        MORADI_DIR
    )

    print(
        "[PASS] Validator initialized"
    )

    # ------------------------------------------------------------------
    # Run real validation
    # ------------------------------------------------------------------

    print()
    print(
        "Running validation on Moradi S1 + S2..."
    )

    print(
        "This may take some time because S1 contains "
        "hundreds of thousands of QTL records."
    )

    try:
        report = (
            validator.validate_core()
        )

    except Exception as exc:
        print(
            "[FAIL] Validation raised "
            f"an exception: {exc}"
        )

        return 1

    print(
        "[PASS] Core validation completed"
    )

    # ------------------------------------------------------------------
    # Validate report structure
    # ------------------------------------------------------------------

    if not validate_report_structure(
        report
    ):
        return 1

    if not validate_expected_sheets(
        report
    ):
        return 1

    # ------------------------------------------------------------------
    # Display S1
    # ------------------------------------------------------------------

    files = report["files"]

    s1 = files[
        "supplementary_s1"
    ]

    print()
    separator("-")

    print(
        "Supplementary Table S1 "
        "— Cis QTL associations"
    )

    separator("-")

    print(
        f"Workbook status: "
        f"{s1.get('status')}"
    )

    for (
        sheet_name,
        sheet_result,
    ) in s1.get(
        "sheets",
        {},
    ).items():

        print_sheet_result(
            sheet_name,
            sheet_result,
        )

    # ------------------------------------------------------------------
    # Display S2
    # ------------------------------------------------------------------

    s2 = files[
        "supplementary_s2"
    ]

    print()
    separator("-")

    print(
        "Supplementary Table S2 "
        "— GWAS-linked cis QTL associations"
    )

    separator("-")

    print(
        f"Workbook status: "
        f"{s2.get('status')}"
    )

    for (
        sheet_name,
        sheet_result,
    ) in s2.get(
        "sheets",
        {},
    ).items():

        print_sheet_result(
            sheet_name,
            sheet_result,
        )

    # ------------------------------------------------------------------
    # Dataset status
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # JSON serialization
    # ------------------------------------------------------------------

    print()

    if not save_report(
        report
    ):
        return 1

    if not verify_saved_report():
        return 1

    # ------------------------------------------------------------------
    # Smoke-test result
    # ------------------------------------------------------------------

    print()
    separator()

    print(
        "SMOKE TEST PASSED"
    )

    separator()

    print()
    print(
        "Note:"
    )

    print(
        "A validation status of FAIL or WARNING does not "
        "automatically mean that the smoke test failed."
    )

    print(
        "The smoke test verifies that the real source data "
        "can be processed and that validation findings are "
        "reported reproducibly."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )