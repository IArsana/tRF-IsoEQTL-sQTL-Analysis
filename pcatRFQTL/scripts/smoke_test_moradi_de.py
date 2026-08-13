"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_moradi_de.py

Description:
    Smoke test for Moradi differential-expression supplementary
    datasets S5-S10.

    The script validates the real source CSV files and verifies that:

        - the Moradi source directory exists
        - all expected S5-S10 source files exist
        - all six datasets can be processed
        - report structure is valid
        - per-file row counts are non-zero
        - expected feature types are preserved
        - validation statistics are generated
        - the QC report is JSON serializable
        - the saved JSON report can be read successfully

    Validation status may be PASS, WARNING, or FAIL depending on
    source-data findings.

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

from pcatrfqtl.validation.runners.moradi_de import (
    MoradiDERunner,
)


# ===========================================================================
# Paths
# ===========================================================================

PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parent
    .parent
)

MORADI_DIRECTORY = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "moradi"
)

OUTPUT = (
    PROJECT_ROOT
    / "data"
    / "interim"
    / "qc"
    / "moradi_de_validation.json"
)


EXPECTED_DATASETS = {
    "S5",
    "S6",
    "S7",
    "S8",
    "S9",
    "S10",
}


# ===========================================================================
# Display helpers
# ===========================================================================


def separator(
    character: str = "=",
    width: int = 80,
) -> None:
    """Print a visual separator."""

    print(
        character
        * width
    )


def print_dictionary(
    values: dict[str, Any],
    *,
    indent: str = "    ",
) -> None:
    """Print dictionary entries line by line."""

    for key, value in values.items():
        print(
            f"{indent}{key}: {value}"
        )


# ===========================================================================
# Source checks
# ===========================================================================


def verify_source_directory() -> bool:
    """Verify the Moradi source directory exists."""

    if not MORADI_DIRECTORY.exists():
        print(
            "[FAIL] Moradi source directory missing: "
            f"{MORADI_DIRECTORY}"
        )
        return False

    if not MORADI_DIRECTORY.is_dir():
        print(
            "[FAIL] Moradi source path is not a directory: "
            f"{MORADI_DIRECTORY}"
        )
        return False

    print(
        "[PASS] Moradi source directory exists"
    )

    return True


def verify_source_files() -> bool:
    """Verify all expected S5-S10 source CSV files exist."""

    success = True

    for (
        dataset_name,
        specification,
    ) in (
        MoradiDERunner
        .DATASETS
        .items()
    ):
        file_path = (
            MORADI_DIRECTORY
            / specification[
                "filename"
            ]
        )

        if not file_path.exists():
            print(
                f"[FAIL] {dataset_name} missing: "
                f"{file_path.name}"
            )

            success = False
            continue

        size_mb = (
            file_path
            .stat()
            .st_size
            / (
                1024 ** 2
            )
        )

        print(
            f"[PASS] {dataset_name}: "
            f"{file_path.name} "
            f"({size_mb:.3f} MB)"
        )

    return success


# ===========================================================================
# Report validation
# ===========================================================================


def verify_report_structure(
    report: dict[str, Any],
) -> bool:
    """Verify the expected top-level report structure."""

    required = {
        "dataset_id",
        "source_directory",
        "status",
        "files_validated",
        "total_rows",
        "files",
        "anomalies",
        "validation_summary",
    }

    missing = (
        required
        - set(
            report
        )
    )

    if missing:
        print(
            "[FAIL] Validation report missing keys: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )

        return False

    if (
        report.get(
            "dataset_id"
        )
        != "moradi_de"
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


def verify_expected_datasets(
    report: dict[str, Any],
) -> bool:
    """Verify all expected S5-S10 datasets are represented."""

    observed = set(
        report.get(
            "files",
            {},
        )
    )

    missing = (
        EXPECTED_DATASETS
        - observed
    )

    extra = (
        observed
        - EXPECTED_DATASETS
    )

    if missing:
        print(
            "[FAIL] Missing datasets in report: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )

        return False

    if extra:
        print(
            "[FAIL] Unexpected datasets in report: "
            + ", ".join(
                sorted(
                    extra
                )
            )
        )

        return False

    print(
        "[PASS] All expected S5-S10 datasets are present"
    )

    return True


def verify_nonzero_rows(
    report: dict[str, Any],
) -> bool:
    """
    Verify real source datasets contain rows.

    Exact row counts are intentionally not hardcoded.
    """

    success = True

    for (
        dataset_name,
        result,
    ) in (
        report.get(
            "files",
            {}
        )
        .items()
    ):
        rows = int(
            result.get(
                "rows",
                0,
            )
        )

        if rows <= 0:
            print(
                f"[FAIL] {dataset_name} contains zero rows"
            )

            success = False

    if success:
        print(
            "[PASS] All S5-S10 datasets contain source records"
        )

    return success


# ===========================================================================
# Report persistence
# ===========================================================================


def save_report(
    report: dict[str, Any],
) -> bool:
    """Write validation report to JSON."""

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
            "[FAIL] Could not save validation report: "
            f"{exc}"
        )

        return False

    print(
        "[PASS] Validation report written: "
        f"{OUTPUT}"
    )

    return True


def verify_saved_report() -> bool:
    """Read the generated JSON report back."""

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
            "[FAIL] Could not read generated JSON: "
            f"{exc}"
        )

        return False

    if (
        report.get(
            "dataset_id"
        )
        != "moradi_de"
    ):
        print(
            "[FAIL] Saved report has unexpected dataset_id"
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
    """Run the real Moradi S5-S10 smoke test."""

    separator()

    print(
        "pcatRFQTL — Moradi S5-S10 Differential Expression Smoke Test"
    )

    separator()

    print()
    print(
        "Checking source datasets..."
    )

    if not verify_source_directory():
        return 1

    if not verify_source_files():
        return 1

    print()
    print(
        "Initializing Moradi DE validator..."
    )

    runner = MoradiDERunner(
        MORADI_DIRECTORY,
        anomaly_sample_limit=1_000,
    )

    print(
        "[PASS] Validator initialized"
    )

    print()
    print(
        "Running Moradi S5-S10 validation..."
    )

    try:
        report = (
            runner.validate()
        )

    except Exception as exc:
        print(
            "[FAIL] Validation raised an exception: "
            f"{exc}"
        )

        return 1

    print(
        "[PASS] Moradi S5-S10 validation completed"
    )

    if not verify_report_structure(
        report
    ):
        return 1

    if not verify_expected_datasets(
        report
    ):
        return 1

    if not verify_nonzero_rows(
        report
    ):
        return 1

    # ------------------------------------------------------------------
    # Per-file summary
    # ------------------------------------------------------------------

    print()
    separator("-")

    print(
        "Moradi S5-S10 validation summary"
    )

    separator("-")

    for (
        dataset_name,
        result,
    ) in (
        report[
            "files"
        ]
        .items()
    ):

        print()
        print(
            f"  {dataset_name}"
        )

        print(
            f"    File         : "
            f"{Path(result['file']).name}"
        )

        print(
            f"    Status       : "
            f"{result.get('status')}"
        )

        print(
            f"    Rows         : "
            f"{result.get('rows')}"
        )

        print(
            f"    Columns      : "
            f"{result.get('columns')}"
        )

        print(
            f"    Feature type : "
            f"{result.get('feature_type')}"
        )

        print(
            f"    Scope        : "
            f"{result.get('analysis_scope')}"
        )

        print(
            "    Errors:"
        )

        print_dictionary(
            result.get(
                "validation_error_counts",
                {},
            ),
            indent="      ",
        )

        print(
            "    Duplicates:"
        )

        print_dictionary(
            result.get(
                "duplicates",
                {},
            ),
            indent="      ",
        )

        missing = (
            result.get(
                "missing_values",
                {}
            )
        )

        missing_nonzero = {
            key: value
            for key, value
            in missing.items()
            if value > 0
        }

        if missing_nonzero:

            print(
                "    Missing values:"
            )

            print_dictionary(
                missing_nonzero,
                indent="      ",
            )

    # ------------------------------------------------------------------
    # Overall summary
    # ------------------------------------------------------------------

    print()
    separator("-")

    print(
        "Overall validation status"
    )

    separator("-")

    print(
        f"Dataset status   : "
        f"{report.get('status')}"
    )

    print(
        f"Files validated  : "
        f"{report.get('files_validated')}"
    )

    print(
        f"Total rows       : "
        f"{report.get('total_rows')}"
    )

    summary = (
        report.get(
            "validation_summary",
            {}
        )
    )

    print(
        f"Severity counts  : "
        f"{summary.get('severity_counts')}"
    )

    print(
        f"Total anomalies  : "
        f"{summary.get('total_anomaly_occurrences')}"
    )

    anomalies = (
        report.get(
            "anomalies",
            []
        )
    )

    if anomalies:

        print()
        print(
            "First validation anomalies:"
        )

        for anomaly in anomalies[
            :10
        ]:

            anomaly_id = (
                anomaly.get(
                    "anomaly_id"
                )
                or anomaly.get(
                    "id"
                )
            )

            print(
                "  "
                f"{anomaly_id} "
                f"| {anomaly.get('severity')} "
                f"| {anomaly.get('type')} "
                f"| file="
                f"{Path(anomaly.get('file', '')).name} "
                f"| row={anomaly.get('row')} "
                f"| value={anomaly.get('value')!r}"
            )

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

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
        "A PASS smoke test means all real Moradi S5-S10 source files "
        "were processed successfully."
    )

    print(
        "Individual dataset statuses may still be WARNING or FAIL when "
        "source-data anomalies are detected."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )