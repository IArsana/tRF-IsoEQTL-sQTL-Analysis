"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_gwas.py

Description:
    Smoke test for the full GWAS Catalog association validation pipeline.

    The script validates the real GWAS Catalog association TSV using
    chunked processing and verifies that:

        - The source file exists.
        - The required schema is present.
        - All records can be processed.
        - Validation statistics are generated.
        - AnomalyFactory integration works.
        - The QC report is JSON serializable.
        - The saved JSON report can be read back successfully.

    Raw GWAS Catalog source data are never modified.

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

from pcatrfqtl.validation.runners.gwas_catalog import (
    GWASCatalogValidator,
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

GWAS_FILE = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "gwas"
    / "gwas-catalog-download-associations-alt-full.tsv"
)

OUTPUT = (
    PROJECT_ROOT
    / "data"
    / "interim"
    / "qc"
    / "gwas_catalog_validation.json"
)


EXPECTED_COLUMNS = {
    "SNPS",
    "CHR_ID",
    "CHR_POS",
    "P-VALUE",
}


# ===========================================================================
# Display
# ===========================================================================


def separator(
    character: str = "=",
    width: int = 80,
) -> None:
    """Print a visual separator."""

    print(
        character * width
    )


def print_dictionary(
    values: dict[str, Any],
    *,
    indent: str = "  ",
) -> None:
    """Print dictionary values line-by-line."""

    for key, value in values.items():
        print(
            f"{indent}{key}: {value}"
        )


# ===========================================================================
# Checks
# ===========================================================================


def verify_source_file() -> bool:
    """Verify the GWAS Catalog TSV exists."""

    if not GWAS_FILE.exists():
        print(
            "[FAIL] GWAS Catalog file missing: "
            f"{GWAS_FILE}"
        )
        return False

    if not GWAS_FILE.is_file():
        print(
            "[FAIL] GWAS Catalog path is not a file: "
            f"{GWAS_FILE}"
        )
        return False

    size_mb = (
        GWAS_FILE.stat().st_size
        / (1024 ** 2)
    )

    print(
        "[PASS] GWAS Catalog source exists"
    )

    print(
        f"       File : {GWAS_FILE.name}"
    )

    print(
        f"       Size : {size_mb:.2f} MB"
    )

    return True


def validate_report_structure(
    report: dict[str, Any],
) -> bool:
    """Validate required QC report keys."""

    required = {
        "dataset_id",
        "file",
        "status",
        "rows",
        "columns",
        "schema",
        "validation_error_counts",
        "missing_values",
        "mapping_statistics",
        "duplicates",
        "anomalies",
        "validation_summary",
    }

    missing = (
        required
        - set(report)
    )

    if missing:
        print(
            "[FAIL] Missing report keys: "
            + ", ".join(
                sorted(missing)
            )
        )
        return False

    if (
        report.get(
            "dataset_id"
        )
        != "gwas_catalog"
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


def validate_expected_schema(
    report: dict[str, Any],
) -> bool:
    """Verify required GWAS columns are present."""

    available = set(
        report.get(
            "schema",
            {},
        ).get(
            "available_columns",
            [],
        )
    )

    missing = (
        EXPECTED_COLUMNS
        - available
    )

    if missing:
        print(
            "[FAIL] Expected GWAS columns missing: "
            + ", ".join(
                sorted(missing)
            )
        )
        return False

    print(
        "[PASS] Required GWAS Catalog schema is present"
    )

    return True


def save_report(
    report: dict[str, Any],
) -> bool:
    """Save QC report as JSON."""

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
            "[FAIL] Could not save QC report: "
            f"{exc}"
        )
        return False

    print(
        "[PASS] Validation report written: "
        f"{OUTPUT}"
    )

    return True


def verify_saved_report() -> bool:
    """Read generated JSON report back."""

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
        != "gwas_catalog"
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
    """Run the real GWAS Catalog smoke test."""

    separator()

    print(
        "pcatRFQTL — GWAS Catalog Smoke Test"
    )

    separator()

    print()
    print(
        "Checking source dataset..."
    )

    if not verify_source_file():
        return 1

    print()
    print(
        "Initializing chunked GWAS Catalog validator..."
    )

    validator = GWASCatalogValidator(
        GWAS_FILE,
        chunk_size=100_000,
        anomaly_sample_limit=1_000,
    )

    print(
        "[PASS] Validator initialized"
    )

    print()
    print(
        "Running GWAS Catalog validation..."
    )

    print(
        "The complete TSV will be processed in "
        "100,000-row chunks."
    )

    try:
        report = (
            validator.validate()
        )

    except Exception as exc:
        print(
            "[FAIL] Validation raised an exception: "
            f"{exc}"
        )
        return 1

    print(
        "[PASS] GWAS Catalog validation completed"
    )

    if not validate_report_structure(
        report
    ):
        return 1

    if not validate_expected_schema(
        report
    ):
        return 1

    # ------------------------------------------------------------------
    # Dataset summary
    # ------------------------------------------------------------------

    print()
    separator("-")

    print(
        "GWAS Catalog validation summary"
    )

    separator("-")

    print(
        f"Status           : "
        f"{report.get('status')}"
    )

    print(
        f"Rows             : "
        f"{report.get('rows')}"
    )

    print(
        f"Columns          : "
        f"{report.get('columns')}"
    )

    print(
        f"Chunks processed : "
        f"{report.get('chunks_processed')}"
    )

    print()
    print(
        "Validation errors:"
    )

    print_dictionary(
        report.get(
            "validation_error_counts",
            {},
        )
    )

    print()
    print(
        "Genomic mapping:"
    )

    print_dictionary(
        report.get(
            "mapping_statistics",
            {},
        )
    )

    print()
    print(
        "Duplicates:"
    )

    print_dictionary(
        report.get(
            "duplicates",
            {},
        )
    )

    print()
    print(
        "Validation summary:"
    )

    summary = report.get(
        "validation_summary",
        {},
    )

    print_dictionary(
        summary
    )

    anomalies = report.get(
        "anomalies",
        [],
    )

    if anomalies:
        print()
        print(
            "First validation anomalies:"
        )

        for anomaly in anomalies[:10]:
            print(
                "  "
                f"{anomaly.get('id') or anomaly.get('anomaly_id')} "
                f"| {anomaly.get('severity')} "
                f"| {anomaly.get('type')} "
                f"| row={anomaly.get('row')} "
                f"| value={anomaly.get('value')!r}"
            )

    # ------------------------------------------------------------------
    # Save report
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
        "A WARNING or FAIL dataset status does not mean "
        "the smoke test failed."
    )

    print(
        "The smoke test confirms that the complete GWAS "
        "Catalog source can be processed reproducibly and "
        "that detected data-quality findings are recorded."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )