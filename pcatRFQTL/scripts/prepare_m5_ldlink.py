"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/prepare_m5_ldlink.py

Description:
    Automated preparation of external LD evidence for M5.3 and M5.4.

    The script:

        1. Loads environment variables from the project-root .env file.
        2. Reads the LDlink API token from LDLINK_TOKEN.
        3. Queries LDproxy sequentially for all M5 lead SNPs.
        4. Preserves raw LDlink responses.
        5. Detects tabular and JSON provider responses.
        6. Handles population-specific unavailable variants without
           terminating the complete batch.
        7. Imports usable LDproxy output into the provider-neutral
           M5 LD schema.
        8. Writes one Parquet file per usable lead SNP and population.
        9. Produces a QC JSON report documenting successes and failures.

    Default analysis:
        Primary exploratory population: EAS

    Optional sensitivity populations:
        EUR
        SAS

    Genome build:
        GRCh37 / hg19

    Screening window:
        +/- 500 kb

    Operational LD thresholds:
        Primary:
            r2 >= 0.80

        Secondary:
            r2 >= 0.50

    Important safeguards:

        - A monoallelic variant in one reference population remains a
          valid lead candidate; it is recorded as unavailable for LD
          analysis in that population.

        - Raw provider error responses are retained for provenance.

        - No dummy or empty LD evidence is fabricated for unavailable
          lead variants.

        - Self-pairs are preserved in imported data for provenance but
          are excluded from downstream LD-mediated bridge eligibility.

        - This script does not perform GWAS matching, Moradi matching,
          colocalization, candidate ranking, or causal inference.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import argparse
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd
from dotenv import load_dotenv

from pcatrfqtl.io.parquet import write_parquet


# ============================================================================
# Project
# ============================================================================

PROJECT_ROOT = (
    Path(__file__)
    .resolve()
    .parents[1]
)


# ============================================================================
# Environment
# ============================================================================

ENV_PATH = (
    PROJECT_ROOT
    / ".env"
)

load_dotenv(
    dotenv_path=ENV_PATH,
    override=False,
)


# ============================================================================
# M5 configuration
# ============================================================================

LDLINK_URL = (
    "https://ldlink.nih.gov/LDlinkRest/ldproxy"
)

LEAD_SNPS = (
    "rs10216902",
    "rs1288100",
    "rs2328376",
)

PRIMARY_POPULATION = "EAS"

SENSITIVITY_POPULATIONS = (
    "EUR",
    "SAS",
)

DEFAULT_POPULATIONS = (
    PRIMARY_POPULATION,
)

GENOME_BUILD = "grch37"

REFERENCE_BUILD = "GRCh37"

REFERENCE_PANEL = (
    "1000_Genomes_Project"
)

PROVIDER = (
    "LDlink_LDproxy"
)

WINDOW_BP = 500_000

PRIMARY_R2_THRESHOLD = 0.80

SECONDARY_R2_THRESHOLD = 0.50

REQUEST_DELAY_SECONDS = 3.0

DEFAULT_TIMEOUT_SECONDS = 180


# ============================================================================
# Generic helpers
# ============================================================================


def first_existing_column(
    dataframe: pd.DataFrame,
    candidates: tuple[str, ...],
) -> str | None:
    """Return first compatible source column."""

    normalized = {
        str(column).strip().lower():
            column
        for column
        in dataframe.columns
    }

    for candidate in candidates:

        resolved = normalized.get(
            candidate.strip().lower()
        )

        if resolved is not None:
            return resolved

    return None


def normalize_rsid(
    series: pd.Series,
) -> pd.Series:
    """Normalize rsID strings conservatively."""

    return (
        series
        .astype("string")
        .str.strip()
        .str.lower()
    )


def classify_provider_failure(
    message: str,
) -> str:
    """Convert provider error text into stable pipeline categories."""

    normalized = (
        message
        .strip()
        .lower()
    )

    if "monoallelic" in normalized:
        return (
            "MONOALLELIC_IN_REFERENCE_POPULATION"
        )

    if (
        "not found"
        in normalized
        or "not present"
        in normalized
    ):
        return (
            "VARIANT_NOT_FOUND_IN_REFERENCE"
        )

    if (
        "invalid token"
        in normalized
        or "unauthorized"
        in normalized
    ):
        return (
            "LDLINK_AUTHENTICATION_ERROR"
        )

    if (
        "rate limit"
        in normalized
        or "too many requests"
        in normalized
    ):
        return (
            "LDLINK_RATE_LIMIT"
        )

    if (
        "biallelic"
        in normalized
        or "bi-allelic"
        in normalized
    ):
        return (
            "VARIANT_NOT_BIALLELIC"
        )

    return (
        "LDLINK_PROVIDER_ERROR"
    )


# ============================================================================
# LDlink URL and download
# ============================================================================


def build_ldproxy_url(
    *,
    snp: str,
    population: str,
    token: str,
) -> str:
    """Construct LDproxy GET URL."""

    parameters = {
        "var":
            snp,

        "pop":
            population,

        "r2_d":
            "r2",

        "window":
            str(
                WINDOW_BP
            ),

        "genome_build":
            GENOME_BUILD,

        "token":
            token,
    }

    return (
        LDLINK_URL
        + "?"
        + urllib.parse.urlencode(
            parameters
        )
    )


def download_ldproxy(
    *,
    snp: str,
    population: str,
    token: str,
    output_path: Path,
    timeout: int = DEFAULT_TIMEOUT_SECONDS,
    insecure: bool = False,
) -> None:
    """
    Download one LDproxy response.

    Requests are performed sequentially by the caller.
    """

    url = build_ldproxy_url(
        snp=snp,
        population=population,
        token=token,
    )

    request = urllib.request.Request(
        url,
        method="GET",
        headers={
            "User-Agent":
                "pcatRFQTL/1.0 "
                "M5-LD-preparation",
        },
    )

    context = None

    if insecure:

        context = (
            ssl._create_unverified_context()
        )

    try:

        with urllib.request.urlopen(
            request,
            timeout=timeout,
            context=context,
        ) as response:

            content = response.read()

    except urllib.error.HTTPError as exc:

        try:

            provider_body = (
                exc.read()
                .decode(
                    "utf-8",
                    errors="replace",
                )
                .strip()
            )

        except Exception:

            provider_body = ""

        message = (
            f"LDlink HTTP error for "
            f"{snp}/{population}: "
            f"{exc.code} {exc.reason}"
        )

        if provider_body:

            message += (
                f". Provider response: "
                f"{provider_body[:500]}"
            )

        raise RuntimeError(
            message
        ) from exc

    except urllib.error.URLError as exc:

        raise RuntimeError(
            "LDlink network error for "
            f"{snp}/{population}: "
            f"{exc.reason}"
        ) from exc

    if not content:

        raise RuntimeError(
            "LDlink returned an empty response "
            f"for {snp}/{population}."
        )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path.write_bytes(
        content
    )


# ============================================================================
# Raw response validation
# ============================================================================


def validate_raw_response(
    path: Path,
) -> dict[str, Any]:
    """
    Validate raw LDlink response.

    Successful LDproxy responses are tab-delimited.

    Provider validation failures may be returned as JSON. JSON failures
    are parsed so the exact provider message can be preserved in QC.
    """

    if not path.exists():

        raise FileNotFoundError(
            f"LDlink response does not exist: {path}"
        )

    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    stripped = text.strip()

    if not stripped:

        raise RuntimeError(
            f"Empty LDlink response: {path}"
        )

    lines = stripped.splitlines()

    first_line = (
        lines[0]
        if lines
        else ""
    )

    # ------------------------------------------------------------------
    # JSON response
    # ------------------------------------------------------------------

    if stripped.startswith(
        (
            "{",
            "[",
        )
    ):

        try:

            payload = json.loads(
                stripped
            )

        except json.JSONDecodeError as exc:

            raise RuntimeError(
                "LDlink returned a JSON-like response "
                "that could not be parsed for "
                f"{path.name}: "
                f"{stripped[:500]}"
            ) from exc

        provider_message = None

        if isinstance(
            payload,
            dict,
        ):

            for key in (
                "error",
                "message",
                "warning",
                "detail",
                "status",
            ):

                value = payload.get(
                    key
                )

                if value:

                    provider_message = (
                        f"{key}: {value}"
                    )

                    break

        if provider_message is None:

            provider_message = json.dumps(
                payload,
                ensure_ascii=False,
            )

        raise RuntimeError(
            "LDlink returned a non-tabular JSON response "
            f"for {path.name}: "
            f"{provider_message}"
        )

    # ------------------------------------------------------------------
    # Plain-text provider errors
    # ------------------------------------------------------------------

    lower = (
        stripped.lower()
    )

    suspicious_terms = (
        "invalid token",
        "error",
        "failed",
        "rate limit",
        "too many requests",
        "not found",
        "not a valid",
        "monoallelic",
    )

    detected_errors = [
        term
        for term
        in suspicious_terms
        if term
        in lower
    ]

    looks_tabular = (
        "\t"
        in first_line
    )

    if (
        detected_errors
        and not looks_tabular
    ):

        raise RuntimeError(
            "LDlink returned an apparent API error "
            f"for {path.name}: "
            f"{stripped[:500]}"
        )

    if not looks_tabular:

        raise RuntimeError(
            "LDlink response is not tab-delimited: "
            f"{path}. "
            "Response preview: "
            f"{stripped[:500]}"
        )

    return {
        "raw_file":
            str(
                path
            ),

        "bytes":
            path.stat().st_size,

        "lines":
            len(
                lines
            ),

        "header":
            first_line,

        "response_format":
            "tsv",
    }


# ============================================================================
# LDproxy importer
# ============================================================================


def import_ldproxy(
    path: Path,
    *,
    lead_rsid: str,
    population: str,
) -> pd.DataFrame:
    """
    Import one usable LDproxy response into the M5 provider-neutral schema.
    """

    source = pd.read_csv(
        path,
        sep="\t",
        dtype="object",
    )

    if source.empty:

        raise RuntimeError(
            f"Parsed LDproxy table is empty: {path}"
        )

    rsid_column = first_existing_column(
        source,
        (
            "RS_Number",
            "RS number",
            "RSNumber",
            "rs_number",
            "rsid",
        ),
    )

    r2_column = first_existing_column(
        source,
        (
            "R2",
            "R^2",
            "R²",
            "r2",
        ),
    )

    dprime_column = first_existing_column(
        source,
        (
            "Dprime",
            "D_Prime",
            "D prime",
            "D'",
            "dprime",
            "d_prime",
        ),
    )

    coordinate_column = first_existing_column(
        source,
        (
            "Coord",
            "Coordinate",
            "coord",
        ),
    )

    distance_column = first_existing_column(
        source,
        (
            "Distance",
            "distance",
        ),
    )

    if rsid_column is None:

        raise RuntimeError(
            "Could not resolve LDproxy rsID column. "
            "Observed columns: "
            f"{list(source.columns)}"
        )

    if r2_column is None:

        raise RuntimeError(
            "Could not resolve LDproxy R2 column. "
            "Observed columns: "
            f"{list(source.columns)}"
        )

    result = pd.DataFrame(
        index=source.index
    )

    result[
        "lead_rsid"
    ] = (
        str(
            lead_rsid
        )
        .strip()
        .lower()
    )

    result[
        "neighbor_rsid"
    ] = normalize_rsid(
        source[
            rsid_column
        ]
    )

    result[
        "r2"
    ] = pd.to_numeric(
        source[
            r2_column
        ],
        errors="coerce",
    )

    if dprime_column is not None:

        result[
            "d_prime"
        ] = pd.to_numeric(
            source[
                dprime_column
            ],
            errors="coerce",
        )

    else:

        result[
            "d_prime"
        ] = pd.Series(
            float("nan"),
            index=source.index,
            dtype="float64",
        )

    if coordinate_column is not None:

        result[
            "provider_coordinate"
        ] = (
            source[
                coordinate_column
            ]
            .astype(
                "string"
            )
        )

    else:

        result[
            "provider_coordinate"
        ] = pd.Series(
            pd.NA,
            index=source.index,
            dtype="string",
        )

    if distance_column is not None:

        result[
            "distance_bp"
        ] = (
            pd.to_numeric(
                source[
                    distance_column
                ],
                errors="coerce",
            )
            .astype(
                "Int64"
            )
        )

    else:

        result[
            "distance_bp"
        ] = pd.Series(
            pd.NA,
            index=source.index,
            dtype="Int64",
        )

    result[
        "population"
    ] = population

    result[
        "reference_panel"
    ] = REFERENCE_PANEL

    result[
        "reference_build"
    ] = REFERENCE_BUILD

    result[
        "provider"
    ] = PROVIDER

    result[
        "source_file"
    ] = path.name

    result[
        "source_row"
    ] = pd.Series(
        range(
            2,
            len(
                result
            )
            + 2,
        ),
        dtype="Int64",
    )

    result[
        "rsid_valid"
    ] = (
        result[
            "neighbor_rsid"
        ]
        .str.match(
            r"^rs\d+$",
            na=False,
        )
        .astype(
            "boolean"
        )
    )

    result[
        "r2_valid"
    ] = (
        result[
            "r2"
        ]
        .between(
            0.0,
            1.0,
            inclusive="both",
        )
        .fillna(
            False
        )
        .astype(
            "boolean"
        )
    )

    result[
        "is_self_pair"
    ] = (
        result[
            "lead_rsid"
        ]
        == result[
            "neighbor_rsid"
        ]
    ).astype(
        "boolean"
    )

    return result


# ============================================================================
# QC
# ============================================================================


def build_success_qc(
    dataframe: pd.DataFrame,
    *,
    snp: str,
    population: str,
    raw_qc: dict[str, Any],
    raw_path: Path,
    parquet_path: Path,
) -> dict[str, Any]:
    """Build QC record for successful LD evidence."""

    usable = (
        dataframe[
            "rsid_valid"
        ]
        .fillna(
            False
        )
        &
        dataframe[
            "r2_valid"
        ]
        .fillna(
            False
        )
    )

    non_self = (
        usable
        &
        ~dataframe[
            "is_self_pair"
        ]
        .fillna(
            False
        )
    )

    moderate = (
        non_self
        &
        (
            dataframe[
                "r2"
            ]
            >= SECONDARY_R2_THRESHOLD
        )
    )

    high = (
        non_self
        &
        (
            dataframe[
                "r2"
            ]
            >= PRIMARY_R2_THRESHOLD
        )
    )

    return {
        "lead_rsid":
            snp,

        "population":
            population,

        "status":
            "SUCCESS",

        "ld_reference_available":
            True,

        "failure_reason":
            None,

        "provider_error":
            None,

        "parquet_written":
            True,

        "raw":
            raw_qc,

        "raw_file":
            str(
                raw_path
            ),

        "output":
            str(
                parquet_path
            ),

        "rows":
            len(
                dataframe
            ),

        "unique_neighbor_rsids":
            int(
                dataframe[
                    "neighbor_rsid"
                ]
                .dropna()
                .nunique()
            ),

        "rsid_valid_rows":
            int(
                dataframe[
                    "rsid_valid"
                ]
                .fillna(
                    False
                )
                .sum()
            ),

        "r2_valid_rows":
            int(
                dataframe[
                    "r2_valid"
                ]
                .fillna(
                    False
                )
                .sum()
            ),

        "self_pair_rows":
            int(
                dataframe[
                    "is_self_pair"
                ]
                .fillna(
                    False
                )
                .sum()
            ),

        "secondary_ld_rows_r2_ge_0_5":
            int(
                moderate.sum()
            ),

        "primary_ld_rows_r2_ge_0_8":
            int(
                high.sum()
            ),

        "maximum_nonself_r2":
            (
                float(
                    dataframe.loc[
                        non_self,
                        "r2",
                    ]
                    .max()
                )
                if non_self.any()
                else None
            ),
    }


def build_failure_qc(
    *,
    snp: str,
    population: str,
    raw_path: Path,
    error_message: str,
) -> dict[str, Any]:
    """Build QC record for provider-unavailable LD evidence."""

    return {
        "lead_rsid":
            snp,

        "population":
            population,

        "status":
            "UNAVAILABLE",

        "ld_reference_available":
            False,

        "failure_reason":
            classify_provider_failure(
                error_message
            ),

        "provider_error":
            error_message,

        "parquet_written":
            False,

        "raw_file":
            (
                str(
                    raw_path
                )
                if raw_path.exists()
                else None
            ),

        "output":
            None,

        "rows":
            0,

        "unique_neighbor_rsids":
            0,

        "rsid_valid_rows":
            0,

        "r2_valid_rows":
            0,

        "self_pair_rows":
            0,

        "secondary_ld_rows_r2_ge_0_5":
            0,

        "primary_ld_rows_r2_ge_0_8":
            0,

        "maximum_nonself_r2":
            None,
    }


# ============================================================================
# Population processing
# ============================================================================


def process_population(
    *,
    population: str,
    token: str,
    force: bool,
    delay: float,
    insecure: bool,
) -> dict[str, Any]:
    """Download, import, and QC all lead SNPs for one population."""

    ld_root = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m5"
        / "ld"
    )

    raw_directory = (
        ld_root
        / "raw"
        / population
    )

    output_directory = (
        ld_root
        / population
    )

    raw_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    records: dict[
        str,
        Any,
    ] = {}

    for index, snp in enumerate(
        LEAD_SNPS
    ):

        print(
            "\n"
            f"[{population}] "
            f"{index + 1}/{len(LEAD_SNPS)} "
            f"{snp}"
        )

        raw_path = (
            raw_directory
            / f"{snp}.txt"
        )

        parquet_path = (
            output_directory
            / f"{snp}.parquet"
        )

        try:

            # ----------------------------------------------------------
            # Raw download / cache
            # ----------------------------------------------------------

            if (
                raw_path.exists()
                and not force
            ):

                print(
                    "  raw response already exists; "
                    "reusing cached file."
                )

            else:

                print(
                    "  querying LDlink LDproxy..."
                )

                download_ldproxy(
                    snp=snp,
                    population=population,
                    token=token,
                    output_path=raw_path,
                    insecure=insecure,
                )

                print(
                    f"  downloaded: {raw_path}"
                )

            # ----------------------------------------------------------
            # Provider response validation
            # ----------------------------------------------------------

            raw_qc = validate_raw_response(
                raw_path
            )

            print(
                "  response lines: "
                f"{raw_qc['lines']}"
            )

            # ----------------------------------------------------------
            # Parse usable response
            # ----------------------------------------------------------

            dataframe = import_ldproxy(
                raw_path,
                lead_rsid=snp,
                population=population,
            )

            # ----------------------------------------------------------
            # Write derived evidence
            # ----------------------------------------------------------

            write_parquet(
                dataframe,
                parquet_path,
                index=False,
            )

            qc = build_success_qc(
                dataframe,
                snp=snp,
                population=population,
                raw_qc=raw_qc,
                raw_path=raw_path,
                parquet_path=parquet_path,
            )

            records[
                snp
            ] = qc

            print(
                "  status: SUCCESS"
            )

            print(
                f"  rows: {qc['rows']}"
            )

            print(
                "  r2 >= 0.5: "
                f"{qc['secondary_ld_rows_r2_ge_0_5']}"
            )

            print(
                "  r2 >= 0.8: "
                f"{qc['primary_ld_rows_r2_ge_0_8']}"
            )

            print(
                f"  parquet: {parquet_path}"
            )

        except (
            RuntimeError,
            ValueError,
        ) as exc:

            error_message = str(
                exc
            )

            qc = build_failure_qc(
                snp=snp,
                population=population,
                raw_path=raw_path,
                error_message=error_message,
            )

            records[
                snp
            ] = qc

            # Do not leave stale derived evidence behind.
            if parquet_path.exists():

                parquet_path.unlink()

            print(
                "  status: UNAVAILABLE"
            )

            print(
                "  reason: "
                f"{qc['failure_reason']}"
            )

            print(
                "  provider/error: "
                f"{error_message}"
            )

            print(
                "  no Parquet evidence written."
            )

            print(
                "  continuing to next lead SNP."
            )

        # --------------------------------------------------------------
        # Sequential API delay
        # --------------------------------------------------------------

        if (
            index
            < len(
                LEAD_SNPS
            )
            - 1
        ):

            time.sleep(
                delay
            )

    # ------------------------------------------------------------------
    # Population-level summary
    # ------------------------------------------------------------------

    successful = [
        snp
        for snp, record
        in records.items()
        if record[
            "status"
        ]
        == "SUCCESS"
    ]

    unavailable = [
        snp
        for snp, record
        in records.items()
        if record[
            "status"
        ]
        == "UNAVAILABLE"
    ]

    failure_counts: dict[
        str,
        int,
    ] = {}

    for record in records.values():

        failure_reason = record.get(
            "failure_reason"
        )

        if failure_reason is None:
            continue

        failure_counts[
            failure_reason
        ] = (
            failure_counts.get(
                failure_reason,
                0,
            )
            + 1
        )

    return {
        "population":
            population,

        "reference_panel":
            REFERENCE_PANEL,

        "reference_build":
            REFERENCE_BUILD,

        "provider":
            PROVIDER,

        "window_bp":
            WINDOW_BP,

        "primary_r2_threshold":
            PRIMARY_R2_THRESHOLD,

        "secondary_r2_threshold":
            SECONDARY_R2_THRESHOLD,

        "summary": {
            "lead_snps_total":
                len(
                    LEAD_SNPS
                ),

            "successful_lead_snps":
                len(
                    successful
                ),

            "unavailable_lead_snps":
                len(
                    unavailable
                ),

            "successful_rsids":
                successful,

            "unavailable_rsids":
                unavailable,

            "failure_reason_counts":
                failure_counts,
        },

        "lead_snps":
            records,
    }


# ============================================================================
# CLI
# ============================================================================


def parse_arguments() -> argparse.Namespace:
    """Parse command-line options."""

    parser = argparse.ArgumentParser(
        description=(
            "Prepare LDlink LDproxy evidence "
            "for M5.3 and M5.4."
        )
    )

    parser.add_argument(
        "--populations",
        nargs="+",
        default=list(
            DEFAULT_POPULATIONS
        ),
        choices=[
            "EAS",
            "EUR",
            "SAS",
        ],
        help=(
            "1000 Genomes continental populations. "
            "Default: EAS."
        ),
    )

    parser.add_argument(
        "--with-sensitivity",
        action="store_true",
        help=(
            "Process EAS primary plus EUR and SAS "
            "sensitivity populations."
        ),
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Re-download raw LDlink responses "
            "even when cached files already exist."
        ),
    )

    parser.add_argument(
        "--delay",
        type=float,
        default=REQUEST_DELAY_SECONDS,
        help=(
            "Delay in seconds between sequential "
            "LDlink requests."
        ),
    )

    parser.add_argument(
        "--insecure",
        action="store_true",
        help=(
            "Disable HTTPS certificate verification. "
            "Use only when local TLS configuration "
            "requires it."
        ),
    )

    return parser.parse_args()


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Run automated M5 LD evidence preparation."""

    args = parse_arguments()

    # ------------------------------------------------------------------
    # Environment validation
    # ------------------------------------------------------------------

    if not ENV_PATH.exists():

        print(
            "WARNING: .env file not found at "
            f"{ENV_PATH}",
            file=sys.stderr,
        )

    token = os.getenv(
        "LDLINK_TOKEN"
    )

    if token is None:

        print(
            "ERROR: LDLINK_TOKEN was not found.",
            file=sys.stderr,
        )

        print(
            f"Expected .env file: {ENV_PATH}",
            file=sys.stderr,
        )

        print(
            "Example:",
            file=sys.stderr,
        )

        print(
            "LDLINK_TOKEN=YOUR_TOKEN",
            file=sys.stderr,
        )

        raise SystemExit(
            1
        )

    token = token.strip()

    if not token:

        print(
            "ERROR: LDLINK_TOKEN is empty.",
            file=sys.stderr,
        )

        raise SystemExit(
            1
        )

    # ------------------------------------------------------------------
    # Population selection
    # ------------------------------------------------------------------

    if args.with_sensitivity:

        populations = [
            PRIMARY_POPULATION,
            *SENSITIVITY_POPULATIONS,
        ]

    else:

        populations = list(
            dict.fromkeys(
                args.populations
            )
        )

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------

    print(
        "=" * 72
    )

    print(
        "M5 LDlink preparation"
    )

    print(
        "=" * 72
    )

    print(
        f"Project root : {PROJECT_ROOT}"
    )

    print(
        f".env path    : {ENV_PATH}"
    )

    print(
        "Token        : loaded "
        "(value hidden)"
    )

    print(
        "Lead SNPs    : "
        + ", ".join(
            LEAD_SNPS
        )
    )

    print(
        "Populations  : "
        + ", ".join(
            populations
        )
    )

    print(
        f"Build        : {REFERENCE_BUILD}"
    )

    print(
        f"Window       : +/- {WINDOW_BP:,} bp"
    )

    print(
        "Primary r2   : "
        f">= {PRIMARY_R2_THRESHOLD}"
    )

    print(
        "Secondary r2 : "
        f">= {SECONDARY_R2_THRESHOLD}"
    )

    # ------------------------------------------------------------------
    # Processing
    # ------------------------------------------------------------------

    reports: dict[
        str,
        Any,
    ] = {}

    for (
        population_index,
        population,
    ) in enumerate(
        populations
    ):

        reports[
            population
        ] = process_population(
            population=population,
            token=token,
            force=args.force,
            delay=args.delay,
            insecure=args.insecure,
        )

        if (
            population_index
            < len(
                populations
            )
            - 1
        ):

            time.sleep(
                args.delay
            )

    # ------------------------------------------------------------------
    # Global QC
    # ------------------------------------------------------------------

    qc_directory = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / "m5"
        / "qc"
    )

    qc_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    qc_path = (
        qc_directory
        / "m5_ldlink_preparation.json"
    )

    successful_total = sum(
        report[
            "summary"
        ][
            "successful_lead_snps"
        ]
        for report
        in reports.values()
    )

    unavailable_total = sum(
        report[
            "summary"
        ][
            "unavailable_lead_snps"
        ]
        for report
        in reports.values()
    )

    report = {
        "stage":
            "M5_external_ld_preparation",

        "provider":
            PROVIDER,

        "reference_panel":
            REFERENCE_PANEL,

        "reference_build":
            REFERENCE_BUILD,

        "window_bp":
            WINDOW_BP,

        "primary_population":
            PRIMARY_POPULATION,

        "sensitivity_populations":
            list(
                SENSITIVITY_POPULATIONS
            ),

        "primary_r2_threshold":
            PRIMARY_R2_THRESHOLD,

        "secondary_r2_threshold":
            SECONDARY_R2_THRESHOLD,

        "ld_calculated_by_pcatrfqtl":
            False,

        "external_ld_imported":
            True,

        "token_source":
            ".env:LDLINK_TOKEN",

        "token_value_stored_in_report":
            False,

        "population_specific_ld_availability":
            True,

        "unavailable_variant_is_not_removed_from_candidate_set":
            True,

        "dummy_ld_evidence_created":
            False,

        "populations_processed":
            populations,

        "summary": {
            "population_runs":
                len(
                    populations
                ),

            "lead_population_queries":
                (
                    len(
                        populations
                    )
                    * len(
                        LEAD_SNPS
                    )
                ),

            "successful_queries":
                successful_total,

            "unavailable_queries":
                unavailable_total,
        },

        "results":
            reports,
    }

    with qc_path.open(
        "w",
        encoding="utf-8",
    ) as handle:

        json.dump(
            report,
            handle,
            indent=2,
            ensure_ascii=False,
        )

    # ------------------------------------------------------------------
    # Final console output
    # ------------------------------------------------------------------

    print(
        "\n"
        + "=" * 72
    )

    print(
        "LD preparation complete."
    )

    print(
        "=" * 72
    )

    print(
        "Successful queries : "
        f"{successful_total}"
    )

    print(
        "Unavailable queries: "
        f"{unavailable_total}"
    )

    print(
        f"QC report          : {qc_path}"
    )

    print(
        "=" * 72
    )


if __name__ == "__main__":
    main()