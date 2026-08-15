"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/probe_m5_3c_3g_ldpair.py

Description:
    Single-pair diagnostic probe for M5.3C.3G.

    This script queries NIH LDlink LDpair for one disease proxy and one
    regulatory proxy across EAS, EUR, and SAS reference populations.

    The purpose is to inspect the actual live LDpair JSON response structure
    before enabling the full disease-proxy × regulatory-proxy pairwise LD
    analysis.

    Default pair:
        Disease proxy:
            rs11786992

        Regulatory proxy:
            rs28413800

    Genome build:
        GRCh37

    Environment:
        LDLINK_TOKEN is loaded from:
            <project_root>/.env

    Security:
        - The LDlink token is never printed.
        - The token is never written to output files.
        - Error text is sanitized before display.

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
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from dotenv import load_dotenv


# ============================================================================
# Project paths
# ============================================================================


ROOT = Path(
    __file__
).resolve().parents[1]

ENV_PATH = (
    ROOT
    / ".env"
)


# ============================================================================
# LDlink configuration
# ============================================================================


LDLINK_LDPAIR_URL = (
    "https://ldlink.nih.gov/LDlinkRest/ldpair"
)

DEFAULT_DISEASE_PROXY = (
    "rs11786992"
)

DEFAULT_REGULATORY_PROXY = (
    "rs28413800"
)

DEFAULT_POPULATIONS = (
    "EAS",
    "EUR",
    "SAS",
)

GENOME_BUILD = (
    "grch37"
)

TIMEOUT_SECONDS = 120

MAX_ATTEMPTS = 3

RETRYABLE_HTTP_CODES = {
    429,
    500,
    502,
    503,
    504,
}


# ============================================================================
# Environment
# ============================================================================


def load_project_environment() -> None:
    """
    Load the project-level .env file.

    The file is resolved relative to the repository root rather than the
    current working directory, so the script can be executed from elsewhere.
    """

    if not ENV_PATH.exists():

        print(
            f"ERROR: .env file not found: {ENV_PATH}",
            file=sys.stderr,
        )

        raise SystemExit(
            1
        )

    loaded = load_dotenv(
        dotenv_path=ENV_PATH,
        override=False,
    )

    if not loaded:

        print(
            f"WARNING: .env could not be loaded: {ENV_PATH}",
            file=sys.stderr,
        )


def get_ldlink_token() -> str:
    """Return LDlink API token loaded from the project environment."""

    token = os.environ.get(
        "LDLINK_TOKEN"
    )

    if token is None:

        print(
            "ERROR: LDLINK_TOKEN was not found in the project .env file.",
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

    return token


# ============================================================================
# Validation helpers
# ============================================================================


def normalize_rsid(
    value: str,
) -> str:
    """Normalize and validate one canonical rsID."""

    normalized = (
        str(
            value
        )
        .strip()
        .lower()
    )

    if not normalized.startswith(
        "rs"
    ):

        raise ValueError(
            f"Invalid rsID: {value}"
        )

    numeric_part = (
        normalized[
            2:
        ]
    )

    if not numeric_part.isdigit():

        raise ValueError(
            f"Invalid rsID: {value}"
        )

    return normalized


def normalize_population(
    population: str,
) -> str:
    """Normalize and validate broad 1000 Genomes population group."""

    normalized = (
        str(
            population
        )
        .strip()
        .upper()
    )

    if normalized not in {
        "EAS",
        "EUR",
        "SAS",
    }:

        raise ValueError(
            "Unsupported population: "
            f"{population}"
        )

    return normalized


def redact_token_from_text(
    text: str,
    token: str,
) -> str:
    """Remove the API token from diagnostic text."""

    if not text:

        return text

    if not token:

        return text

    return text.replace(
        token,
        "[REDACTED]",
    )


# ============================================================================
# HTTP helpers
# ============================================================================


def build_ldpair_request(
    *,
    variant_a: str,
    variant_b: str,
    population: str,
    token: str,
) -> urllib.request.Request:
    """Build one LDlink LDpair POST request."""

    variant_a = normalize_rsid(
        variant_a
    )

    variant_b = normalize_rsid(
        variant_b
    )

    population = normalize_population(
        population
    )

    url = (
        LDLINK_LDPAIR_URL
        + "?"
        + urllib.parse.urlencode(
            {
                "token":
                    token,
            }
        )
    )

    payload = {
        "snp_pairs": [
            [
                variant_a,
                variant_b,
            ]
        ],
        "pop":
            population,
        "genome_build":
            GENOME_BUILD,
        "json_out":
            True,
    }

    body = json.dumps(
        payload
    ).encode(
        "utf-8"
    )

    return urllib.request.Request(
        url=url,
        data=body,
        method="POST",
        headers={
            "Content-Type":
                "application/json",
            "Accept":
                "application/json",
            "User-Agent":
                "pcatRFQTL/1.0",
        },
    )


def query_ldpair(
    *,
    variant_a: str,
    variant_b: str,
    population: str,
    token: str,
) -> tuple[Any, int]:
    """
    Query one disease-proxy ↔ regulatory-proxy pair.

    Returns:
        parsed response
        number of HTTP attempts
    """

    last_error: Exception | None = None

    for attempt in range(
        1,
        MAX_ATTEMPTS + 1,
    ):

        request = build_ldpair_request(
            variant_a=variant_a,
            variant_b=variant_b,
            population=population,
            token=token,
        )

        try:

            with urllib.request.urlopen(
                request,
                timeout=TIMEOUT_SECONDS,
            ) as response:

                raw_bytes = (
                    response.read()
                )

                status_code = (
                    response.status
                )

            raw_text = raw_bytes.decode(
                "utf-8",
                errors="replace",
            )

            if status_code != 200:

                raise RuntimeError(
                    "Unexpected HTTP status: "
                    f"{status_code}"
                )

            try:

                parsed = json.loads(
                    raw_text
                )

            except json.JSONDecodeError as exc:

                safe_raw = redact_token_from_text(
                    raw_text,
                    token,
                )

                raise RuntimeError(
                    "LDlink returned non-JSON output:\n"
                    f"{safe_raw}"
                ) from exc

            return (
                parsed,
                attempt,
            )

        except urllib.error.HTTPError as exc:

            last_error = exc

            try:

                response_body = (
                    exc.read()
                    .decode(
                        "utf-8",
                        errors="replace",
                    )
                )

            except Exception:

                response_body = ""

            response_body = (
                redact_token_from_text(
                    response_body,
                    token,
                )
            )

            print(
                f"HTTP error {exc.code} "
                f"on attempt {attempt}."
            )

            if response_body:

                print(
                    "Provider response:"
                )

                print(
                    response_body
                )

            if (
                exc.code
                not in RETRYABLE_HTTP_CODES
            ):

                raise RuntimeError(
                    "LDlink returned a non-retryable "
                    f"HTTP {exc.code} response."
                ) from exc

        except (
            urllib.error.URLError,
            TimeoutError,
            ConnectionError,
        ) as exc:

            last_error = exc

            safe_error = (
                redact_token_from_text(
                    str(
                        exc
                    ),
                    token,
                )
            )

            print(
                "Transient network error "
                f"on attempt {attempt}: "
                f"{safe_error}"
            )

        if attempt < MAX_ATTEMPTS:

            delay = (
                2 ** (
                    attempt - 1
                )
            )

            print(
                f"Retrying in {delay} second(s)..."
            )

            time.sleep(
                delay
            )

    safe_last_error = (
        redact_token_from_text(
            str(
                last_error
            ),
            token,
        )
    )

    raise RuntimeError(
        "LDlink request failed after "
        f"{MAX_ATTEMPTS} attempts: "
        f"{safe_last_error}"
    )


# ============================================================================
# Response diagnostics
# ============================================================================


def inspect_response(
    response: Any,
) -> None:
    """Print the structure and parsed JSON response."""

    print(
        "Python response type: "
        f"{type(response).__name__}"
    )

    if isinstance(
        response,
        list,
    ):

        print(
            "Top-level list items: "
            f"{len(response)}"
        )

        if response:

            first = response[
                0
            ]

            print(
                "First item type: "
                f"{type(first).__name__}"
            )

            if isinstance(
                first,
                dict,
            ):

                print(
                    "First item keys:"
                )

                for key in first.keys():

                    print(
                        f"  - {key}"
                    )

    elif isinstance(
        response,
        dict,
    ):

        print(
            "Top-level keys:"
        )

        for key in response.keys():

            print(
                f"  - {key}"
            )

    else:

        print(
            "Response is neither a "
            "dictionary nor a list."
        )

    print()

    print(
        "Raw parsed JSON:"
    )

    print(
        json.dumps(
            response,
            indent=2,
            ensure_ascii=False,
        )
    )


# ============================================================================
# CLI
# ============================================================================


def build_parser() -> argparse.ArgumentParser:
    """Build command-line interface."""

    parser = argparse.ArgumentParser(
        description=(
            "Probe NIH LDlink LDpair for one "
            "disease-regulatory variant pair."
        )
    )

    parser.add_argument(
        "--disease",
        default=DEFAULT_DISEASE_PROXY,
        help=(
            "Disease proxy rsID. "
            f"Default: {DEFAULT_DISEASE_PROXY}"
        ),
    )

    parser.add_argument(
        "--regulatory",
        default=DEFAULT_REGULATORY_PROXY,
        help=(
            "Regulatory proxy rsID. "
            f"Default: {DEFAULT_REGULATORY_PROXY}"
        ),
    )

    parser.add_argument(
        "--population",
        default=None,
        choices=[
            "EAS",
            "EUR",
            "SAS",
        ],
        help=(
            "Query only one population. "
            "If omitted, EAS, EUR, and SAS "
            "are queried sequentially."
        ),
    )

    return parser


# ============================================================================
# Main
# ============================================================================


def main() -> None:
    """Run the M5.3C.3G LDpair diagnostic probe."""

    load_project_environment()

    token = get_ldlink_token()

    parser = build_parser()

    args = parser.parse_args()

    disease_proxy = normalize_rsid(
        args.disease
    )

    regulatory_proxy = normalize_rsid(
        args.regulatory
    )

    if args.population is not None:

        populations = (
            normalize_population(
                args.population
            ),
        )

    else:

        populations = (
            DEFAULT_POPULATIONS
        )

    print()

    print(
        "=" * 72
    )

    print(
        "M5.3C.3G LDLINK LDPAIR LIVE PROBE"
    )

    print(
        "=" * 72
    )

    print(
        f"Project root:       {ROOT}"
    )

    print(
        f".env:              {ENV_PATH}"
    )

    print(
        f"Disease proxy:     {disease_proxy}"
    )

    print(
        f"Regulatory proxy:  {regulatory_proxy}"
    )

    print(
        f"Genome build:      {GENOME_BUILD}"
    )

    print(
        "Populations:       "
        + ", ".join(
            populations
        )
    )

    print(
        "LDLINK_TOKEN:      [LOADED / REDACTED]"
    )

    print(
        "=" * 72
    )

    for index, population in enumerate(
        populations,
        start=1,
    ):

        if index > 1:

            # Requests are deliberately sequential.
            time.sleep(
                1
            )

        print()

        print(
            "-" * 72
        )

        print(
            f"Population: {population}"
        )

        print(
            "-" * 72
        )

        try:

            response, attempts = query_ldpair(
                variant_a=disease_proxy,
                variant_b=regulatory_proxy,
                population=population,
                token=token,
            )

        except Exception as exc:

            safe_error = redact_token_from_text(
                str(
                    exc
                ),
                token,
            )

            print(
                f"FAILED: {safe_error}"
            )

            continue

        print(
            f"HTTP attempts: {attempts}"
        )

        inspect_response(
            response
        )

    print()

    print(
        "=" * 72
    )

    print(
        "Probe complete."
    )

    print(
        "=" * 72
    )


if __name__ == "__main__":
    main()