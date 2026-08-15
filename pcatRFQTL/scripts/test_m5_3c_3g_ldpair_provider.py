"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/test_m5_3c_3g_ldpair_provider.py

Description:
    Direct smoke test for the M5.3C.3G NIH LDlink LDpair provider.

    This script:
        - Loads LDLINK_TOKEN from the project-level .env file.
        - Instantiates LDlinkPairwiseProvider.
        - Queries one known disease-proxy ↔ regulatory-proxy pair.
        - Prints the parsed PairwiseLDRecord.
        - Performs basic assertions against the known live response structure.

    Default test:
        Disease proxy:
            rs11786992

        Regulatory proxy:
            rs28413800

        Population:
            EAS

        Genome build:
            GRCh37

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

from dotenv import load_dotenv

from pcatrfqtl.analysis.m5.ldlink_pairwise_provider import (
    LDlinkPairwiseProvider,
)


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
# Test configuration
# ============================================================================


DISEASE_PROXY = (
    "rs11786992"
)

REGULATORY_PROXY = (
    "rs28413800"
)

POPULATION = (
    "EAS"
)


# ============================================================================
# Environment
# ============================================================================


def load_environment() -> str:
    """Load and return LDLINK_TOKEN from the project .env file."""

    if not ENV_PATH.exists():

        print(
            f"ERROR: .env file not found: {ENV_PATH}",
            file=sys.stderr,
        )

        raise SystemExit(
            1
        )

    load_dotenv(
        dotenv_path=ENV_PATH,
        override=False,
    )

    token = os.environ.get(
        "LDLINK_TOKEN"
    )

    if token is None:

        print(
            "ERROR: LDLINK_TOKEN was not found in .env.",
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
# Main
# ============================================================================


def main() -> None:
    """Run direct LDlink provider smoke test."""

    token = load_environment()

    print()
    print(
        "=" * 72
    )

    print(
        "M5.3C.3G LDLINK PAIRWISE PROVIDER TEST"
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
        f"Disease proxy:     {DISEASE_PROXY}"
    )

    print(
        f"Regulatory proxy:  {REGULATORY_PROXY}"
    )

    print(
        f"Population:        {POPULATION}"
    )

    print(
        "LDLINK_TOKEN:      [LOADED / REDACTED]"
    )

    print(
        "=" * 72
    )

    provider = LDlinkPairwiseProvider(
        token=token,
    )

    print()
    print(
        "Querying LDlink..."
    )

    result = provider.query_pair(
        variant_a=DISEASE_PROXY,
        variant_b=REGULATORY_PROXY,
        population=POPULATION,
    )

    print()
    print(
        "Parsed PairwiseLDRecord:"
    )

    print(
        result
    )

    print()
    print(
        "Parsed fields:"
    )

    print(
        f"  variant_a : {result.variant_a}"
    )

    print(
        f"  variant_b : {result.variant_b}"
    )

    print(
        f"  population: {result.population}"
    )

    print(
        f"  status    : {result.status}"
    )

    print(
        f"  r2        : {result.r2}"
    )

    print(
        f"  d_prime   : {result.d_prime}"
    )

    print(
        f"  provider  : {result.provider}"
    )

    print(
        f"  reason    : {result.reason}"
    )

    # ======================================================================
    # Basic assertions
    # ======================================================================

    assert (
        result.variant_a
        == DISEASE_PROXY
    ), (
        "Unexpected variant_a: "
        f"{result.variant_a}"
    )

    assert (
        result.variant_b
        == REGULATORY_PROXY
    ), (
        "Unexpected variant_b: "
        f"{result.variant_b}"
    )

    assert (
        result.population
        == POPULATION
    ), (
        "Unexpected population: "
        f"{result.population}"
    )

    assert (
        result.status
        == "SUCCESS"
    ), (
        "Expected SUCCESS, received "
        f"{result.status}: "
        f"{result.reason}"
    )

    assert (
        result.r2 is not None
    ), (
        "Expected numeric r2."
    )

    assert (
        result.d_prime is not None
    ), (
        "Expected numeric D'."
    )

    print()
    print(
        "Assertions: PASS"
    )

    print(
        "=" * 72
    )


if __name__ == "__main__":
    main()