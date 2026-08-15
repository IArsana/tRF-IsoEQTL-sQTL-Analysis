"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/test_m5_3c_3g_pairwise_cache.py

Description:
    Persistent-cache smoke test for M5.3C.3G.

    The script queries the same LDpair twice.

    Expected behavior:
        Query 1:
            remote provider if cache is initially empty.

        Query 2:
            persistent SQLite cache hit.

    Known test pair:
        rs11786992 ↔ rs28413800
        EAS

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
from pathlib import Path

from dotenv import load_dotenv

from pcatrfqtl.analysis.m5.cached_pairwise_ld_provider import (
    CachedPairwiseLDProvider,
)
from pcatrfqtl.analysis.m5.ldlink_pairwise_provider import (
    LDlinkPairwiseProvider,
    PROVIDER_NAME,
)
from pcatrfqtl.analysis.m5.pairwise_ld_cache import (
    PairwiseLDCache,
)


ROOT = Path(
    __file__
).resolve().parents[1]

ENV_PATH = (
    ROOT
    / ".env"
)

CACHE_PATH = (
    ROOT
    / "data/interim/m5/pairwise_ld/"
    "pairwise_ld_cache.sqlite"
)


DISEASE_PROXY = "rs11786992"

REGULATORY_PROXY = "rs28413800"

POPULATION = "EAS"


def main() -> None:
    """Run persistent pairwise LD cache smoke test."""

    load_dotenv(
        ENV_PATH,
        override=False,
    )

    token = os.environ.get(
        "LDLINK_TOKEN"
    )

    if not token:

        raise RuntimeError(
            "LDLINK_TOKEN was not found."
        )

    print()
    print("=" * 72)
    print("M5.3C.3G PAIRWISE LD PERSISTENT CACHE TEST")
    print("=" * 72)

    print(
        f"Pair:       {DISEASE_PROXY} ↔ {REGULATORY_PROXY}"
    )

    print(
        f"Population: {POPULATION}"
    )

    print(
        f"Cache:      {CACHE_PATH}"
    )

    print()

    with PairwiseLDCache(
        CACHE_PATH
    ) as cache:

        provider = LDlinkPairwiseProvider(
            token=token,
        )

        cached_provider = CachedPairwiseLDProvider(
            provider=provider,
            cache=cache,
            provider_name=PROVIDER_NAME,
        )

        before = cache.count()

        print(
            f"Cache records before: {before}"
        )

        print()
        print(
            "Query #1"
        )

        result_1 = cached_provider.query_pair(
            variant_a=DISEASE_PROXY,
            variant_b=REGULATORY_PROXY,
            population=POPULATION,
        )

        after_first = cache.count()

        print(
            result_1
        )

        print(
            f"Cache records after query #1: "
            f"{after_first}"
        )

        print()
        print(
            "Query #2"
        )

        result_2 = cached_provider.query_pair(
            variant_a=DISEASE_PROXY,
            variant_b=REGULATORY_PROXY,
            population=POPULATION,
        )

        after_second = cache.count()

        print(
            result_2
        )

        print(
            f"Cache records after query #2: "
            f"{after_second}"
        )

        assert (
            result_1.status
            == "SUCCESS"
        )

        assert (
            result_2.status
            == "SUCCESS"
        )

        assert (
            result_1.r2
            == result_2.r2
        )

        assert (
            result_1.d_prime
            == result_2.d_prime
        )

        assert (
            after_first
            == after_second
        ), (
            "Second identical query unexpectedly "
            "created another cache row."
        )

        cached_row = cache.get(
            provider=PROVIDER_NAME,
            population=POPULATION,
            variant_a=DISEASE_PROXY,
            variant_b=REGULATORY_PROXY,
        )

        assert (
            cached_row is not None
        )

        print()
        print(
            "Cached record:"
        )

        print(
            cached_row
        )

        print()
        print(
            "Cache by status:"
        )

        print(
            cache.count_by_status()
        )

        print()
        print(
            "Cache by population:"
        )

        print(
            cache.count_by_population()
        )

    print()
    print(
        "Assertions: PASS"
    )

    print("=" * 72)


if __name__ == "__main__":
    main()