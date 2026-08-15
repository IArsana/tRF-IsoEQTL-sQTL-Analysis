"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/test_m5_3c_3g_batch_cache.py

Description:
    Batch-aware cache smoke test for M5.3C.3G.

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

load_dotenv(
    ROOT / ".env",
    override=False,
)

CACHE_PATH = (
    ROOT
    / "data/interim/m5/pairwise_ld/"
    "pairwise_ld_cache.sqlite"
)

PAIRS = [
    (
        "rs11786992",
        "rs28413800",
    ),
    (
        "rs576286012",
        "rs28413800",
    ),
    (
        "rs67763258",
        "rs28413800",
    ),
]

POPULATION = "EAS"


def main() -> None:
    """Run batch-aware cache test."""

    token = os.environ.get(
        "LDLINK_TOKEN"
    )

    if not token:

        raise RuntimeError(
            "LDLINK_TOKEN is missing."
        )

    print()
    print("=" * 72)
    print("M5.3C.3G BATCH-AWARE CACHE TEST")
    print("=" * 72)

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
            f"Cache before: {before}"
        )

        results = cached_provider.query_pairs(
            pairs=PAIRS,
            population=POPULATION,
        )

        print()
        print("Results:")

        for result in results:

            print(
                result
            )

        after = cache.count()

        print()
        print(
            f"Cache after: {after}"
        )

        print()
        print(
            "Diagnostics:"
        )

        print(
            cached_provider.diagnostics()
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
            "Cache by status:"
        )

        print(
            cache.count_by_status()
        )

        assert (
            len(
                results
            )
            == len(
                PAIRS
            )
        )

        assert all(
            result.population
            == POPULATION
            for result in results
        )

    print()
    print("Assertions: PASS")
    print("=" * 72)


if __name__ == "__main__":
    main()