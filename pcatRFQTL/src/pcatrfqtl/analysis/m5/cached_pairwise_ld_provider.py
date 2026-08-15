"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/cached_pairwise_ld_provider.py

Description:
    Batch-aware persistent caching wrapper for M5.3C.3G pairwise LD queries.

    The wrapper sits between the analysis runner and a remote pairwise LD
    provider.

    Query behavior:

        CACHE HIT
            Return the previously persisted scientific result without
            contacting the remote provider.

        CACHE MISS
            Collect missing pairs and query them through the underlying
            provider's native query_pairs() method.

            For LDlink, the underlying provider internally batches requests
            according to its API limit.

    Terminal states persisted:
        SUCCESS
        VARIANT_UNAVAILABLE
        MONOALLELIC
        PAIR_UNAVAILABLE

    Non-terminal states are not persisted:
        RATE_LIMITED
        REMOTE_ERROR
        TIMEOUT
        PARSE_ERROR
        PAIR_MISMATCH

    Cache identity:
        provider
        population
        unordered canonical rsID pair

    Output guarantees:
        - Result cardinality equals request cardinality.
        - Result ordering follows the caller's requested pair order.
        - Caller-requested pair orientation is preserved.
        - Duplicate requested pairs are supported.
        - Cached unavailable results remain unavailable, never r² = 0.

    Diagnostics:
        The wrapper records cumulative counters for:
            requested_pairs
            cache_hits
            cache_misses
            remote_pairs
            remote_calls
            newly_cached_results
            nonterminal_results

    Scientific safeguards:
        - Low LD and unavailable LD remain distinct.
        - Cache reuse does not alter scientific values.
        - Pairwise LD does not establish colocalization or causality.
        - Transient provider failures are allowed to retry on later runs.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from pcatrfqtl.analysis.m5.pairwise_ld_cache import (
    PairwiseLDCache,
    canonical_pair,
    normalize_population,
    normalize_rsid,
)
from pcatrfqtl.analysis.m5.pairwise_ld_provider import (
    PairwiseLDProvider,
    PairwiseLDRecord,
)


# ============================================================================
# Cache policy
# ============================================================================


TERMINAL_CACHE_STATUSES = {
    "SUCCESS",
    "VARIANT_UNAVAILABLE",
    "MONOALLELIC",
    "PAIR_UNAVAILABLE",
}


NON_TERMINAL_STATUSES = {
    "RATE_LIMITED",
    "REMOTE_ERROR",
    "TIMEOUT",
    "PARSE_ERROR",
    "PAIR_MISMATCH",
}


# ============================================================================
# Native batch-provider protocol
# ============================================================================


@runtime_checkable
class NativeBatchPairwiseLDProvider(Protocol):
    """
    Optional protocol for providers supporting native multi-pair queries.

    LDlinkPairwiseProvider implements this through query_pairs().
    """

    def query_pairs(
        self,
        *,
        pairs: list[tuple[str, str]],
        population: str,
    ) -> list[PairwiseLDRecord]:
        """Query multiple pairwise LD combinations."""
        ...


# ============================================================================
# Diagnostics
# ============================================================================


@dataclass
class CachedProviderStats:
    """Cumulative cached-provider diagnostics."""

    requested_pairs: int = 0

    cache_hits: int = 0

    cache_misses: int = 0

    remote_pairs: int = 0

    remote_calls: int = 0

    newly_cached_results: int = 0

    nonterminal_results: int = 0

    def as_dict(
        self,
    ) -> dict[str, int]:
        """Return diagnostics as plain dictionary."""

        return {
            "requested_pairs":
                int(
                    self.requested_pairs
                ),

            "cache_hits":
                int(
                    self.cache_hits
                ),

            "cache_misses":
                int(
                    self.cache_misses
                ),

            "remote_pairs":
                int(
                    self.remote_pairs
                ),

            "remote_calls":
                int(
                    self.remote_calls
                ),

            "newly_cached_results":
                int(
                    self.newly_cached_results
                ),

            "nonterminal_results":
                int(
                    self.nonterminal_results
                ),
        }

    def reset(
        self,
    ) -> None:
        """Reset all diagnostics to zero."""

        self.requested_pairs = 0

        self.cache_hits = 0

        self.cache_misses = 0

        self.remote_pairs = 0

        self.remote_calls = 0

        self.newly_cached_results = 0

        self.nonterminal_results = 0


# ============================================================================
# Internal helpers
# ============================================================================


def _normalize_status(
    status: str,
) -> str:
    """Normalize result status."""

    return (
        str(status)
        .strip()
        .upper()
    )


def _request_key(
    *,
    variant_a: str,
    variant_b: str,
) -> tuple[str, str]:
    """
    Return canonical unordered pair key.

    Used only for cache/dedup identity.
    """

    return canonical_pair(
        variant_a,
        variant_b,
    )


def _caller_oriented_record(
    *,
    record: PairwiseLDRecord,
    requested_variant_a: str,
    requested_variant_b: str,
    population: str,
) -> PairwiseLDRecord:
    """
    Rebuild one record in caller-requested pair orientation.

    Pairwise LD is symmetric, but retaining caller orientation simplifies
    downstream disease_proxy ↔ regulatory_proxy provenance.
    """

    return PairwiseLDRecord(
        variant_a=requested_variant_a,
        variant_b=requested_variant_b,
        population=population,
        r2=record.r2,
        d_prime=record.d_prime,
        status=record.status,
        provider=record.provider,
        reason=record.reason,
    )


# ============================================================================
# Cached provider
# ============================================================================


class CachedPairwiseLDProvider:
    """Batch-aware persistent cache wrapper for pairwise LD providers."""

    def __init__(
        self,
        *,
        provider: PairwiseLDProvider,
        cache: PairwiseLDCache,
        provider_name: str,
    ) -> None:

        self.provider = provider

        self.cache = cache

        self.provider_name = (
            str(provider_name)
            .strip()
        )

        if not self.provider_name:

            raise ValueError(
                "provider_name cannot be empty."
            )

        self.stats = CachedProviderStats()

    # ======================================================================
    # Cache policy
    # ======================================================================

    def _should_cache(
        self,
        result: PairwiseLDRecord,
    ) -> bool:
        """Return whether one provider result is terminal."""

        status = _normalize_status(
            result.status
        )

        return (
            status
            in TERMINAL_CACHE_STATUSES
        )

    # ======================================================================
    # Cache retrieval
    # ======================================================================

    def _from_cache(
        self,
        *,
        variant_a: str,
        variant_b: str,
        population: str,
    ) -> PairwiseLDRecord | None:
        """Return caller-oriented record from persistent cache."""

        normalized_a = normalize_rsid(
            variant_a
        )

        normalized_b = normalize_rsid(
            variant_b
        )

        normalized_population = (
            normalize_population(
                population
            )
        )

        cached = self.cache.get(
            provider=self.provider_name,
            population=normalized_population,
            variant_a=normalized_a,
            variant_b=normalized_b,
        )

        if cached is None:

            return None

        return PairwiseLDRecord(
            variant_a=normalized_a,
            variant_b=normalized_b,
            population=normalized_population,
            r2=cached.r2,
            d_prime=cached.d_prime,
            status=cached.status,
            provider=cached.provider,
            reason=cached.reason,
        )

    # ======================================================================
    # Cache persistence
    # ======================================================================

    def _persist_result(
        self,
        result: PairwiseLDRecord,
    ) -> bool:
        """
        Persist one terminal provider result.

        Returns:
            True when result was written to cache.
        """

        if not self._should_cache(
            result
        ):

            self.stats.nonterminal_results += 1

            return False

        existing = self.cache.get(
            provider=self.provider_name,
            population=result.population,
            variant_a=result.variant_a,
            variant_b=result.variant_b,
        )

        self.cache.put(
            provider=self.provider_name,
            population=result.population,
            variant_a=result.variant_a,
            variant_b=result.variant_b,
            r2=result.r2,
            d_prime=result.d_prime,
            status=result.status,
            reason=result.reason,
        )

        if existing is None:

            self.stats.newly_cached_results += 1

        return True

    # ======================================================================
    # Native remote batching
    # ======================================================================

    def _query_remote_pairs(
        self,
        *,
        pairs: list[tuple[str, str]],
        population: str,
    ) -> list[PairwiseLDRecord]:
        """
        Query cache misses through native provider batching when available.

        Falls back to individual query_pair() calls for generic providers.
        """

        if not pairs:

            return []

        self.stats.remote_pairs += int(
            len(
                pairs
            )
        )

        normalized_population = (
            normalize_population(
                population
            )
        )

        # ------------------------------------------------------------------
        # Native batch provider.
        # ------------------------------------------------------------------

        if isinstance(
            self.provider,
            NativeBatchPairwiseLDProvider,
        ):

            self.stats.remote_calls += 1

            records = self.provider.query_pairs(
                pairs=pairs,
                population=normalized_population,
            )

        else:

            # --------------------------------------------------------------
            # Generic correctness-first fallback.
            # --------------------------------------------------------------

            records = []

            for variant_a, variant_b in pairs:

                self.stats.remote_calls += 1

                records.append(
                    self.provider.query_pair(
                        variant_a=variant_a,
                        variant_b=variant_b,
                        population=normalized_population,
                    )
                )

        if len(
            records
        ) != len(
            pairs
        ):

            raise RuntimeError(
                "Remote pairwise provider cardinality mismatch: "
                f"requested={len(pairs)} "
                f"returned={len(records)}"
            )

        return records

    # ======================================================================
    # Single-pair API
    # ======================================================================

    def query_pair(
        self,
        *,
        variant_a: str,
        variant_b: str,
        population: str,
    ) -> PairwiseLDRecord:
        """
        Query one pair using the same batch-aware cache machinery.

        This keeps single-pair and multi-pair behavior consistent.
        """

        records = self.query_pairs(
            pairs=[
                (
                    variant_a,
                    variant_b,
                )
            ],
            population=population,
        )

        if len(
            records
        ) != 1:

            raise RuntimeError(
                "Single cached pairwise query did not return "
                "exactly one result."
            )

        return records[
            0
        ]

    # ======================================================================
    # Batch API
    # ======================================================================

    def query_pairs(
        self,
        *,
        pairs: list[tuple[str, str]],
        population: str,
    ) -> list[PairwiseLDRecord]:
        """
        Query multiple pairwise LD combinations with cache-aware batching.

        Algorithm:
            1. Normalize all caller pairs.
            2. Resolve each unique unordered pair against SQLite.
            3. Return cache hits immediately.
            4. Collect unique cache misses.
            5. Query only those misses using native provider batching.
            6. Persist terminal remote results.
            7. Reconstruct output in the exact caller-request order.
        """

        normalized_population = (
            normalize_population(
                population
            )
        )

        if not pairs:

            return []

        # ------------------------------------------------------------------
        # Normalize caller input while preserving order.
        # ------------------------------------------------------------------

        normalized_requests: list[
            tuple[str, str]
        ] = []

        for variant_a, variant_b in pairs:

            normalized_a = normalize_rsid(
                variant_a
            )

            normalized_b = normalize_rsid(
                variant_b
            )

            if normalized_a == normalized_b:

                raise ValueError(
                    "Pairwise LD query received self-pair: "
                    f"{normalized_a}"
                )

            normalized_requests.append(
                (
                    normalized_a,
                    normalized_b,
                )
            )

        self.stats.requested_pairs += int(
            len(
                normalized_requests
            )
        )

        # ------------------------------------------------------------------
        # One scientific result per canonical unordered pair.
        #
        # Several caller rows may request the same pair. We only perform one
        # cache lookup / remote query for each unique pair.
        # ------------------------------------------------------------------

        result_by_key: dict[
            tuple[str, str],
            PairwiseLDRecord,
        ] = {}

        first_orientation_by_key: dict[
            tuple[str, str],
            tuple[str, str],
        ] = {}

        unique_keys_in_order: list[
            tuple[str, str]
        ] = []

        for variant_a, variant_b in normalized_requests:

            key = _request_key(
                variant_a=variant_a,
                variant_b=variant_b,
            )

            if key not in first_orientation_by_key:

                first_orientation_by_key[
                    key
                ] = (
                    variant_a,
                    variant_b,
                )

                unique_keys_in_order.append(
                    key
                )

        # ------------------------------------------------------------------
        # Cache lookup.
        # ------------------------------------------------------------------

        missing_keys: list[
            tuple[str, str]
        ] = []

        for key in unique_keys_in_order:

            variant_a, variant_b = (
                first_orientation_by_key[
                    key
                ]
            )

            cached = self._from_cache(
                variant_a=variant_a,
                variant_b=variant_b,
                population=normalized_population,
            )

            if cached is not None:

                self.stats.cache_hits += 1

                result_by_key[
                    key
                ] = cached

            else:

                self.stats.cache_misses += 1

                missing_keys.append(
                    key
                )

        # ------------------------------------------------------------------
        # Remote query for misses only.
        # ------------------------------------------------------------------

        if missing_keys:

            remote_pairs = [
                first_orientation_by_key[
                    key
                ]
                for key in missing_keys
            ]

            remote_records = (
                self._query_remote_pairs(
                    pairs=remote_pairs,
                    population=normalized_population,
                )
            )

            for key, (
                requested_pair
            ), remote_record in zip(
                missing_keys,
                remote_pairs,
                remote_records,
                strict=True,
            ):

                requested_a, requested_b = (
                    requested_pair
                )

                # ----------------------------------------------------------
                # Validate provider returned pair identity.
                # ----------------------------------------------------------

                returned_key = _request_key(
                    variant_a=remote_record.variant_a,
                    variant_b=remote_record.variant_b,
                )

                if returned_key != key:

                    mismatch = PairwiseLDRecord(
                        variant_a=requested_a,
                        variant_b=requested_b,
                        population=normalized_population,
                        r2=None,
                        d_prime=None,
                        status="PAIR_MISMATCH",
                        provider=remote_record.provider,
                        reason=(
                            "Underlying provider returned pair "
                            f"{remote_record.variant_a} ↔ "
                            f"{remote_record.variant_b} "
                            "for requested pair "
                            f"{requested_a} ↔ {requested_b}."
                        ),
                    )

                    result_by_key[
                        key
                    ] = mismatch

                    self.stats.nonterminal_results += 1

                    continue

                oriented = _caller_oriented_record(
                    record=remote_record,
                    requested_variant_a=requested_a,
                    requested_variant_b=requested_b,
                    population=normalized_population,
                )

                result_by_key[
                    key
                ] = oriented

                self._persist_result(
                    oriented
                )

        # ------------------------------------------------------------------
        # Reconstruct exact caller request order.
        # ------------------------------------------------------------------

        output: list[
            PairwiseLDRecord
        ] = []

        for requested_a, requested_b in normalized_requests:

            key = _request_key(
                variant_a=requested_a,
                variant_b=requested_b,
            )

            scientific_result = result_by_key.get(
                key
            )

            if scientific_result is None:

                raise RuntimeError(
                    "Internal cached provider error: "
                    "missing reconstructed result for "
                    f"{requested_a} ↔ {requested_b}."
                )

            output.append(
                _caller_oriented_record(
                    record=scientific_result,
                    requested_variant_a=requested_a,
                    requested_variant_b=requested_b,
                    population=normalized_population,
                )
            )

        if len(
            output
        ) != len(
            normalized_requests
        ):

            raise RuntimeError(
                "Cached provider output cardinality mismatch: "
                f"requested={len(normalized_requests)} "
                f"returned={len(output)}"
            )

        return output

    # ======================================================================
    # Diagnostics
    # ======================================================================

    def diagnostics(
        self,
    ) -> dict[str, int]:
        """Return cumulative cache/provider diagnostics."""

        return self.stats.as_dict()

    def reset_diagnostics(
        self,
    ) -> None:
        """Reset cumulative diagnostics."""

        self.stats.reset()