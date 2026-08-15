"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/ldlink_pairwise_provider.py

Description:
    NIH LDlink LDpair provider for M5.3C.3G disease-proxy ↔
    regulatory-proxy pairwise LD analysis.

    This provider queries pairwise linkage disequilibrium between canonical
    rsIDs using the NIH LDlink LDpair REST API.

    Confirmed live LDlink LDpair JSON structure:
        [
            {
                "corr_alleles": [...],
                "haplotypes": {...},
                "pair": ["rs...", "rs..."],
                "request": "...",
                "snp1": {...},
                "snp2": {...},
                "statistics": {
                    "chisq": "...",
                    "d_prime": "...",
                    "p": "...",
                    "r2": "..."
                },
                "two_by_two": {...}
            }
        ]

    Reference configuration:
        Provider:
            NIH LDlink LDpair

        Reference panel:
            1000 Genomes Project

        Genome build:
            GRCh37

        Populations:
            EAS
            EUR
            SAS

    Request policy:
        - API token is loaded from LDLINK_TOKEN.
        - Requests are executed sequentially.
        - Batch POST requests contain at most 10 SNP pairs.
        - Transient HTTP/network failures use bounded retry.
        - Biological unavailability is kept distinct from low LD.
        - Missing pairwise LD is never converted to r² = 0.

    Result states:
        SUCCESS
            Numeric r² was returned.

        VARIANT_UNAVAILABLE
            One or more requested variants were unavailable in the
            reference panel.

        MONOALLELIC
            Variant could not be evaluated because it was monoallelic.

        PAIR_UNAVAILABLE
            LDpair did not provide usable pairwise statistics.

        PAIR_MISMATCH
            Returned variant identity did not match the requested pair.

        PARSE_ERROR
            Response structure could not be interpreted safely.

    Scientific safeguards:
        - Physical proximity is not interpreted as LD.
        - Pairwise LD does not establish causality.
        - Pairwise LD does not establish colocalization.
        - Cross-population consistency is sensitivity evidence rather than
          independent biological replication.
        - No fine-mapping is performed.
        - No causal inference is performed.

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
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Any, Iterable

from pcatrfqtl.analysis.m5.pairwise_ld_provider import (
    PairwiseLDRecord,
)


# ============================================================================
# Provider configuration
# ============================================================================


LDLINK_LDPAIR_URL = (
    "https://ldlink.nih.gov/LDlinkRest/ldpair"
)

PROVIDER_NAME = (
    "LDlink_LDpair"
)

REFERENCE_PANEL = (
    "1000_Genomes_Project"
)

GENOME_BUILD = (
    "grch37"
)

SUPPORTED_POPULATIONS = {
    "EAS",
    "EUR",
    "SAS",
}

MAX_BATCH_SIZE = 10

DEFAULT_TIMEOUT_SECONDS = 120.0

DEFAULT_MAX_ATTEMPTS = 3

RETRYABLE_HTTP_CODES = {
    429,
    500,
    502,
    503,
    504,
}


# ============================================================================
# Models
# ============================================================================


@dataclass(frozen=True)
class LDPairBatchResult:
    """Result metadata for one LDlink LDpair batch request."""

    records: list[PairwiseLDRecord]

    request_pairs: int

    returned_records: int

    population: str

    attempts: int


# ============================================================================
# Scalar helpers
# ============================================================================


def _normalize_rsid(
    value: str,
) -> str:
    """Normalize and validate a canonical rsID."""

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

    numeric_part = normalized[
        2:
    ]

    if not numeric_part.isdigit():

        raise ValueError(
            f"Invalid rsID: {value}"
        )

    return normalized


def _normalize_population(
    population: str,
) -> str:
    """Normalize and validate a supported population label."""

    normalized = (
        str(
            population
        )
        .strip()
        .upper()
    )

    if normalized not in SUPPORTED_POPULATIONS:

        raise ValueError(
            "Unsupported LDlink population: "
            f"{population}. "
            f"Supported={sorted(SUPPORTED_POPULATIONS)}"
        )

    return normalized


def _safe_float(
    value: Any,
) -> float | None:
    """Convert provider scalar to finite float when possible."""

    if value is None:
        return None

    text = (
        str(
            value
        )
        .strip()
    )

    if not text:

        return None

    if text.upper() in {
        "NA",
        "N/A",
        "NONE",
        "NULL",
        ".",
        "NAN",
    }:

        return None

    try:

        result = float(
            text
        )

    except (
        TypeError,
        ValueError,
    ):

        return None

    return result


def _chunks(
    values: list[
        tuple[str, str]
    ],
    size: int,
) -> Iterable[
    list[
        tuple[str, str]
    ]
]:
    """Yield fixed-size pair batches."""

    if size <= 0:

        raise ValueError(
            "Chunk size must be > 0."
        )

    for index in range(
        0,
        len(values),
        size,
    ):

        yield values[
            index:
            index + size
        ]


def _redact_token(
    text: str,
    token: str,
) -> str:
    """Remove API token from error text."""

    if not text:
        return text

    if not token:
        return text

    return text.replace(
        token,
        "[REDACTED]",
    )


# ============================================================================
# Response helpers
# ============================================================================


def _extract_provider_message(
    item: dict[str, Any],
) -> str | None:
    """Extract warning/error information from an LDlink item."""

    for key in (
        "error",
        "Error",
        "warning",
        "Warning",
        "message",
        "Message",
    ):

        value = item.get(
            key
        )

        if value is not None:

            text = (
                str(
                    value
                )
                .strip()
            )

            if text:

                return text

    corr_alleles = item.get(
        "corr_alleles"
    )

    if isinstance(
        corr_alleles,
        list,
    ):

        values = [
            str(
                value
            ).strip()
            for value in corr_alleles
            if str(
                value
            ).strip()
        ]

        if values:

            return " | ".join(
                values
            )

    return None


def _classify_unavailable_message(
    message: str | None,
) -> str:
    """Classify a non-numeric LDlink response."""

    if message is None:

        return "PAIR_UNAVAILABLE"

    normalized = (
        message
        .strip()
        .lower()
    )

    if "monoallelic" in normalized:

        return "MONOALLELIC"

    unavailable_patterns = (
        "not found",
        "not in 1000 genomes",
        "not in the 1000 genomes",
        "not found in",
        "does not match",
        "not available",
        "variant unavailable",
    )

    if any(
        pattern in normalized
        for pattern in unavailable_patterns
    ):

        return "VARIANT_UNAVAILABLE"

    return "PAIR_UNAVAILABLE"


# ============================================================================
# LDlink provider
# ============================================================================


class LDlinkPairwiseProvider:
    """NIH LDlink LDpair implementation of the pairwise LD provider."""

    def __init__(
        self,
        *,
        token: str | None = None,
        timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
        max_attempts: int = DEFAULT_MAX_ATTEMPTS,
        inter_batch_delay_seconds: float = 1.0,
    ) -> None:

        resolved_token = (
            token
            or os.environ.get(
                "LDLINK_TOKEN"
            )
        )

        if resolved_token is None:

            raise RuntimeError(
                "LDLINK_TOKEN is not configured."
            )

        resolved_token = (
            resolved_token
            .strip()
        )

        if not resolved_token:

            raise RuntimeError(
                "LDLINK_TOKEN is empty."
            )

        self.token = resolved_token

        self.timeout_seconds = float(
            timeout_seconds
        )

        self.max_attempts = int(
            max_attempts
        )

        self.inter_batch_delay_seconds = float(
            inter_batch_delay_seconds
        )

        if self.timeout_seconds <= 0:

            raise ValueError(
                "timeout_seconds must be > 0."
            )

        if self.max_attempts < 1:

            raise ValueError(
                "max_attempts must be >= 1."
            )

        if self.inter_batch_delay_seconds < 0:

            raise ValueError(
                "inter_batch_delay_seconds must be >= 0."
            )

    # ======================================================================
    # Request construction
    # ======================================================================

    def _build_request(
        self,
        *,
        pairs: list[
            tuple[str, str]
        ],
        population: str,
    ) -> urllib.request.Request:
        """Construct one LDlink LDpair POST request."""

        if not pairs:

            raise ValueError(
                "Cannot build LDpair request with zero variant pairs."
            )

        if len(
            pairs
        ) > MAX_BATCH_SIZE:

            raise ValueError(
                "LDpair batch exceeds maximum size: "
                f"{len(pairs)} > {MAX_BATCH_SIZE}"
            )

        population = _normalize_population(
            population
        )

        normalized_pairs = [
            (
                _normalize_rsid(
                    variant_a
                ),
                _normalize_rsid(
                    variant_b
                ),
            )
            for (
                variant_a,
                variant_b
            )
            in pairs
        ]

        url = (
            LDLINK_LDPAIR_URL
            + "?"
            + urllib.parse.urlencode(
                {
                    "token":
                        self.token,
                }
            )
        )

        payload = {
            "snp_pairs": [
                [
                    variant_a,
                    variant_b,
                ]
                for (
                    variant_a,
                    variant_b
                )
                in normalized_pairs
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

    # ======================================================================
    # HTTP
    # ======================================================================

    def _post_batch(
        self,
        *,
        pairs: list[
            tuple[str, str]
        ],
        population: str,
    ) -> tuple[
        Any,
        int,
    ]:
        """
        Execute one LDpair batch with bounded retries.

        Returns:
            Parsed JSON response.
            Number of attempts used.
        """

        last_error: Exception | None = None

        for attempt in range(
            1,
            self.max_attempts + 1,
        ):

            request = self._build_request(
                pairs=pairs,
                population=population,
            )

            try:

                with urllib.request.urlopen(
                    request,
                    timeout=self.timeout_seconds,
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
                        "Unexpected LDlink HTTP status: "
                        f"{status_code}"
                    )

                try:

                    parsed = json.loads(
                        raw_text
                    )

                except json.JSONDecodeError as exc:

                    safe_raw = _redact_token(
                        raw_text,
                        self.token,
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
                    _redact_token(
                        response_body,
                        self.token,
                    )
                )

                if (
                    exc.code
                    not in RETRYABLE_HTTP_CODES
                ):

                    raise RuntimeError(
                        "LDlink returned non-retryable "
                        f"HTTP {exc.code}: "
                        f"{response_body}"
                    ) from exc

            except (
                urllib.error.URLError,
                TimeoutError,
                ConnectionError,
            ) as exc:

                last_error = exc

            if attempt < self.max_attempts:

                delay = (
                    2 ** (
                        attempt - 1
                    )
                )

                time.sleep(
                    delay
                )

        safe_error = _redact_token(
            str(
                last_error
            ),
            self.token,
        )

        raise RuntimeError(
            "LDlink LDpair request failed "
            f"after {self.max_attempts} attempts: "
            f"{safe_error}"
        )

    # ======================================================================
    # Item parser
    # ======================================================================

    def _parse_ldpair_item(
        self,
        *,
        item: dict[str, Any],
        requested_variant_a: str,
        requested_variant_b: str,
        population: str,
    ) -> PairwiseLDRecord:
        """Parse one confirmed LDlink LDpair response object."""

        requested_variant_a = _normalize_rsid(
            requested_variant_a
        )

        requested_variant_b = _normalize_rsid(
            requested_variant_b
        )

        population = _normalize_population(
            population
        )

        # ------------------------------------------------------------------
        # Pair identity validation
        # ------------------------------------------------------------------

        pair = item.get(
            "pair"
        )

        if not isinstance(
            pair,
            list,
        ):

            return PairwiseLDRecord(
                variant_a=requested_variant_a,
                variant_b=requested_variant_b,
                population=population,
                r2=None,
                d_prime=None,
                status="PARSE_ERROR",
                provider=PROVIDER_NAME,
                reason=(
                    "Missing or malformed LDlink 'pair' field."
                ),
            )

        if len(
            pair
        ) != 2:

            return PairwiseLDRecord(
                variant_a=requested_variant_a,
                variant_b=requested_variant_b,
                population=population,
                r2=None,
                d_prime=None,
                status="PARSE_ERROR",
                provider=PROVIDER_NAME,
                reason=(
                    "LDlink 'pair' field did not contain "
                    "exactly two variant identifiers."
                ),
            )

        try:

            returned_pair = {
                _normalize_rsid(
                    pair[
                        0
                    ]
                ),
                _normalize_rsid(
                    pair[
                        1
                    ]
                ),
            }

        except ValueError:

            return PairwiseLDRecord(
                variant_a=requested_variant_a,
                variant_b=requested_variant_b,
                population=population,
                r2=None,
                d_prime=None,
                status="PARSE_ERROR",
                provider=PROVIDER_NAME,
                reason=(
                    "LDlink returned invalid rsIDs in the "
                    "'pair' field."
                ),
            )

        requested_pair = {
            requested_variant_a,
            requested_variant_b,
        }

        if returned_pair != requested_pair:

            return PairwiseLDRecord(
                variant_a=requested_variant_a,
                variant_b=requested_variant_b,
                population=population,
                r2=None,
                d_prime=None,
                status="PAIR_MISMATCH",
                provider=PROVIDER_NAME,
                reason=(
                    "Requested pair "
                    f"{sorted(requested_pair)} "
                    "but LDlink returned "
                    f"{sorted(returned_pair)}."
                ),
            )

        # ------------------------------------------------------------------
        # Statistics
        # ------------------------------------------------------------------

        statistics = item.get(
            "statistics"
        )

        provider_message = (
            _extract_provider_message(
                item
            )
        )

        if not isinstance(
            statistics,
            dict,
        ):

            status = (
                _classify_unavailable_message(
                    provider_message
                )
            )

            return PairwiseLDRecord(
                variant_a=requested_variant_a,
                variant_b=requested_variant_b,
                population=population,
                r2=None,
                d_prime=None,
                status=status,
                provider=PROVIDER_NAME,
                reason=(
                    provider_message
                    or
                    "LDlink response contained no "
                    "'statistics' object."
                ),
            )

        r2 = _safe_float(
            statistics.get(
                "r2"
            )
        )

        d_prime = _safe_float(
            statistics.get(
                "d_prime"
            )
        )

        # ------------------------------------------------------------------
        # Successful LD result
        #
        # r² is mandatory for a successful pairwise LD result.
        # D' may still be absent without invalidating r².
        # ------------------------------------------------------------------

        if r2 is not None:

            return PairwiseLDRecord(
                variant_a=requested_variant_a,
                variant_b=requested_variant_b,
                population=population,
                r2=r2,
                d_prime=d_prime,
                status="SUCCESS",
                provider=PROVIDER_NAME,
                reason=None,
            )

        # ------------------------------------------------------------------
        # Statistics object present, but no numeric r²
        # ------------------------------------------------------------------

        status = (
            _classify_unavailable_message(
                provider_message
            )
        )

        return PairwiseLDRecord(
            variant_a=requested_variant_a,
            variant_b=requested_variant_b,
            population=population,
            r2=None,
            d_prime=d_prime,
            status=status,
            provider=PROVIDER_NAME,
            reason=(
                provider_message
                or
                "LDlink returned no numeric r2."
            ),
        )

    # ======================================================================
    # Batch parser
    # ======================================================================

    def _parse_batch_response(
        self,
        *,
        response: Any,
        requested_pairs: list[
            tuple[str, str]
        ],
        population: str,
    ) -> list[
        PairwiseLDRecord
    ]:
        """Parse a confirmed LDlink LDpair top-level list response."""

        population = _normalize_population(
            population
        )

        normalized_pairs = [
            (
                _normalize_rsid(
                    variant_a
                ),
                _normalize_rsid(
                    variant_b
                ),
            )
            for (
                variant_a,
                variant_b
            )
            in requested_pairs
        ]

        if not isinstance(
            response,
            list,
        ):

            return [
                PairwiseLDRecord(
                    variant_a=variant_a,
                    variant_b=variant_b,
                    population=population,
                    r2=None,
                    d_prime=None,
                    status="PARSE_ERROR",
                    provider=PROVIDER_NAME,
                    reason=(
                        "Expected top-level LDlink response "
                        "to be a JSON list, received "
                        f"{type(response).__name__}."
                    ),
                )
                for (
                    variant_a,
                    variant_b
                )
                in normalized_pairs
            ]

        records: list[
            PairwiseLDRecord
        ] = []

        for index, (
            variant_a,
            variant_b,
        ) in enumerate(
            normalized_pairs
        ):

            if index >= len(
                response
            ):

                records.append(
                    PairwiseLDRecord(
                        variant_a=variant_a,
                        variant_b=variant_b,
                        population=population,
                        r2=None,
                        d_prime=None,
                        status="PAIR_UNAVAILABLE",
                        provider=PROVIDER_NAME,
                        reason=(
                            "LDlink returned fewer response "
                            "items than requested pairs."
                        ),
                    )
                )

                continue

            item = response[
                index
            ]

            if not isinstance(
                item,
                dict,
            ):

                records.append(
                    PairwiseLDRecord(
                        variant_a=variant_a,
                        variant_b=variant_b,
                        population=population,
                        r2=None,
                        d_prime=None,
                        status="PARSE_ERROR",
                        provider=PROVIDER_NAME,
                        reason=(
                            "LDlink response item was not "
                            "a JSON object."
                        ),
                    )
                )

                continue

            record = (
                self._parse_ldpair_item(
                    item=item,
                    requested_variant_a=variant_a,
                    requested_variant_b=variant_b,
                    population=population,
                )
            )

            records.append(
                record
            )

        return records

    # ======================================================================
    # Public batch API
    # ======================================================================

    def query_pairs(
        self,
        *,
        pairs: list[
            tuple[str, str]
        ],
        population: str,
    ) -> list[
        PairwiseLDRecord
    ]:
        """
        Query one or more pairwise LD combinations.

        Requests are split into sequential batches of at most 10 pairs.
        """

        population = _normalize_population(
            population
        )

        if not pairs:

            return []

        normalized_pairs = [
            (
                _normalize_rsid(
                    variant_a
                ),
                _normalize_rsid(
                    variant_b
                ),
            )
            for (
                variant_a,
                variant_b
            )
            in pairs
        ]

        all_records: list[
            PairwiseLDRecord
        ] = []

        batches = list(
            _chunks(
                normalized_pairs,
                MAX_BATCH_SIZE,
            )
        )

        for batch_index, batch in enumerate(
            batches
        ):

            if (
                batch_index > 0
                and
                self.inter_batch_delay_seconds
                > 0
            ):

                time.sleep(
                    self.inter_batch_delay_seconds
                )

            response, _ = (
                self._post_batch(
                    pairs=batch,
                    population=population,
                )
            )

            records = (
                self._parse_batch_response(
                    response=response,
                    requested_pairs=batch,
                    population=population,
                )
            )

            all_records.extend(
                records
            )

        if len(
            all_records
        ) != len(
            normalized_pairs
        ):

            raise RuntimeError(
                "Pairwise LD result cardinality mismatch: "
                f"requested={len(normalized_pairs)} "
                f"returned={len(all_records)}"
            )

        return all_records

    # ======================================================================
    # Public single-pair API
    # ======================================================================

    def query_pair(
        self,
        *,
        variant_a: str,
        variant_b: str,
        population: str,
    ) -> PairwiseLDRecord:
        """Query one pair, satisfying PairwiseLDProvider protocol."""

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
                "Single-pair LDlink query did not return exactly "
                "one PairwiseLDRecord."
            )

        return records[
            0
        ]