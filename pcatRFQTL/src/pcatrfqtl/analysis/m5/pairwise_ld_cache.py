"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/pairwise_ld_cache.py

Description:
    Persistent SQLite cache for M5.3C.3G pairwise LD queries.

    The cache stores one result per:

        provider
        population
        canonical unordered rsID pair

    Pair identity is canonicalized so:

        rsA ↔ rsB

    and:

        rsB ↔ rsA

    resolve to the same cache entry.

    Terminal biological/provider states may be persisted, including:

        SUCCESS
        VARIANT_UNAVAILABLE
        MONOALLELIC
        PAIR_UNAVAILABLE

    Transient or parser-related states should normally not be persisted:

        RATE_LIMITED
        REMOTE_ERROR
        TIMEOUT
        PARSE_ERROR
        PAIR_MISMATCH

    Schema migration:
        Schema version 2 introduces the updated_at column.

        Existing schema-version-1 caches are migrated automatically.
        Existing cached scientific results are preserved.

    Scientific safeguards:
        - Missing LD is never represented as r² = 0.
        - Variant unavailability is distinct from low LD.
        - Cache retrieval does not modify the scientific result.
        - Pair orientation is canonicalized only for cache identity.
        - Returned result orientation is handled by the provider wrapper.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path


# ============================================================================
# Constants
# ============================================================================


CACHE_SCHEMA_VERSION = 2


# ============================================================================
# Models
# ============================================================================


@dataclass(frozen=True)
class CachedPairwiseLD:
    """One persisted pairwise LD result."""

    provider: str

    population: str

    variant_a: str

    variant_b: str

    r2: float | None

    d_prime: float | None

    status: str

    reason: str | None

    created_at: str

    updated_at: str


# ============================================================================
# Normalization helpers
# ============================================================================


def normalize_rsid(
    value: str,
) -> str:
    """Normalize and validate one canonical rsID."""

    normalized = (
        str(value)
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


def normalize_population(
    value: str,
) -> str:
    """Normalize population label."""

    normalized = (
        str(value)
        .strip()
        .upper()
    )

    if not normalized:

        raise ValueError(
            "Population cannot be empty."
        )

    return normalized


def normalize_provider(
    value: str,
) -> str:
    """Normalize provider name for cache identity."""

    normalized = (
        str(value)
        .strip()
    )

    if not normalized:

        raise ValueError(
            "Provider name cannot be empty."
        )

    return normalized


def normalize_status(
    value: str,
) -> str:
    """Normalize provider/cache status."""

    normalized = (
        str(value)
        .strip()
        .upper()
    )

    if not normalized:

        raise ValueError(
            "Cache status cannot be empty."
        )

    return normalized


def canonical_pair(
    variant_a: str,
    variant_b: str,
) -> tuple[str, str]:
    """
    Return deterministic unordered pair identity.

    Pairwise LD is symmetric for cache identity. This does not modify the
    requested orientation returned to the calling analysis layer.
    """

    a = normalize_rsid(
        variant_a
    )

    b = normalize_rsid(
        variant_b
    )

    if a == b:

        raise ValueError(
            "Pairwise LD cache received a self-pair: "
            f"{a}"
        )

    first, second = sorted(
        (
            a,
            b,
        )
    )

    return (
        first,
        second,
    )


# ============================================================================
# Cache
# ============================================================================


class PairwiseLDCache:
    """SQLite-backed persistent cache for pairwise LD results."""

    def __init__(
        self,
        path: str | Path,
    ) -> None:

        self.path = Path(
            path
        )

        self.path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.connection = sqlite3.connect(
            self.path
        )

        self.connection.row_factory = (
            sqlite3.Row
        )

        self._initialize()

    # ======================================================================
    # Schema initialization / migration
    # ======================================================================

    def _initialize(
        self,
    ) -> None:
        """
        Create or migrate the cache schema.

        Schema history:

            version 1
                Initial pairwise_ld table.

            version 2
                Adds updated_at while preserving existing cached records.
        """

        self.connection.execute(
            """
            CREATE TABLE IF NOT EXISTS metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            )
            """
        )

        # ------------------------------------------------------------------
        # Base table.
        #
        # updated_at is deliberately not relied upon here for migration,
        # because CREATE TABLE IF NOT EXISTS does not alter an existing
        # version-1 table.
        # ------------------------------------------------------------------

        self.connection.execute(
            """
            CREATE TABLE IF NOT EXISTS pairwise_ld (
                provider TEXT NOT NULL,
                population TEXT NOT NULL,
                variant_a TEXT NOT NULL,
                variant_b TEXT NOT NULL,

                r2 REAL,
                d_prime REAL,

                status TEXT NOT NULL,
                reason TEXT,

                created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,

                PRIMARY KEY (
                    provider,
                    population,
                    variant_a,
                    variant_b
                )
            )
            """
        )

        # ------------------------------------------------------------------
        # Inspect current table schema.
        # ------------------------------------------------------------------

        cursor = self.connection.execute(
            """
            PRAGMA table_info(pairwise_ld)
            """
        )

        existing_columns = {
            str(
                row[
                    "name"
                ]
            )
            for row in cursor.fetchall()
        }

        # ------------------------------------------------------------------
        # Migration to schema version 2.
        #
        # SQLite ALTER TABLE ADD COLUMN cannot safely add a dynamic
        # CURRENT_TIMESTAMP default on older SQLite builds, so the column is
        # introduced nullable and populated explicitly.
        # ------------------------------------------------------------------

        if "updated_at" not in existing_columns:

            self.connection.execute(
                """
                ALTER TABLE pairwise_ld
                ADD COLUMN updated_at TEXT
                """
            )

            self.connection.execute(
                """
                UPDATE pairwise_ld
                SET updated_at = created_at
                WHERE updated_at IS NULL
                """
            )

        # ------------------------------------------------------------------
        # Defensive migration for databases where created_at might somehow
        # predate the current schema.
        # ------------------------------------------------------------------

        cursor = self.connection.execute(
            """
            PRAGMA table_info(pairwise_ld)
            """
        )

        migrated_columns = {
            str(
                row[
                    "name"
                ]
            )
            for row in cursor.fetchall()
        }

        required_columns = {
            "provider",
            "population",
            "variant_a",
            "variant_b",
            "r2",
            "d_prime",
            "status",
            "reason",
            "created_at",
            "updated_at",
        }

        missing_columns = (
            required_columns
            - migrated_columns
        )

        if missing_columns:

            raise RuntimeError(
                "Pairwise LD cache schema is missing required columns: "
                f"{sorted(missing_columns)}"
            )

        # ------------------------------------------------------------------
        # Indexes.
        # ------------------------------------------------------------------

        self.connection.execute(
            """
            CREATE INDEX IF NOT EXISTS
                idx_pairwise_ld_population
            ON pairwise_ld (
                population
            )
            """
        )

        self.connection.execute(
            """
            CREATE INDEX IF NOT EXISTS
                idx_pairwise_ld_status
            ON pairwise_ld (
                status
            )
            """
        )

        self.connection.execute(
            """
            CREATE INDEX IF NOT EXISTS
                idx_pairwise_ld_provider
            ON pairwise_ld (
                provider
            )
            """
        )

        # ------------------------------------------------------------------
        # Record current schema version.
        #
        # DO UPDATE is intentional. Version-1 caches must be upgraded to 2.
        # ------------------------------------------------------------------

        self.connection.execute(
            """
            INSERT INTO metadata (
                key,
                value
            )
            VALUES (
                'schema_version',
                ?
            )

            ON CONFLICT(key)
            DO UPDATE SET
                value = excluded.value
            """,
            (
                str(
                    CACHE_SCHEMA_VERSION
                ),
            ),
        )

        self.connection.commit()

        # ------------------------------------------------------------------
        # Final schema-version verification.
        # ------------------------------------------------------------------

        schema_version = self.get_metadata(
            "schema_version"
        )

        if schema_version != str(
            CACHE_SCHEMA_VERSION
        ):

            raise RuntimeError(
                "Unsupported pairwise LD cache schema version: "
                f"{schema_version}. "
                f"Expected={CACHE_SCHEMA_VERSION}"
            )

    # ======================================================================
    # Metadata
    # ======================================================================

    def get_metadata(
        self,
        key: str,
    ) -> str | None:
        """Read one metadata value."""

        cursor = self.connection.execute(
            """
            SELECT value
            FROM metadata
            WHERE key = ?
            """,
            (
                key,
            ),
        )

        row = cursor.fetchone()

        if row is None:

            return None

        return str(
            row[
                "value"
            ]
        )

    # ======================================================================
    # Retrieval
    # ======================================================================

    def get(
        self,
        *,
        provider: str,
        population: str,
        variant_a: str,
        variant_b: str,
    ) -> CachedPairwiseLD | None:
        """Retrieve one cached pairwise LD result."""

        provider = normalize_provider(
            provider
        )

        population = normalize_population(
            population
        )

        canonical_a, canonical_b = (
            canonical_pair(
                variant_a,
                variant_b,
            )
        )

        cursor = self.connection.execute(
            """
            SELECT
                provider,
                population,
                variant_a,
                variant_b,
                r2,
                d_prime,
                status,
                reason,
                created_at,
                updated_at
            FROM pairwise_ld
            WHERE
                provider = ?
                AND population = ?
                AND variant_a = ?
                AND variant_b = ?
            """,
            (
                provider,
                population,
                canonical_a,
                canonical_b,
            ),
        )

        row = cursor.fetchone()

        if row is None:

            return None

        created_at = (
            ""
            if row[
                "created_at"
            ] is None
            else str(
                row[
                    "created_at"
                ]
            )
        )

        updated_at = (
            created_at
            if row[
                "updated_at"
            ] is None
            else str(
                row[
                    "updated_at"
                ]
            )
        )

        return CachedPairwiseLD(
            provider=str(
                row[
                    "provider"
                ]
            ),

            population=str(
                row[
                    "population"
                ]
            ),

            variant_a=str(
                row[
                    "variant_a"
                ]
            ),

            variant_b=str(
                row[
                    "variant_b"
                ]
            ),

            r2=(
                None
                if row[
                    "r2"
                ] is None
                else float(
                    row[
                        "r2"
                    ]
                )
            ),

            d_prime=(
                None
                if row[
                    "d_prime"
                ] is None
                else float(
                    row[
                        "d_prime"
                    ]
                )
            ),

            status=str(
                row[
                    "status"
                ]
            ),

            reason=(
                None
                if row[
                    "reason"
                ] is None
                else str(
                    row[
                        "reason"
                    ]
                )
            ),

            created_at=created_at,

            updated_at=updated_at,
        )

    def contains(
        self,
        *,
        provider: str,
        population: str,
        variant_a: str,
        variant_b: str,
    ) -> bool:
        """Return whether one pairwise result exists in the cache."""

        return (
            self.get(
                provider=provider,
                population=population,
                variant_a=variant_a,
                variant_b=variant_b,
            )
            is not None
        )

    # ======================================================================
    # Write
    # ======================================================================

    def put(
        self,
        *,
        provider: str,
        population: str,
        variant_a: str,
        variant_b: str,
        r2: float | None,
        d_prime: float | None,
        status: str,
        reason: str | None,
    ) -> None:
        """Insert or update one pairwise LD result."""

        provider = normalize_provider(
            provider
        )

        population = normalize_population(
            population
        )

        canonical_a, canonical_b = (
            canonical_pair(
                variant_a,
                variant_b,
            )
        )

        status = normalize_status(
            status
        )

        self.connection.execute(
            """
            INSERT INTO pairwise_ld (
                provider,
                population,
                variant_a,
                variant_b,
                r2,
                d_prime,
                status,
                reason,
                updated_at
            )
            VALUES (
                ?, ?, ?, ?, ?, ?, ?, ?,
                CURRENT_TIMESTAMP
            )

            ON CONFLICT (
                provider,
                population,
                variant_a,
                variant_b
            )
            DO UPDATE SET
                r2 = excluded.r2,
                d_prime = excluded.d_prime,
                status = excluded.status,
                reason = excluded.reason,
                updated_at = CURRENT_TIMESTAMP
            """,
            (
                provider,
                population,
                canonical_a,
                canonical_b,
                r2,
                d_prime,
                status,
                reason,
            ),
        )

        self.connection.commit()

    # ======================================================================
    # Delete
    # ======================================================================

    def delete(
        self,
        *,
        provider: str,
        population: str,
        variant_a: str,
        variant_b: str,
    ) -> bool:
        """Delete one cached pair and report whether a row was removed."""

        provider = normalize_provider(
            provider
        )

        population = normalize_population(
            population
        )

        canonical_a, canonical_b = (
            canonical_pair(
                variant_a,
                variant_b,
            )
        )

        cursor = self.connection.execute(
            """
            DELETE FROM pairwise_ld
            WHERE
                provider = ?
                AND population = ?
                AND variant_a = ?
                AND variant_b = ?
            """,
            (
                provider,
                population,
                canonical_a,
                canonical_b,
            ),
        )

        self.connection.commit()

        return bool(
            cursor.rowcount
        )

    # ======================================================================
    # Counts
    # ======================================================================

    def count(
        self,
    ) -> int:
        """Return total number of cached pair-population results."""

        cursor = self.connection.execute(
            """
            SELECT COUNT(*) AS n
            FROM pairwise_ld
            """
        )

        row = cursor.fetchone()

        return int(
            row[
                "n"
            ]
        )

    def count_by_status(
        self,
    ) -> dict[str, int]:
        """Return cache counts grouped by result status."""

        cursor = self.connection.execute(
            """
            SELECT
                status,
                COUNT(*) AS n
            FROM pairwise_ld
            GROUP BY status
            ORDER BY status
            """
        )

        return {
            str(
                row[
                    "status"
                ]
            ):
                int(
                    row[
                        "n"
                    ]
                )
            for row in cursor.fetchall()
        }

    def count_by_population(
        self,
    ) -> dict[str, int]:
        """Return cache counts grouped by reference population."""

        cursor = self.connection.execute(
            """
            SELECT
                population,
                COUNT(*) AS n
            FROM pairwise_ld
            GROUP BY population
            ORDER BY population
            """
        )

        return {
            str(
                row[
                    "population"
                ]
            ):
                int(
                    row[
                        "n"
                    ]
                )
            for row in cursor.fetchall()
        }

    def count_by_provider(
        self,
    ) -> dict[str, int]:
        """Return cache counts grouped by provider."""

        cursor = self.connection.execute(
            """
            SELECT
                provider,
                COUNT(*) AS n
            FROM pairwise_ld
            GROUP BY provider
            ORDER BY provider
            """
        )

        return {
            str(
                row[
                    "provider"
                ]
            ):
                int(
                    row[
                        "n"
                    ]
                )
            for row in cursor.fetchall()
        }

    # ======================================================================
    # Inspection
    # ======================================================================

    def schema_version(
        self,
    ) -> int:
        """Return current cache schema version."""

        value = self.get_metadata(
            "schema_version"
        )

        if value is None:

            raise RuntimeError(
                "Cache schema_version metadata is missing."
            )

        return int(
            value
        )

    def table_columns(
        self,
    ) -> list[str]:
        """Return current pairwise_ld table columns."""

        cursor = self.connection.execute(
            """
            PRAGMA table_info(pairwise_ld)
            """
        )

        return [
            str(
                row[
                    "name"
                ]
            )
            for row in cursor.fetchall()
        ]

    # ======================================================================
    # Lifecycle
    # ======================================================================

    def close(
        self,
    ) -> None:
        """Close the SQLite connection."""

        self.connection.close()

    def __enter__(
        self,
    ) -> "PairwiseLDCache":
        """Enter context manager."""

        return self

    def __exit__(
        self,
        exc_type,
        exc_value,
        traceback,
    ) -> None:
        """Exit context manager."""

        self.close()