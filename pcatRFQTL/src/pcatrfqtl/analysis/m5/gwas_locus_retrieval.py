"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_locus_retrieval.py

Description:
    M5.3C.3C utilities for indexed locus-level retrieval from harmonised
    GWAS Catalog summary statistics using pysam/HTSlib.

    Remote harmonised GWAS files are accessed through their BGZF/Tabix
    indexes. One remote Tabix session can be reused to retrieve multiple
    candidate loci from the same GWAS study.

    Header handling:
        GWAS Catalog harmonised summary-statistics files store the TSV
        column names in the first decompressed line of the BGZF file.

        In practice, pysam.TabixFile.header may not expose this line because
        the column header is not necessarily represented as Tabix metadata.

        Therefore:
            1. TabixFile is opened for indexed locus retrieval.
            2. Only the first decompressed line of the remote BGZF file is
               independently read to resolve the TSV columns.
            3. The connection used for header retrieval is immediately closed.
            4. Genomic locus retrieval itself remains indexed through pysam.

    Coordinate conventions:
        Biological locus windows:
            1-based inclusive.

        pysam.TabixFile.fetch():
            0-based half-open.

        Conversion:
            fetch_start = window_start_1based - 1
            fetch_end   = window_end_1based

    Important safeguards:
        - No external tabix executable is required.
        - No genome-wide streaming fallback is performed.
        - Whole GWAS files are not intentionally downloaded.
        - Only the first decompressed line is read for header resolution.
        - Remote locus access remains Tabix-indexed.
        - Actual index contigs are used to resolve chromosome notation.
        - Empty loci are retained as valid retrieval outcomes.
        - Lead absence does not imply absence of regional disease signal.
        - Physical proximity is not interpreted as linkage disequilibrium.
        - No colocalization is performed.
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

import gzip
import re
import urllib.request
from dataclasses import dataclass
from typing import Any, Iterable

import pandas as pd
import pysam


# ============================================================================
# Exceptions
# ============================================================================


class RemoteTabixError(
    RuntimeError
):
    """Base exception for remote indexed GWAS retrieval failures."""


class RemoteTabixOpenError(
    RemoteTabixError
):
    """Raised when a remote Tabix-indexed resource cannot be opened."""


class ContigResolutionError(
    RemoteTabixError
):
    """Raised when a requested chromosome cannot be mapped to an index."""


class HeaderResolutionError(
    RemoteTabixError
):
    """Raised when the harmonised GWAS column header cannot be resolved."""


# ============================================================================
# Data models
# ============================================================================


@dataclass(frozen=True)
class RemoteTabixProbe:
    """Result of testing a remote Tabix-indexed resource."""

    url: str

    accessible: bool

    contigs: tuple[str, ...]

    reason: str | None


@dataclass(frozen=True)
class LocusQueryMetadata:
    """Coordinate provenance for one indexed genomic locus query."""

    requested_chromosome: str

    resolved_contig: str

    window_start_1based: int

    window_end_1based: int

    fetch_start_0based: int

    fetch_end_0based: int


# ============================================================================
# Expected GWAS column aliases
# ============================================================================


VARIANT_ID_ALIASES = (
    "rsid",
    "variant_id",
    "hm_variant_id",
    "rs_id",
    "snp",
)


CHROMOSOME_ALIASES = (
    "chromosome",
    "hm_chromosome",
    "chrom",
    "chr",
)


POSITION_ALIASES = (
    "base_pair_location",
    "hm_base_pair_location",
    "position",
    "pos",
    "bp",
)


KNOWN_GWAS_HEADER_FIELDS = {
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "beta",
    "standard_error",
    "effect_allele_frequency",
    "p_value",
    "rsid",
    "reference_allele",
    "hm_coordinate_conversion",
    "hm_code",
    "variant_id",
}


# ============================================================================
# Generic normalization helpers
# ============================================================================


def _normalize_header(
    value: Any,
) -> str:
    """Normalize one GWAS summary-statistics column name."""

    text = (
        str(
            value
        )
        .replace(
            "\ufeff",
            "",
        )
        .strip()
        .lstrip(
            "#"
        )
    )

    text = re.sub(
        r"[\s\-]+",
        "_",
        text,
    )

    return text.lower()


def _resolve_column(
    columns: Iterable[str],
    aliases: tuple[str, ...],
) -> str | None:
    """Resolve one DataFrame column using normalized aliases."""

    lookup = {
        _normalize_header(
            column
        ):
            column
        for column
        in columns
    }

    for alias in aliases:

        normalized_alias = (
            _normalize_header(
                alias
            )
        )

        if normalized_alias in lookup:

            return lookup[
                normalized_alias
            ]

    return None


def _normalize_chromosome(
    chromosome: Any,
) -> str:
    """Normalize chromosome notation independently of contig naming style."""

    value = (
        str(
            chromosome
        )
        .strip()
    )

    value = re.sub(
        r"^chr",
        "",
        value,
        flags=re.IGNORECASE,
    )

    value = value.upper()

    if value in {
        "M",
        "MT",
    }:

        return "MT"

    return value


# ============================================================================
# Contig handling
# ============================================================================


def _resolve_contig(
    chromosome: str,
    contigs: Iterable[str],
) -> str:
    """
    Resolve one chromosome against actual Tabix contigs.

    Handles common equivalents such as:
        8    <-> chr8
        X    <-> chrX
        MT   <-> chrM / chrMT
    """

    available = tuple(
        str(
            contig
        )
        for contig
        in contigs
    )

    if not available:

        raise ContigResolutionError(
            "Remote Tabix index contains no contigs."
        )

    requested = (
        _normalize_chromosome(
            chromosome
        )
    )

    # ------------------------------------------------------------------
    # Normalized direct match
    # ------------------------------------------------------------------

    for contig in available:

        if (
            _normalize_chromosome(
                contig
            )
            == requested
        ):

            return contig

    # ------------------------------------------------------------------
    # Explicit mitochondrial fallback
    # ------------------------------------------------------------------

    if requested == "MT":

        for candidate in (
            "MT",
            "M",
            "chrMT",
            "chrM",
        ):

            if candidate in available:

                return candidate

    raise ContigResolutionError(
        "Could not resolve requested chromosome "
        f"{chromosome!r} against remote Tabix contigs."
    )


# ============================================================================
# Remote column-header retrieval
# ============================================================================


def _read_remote_column_header(
    url: str,
) -> tuple[str, ...]:
    """
    Read only the first decompressed line of one remote harmonised GWAS file.

    GWAS Catalog harmonised summary statistics store the TSV column names in
    the first decompressed line.

    This function intentionally reads only that first line and immediately
    closes the HTTP stream. It is not a genome-wide streaming fallback.
    """

    request = urllib.request.Request(
        url,
        headers={
            "Accept":
                "*/*",

            "User-Agent":
                (
                    "pcatRFQTL/1.0 "
                    "M5.3C.3-GWAS-header-resolution"
                ),
        },
    )

    try:

        with urllib.request.urlopen(
            request,
            timeout=60,
        ) as response:

            with gzip.GzipFile(
                fileobj=response,
            ) as gzip_file:

                raw_line = (
                    gzip_file.readline()
                )

    except Exception as exc:

        raise HeaderResolutionError(
            "Could not retrieve harmonised GWAS column header from "
            f"{url}: {type(exc).__name__}: {exc}"
        ) from exc

    if not raw_line:

        raise HeaderResolutionError(
            "Remote harmonised GWAS file returned an empty first line."
        )

    try:

        line = (
            raw_line
            .decode(
                "utf-8",
                errors="strict",
            )
            .rstrip(
                "\r\n"
            )
            .lstrip(
                "\ufeff"
            )
        )

    except UnicodeDecodeError as exc:

        raise HeaderResolutionError(
            "Remote harmonised GWAS column header is not valid UTF-8."
        ) from exc

    if "\t" not in line:

        raise HeaderResolutionError(
            "First decompressed line of remote harmonised GWAS file "
            "is not tab-delimited."
        )

    columns = tuple(
        field.strip()
        for field
        in line.split(
            "\t"
        )
    )

    if len(
        columns
    ) < 2:

        raise HeaderResolutionError(
            "Remote harmonised GWAS header contains fewer than two columns."
        )

    if any(
        not column
        for column
        in columns
    ):

        raise HeaderResolutionError(
            "Remote harmonised GWAS header contains an empty column name."
        )

    # ------------------------------------------------------------------
    # Duplicate-column safeguard
    # ------------------------------------------------------------------

    normalized_columns = [
        _normalize_header(
            column
        )
        for column
        in columns
    ]

    if len(
        normalized_columns
    ) != len(
        set(
            normalized_columns
        )
    ):

        raise HeaderResolutionError(
            "Remote harmonised GWAS header contains duplicated "
            "normalized column names."
        )

    # ------------------------------------------------------------------
    # Basic semantic safeguard
    #
    # We avoid accepting an arbitrary first data row as a header.
    # ------------------------------------------------------------------

    recognized_fields = (
        set(
            normalized_columns
        )
        & KNOWN_GWAS_HEADER_FIELDS
    )

    if len(
        recognized_fields
    ) < 2:

        raise HeaderResolutionError(
            "First decompressed line does not resemble a harmonised "
            "GWAS summary-statistics header."
        )

    return columns


# ============================================================================
# Remote capability probe
# ============================================================================


def probe_remote_tabix(
    url: str,
) -> RemoteTabixProbe:
    """
    Test whether pysam/HTSlib can open one remote Tabix-indexed file.

    No genomic locus is fetched.
    """

    try:

        with pysam.TabixFile(
            url
        ) as tabix_file:

            contigs = tuple(
                str(
                    value
                )
                for value
                in tabix_file.contigs
            )

    except Exception as exc:

        return RemoteTabixProbe(
            url=url,
            accessible=False,
            contigs=(),
            reason=(
                f"{type(exc).__name__}: {exc}"
            ),
        )

    return RemoteTabixProbe(
        url=url,
        accessible=True,
        contigs=contigs,
        reason=None,
    )


# ============================================================================
# Remote Tabix session
# ============================================================================


class PysamRemoteTabix:
    """
    Reusable indexed remote-access session for one harmonised GWAS study.

    One instance should normally be opened once per GWAS study and reused
    for all candidate loci belonging to that study.
    """

    def __init__(
        self,
        url: str,
    ) -> None:

        self.url = str(
            url
        )

        self._tabix: pysam.TabixFile | None = None

        self._contigs: tuple[
            str,
            ...
        ] = ()

        self._header_lines: tuple[
            str,
            ...
        ] = ()

        self._columns: tuple[
            str,
            ...
        ] = ()

    # ------------------------------------------------------------------
    # Context manager
    # ------------------------------------------------------------------

    def __enter__(
        self,
    ) -> PysamRemoteTabix:
        """Open the remote Tabix resource and resolve its schema."""

        try:

            self._tabix = (
                pysam.TabixFile(
                    self.url
                )
            )

        except Exception as exc:

            raise RemoteTabixOpenError(
                "Could not open remote Tabix GWAS resource "
                f"{self.url}: {type(exc).__name__}: {exc}"
            ) from exc

        try:

            # ----------------------------------------------------------
            # Indexed contigs
            # ----------------------------------------------------------

            self._contigs = tuple(
                str(
                    value
                )
                for value
                in self._tabix.contigs
            )

            if not self._contigs:

                raise RemoteTabixOpenError(
                    "Remote Tabix resource opened but index contains "
                    "zero contigs."
                )

            # ----------------------------------------------------------
            # Tabix metadata header, retained only for provenance/debug.
            #
            # It must NOT be relied upon for the actual TSV columns.
            # ----------------------------------------------------------

            self._header_lines = tuple(
                str(
                    line
                )
                for line
                in self._tabix.header
            )

            # ----------------------------------------------------------
            # Actual TSV column header.
            #
            # Read only the first decompressed BGZF line.
            # ----------------------------------------------------------

            self._columns = (
                _read_remote_column_header(
                    self.url
                )
            )

        except Exception:

            self.close()

            raise

        return self

    def __exit__(
        self,
        exc_type,
        exc_value,
        traceback,
    ) -> None:
        """Close the underlying remote Tabix resource."""

        self.close()

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def is_open(
        self,
    ) -> bool:
        """Return whether the indexed remote handle is open."""

        return (
            self._tabix
            is not None
        )

    @property
    def contigs(
        self,
    ) -> tuple[str, ...]:
        """Return contigs exposed by the Tabix index."""

        return self._contigs

    @property
    def header_lines(
        self,
    ) -> tuple[str, ...]:
        """Return metadata/header lines exposed by Tabix itself."""

        return self._header_lines

    @property
    def columns(
        self,
    ) -> tuple[str, ...]:
        """Return actual harmonised GWAS TSV column names."""

        return self._columns

    # ------------------------------------------------------------------
    # Close
    # ------------------------------------------------------------------

    def close(
        self,
    ) -> None:
        """Close the underlying pysam Tabix handle."""

        if self._tabix is not None:

            try:

                self._tabix.close()

            finally:

                self._tabix = None

    # ------------------------------------------------------------------
    # Contig resolver
    # ------------------------------------------------------------------

    def resolve_contig(
        self,
        chromosome: str,
    ) -> str:
        """Resolve chromosome notation against this study's Tabix index."""

        if not self.is_open:

            raise RemoteTabixError(
                "Remote Tabix session is not open."
            )

        return _resolve_contig(
            chromosome,
            self._contigs,
        )

    # ------------------------------------------------------------------
    # Fetch
    # ------------------------------------------------------------------

    def fetch_dataframe(
        self,
        *,
        chromosome: str,
        window_start_1based: int,
        window_end_1based: int,
    ) -> tuple[
        pd.DataFrame,
        LocusQueryMetadata,
    ]:
        """
        Retrieve one indexed GWAS locus as a pandas DataFrame.

        Input locus coordinates are biological 1-based inclusive.

        pysam fetch coordinates are 0-based half-open:

            fetch_start = window_start_1based - 1
            fetch_end   = window_end_1based
        """

        if self._tabix is None:

            raise RemoteTabixError(
                "Remote Tabix session is not open."
            )

        if not self._columns:

            raise HeaderResolutionError(
                "Harmonised GWAS column header has not been initialized."
            )

        start_1based = int(
            window_start_1based
        )

        end_1based = int(
            window_end_1based
        )

        if start_1based < 1:

            raise ValueError(
                "window_start_1based must be >= 1."
            )

        if end_1based < start_1based:

            raise ValueError(
                "window_end_1based must be >= window_start_1based."
            )

        # --------------------------------------------------------------
        # Resolve actual indexed contig.
        # --------------------------------------------------------------

        resolved_contig = (
            self.resolve_contig(
                chromosome
            )
        )

        # --------------------------------------------------------------
        # Explicit coordinate-system conversion.
        # --------------------------------------------------------------

        fetch_start_0based = (
            start_1based
            - 1
        )

        fetch_end_0based = (
            end_1based
        )

        metadata = (
            LocusQueryMetadata(
                requested_chromosome=str(
                    chromosome
                ),
                resolved_contig=resolved_contig,
                window_start_1based=start_1based,
                window_end_1based=end_1based,
                fetch_start_0based=fetch_start_0based,
                fetch_end_0based=fetch_end_0based,
            )
        )

        # --------------------------------------------------------------
        # Indexed remote query.
        # --------------------------------------------------------------

        try:

            iterator = (
                self._tabix.fetch(
                    reference=resolved_contig,
                    start=fetch_start_0based,
                    end=fetch_end_0based,
                )
            )

            split_rows: list[
                list[str]
            ] = []

            expected_field_count = (
                len(
                    self._columns
                )
            )

            malformed_row_count = 0

            for raw_row in iterator:

                text = str(
                    raw_row
                ).rstrip(
                    "\r\n"
                )

                if not text:
                    continue

                fields = text.split(
                    "\t"
                )

                if len(
                    fields
                ) != expected_field_count:

                    malformed_row_count += 1

                    continue

                split_rows.append(
                    fields
                )

        except Exception as exc:

            raise RemoteTabixError(
                "Remote Tabix locus query failed for "
                f"{resolved_contig}:"
                f"{start_1based}-{end_1based}: "
                f"{type(exc).__name__}: {exc}"
            ) from exc

        # --------------------------------------------------------------
        # We should not silently discard malformed rows.
        # --------------------------------------------------------------

        if malformed_row_count:

            raise HeaderResolutionError(
                "Fetched harmonised GWAS rows do not match resolved "
                "header width. "
                f"header_columns={len(self._columns)}, "
                f"malformed_rows={malformed_row_count}"
            )

        # --------------------------------------------------------------
        # Empty locus is a valid retrieval outcome.
        # Preserve schema.
        # --------------------------------------------------------------

        dataframe = pd.DataFrame(
            split_rows,
            columns=list(
                self._columns
            ),
        )

        return (
            dataframe,
            metadata,
        )


# ============================================================================
# Lead-presence helper
# ============================================================================


def find_lead_in_dataframe(
    dataframe: pd.DataFrame,
    lead_rsid: str,
) -> bool:
    """
    Determine whether one retrieved GWAS locus directly contains the lead rsID.

    This result is descriptive only.

    False does NOT mean:
        - no disease association exists in the region,
        - no proxy is present,
        - the candidate locus is unsupported,
        - the lead has no LD relationships.
    """

    if dataframe.empty:

        return False

    variant_column = (
        _resolve_column(
            dataframe.columns,
            VARIANT_ID_ALIASES,
        )
    )

    if variant_column is None:

        return False

    target = (
        str(
            lead_rsid
        )
        .strip()
        .lower()
    )

    values = (
        dataframe[
            variant_column
        ]
        .astype(
            "string"
        )
        .str.strip()
        .str.lower()
    )

    return bool(
        values
        .eq(
            target
        )
        .fillna(
            False
        )
        .any()
    )