"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/retrieve_gwas_loci.py

Description:
    Runner for M5.3C.3B-C harmonised GWAS locus-level summary-statistics
    retrieval.

    Workflow:
        M5.3C.3B
            Resolve one canonical GRCh38 coordinate for each candidate
            lead rsID using the configured multi-provider resolver.

        M5.3C.3C
            For each selected GWAS study:
                - use harmonised GWAS Catalog summary statistics,
                - open the remote BGZF/Tabix file using pysam/HTSlib,
                - reuse one remote Tabix session for all candidate loci,
                - retrieve the configured genomic window around each lead,
                - preserve one output Parquet file per successfully queried
                  study-lead locus.

    Coordinate policy:
        - Lead coordinates must explicitly resolve to GRCh38.
        - Candidate GRCh37 coordinates are not reused for harmonised queries.
        - No silent liftover is performed.
        - Biological locus windows are 1-based inclusive.
        - pysam Tabix queries use 0-based half-open coordinates internally.
        - Coordinate conversion is explicit and stored in the retrieval
          manifest.

    Retrieval policy:
        - pysam/HTSlib is the remote indexed retrieval backend.
        - External tabix executable is not used.
        - Genome-wide streaming fallback is disabled.
        - Whole GWAS files are not intentionally downloaded.
        - Studies without a discovered Tabix index are explicitly deferred.
        - Remote Tabix open/query failures are explicitly recorded.
        - Study-specific presence of the lead rsID is not required for
          regional locus retrieval.
        - Lead absence is not interpreted as absence of regional disease
          association.
        - Physical proximity is not interpreted as LD.
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

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq
import yaml

from pcatrfqtl.analysis.m5.gwas_lead_coordinate import (
    CanonicalLeadCoordinate,
    resolve_lead_coordinate,
)
from pcatrfqtl.analysis.m5.gwas_locus_retrieval import (
    ContigResolutionError,
    HeaderResolutionError,
    PysamRemoteTabix,
    RemoteTabixError,
    RemoteTabixOpenError,
    find_lead_in_dataframe,
)
from pcatrfqtl.io.parquet import (
    write_parquet,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


logger = get_logger(
    __name__
)


# ============================================================================
# Parquet compatibility reader
# ============================================================================


def _read_parquet_compat(
    path: Path,
) -> pd.DataFrame:
    """
    Read an M5 Parquet artifact while ignoring incompatible pandas metadata.

    Some earlier M5 artifacts contain Arrow-native list values. Using the
    physical Arrow schema while ignoring pandas reconstruction metadata avoids
    compatibility problems in certain pandas/PyArrow versions.
    """

    if not path.exists():

        raise FileNotFoundError(
            f"Parquet input not found: {path}"
        )

    table = pq.read_table(
        path,
        use_pandas_metadata=False,
    )

    dataframe = table.to_pandas(
        ignore_metadata=True,
    )

    if dataframe.empty:

        raise RuntimeError(
            f"Parquet input contains zero rows: {path}"
        )

    return dataframe


# ============================================================================
# Runner
# ============================================================================


class M53C3LocusRetrievalRunner:
    """
    Execute M5.3C.3B-C canonical lead-coordinate resolution and indexed
    GWAS locus retrieval.
    """

    def __init__(
        self,
        *,
        config_path: str | Path,
        remote_manifest_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.config_path = Path(
            config_path
        )

        self.remote_manifest_path = Path(
            remote_manifest_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    # ----------------------------------------------------------------------
    # Paths
    # ----------------------------------------------------------------------

    @property
    def coordinate_output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "gwas_lead_harmonised_coordinates.parquet"
        )

    @property
    def retrieval_manifest_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "gwas_locus_retrieval_manifest.parquet"
        )

    @property
    def locus_root(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "loci"
        )

    @property
    def qc_path(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_3c_3bc_gwas_locus_retrieval.json"
        )

    # ----------------------------------------------------------------------
    # Configuration
    # ----------------------------------------------------------------------

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load and minimally validate M5.3C.3 configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M5.3C.3 configuration not found: {self.config_path}"
            )

        with self.config_path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            config = yaml.safe_load(
                handle
            )

        if not isinstance(
            config,
            dict,
        ):

            raise ValueError(
                "M5.3C.3 configuration must be a YAML mapping."
            )

        if (
            "candidate_leads"
            not in config
        ):

            raise ValueError(
                "M5.3C.3 configuration is missing candidate_leads."
            )

        if (
            "locus_window_bp"
            not in config
        ):

            raise ValueError(
                "M5.3C.3 configuration is missing locus_window_bp."
            )

        return config

    # ----------------------------------------------------------------------
    # Coordinate resolution
    # ----------------------------------------------------------------------

    @staticmethod
    def _resolve_coordinates(
        leads: list[str],
    ) -> dict[
        str,
        CanonicalLeadCoordinate,
    ]:
        """
        Resolve one canonical GRCh38 coordinate per lead rsID.

        Resolution is study-independent and performed once per candidate lead.
        """

        coordinates: dict[
            str,
            CanonicalLeadCoordinate,
        ] = {}

        for lead in leads:

            logger.info(
                "Resolving canonical coordinate for %s.",
                lead,
            )

            coordinate = resolve_lead_coordinate(
                lead
            )

            coordinates[
                lead
            ] = coordinate

            if coordinate.found:

                logger.info(
                    "%s resolved to %s:%s (%s) "
                    "provider=%s method=%s fallback=%s attempts=%d.",
                    lead,
                    coordinate.chromosome,
                    coordinate.base_pair_location,
                    coordinate.genome_assembly,
                    coordinate.provider,
                    coordinate.lookup_method,
                    coordinate.fallback_used,
                    coordinate.attempt_count,
                )

            else:

                logger.warning(
                    "%s coordinate unresolved after %d attempt(s): %s",
                    lead,
                    coordinate.attempt_count,
                    coordinate.reason,
                )

        return coordinates

    # ----------------------------------------------------------------------
    # Record builders
    # ----------------------------------------------------------------------

    @staticmethod
    def _coordinate_unresolved_record(
        *,
        accession: str,
        lead: str,
        coordinate: CanonicalLeadCoordinate,
    ) -> dict[str, Any]:
        """Build one unresolved-coordinate retrieval record."""

        return {
            "study_accession":
                accession,

            "lead_rsid":
                lead,

            "chromosome":
                coordinate.chromosome,

            "resolved_contig":
                None,

            "lead_position":
                coordinate.base_pair_location,

            "genome_assembly":
                coordinate.genome_assembly,

            "coordinate_provider":
                coordinate.provider,

            "coordinate_lookup_method":
                coordinate.lookup_method,

            "coordinate_fallback_used":
                coordinate.fallback_used,

            "coordinate_attempt_count":
                coordinate.attempt_count,

            "window_start":
                None,

            "window_end":
                None,

            "fetch_start_0based":
                None,

            "fetch_end_0based":
                None,

            "variant_rows":
                0,

            "lead_present_in_locus":
                False,

            "retrieval_method":
                "NONE",

            "status":
                "LEAD_COORDINATE_UNRESOLVED",

            "reason":
                coordinate.reason,

            "output_path":
                None,
        }

    @staticmethod
    def _no_index_record(
        *,
        accession: str,
        lead: str,
        coordinate: CanonicalLeadCoordinate,
        window_bp: int,
    ) -> dict[str, Any]:
        """Build one intentionally deferred record for a study without index."""

        if (
            coordinate.found
            and coordinate.base_pair_location
            is not None
        ):

            start = max(
                1,
                int(
                    coordinate.base_pair_location
                )
                - window_bp,
            )

            end = (
                int(
                    coordinate.base_pair_location
                )
                + window_bp
            )

        else:

            start = None
            end = None

        return {
            "study_accession":
                accession,

            "lead_rsid":
                lead,

            "chromosome":
                coordinate.chromosome,

            "resolved_contig":
                None,

            "lead_position":
                coordinate.base_pair_location,

            "genome_assembly":
                coordinate.genome_assembly,

            "coordinate_provider":
                coordinate.provider,

            "coordinate_lookup_method":
                coordinate.lookup_method,

            "coordinate_fallback_used":
                coordinate.fallback_used,

            "coordinate_attempt_count":
                coordinate.attempt_count,

            "window_start":
                start,

            "window_end":
                end,

            "fetch_start_0based":
                (
                    start
                    - 1
                    if start is not None
                    else None
                ),

            "fetch_end_0based":
                end,

            "variant_rows":
                0,

            "lead_present_in_locus":
                False,

            "retrieval_method":
                "NONE",

            "status":
                "LOCUS_RETRIEVAL_DEFERRED_NO_INDEX",

            "reason":
                (
                    "Harmonised summary statistics are available, "
                    "but no Tabix index was discovered. Genome-wide "
                    "streaming fallback is intentionally disabled."
                ),

            "output_path":
                None,
        }

    @staticmethod
    def _study_failure_record(
        *,
        accession: str,
        lead: str,
        coordinate: CanonicalLeadCoordinate,
        window_bp: int,
        status: str,
        reason: str,
        retrieval_method: str,
    ) -> dict[str, Any]:
        """Build one indexed study/locus failure record."""

        if (
            coordinate.found
            and coordinate.base_pair_location
            is not None
        ):

            start = max(
                1,
                int(
                    coordinate.base_pair_location
                )
                - window_bp,
            )

            end = (
                int(
                    coordinate.base_pair_location
                )
                + window_bp
            )

        else:

            start = None
            end = None

        return {
            "study_accession":
                accession,

            "lead_rsid":
                lead,

            "chromosome":
                coordinate.chromosome,

            "resolved_contig":
                None,

            "lead_position":
                coordinate.base_pair_location,

            "genome_assembly":
                coordinate.genome_assembly,

            "coordinate_provider":
                coordinate.provider,

            "coordinate_lookup_method":
                coordinate.lookup_method,

            "coordinate_fallback_used":
                coordinate.fallback_used,

            "coordinate_attempt_count":
                coordinate.attempt_count,

            "window_start":
                start,

            "window_end":
                end,

            "fetch_start_0based":
                (
                    start
                    - 1
                    if start is not None
                    else None
                ),

            "fetch_end_0based":
                end,

            "variant_rows":
                0,

            "lead_present_in_locus":
                False,

            "retrieval_method":
                retrieval_method,

            "status":
                status,

            "reason":
                reason,

            "output_path":
                None,
        }

    # ----------------------------------------------------------------------
    # Main
    # ----------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M5.3C.3B-C."""

        # ==================================================================
        # Configuration
        # ==================================================================

        config = self._load_config()

        leads = [
            str(
                lead
            )
            .strip()
            .lower()
            for lead
            in config[
                "candidate_leads"
            ]
        ]

        if not leads:

            raise RuntimeError(
                "M5.3C.3 configuration contains zero candidate leads."
            )

        if len(
            set(
                leads
            )
        ) != len(
            leads
        ):

            raise ValueError(
                "M5.3C.3 candidate_leads contains duplicated rsIDs."
            )

        window_bp = int(
            config[
                "locus_window_bp"
            ]
        )

        if window_bp <= 0:

            raise ValueError(
                "locus_window_bp must be greater than zero."
            )

        # ==================================================================
        # Remote discovery manifest
        # ==================================================================

        remote = _read_parquet_compat(
            self.remote_manifest_path
        )

        required_remote_columns = {
            "study_accession",
            "harmonised_sumstats_url",
            "harmonised_available",
            "tabix_available",
            "retrieval_status",
        }

        missing_remote_columns = (
            required_remote_columns
            - set(
                remote.columns
            )
        )

        if missing_remote_columns:

            raise ValueError(
                "M5.3C.3A remote manifest is missing required columns: "
                f"{sorted(missing_remote_columns)}"
            )

        duplicated_studies = int(
            remote[
                "study_accession"
            ]
            .astype("string")
            .str.upper()
            .duplicated()
            .sum()
        )

        if duplicated_studies:

            raise RuntimeError(
                "M5.3C.3A remote manifest contains duplicated study rows: "
                f"{duplicated_studies}"
            )

        # ==================================================================
        # Output directories
        # ==================================================================

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.locus_root.mkdir(
            parents=True,
            exist_ok=True,
        )

        # ==================================================================
        # M5.3C.3B
        # Resolve coordinates once per lead.
        # ==================================================================

        coordinates = self._resolve_coordinates(
            leads
        )

        coordinate_records: list[
            dict[str, Any]
        ] = []

        for lead in leads:

            coordinate = coordinates[
                lead
            ]

            coordinate_records.append(
                {
                    "lead_rsid":
                        lead,

                    "found":
                        coordinate.found,

                    "chromosome":
                        coordinate.chromosome,

                    "base_pair_location":
                        coordinate.base_pair_location,

                    "genome_assembly":
                        coordinate.genome_assembly,

                    "provider":
                        coordinate.provider,

                    "lookup_method":
                        coordinate.lookup_method,

                    "fallback_used":
                        coordinate.fallback_used,

                    "attempt_count":
                        coordinate.attempt_count,

                    "source_url":
                        coordinate.source_url,

                    "reason":
                        coordinate.reason,
                }
            )

        coordinates_df = pd.DataFrame(
            coordinate_records
        )

        write_parquet(
            coordinates_df,
            self.coordinate_output_path,
            index=False,
        )

        # ==================================================================
        # M5.3C.3C
        # Remote indexed locus retrieval.
        # ==================================================================

        retrieval_records: list[
            dict[str, Any]
        ] = []

        remote_tabix_opened_studies = 0
        remote_tabix_open_failed_studies = 0

        for _, study in remote.iterrows():

            accession = (
                str(
                    study[
                        "study_accession"
                    ]
                )
                .strip()
                .upper()
            )

            logger.info(
                "Processing %s.",
                accession,
            )

            harmonised_available = bool(
                study.get(
                    "harmonised_available",
                    False,
                )
            )

            indexed = bool(
                study.get(
                    "tabix_available",
                    False,
                )
            )

            sumstats_url_value = study.get(
                "harmonised_sumstats_url"
            )

            sumstats_url = (
                None
                if pd.isna(
                    sumstats_url_value
                )
                else str(
                    sumstats_url_value
                )
            )

            # --------------------------------------------------------------
            # No harmonised summary statistics
            # --------------------------------------------------------------

            if (
                not harmonised_available
                or sumstats_url is None
            ):

                logger.warning(
                    "%s has no usable harmonised summary statistics.",
                    accession,
                )

                for lead in leads:

                    coordinate = coordinates[
                        lead
                    ]

                    retrieval_records.append(
                        self._study_failure_record(
                            accession=accession,
                            lead=lead,
                            coordinate=coordinate,
                            window_bp=window_bp,
                            status="NO_HARMONISED_SUMSTATS",
                            reason=(
                                "No usable harmonised summary-statistics "
                                "resource was discovered."
                            ),
                            retrieval_method="NONE",
                        )
                    )

                continue

            # --------------------------------------------------------------
            # Harmonised summary statistics but no index
            # --------------------------------------------------------------

            if not indexed:

                logger.warning(
                    "%s has harmonised summary statistics but no "
                    "Tabix index; retrieval deferred.",
                    accession,
                )

                for lead in leads:

                    coordinate = coordinates[
                        lead
                    ]

                    if not coordinate.found:

                        retrieval_records.append(
                            self._coordinate_unresolved_record(
                                accession=accession,
                                lead=lead,
                                coordinate=coordinate,
                            )
                        )

                    else:

                        retrieval_records.append(
                            self._no_index_record(
                                accession=accession,
                                lead=lead,
                                coordinate=coordinate,
                                window_bp=window_bp,
                            )
                        )

                continue

            # --------------------------------------------------------------
            # Open one remote pysam Tabix session per study.
            # --------------------------------------------------------------

            session = PysamRemoteTabix(
                sumstats_url
            )

            try:

                session.__enter__()

            except RemoteTabixOpenError as exc:

                remote_tabix_open_failed_studies += 1

                logger.warning(
                    "%s remote Tabix open failed: %s",
                    accession,
                    exc,
                )

                for lead in leads:

                    coordinate = coordinates[
                        lead
                    ]

                    if not coordinate.found:

                        retrieval_records.append(
                            self._coordinate_unresolved_record(
                                accession=accession,
                                lead=lead,
                                coordinate=coordinate,
                            )
                        )

                    else:

                        retrieval_records.append(
                            self._study_failure_record(
                                accession=accession,
                                lead=lead,
                                coordinate=coordinate,
                                window_bp=window_bp,
                                status="REMOTE_TABIX_OPEN_FAILED",
                                reason=str(
                                    exc
                                ),
                                retrieval_method="PYSAM_REMOTE_TABIX",
                            )
                        )

                continue

            remote_tabix_opened_studies += 1

            logger.info(
                "%s remote Tabix opened successfully; contigs=%d.",
                accession,
                len(
                    session.contigs
                ),
            )

            # --------------------------------------------------------------
            # Query three candidate loci through the same remote handle.
            # --------------------------------------------------------------

            try:

                for lead in leads:

                    coordinate = coordinates[
                        lead
                    ]

                    # ------------------------------------------------------
                    # Unresolved canonical coordinate
                    # ------------------------------------------------------

                    if (
                        not coordinate.found
                        or coordinate.chromosome
                        is None
                        or coordinate.base_pair_location
                        is None
                    ):

                        retrieval_records.append(
                            self._coordinate_unresolved_record(
                                accession=accession,
                                lead=lead,
                                coordinate=coordinate,
                            )
                        )

                        continue

                    # ------------------------------------------------------
                    # Biological window: 1-based inclusive
                    # ------------------------------------------------------

                    window_start = max(
                        1,
                        int(
                            coordinate.base_pair_location
                        )
                        - window_bp,
                    )

                    window_end = (
                        int(
                            coordinate.base_pair_location
                        )
                        + window_bp
                    )

                    logger.info(
                        "%s %s: querying %s:%d-%d (%s).",
                        accession,
                        lead,
                        coordinate.chromosome,
                        window_start,
                        window_end,
                        coordinate.genome_assembly,
                    )

                    # ------------------------------------------------------
                    # Indexed remote query
                    # ------------------------------------------------------

                    try:

                        (
                            dataframe,
                            query_metadata,
                        ) = session.fetch_dataframe(
                            chromosome=str(
                                coordinate.chromosome
                            ),
                            window_start_1based=window_start,
                            window_end_1based=window_end,
                        )

                    except ContigResolutionError as exc:

                        logger.warning(
                            "%s %s contig resolution failed: %s",
                            accession,
                            lead,
                            exc,
                        )

                        retrieval_records.append(
                            self._study_failure_record(
                                accession=accession,
                                lead=lead,
                                coordinate=coordinate,
                                window_bp=window_bp,
                                status="CONTIG_RESOLUTION_FAILED",
                                reason=str(
                                    exc
                                ),
                                retrieval_method="PYSAM_REMOTE_TABIX",
                            )
                        )

                        continue

                    except HeaderResolutionError as exc:

                        logger.warning(
                            "%s %s harmonised header invalid: %s",
                            accession,
                            lead,
                            exc,
                        )

                        retrieval_records.append(
                            self._study_failure_record(
                                accession=accession,
                                lead=lead,
                                coordinate=coordinate,
                                window_bp=window_bp,
                                status="HARMONISED_HEADER_INVALID",
                                reason=str(
                                    exc
                                ),
                                retrieval_method="PYSAM_REMOTE_TABIX",
                            )
                        )

                        continue

                    except RemoteTabixError as exc:

                        logger.warning(
                            "%s %s remote Tabix query failed: %s",
                            accession,
                            lead,
                            exc,
                        )

                        retrieval_records.append(
                            self._study_failure_record(
                                accession=accession,
                                lead=lead,
                                coordinate=coordinate,
                                window_bp=window_bp,
                                status="REMOTE_TABIX_QUERY_FAILED",
                                reason=str(
                                    exc
                                ),
                                retrieval_method="PYSAM_REMOTE_TABIX",
                            )
                        )

                        continue

                    except Exception as exc:

                        logger.warning(
                            "%s %s unexpected indexed query failure: %s",
                            accession,
                            lead,
                            exc,
                        )

                        retrieval_records.append(
                            self._study_failure_record(
                                accession=accession,
                                lead=lead,
                                coordinate=coordinate,
                                window_bp=window_bp,
                                status="REMOTE_TABIX_QUERY_FAILED",
                                reason=(
                                    f"{type(exc).__name__}: {exc}"
                                ),
                                retrieval_method="PYSAM_REMOTE_TABIX",
                            )
                        )

                        continue

                    # ------------------------------------------------------
                    # Persist locus
                    # ------------------------------------------------------

                    study_directory = (
                        self.locus_root
                        / accession
                    )

                    study_directory.mkdir(
                        parents=True,
                        exist_ok=True,
                    )

                    output_path = (
                        study_directory
                        / f"{lead}.parquet"
                    )

                    write_parquet(
                        dataframe,
                        output_path,
                        index=False,
                    )

                    lead_present = (
                        find_lead_in_dataframe(
                            dataframe,
                            lead,
                        )
                    )

                    status = (
                        "LOCUS_RETRIEVED"
                        if not dataframe.empty
                        else "LOCUS_EMPTY"
                    )

                    retrieval_records.append(
                        {
                            "study_accession":
                                accession,

                            "lead_rsid":
                                lead,

                            "chromosome":
                                coordinate.chromosome,

                            "resolved_contig":
                                query_metadata.resolved_contig,

                            "lead_position":
                                coordinate.base_pair_location,

                            "genome_assembly":
                                coordinate.genome_assembly,

                            "coordinate_provider":
                                coordinate.provider,

                            "coordinate_lookup_method":
                                coordinate.lookup_method,

                            "coordinate_fallback_used":
                                coordinate.fallback_used,

                            "coordinate_attempt_count":
                                coordinate.attempt_count,

                            "window_start":
                                query_metadata.window_start_1based,

                            "window_end":
                                query_metadata.window_end_1based,

                            "fetch_start_0based":
                                query_metadata.fetch_start_0based,

                            "fetch_end_0based":
                                query_metadata.fetch_end_0based,

                            "variant_rows":
                                int(
                                    len(
                                        dataframe
                                    )
                                ),

                            "lead_present_in_locus":
                                bool(
                                    lead_present
                                ),

                            "retrieval_method":
                                "PYSAM_REMOTE_TABIX",

                            "status":
                                status,

                            "reason":
                                None,

                            "output_path":
                                str(
                                    output_path
                                ),
                        }
                    )

                    logger.info(
                        "%s %s: %d rows retrieved; "
                        "contig=%s; lead present=%s.",
                        accession,
                        lead,
                        len(
                            dataframe
                        ),
                        query_metadata.resolved_contig,
                        lead_present,
                    )

            finally:

                session.close()

        # ==================================================================
        # Retrieval manifest
        # ==================================================================

        retrieval_df = pd.DataFrame(
            retrieval_records
        )

        expected_pairs = (
            len(
                remote
            )
            * len(
                leads
            )
        )

        if len(
            retrieval_df
        ) != expected_pairs:

            raise RuntimeError(
                "M5.3C.3B-C retrieval cardinality mismatch: "
                f"expected={expected_pairs}, "
                f"observed={len(retrieval_df)}"
            )

        write_parquet(
            retrieval_df,
            self.retrieval_manifest_path,
            index=False,
        )

        # ==================================================================
        # QC counts
        # ==================================================================

        coordinates_resolved = int(
            coordinates_df[
                "found"
            ]
            .fillna(
                False
            )
            .sum()
        )

        fallback_coordinates = int(
            coordinates_df[
                "fallback_used"
            ]
            .fillna(
                False
            )
            .sum()
        )

        indexed_studies = int(
            remote[
                "tabix_available"
            ]
            .fillna(
                False
            )
            .sum()
        )

        indexed_pairs = (
            indexed_studies
            * len(
                leads
            )
        )

        loci_retrieved = int(
            retrieval_df[
                "status"
            ]
            .eq(
                "LOCUS_RETRIEVED"
            )
            .sum()
        )

        loci_empty = int(
            retrieval_df[
                "status"
            ]
            .eq(
                "LOCUS_EMPTY"
            )
            .sum()
        )

        deferred_no_index = int(
            retrieval_df[
                "status"
            ]
            .eq(
                "LOCUS_RETRIEVAL_DEFERRED_NO_INDEX"
            )
            .sum()
        )

        remote_open_failed_pairs = int(
            retrieval_df[
                "status"
            ]
            .eq(
                "REMOTE_TABIX_OPEN_FAILED"
            )
            .sum()
        )

        remote_query_failed = int(
            retrieval_df[
                "status"
            ]
            .eq(
                "REMOTE_TABIX_QUERY_FAILED"
            )
            .sum()
        )

        contig_failed = int(
            retrieval_df[
                "status"
            ]
            .eq(
                "CONTIG_RESOLUTION_FAILED"
            )
            .sum()
        )

        header_failed = int(
            retrieval_df[
                "status"
            ]
            .eq(
                "HARMONISED_HEADER_INVALID"
            )
            .sum()
        )

        unresolved_pairs = int(
            retrieval_df[
                "status"
            ]
            .eq(
                "LEAD_COORDINATE_UNRESOLVED"
            )
            .sum()
        )

        lead_present = int(
            retrieval_df[
                "lead_present_in_locus"
            ]
            .fillna(
                False
            )
            .sum()
        )

        total_variant_rows = int(
            pd.to_numeric(
                retrieval_df[
                    "variant_rows"
                ],
                errors="coerce",
            )
            .fillna(
                0
            )
            .sum()
        )

        # ==================================================================
        # Distributions
        # ==================================================================

        status_counts = (
            retrieval_df[
                "status"
            ]
            .value_counts(
                dropna=False
            )
            .to_dict()
        )

        method_counts = (
            retrieval_df[
                "retrieval_method"
            ]
            .value_counts(
                dropna=False
            )
            .to_dict()
        )

        coordinate_provider_counts = (
            coordinates_df[
                "provider"
            ]
            .fillna(
                "UNRESOLVED"
            )
            .value_counts(
                dropna=False
            )
            .to_dict()
        )

        coordinate_method_counts = (
            coordinates_df[
                "lookup_method"
            ]
            .fillna(
                "UNRESOLVED"
            )
            .value_counts(
                dropna=False
            )
            .to_dict()
        )

        # ==================================================================
        # QC report
        # ==================================================================

        report = {
            "milestone":
                "M5.3C.3B-C",

            "stage":
                "canonical_lead_resolution_and_pysam_indexed_locus_retrieval",

            "policy": {
                "candidate_grch37_coordinates_used_for_query":
                    False,

                "silent_liftover_performed":
                    False,

                "canonical_coordinate_resolved_once_per_lead":
                    True,

                "coordinate_resolution_multi_provider":
                    True,

                "coordinate_provider_fallback_allowed":
                    True,

                "coordinate_target_assembly":
                    "GRCh38",

                "study_specific_lead_presence_required_for_locus_query":
                    False,

                "harmonised_sumstats_used":
                    True,

                "remote_index_backend":
                    "pysam_htslib",

                "external_tabix_executable_used":
                    False,

                "remote_bgzf_random_access_used":
                    True,

                "genome_wide_sumstats_downloaded":
                    False,

                "genome_wide_streaming_fallback_allowed":
                    False,

                "query_failure_triggers_streaming":
                    False,

                "no_index_study_deferred":
                    True,

                "biological_window_coordinate_system":
                    "1_based_inclusive",

                "pysam_fetch_coordinate_system":
                    "0_based_half_open",

                "coordinate_conversion_explicit":
                    True,

                "physical_proximity_interpreted_as_ld":
                    False,

                "colocalization_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "studies":
                    len(
                        remote
                    ),

                "candidate_leads":
                    len(
                        leads
                    ),

                "expected_study_lead_pairs":
                    expected_pairs,

                "canonical_lead_coordinates_resolved":
                    coordinates_resolved,

                "canonical_lead_coordinates_unresolved":
                    len(
                        leads
                    )
                    - coordinates_resolved,

                "canonical_coordinates_using_fallback":
                    fallback_coordinates,

                "indexed_studies":
                    indexed_studies,

                "indexed_study_lead_pairs":
                    indexed_pairs,

                "remote_tabix_opened_studies":
                    remote_tabix_opened_studies,

                "remote_tabix_open_failed_studies":
                    remote_tabix_open_failed_studies,

                "loci_retrieved":
                    loci_retrieved,

                "loci_empty":
                    loci_empty,

                "lead_coordinate_unresolved_pairs":
                    unresolved_pairs,

                "remote_tabix_open_failed_pairs":
                    remote_open_failed_pairs,

                "remote_tabix_query_failed":
                    remote_query_failed,

                "contig_resolution_failed":
                    contig_failed,

                "harmonised_header_invalid":
                    header_failed,

                "deferred_no_index":
                    deferred_no_index,

                "lead_present_in_retrieved_locus":
                    lead_present,

                "total_retrieved_variant_rows":
                    total_variant_rows,
            },

            "coordinate_provider_counts": {
                str(
                    key
                ):
                    int(
                        value
                    )
                for key, value
                in coordinate_provider_counts.items()
            },

            "coordinate_lookup_method_counts": {
                str(
                    key
                ):
                    int(
                        value
                    )
                for key, value
                in coordinate_method_counts.items()
            },

            "retrieval_status_counts": {
                str(
                    key
                ):
                    int(
                        value
                    )
                for key, value
                in status_counts.items()
            },

            "retrieval_method_counts": {
                str(
                    key
                ):
                    int(
                        value
                    )
                for key, value
                in method_counts.items()
            },

            "coordinate_records":
                coordinates_df[
                    [
                        "lead_rsid",
                        "chromosome",
                        "base_pair_location",
                        "genome_assembly",
                        "provider",
                        "lookup_method",
                        "fallback_used",
                        "attempt_count",
                    ]
                ]
                .astype("string")
                .to_dict(
                    orient="records"
                ),

            "unresolved_leads":
                coordinates_df.loc[
                    ~coordinates_df[
                        "found"
                    ]
                    .fillna(
                        False
                    ),
                    "lead_rsid",
                ]
                .astype(str)
                .tolist(),

            "failed_pairs":
                retrieval_df.loc[
                    retrieval_df[
                        "status"
                    ]
                    .isin(
                        [
                            "LEAD_COORDINATE_UNRESOLVED",
                            "REMOTE_TABIX_OPEN_FAILED",
                            "REMOTE_TABIX_QUERY_FAILED",
                            "CONTIG_RESOLUTION_FAILED",
                            "HARMONISED_HEADER_INVALID",
                        ]
                    ),
                    [
                        "study_accession",
                        "lead_rsid",
                        "status",
                        "reason",
                    ],
                ]
                .astype("string")
                .to_dict(
                    orient="records"
                ),

            "deferred_no_index_pairs":
                retrieval_df.loc[
                    retrieval_df[
                        "status"
                    ]
                    .eq(
                        "LOCUS_RETRIEVAL_DEFERRED_NO_INDEX"
                    ),
                    [
                        "study_accession",
                        "lead_rsid",
                    ],
                ]
                .astype(str)
                .to_dict(
                    orient="records"
                ),

            "coordinate_output":
                str(
                    self.coordinate_output_path
                ),

            "retrieval_manifest":
                str(
                    self.retrieval_manifest_path
                ),

            "locus_output_root":
                str(
                    self.locus_root
                ),
        }

        # ==================================================================
        # Write QC
        # ==================================================================

        with self.qc_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        # ==================================================================
        # Logging
        # ==================================================================

        logger.info(
            "M5.3C.3B-C complete."
        )

        logger.info(
            "Canonical coordinates resolved: %d/%d | fallback used: %d",
            coordinates_resolved,
            len(
                leads
            ),
            fallback_coordinates,
        )

        logger.info(
            "Remote indexed studies opened: %d/%d.",
            remote_tabix_opened_studies,
            indexed_studies,
        )

        logger.info(
            "Loci retrieved: %d | Empty: %d | Deferred no index: %d",
            loci_retrieved,
            loci_empty,
            deferred_no_index,
        )

        logger.info(
            "Failures — open: %d | query: %d | "
            "contig: %d | header: %d | unresolved coordinate pairs: %d",
            remote_tabix_open_failed_studies,
            remote_query_failed,
            contig_failed,
            header_failed,
            unresolved_pairs,
        )

        logger.info(
            "Lead variants directly present in retrieved loci: %d",
            lead_present,
        )

        logger.info(
            "Total locus rows retrieved: %d",
            total_variant_rows,
        )

        logger.info(
            "Coordinate manifest: %s",
            self.coordinate_output_path,
        )

        logger.info(
            "Retrieval manifest: %s",
            self.retrieval_manifest_path,
        )

        logger.info(
            "Locus root: %s",
            self.locus_root,
        )

        logger.info(
            "QC output: %s",
            self.qc_path,
        )

        return report