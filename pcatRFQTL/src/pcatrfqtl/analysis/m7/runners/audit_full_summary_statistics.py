"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/audit_full_summary_statistics.py

Description:
    Runner for M7.5B Full Summary Statistics Audit.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin
from urllib.request import Request, urlopen

import pandas as pd
import yaml

from pcatrfqtl.analysis.m7.full_summary_statistics_audit import (
    assess_full_summary_statistics_audit,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


STUDY_COLUMNS = [
    "accession_id",
    "pubmed_id",
    "disease_trait",
    "efo_ids",
    "efo_traits",
    "initial_sample_size",
    "replication_sample_size",
    "discovery_ancestry",
    "replication_ancestry",
    "cohort",
    "platforms",
    "genotyping_technologies",
    "array_manufacturer",
    "imputed",
    "gxe",
    "gxg",
    "terms_of_license",
    "full_summary_stats_available",
    "full_summary_stats",
    "direct_prostate_phenotype",
    "phenotype_screening_status",
    "eligible_for_exact_lookup",
    "study_status",
    "phenotype_compatibility_verified",
    "ancestry_compatibility_verified",
    "study_independence_verified",
]


FILE_COLUMNS = [
    "accession_id",
    "source_directory_url",
    "file_url",
    "filename",
    "file_representation",
    "compression",
    "harmonised",
    "candidate_files_discovered",
    "formatted_files_available",
    "selected_for_scan",
    "scan_success",
    "eof_reached",
    "scan_error",
    "source_file_sha256",
    "candidate_hits",
    "file_status",
]


HIT_COLUMNS = [
    "accession_id",
    "rs_id",
    "source_file_url",
    "exact_rsid_match",
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "effect_allele_frequency",
    "beta",
    "odds_ratio",
    "standard_error",
    "p_value",
    "neg_log_10_p_value",
]


@dataclass(frozen=True)
class HttpResult:
    """Simple HTTP response container."""

    success: bool
    status_code: int | None
    payload: bytes | None
    error_type: str | None
    error_message: str | None


class HashingReader(io.RawIOBase):
    """
    Read-through SHA-256 wrapper.

    Bytes are hashed exactly as transferred from the source without loading
    the complete file into memory.
    """

    def __init__(
        self,
        raw: Any,
        hasher: Any,
    ) -> None:

        super().__init__()

        self.raw = raw
        self.hasher = hasher

    def readable(
        self,
    ) -> bool:

        return True

    def readinto(
        self,
        buffer: bytearray,
    ) -> int:

        data = self.raw.read(
            len(buffer)
        )

        if not data:

            return 0

        self.hasher.update(
            data
        )

        length = len(data)

        buffer[
            :length
        ] = data

        return length


def _json_safe_value(
    value: Any,
) -> Any:

    if value is None:

        return None

    if isinstance(
        value,
        (
            str,
            int,
            float,
            bool,
            list,
            dict,
        ),
    ):

        return value

    if hasattr(
        value,
        "item",
    ):

        return value.item()

    return str(value)


def _records(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:

    output: list[dict[str, Any]] = []

    for _, row in dataframe.iterrows():

        item: dict[str, Any] = {}

        for key, value in row.to_dict().items():

            try:

                if (
                    not isinstance(
                        value,
                        (
                            list,
                            dict,
                        ),
                    )
                    and
                    pd.isna(value)
                ):

                    value = None

            except (
                TypeError,
                ValueError,
            ):
                pass

            item[str(key)] = _json_safe_value(
                value
            )

        output.append(item)

    return output


class M75BFullSummaryStatisticsAuditRunner:
    """Execute M7.5B."""

    def __init__(
        self,
        *,
        project_root: str | Path,
        config_path: str | Path,
    ) -> None:

        self.project_root = Path(
            project_root
        ).resolve()

        self.config_path = Path(
            config_path
        ).resolve()

    # ======================================================================
    # Config
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:

        with self.config_path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            config = yaml.safe_load(handle)

        if not isinstance(
            config,
            dict,
        ):

            raise RuntimeError(
                "M7.5B config must be a mapping."
            )

        if config.get(
            "milestone"
        ) != "M7.5B":

            raise RuntimeError(
                "Unexpected milestone."
            )

        if config.get(
            "stage"
        ) != "full_summary_statistics_audit":

            raise RuntimeError(
                "Unexpected stage."
            )

        return config

    # ======================================================================
    # Upstream
    # ======================================================================

    def _validate_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> dict[str, Any]:

        spec = config[
            "upstream_qc"
        ][
            "m7_5"
        ]

        path = (
            self.project_root
            /
            spec["relative_path"]
        )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            payload = json.load(handle)

        if payload.get(
            "milestone"
        ) != spec[
            "expected_milestone"
        ]:

            raise RuntimeError(
                "M7.5 milestone mismatch."
            )

        if payload.get(
            "stage"
        ) != spec[
            "expected_stage"
        ]:

            raise RuntimeError(
                "M7.5 stage mismatch."
            )

        summary = payload.get(
            "summary"
        )

        if not isinstance(
            summary,
            dict,
        ):

            raise RuntimeError(
                "M7.5 summary missing."
            )

        status_field = str(
            spec[
                "status_field"
            ]
        )

        if summary.get(
            status_field
        ) != spec[
            "expected_status"
        ]:

            raise RuntimeError(
                "M7.5 status mismatch."
            )

        for field, expected in spec[
            "required_summary_state"
        ].items():

            if summary.get(
                field
            ) != expected:

                raise RuntimeError(
                    f"M7.5 state mismatch: {field}."
                )

        logger.info(
            "Validated locked M7.5 | %s=%s.",
            status_field,
            summary[status_field],
        )

        return {
            "path":
                str(path),

            "status":
                summary[status_field],
        }

    # ======================================================================
    # Small HTTP requests
    # ======================================================================

    def _request_bytes(
        self,
        *,
        url: str,
        timeout: int,
        user_agent: str,
        retries: int = 3,
        backoff: int = 3,
    ) -> HttpResult:

        last_status = None
        last_type = None
        last_message = None

        for attempt in range(
            1,
            retries + 1,
        ):

            request = Request(
                url,
                headers={
                    "User-Agent":
                        user_agent,

                    "Accept":
                        "*/*",
                },
            )

            try:

                with urlopen(
                    request,
                    timeout=timeout,
                ) as response:

                    return HttpResult(
                        success=True,
                        status_code=int(
                            response.status
                        ),
                        payload=response.read(),
                        error_type=None,
                        error_message=None,
                    )

            except HTTPError as exc:

                last_status = int(
                    exc.code
                )

                last_type = "HTTP_ERROR"

                last_message = (
                    f"HTTP {exc.code}: {exc.reason}"
                )

            except URLError as exc:

                last_type = "URL_ERROR"
                last_message = str(
                    exc.reason
                )

            if attempt < retries:

                time.sleep(
                    backoff * attempt
                )

        return HttpResult(
            success=False,
            status_code=last_status,
            payload=None,
            error_type=last_type,
            error_message=last_message,
        )

    # ======================================================================
    # Study API parsing
    # ======================================================================

    @staticmethod
    def _extract_page(
        payload: dict[str, Any],
    ) -> dict[str, int]:

        page = payload.get("page")

        if not isinstance(
            page,
            dict,
        ):

            raise RuntimeError(
                "Study response missing page metadata."
            )

        total_elements = page.get(
            "total_elements",
            page.get("totalElements"),
        )

        total_pages = page.get(
            "total_pages",
            page.get("totalPages"),
        )

        if total_elements is None:

            raise RuntimeError(
                "Missing totalElements."
            )

        if total_pages is None:

            raise RuntimeError(
                "Missing totalPages."
            )

        return {
            "size":
                int(
                    page.get(
                        "size",
                        0,
                    )
                ),

            "total_elements":
                int(total_elements),

            "total_pages":
                int(total_pages),

            "number":
                int(
                    page.get(
                        "number",
                        0,
                    )
                ),
        }

    @staticmethod
    def _extract_studies(
        payload: dict[str, Any],
    ) -> list[dict[str, Any]]:

        if "content" in payload:

            content = payload["content"]

            if not isinstance(
                content,
                list,
            ):

                raise RuntimeError(
                    "Study content must be a list."
                )

            return [
                row
                for row in content
                if isinstance(
                    row,
                    dict,
                )
            ]

        embedded = payload.get(
            "_embedded"
        )

        if isinstance(
            embedded,
            dict,
        ):

            studies = embedded.get(
                "studies"
            )

            if isinstance(
                studies,
                list,
            ):

                return [
                    row
                    for row in studies
                    if isinstance(
                        row,
                        dict,
                    )
                ]

        page = (
            M75BFullSummaryStatisticsAuditRunner
            ._extract_page(payload)
        )

        if page[
            "total_elements"
        ] == 0:

            return []

        raise RuntimeError(
            "Study collection missing."
        )

    # ======================================================================
    # Immutable study snapshots
    # ======================================================================

    def _study_snapshot_path(
        self,
        *,
        page_number: int,
        config: dict[str, Any],
    ) -> Path:

        directory = (
            self.project_root
            /
            config[
                "study_discovery"
            ][
                "raw_snapshot_directory"
            ]
        )

        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        return (
            directory
            /
            f"page_{page_number:03d}.json"
        )

    def _load_or_fetch_study_page(
        self,
        *,
        page_number: int,
        config: dict[str, Any],
    ) -> dict[str, Any]:

        spec = config[
            "study_discovery"
        ]

        snapshot = self._study_snapshot_path(
            page_number=page_number,
            config=config,
        )

        if (
            bool(
                spec[
                    "reuse_existing_raw_snapshots"
                ]
            )
            and
            snapshot.exists()
        ):

            with snapshot.open(
                "r",
                encoding="utf-8",
            ) as handle:

                payload = json.load(handle)

            return payload

        query = {
            spec["trait_parameter"]:
                spec["trait_query"],

            spec["show_child_traits_parameter"]:
                str(
                    bool(
                        spec[
                            "show_child_traits"
                        ]
                    )
                ).lower(),

            spec["page_parameter"]:
                page_number,

            spec["size_parameter"]:
                int(
                    spec[
                        "page_size"
                    ]
                ),
        }

        url = (
            str(
                spec["base_url"]
            ).rstrip("/")
            +
            str(
                spec["endpoint"]
            )
            +
            "?"
            +
            urlencode(query)
        )

        response = self._request_bytes(
            url=url,
            timeout=int(
                spec[
                    "timeout_seconds"
                ]
            ),
            user_agent=str(
                spec[
                    "user_agent"
                ]
            ),
            retries=int(
                spec[
                    "max_retries"
                ]
            ),
            backoff=int(
                spec[
                    "retry_backoff_seconds"
                ]
            ),
        )

        if (
            not response.success
            or
            response.payload is None
        ):

            raise RuntimeError(
                "Study discovery request failed."
            )

        payload = json.loads(
            response.payload.decode(
                "utf-8"
            )
        )

        if snapshot.exists():

            raise RuntimeError(
                f"Refusing overwrite: {snapshot}"
            )

        snapshot.write_text(
            json.dumps(
                payload,
                indent=2,
                ensure_ascii=False,
                sort_keys=True,
            ),
            encoding="utf-8",
        )

        return payload

    def _discover_studies(
        self,
        *,
        config: dict[str, Any],
    ) -> list[dict[str, Any]]:

        spec = config[
            "study_discovery"
        ]

        page_number = int(
            spec["start_page"]
        )

        records: list[dict[str, Any]] = []

        expected_total = None
        expected_pages = None

        while True:

            payload = self._load_or_fetch_study_page(
                page_number=page_number,
                config=config,
            )

            page = self._extract_page(
                payload
            )

            studies = self._extract_studies(
                payload
            )

            if expected_total is None:

                expected_total = page[
                    "total_elements"
                ]

                expected_pages = page[
                    "total_pages"
                ]

            else:

                if page[
                    "total_elements"
                ] != expected_total:

                    raise RuntimeError(
                        "totalElements changed."
                    )

                if page[
                    "total_pages"
                ] != expected_pages:

                    raise RuntimeError(
                        "totalPages changed."
                    )

            records.extend(studies)

            if expected_pages == 0:

                break

            if (
                page_number + 1
                >=
                expected_pages
            ):

                break

            page_number += 1

            rate = float(
                spec[
                    "client_queries_per_second"
                ]
            )

            if rate > 0:

                time.sleep(
                    1.0 / rate
                )

        unique: dict[str, dict[str, Any]] = {}

        for record in records:

            accession = record.get(
                "accession_id"
            )

            if accession:

                unique[
                    str(accession)
                ] = record

        studies = list(
            unique.values()
        )

        if (
            expected_total is not None
            and
            len(studies)
            != expected_total
        ):

            raise RuntimeError(
                "Study count mismatch."
            )

        return studies

    # ======================================================================
    # Study metadata
    # ======================================================================

    @staticmethod
    def _normalize_list(
        value: Any,
    ) -> list[str]:

        if not isinstance(
            value,
            list,
        ):

            return []

        return list(
            dict.fromkeys(
                str(item).strip()
                for item in value
                if item is not None
                and str(item).strip()
            )
        )

    def _classify_phenotype_eligibility(
        self,
        *,
        disease_trait: Any,
        config: dict[str, Any],
    ) -> tuple[bool, str]:

        spec = config[
            "phenotype_eligibility"
        ]

        statuses = config[
            "statuses"
        ][
            "phenotype"
        ]

        trait = (
            str(disease_trait).strip()
            if disease_trait is not None
            else ""
        )

        trait_lower = trait.lower()

        for fragment in spec[
            "excluded_trait_fragments"
        ]:

            if (
                str(fragment).lower()
                in
                trait_lower
            ):

                return (
                    False,
                    statuses["excluded"],
                )

        accepted_exact = {
            str(value).lower()
            for value in spec[
                "accepted_exact_traits"
            ]
        }

        if trait_lower in accepted_exact:

            return (
                True,
                statuses["direct"],
            )

        for prefix in spec[
            "accepted_trait_prefixes"
        ]:

            if trait_lower.startswith(
                str(prefix).lower()
            ):

                return (
                    True,
                    statuses["direct"],
                )

        return (
            False,
            statuses["unresolved"],
        )

    def _build_study_inventory(
        self,
        *,
        studies: list[dict[str, Any]],
        config: dict[str, Any],
    ) -> pd.DataFrame:

        rows: list[dict[str, Any]] = []

        study_statuses = config[
            "statuses"
        ][
            "study"
        ]

        phenotype_statuses = config[
            "statuses"
        ][
            "phenotype"
        ]

        for record in studies:

            accession = record.get(
                "accession_id"
            )

            if not accession:

                continue

            efo_ids: list[str] = []
            efo_traits: list[str] = []

            for item in record.get(
                "efo_traits",
                [],
            ):

                if not isinstance(
                    item,
                    dict,
                ):

                    continue

                if item.get(
                    "efo_id"
                ):

                    efo_ids.append(
                        str(
                            item["efo_id"]
                        )
                    )

                if item.get(
                    "efo_trait"
                ):

                    efo_traits.append(
                        str(
                            item["efo_trait"]
                        )
                    )

            full_available = (
                record.get(
                    "full_summary_stats_available"
                )
                is True
            )

            (
                direct_phenotype,
                phenotype_status,
            ) = self._classify_phenotype_eligibility(
                disease_trait=record.get(
                    "disease_trait"
                ),
                config=config,
            )

            eligible = bool(
                full_available
                and
                direct_phenotype
            )

            if eligible:

                study_status = study_statuses[
                    "eligible"
                ]

            elif not full_available:

                study_status = study_statuses[
                    "no_full_sumstats"
                ]

            elif (
                phenotype_status
                ==
                phenotype_statuses[
                    "excluded"
                ]
            ):

                study_status = study_statuses[
                    "phenotype_excluded"
                ]

            else:

                study_status = study_statuses[
                    "phenotype_unresolved"
                ]

            rows.append(
                {
                    "accession_id":
                        str(accession),

                    "pubmed_id":
                        record.get(
                            "pubmed_id"
                        ),

                    "disease_trait":
                        record.get(
                            "disease_trait"
                        ),

                    "efo_ids":
                        list(
                            dict.fromkeys(
                                efo_ids
                            )
                        ),

                    "efo_traits":
                        list(
                            dict.fromkeys(
                                efo_traits
                            )
                        ),

                    "initial_sample_size":
                        record.get(
                            "initial_sample_size"
                        ),

                    "replication_sample_size":
                        record.get(
                            "replication_sample_size"
                        ),

                    "discovery_ancestry":
                        self._normalize_list(
                            record.get(
                                "discovery_ancestry"
                            )
                        ),

                    "replication_ancestry":
                        self._normalize_list(
                            record.get(
                                "replication_ancestry"
                            )
                        ),

                    "cohort":
                        self._normalize_list(
                            record.get(
                                "cohort"
                            )
                        ),

                    "platforms":
                        record.get(
                            "platforms"
                        ),

                    "genotyping_technologies":
                        self._normalize_list(
                            record.get(
                                "genotyping_technologies"
                            )
                        ),

                    "array_manufacturer":
                        self._normalize_list(
                            record.get(
                                "array_manufacturer"
                            )
                        ),

                    "imputed":
                        record.get(
                            "imputed"
                        ),

                    "gxe":
                        record.get(
                            "gxe"
                        ),

                    "gxg":
                        record.get(
                            "gxg"
                        ),

                    "terms_of_license":
                        record.get(
                            "terms_of_license"
                        ),

                    "full_summary_stats_available":
                        full_available,

                    "full_summary_stats":
                        record.get(
                            "full_summary_stats"
                        ),

                    "direct_prostate_phenotype":
                        direct_phenotype,

                    "phenotype_screening_status":
                        phenotype_status,

                    "eligible_for_exact_lookup":
                        eligible,

                    "study_status":
                        study_status,

                    "phenotype_compatibility_verified":
                        False,

                    "ancestry_compatibility_verified":
                        False,

                    "study_independence_verified":
                        False,
                }
            )

        return (
            pd.DataFrame(
                rows,
                columns=STUDY_COLUMNS,
            )
            .drop_duplicates(
                subset=[
                    "accession_id",
                ]
            )
            .sort_values(
                "accession_id",
                kind="stable",
            )
            .reset_index(
                drop=True
            )
        )

    # ======================================================================
    # File discovery
    # ======================================================================

    def _list_directory(
        self,
        *,
        url: str,
        config: dict[str, Any],
    ) -> list[str]:

        source = config[
            "summary_statistics"
        ][
            "public_collection"
        ]

        response = self._request_bytes(
            url=url,
            timeout=int(
                source[
                    "timeout_seconds"
                ]
            ),
            user_agent=str(
                source[
                    "user_agent"
                ]
            ),
        )

        if (
            not response.success
            or
            response.payload is None
        ):

            return []

        html = response.payload.decode(
            "utf-8",
            errors="replace",
        )

        return re.findall(
            r'href=["\']([^"\']+)["\']',
            html,
            flags=re.IGNORECASE,
        )

    @staticmethod
    def _clean_directory_files(
        values: list[str],
    ) -> list[str]:

        return [
            value
            for value in values
            if value
            and
            not value.startswith("?")
            and
            not value.startswith("/")
            and
            not value.endswith("/")
        ]

    def _discover_file(
        self,
        *,
        accession: str,
        source_directory_url: str,
        config: dict[str, Any],
    ) -> dict[str, Any]:

        selection = config[
            "summary_statistics"
        ][
            "file_selection"
        ]

        statuses = config[
            "statuses"
        ][
            "file"
        ]

        root_url = (
            source_directory_url.rstrip("/")
            +
            "/"
        )

        root_links = self._list_directory(
            url=root_url,
            config=config,
        )

        root_files = self._clean_directory_files(
            root_links
        )

        harmonised_dir = (
            str(
                selection[
                    "harmonised_subdirectory"
                ]
            ).rstrip("/")
            +
            "/"
        )

        harmonised_files: list[str] = []

        # ------------------------------------------------------------------
        # Priority 1: *.h.tsv.gz
        # ------------------------------------------------------------------

        if harmonised_dir in root_links:

            harmonised_url = urljoin(
                root_url,
                harmonised_dir,
            )

            harmonised_links = self._list_directory(
                url=harmonised_url,
                config=config,
            )

            harmonised_files = (
                self._clean_directory_files(
                    harmonised_links
                )
            )

            h_files = sorted(
                filename
                for filename in harmonised_files
                if filename.endswith(
                    str(
                        selection[
                            "harmonised_preferred_suffix"
                        ]
                    )
                )
                and
                "-meta.yaml"
                not in
                filename.lower()
            )

            formatted_files = sorted(
                filename
                for filename in harmonised_files
                if filename.endswith(
                    str(
                        selection[
                            "formatted_suffix"
                        ]
                    )
                )
            )

            if h_files:

                filename = h_files[0]

                return {
                    "accession_id":
                        accession,

                    "source_directory_url":
                        harmonised_url,

                    "file_url":
                        urljoin(
                            harmonised_url,
                            filename,
                        ),

                    "filename":
                        filename,

                    "file_representation":
                        "HARMONISED",

                    "compression":
                        "GZIP",

                    "harmonised":
                        True,

                    "candidate_files_discovered":
                        len(
                            h_files
                        ),

                    "formatted_files_available":
                        len(
                            formatted_files
                        ),

                    "selected_for_scan":
                        True,

                    "scan_success":
                        False,

                    "eof_reached":
                        False,

                    "scan_error":
                        None,

                    "source_file_sha256":
                        None,

                    "candidate_hits":
                        0,

                    "file_status":
                        statuses[
                            "harmonised"
                        ],
                }

        # ------------------------------------------------------------------
        # Priority 2: root *.tsv.gz
        # ------------------------------------------------------------------

        root_tsv_gz = sorted(
            filename
            for filename in root_files
            if filename.lower().endswith(
                ".tsv.gz"
            )
            and
            "-meta.yaml"
            not in
            filename.lower()
        )

        if root_tsv_gz:

            filename = root_tsv_gz[0]

            return {
                "accession_id":
                    accession,

                "source_directory_url":
                    root_url,

                "file_url":
                    urljoin(
                        root_url,
                        filename,
                    ),

                "filename":
                    filename,

                "file_representation":
                    "ORIGINAL",

                "compression":
                    "GZIP",

                "harmonised":
                    False,

                "candidate_files_discovered":
                    len(
                        root_tsv_gz
                    ),

                "formatted_files_available":
                    0,

                "selected_for_scan":
                    True,

                "scan_success":
                    False,

                "eof_reached":
                    False,

                "scan_error":
                    None,

                "source_file_sha256":
                    None,

                "candidate_hits":
                    0,

                "file_status":
                    statuses[
                        "root_compressed"
                    ],
            }

        # ------------------------------------------------------------------
        # Priority 3: root *.txt.gz
        # ------------------------------------------------------------------

        root_txt_gz = sorted(
            filename
            for filename in root_files
            if filename.lower().endswith(
                ".txt.gz"
            )
        )

        if root_txt_gz:

            filename = root_txt_gz[0]

            return {
                "accession_id": accession,
                "source_directory_url": root_url,
                "file_url": urljoin(root_url, filename),
                "filename": filename,
                "file_representation": "ORIGINAL",
                "compression": "GZIP",
                "harmonised": False,
                "candidate_files_discovered": len(root_txt_gz),
                "formatted_files_available": 0,
                "selected_for_scan": True,
                "scan_success": False,
                "eof_reached": False,
                "scan_error": None,
                "source_file_sha256": None,
                "candidate_hits": 0,
                "file_status": statuses["root_compressed"],
            }

        # ------------------------------------------------------------------
        # Priority 4: root *.tsv
        # ------------------------------------------------------------------

        root_tsv = sorted(
            filename
            for filename in root_files
            if filename.lower().endswith(
                ".tsv"
            )
            and
            not filename.lower().endswith(
                ".tsv.gz"
            )
        )

        if root_tsv:

            filename = root_tsv[0]

            return {
                "accession_id": accession,
                "source_directory_url": root_url,
                "file_url": urljoin(root_url, filename),
                "filename": filename,
                "file_representation": "ORIGINAL",
                "compression": "NONE",
                "harmonised": False,
                "candidate_files_discovered": len(root_tsv),
                "formatted_files_available": 0,
                "selected_for_scan": True,
                "scan_success": False,
                "eof_reached": False,
                "scan_error": None,
                "source_file_sha256": None,
                "candidate_hits": 0,
                "file_status": statuses["root_uncompressed"],
            }

        # ------------------------------------------------------------------
        # Priority 5: root *.txt
        # ------------------------------------------------------------------

        root_txt = sorted(
            filename
            for filename in root_files
            if filename.lower().endswith(
                ".txt"
            )
            and
            "readme"
            not in
            filename.lower()
        )

        if root_txt:

            filename = root_txt[0]

            return {
                "accession_id": accession,
                "source_directory_url": root_url,
                "file_url": urljoin(root_url, filename),
                "filename": filename,
                "file_representation": "ORIGINAL",
                "compression": "NONE",
                "harmonised": False,
                "candidate_files_discovered": len(root_txt),
                "formatted_files_available": 0,
                "selected_for_scan": True,
                "scan_success": False,
                "eof_reached": False,
                "scan_error": None,
                "source_file_sha256": None,
                "candidate_hits": 0,
                "file_status": statuses["root_uncompressed"],
            }

        return {
            "accession_id": accession,
            "source_directory_url": root_url,
            "file_url": None,
            "filename": None,
            "file_representation": None,
            "compression": None,
            "harmonised": False,
            "candidate_files_discovered": 0,
            "formatted_files_available": 0,
            "selected_for_scan": False,
            "scan_success": False,
            "eof_reached": False,
            "scan_error": "SUMMARY_STATISTICS_FILE_LAYOUT_UNRESOLVED",
            "source_file_sha256": None,
            "candidate_hits": 0,
            "file_status": statuses["unavailable"],
        }

    # ======================================================================
    # Column resolution
    # ======================================================================

    @staticmethod
    def _resolve_column(
        *,
        header: list[str],
        aliases: list[str],
    ) -> str | None:

        exact = {
            value:
                value
            for value in header
        }

        lower = {
            value.lower():
                value
            for value in header
        }

        for alias in aliases:

            if alias in exact:

                return exact[alias]

            if alias.lower() in lower:

                return lower[
                    alias.lower()
                ]

        return None

    # ======================================================================
    # Exact-rsID streaming scan
    # ======================================================================

    def _scan_file(
        self,
        *,
        accession: str,
        file_url: str,
        compression: str,
        config: dict[str, Any],
    ) -> tuple[
        list[dict[str, Any]],
        str,
        bool,
    ]:

        targets = {
            str(
                row["rsid"]
            )
            for row in config[
                "candidates"
            ]
        }

        aliases = config[
            "summary_statistics"
        ][
            "column_aliases"
        ]

        source = config[
            "summary_statistics"
        ][
            "public_collection"
        ]

        request = Request(
            file_url,
            headers={
                "User-Agent":
                    str(
                        source[
                            "user_agent"
                        ]
                    ),

                "Accept":
                    "*/*",
            },
        )

        hasher = hashlib.sha256()

        hits: list[dict[str, Any]] = []

        eof_reached = False

        with urlopen(
            request,
            timeout=int(
                source[
                    "timeout_seconds"
                ]
            ),
        ) as response:

            hashing_reader = HashingReader(
                response,
                hasher,
            )

            buffered = io.BufferedReader(
                hashing_reader,
                buffer_size=1024 * 1024,
            )

            if compression == "GZIP":

                binary_stream: Any = gzip.GzipFile(
                    fileobj=buffered,
                    mode="rb",
                )

            elif compression == "NONE":

                binary_stream = buffered

            else:

                raise RuntimeError(
                    f"Unsupported compression: {compression}"
                )

            with binary_stream:

                with io.TextIOWrapper(
                    binary_stream,
                    encoding="utf-8",
                    errors="replace",
                    newline="",
                ) as text:

                    reader = csv.DictReader(
                        text,
                        delimiter="\t",
                    )

                    if reader.fieldnames is None:

                        raise RuntimeError(
                            "Summary-statistics file has no header."
                        )

                    resolved = {
                        standard_name:
                            self._resolve_column(
                                header=reader.fieldnames,
                                aliases=alias_values,
                            )
                        for (
                            standard_name,
                            alias_values,
                        )
                        in aliases.items()
                    }

                    rs_column = resolved.get(
                        "rs_id"
                    )

                    if rs_column is None:

                        raise RuntimeError(
                            "No defensible rsID column found."
                        )

                    for row in reader:

                        observed_rsid = str(
                            row.get(
                                rs_column,
                                "",
                            )
                        ).strip()

                        if observed_rsid not in targets:

                            continue

                        hit: dict[str, Any] = {
                            "accession_id":
                                accession,

                            "rs_id":
                                observed_rsid,

                            "source_file_url":
                                file_url,

                            "exact_rsid_match":
                                True,
                        }

                        for field in (
                            "chromosome",
                            "base_pair_location",
                            "effect_allele",
                            "other_allele",
                            "effect_allele_frequency",
                            "beta",
                            "odds_ratio",
                            "standard_error",
                            "p_value",
                            "neg_log_10_p_value",
                        ):

                            source_column = resolved.get(
                                field
                            )

                            hit[field] = (
                                row.get(
                                    source_column
                                )
                                if source_column
                                is not None
                                else None
                            )

                        hits.append(hit)

                    eof_reached = True

        if not eof_reached:

            raise RuntimeError(
                "Source file did not reach EOF."
            )

        return (
            hits,
            hasher.hexdigest(),
            eof_reached,
        )

    # ======================================================================
    # Audit files
    # ======================================================================

    def _audit_files(
        self,
        *,
        study_inventory: pd.DataFrame,
        config: dict[str, Any],
    ) -> tuple[pd.DataFrame, pd.DataFrame]:

        manifests: list[dict[str, Any]] = []
        hits: list[dict[str, Any]] = []

        eligible = study_inventory.loc[
            study_inventory[
                "eligible_for_exact_lookup"
            ]
            .fillna(False)
            .astype(bool)
        ]

        logger.info(
            "Eligible studies=%d.",
            len(eligible),
        )

        for index, (_, study) in enumerate(
            eligible.iterrows(),
            start=1,
        ):

            accession = str(
                study[
                    "accession_id"
                ]
            )

            source_directory = study.get(
                "full_summary_stats"
            )

            logger.info(
                "Study %d/%d | %s.",
                index,
                len(eligible),
                accession,
            )

            if (
                source_directory is None
                or
                str(
                    source_directory
                ).strip().upper()
                in {
                    "",
                    "NA",
                    "NONE",
                }
            ):

                manifest = {
                    "accession_id": accession,
                    "source_directory_url": source_directory,
                    "file_url": None,
                    "filename": None,
                    "file_representation": None,
                    "compression": None,
                    "harmonised": False,
                    "candidate_files_discovered": 0,
                    "formatted_files_available": 0,
                    "selected_for_scan": False,
                    "scan_success": False,
                    "eof_reached": False,
                    "scan_error": "FULL_SUMMARY_STATS_ROUTE_MISSING",
                    "source_file_sha256": None,
                    "candidate_hits": 0,
                    "file_status":
                        config[
                            "statuses"
                        ][
                            "file"
                        ][
                            "unavailable"
                        ],
                }

                manifests.append(
                    manifest
                )

                continue

            manifest = self._discover_file(
                accession=accession,
                source_directory_url=str(
                    source_directory
                ),
                config=config,
            )

            if not manifest[
                "selected_for_scan"
            ]:

                manifests.append(
                    manifest
                )

                continue

            try:

                (
                    file_hits,
                    digest,
                    eof_reached,
                ) = self._scan_file(
                    accession=accession,
                    file_url=str(
                        manifest[
                            "file_url"
                        ]
                    ),
                    compression=str(
                        manifest[
                            "compression"
                        ]
                    ),
                    config=config,
                )

                manifest[
                    "scan_success"
                ] = True

                manifest[
                    "eof_reached"
                ] = eof_reached

                manifest[
                    "scan_error"
                ] = None

                manifest[
                    "source_file_sha256"
                ] = digest

                manifest[
                    "candidate_hits"
                ] = len(
                    file_hits
                )

                manifest[
                    "file_status"
                ] = config[
                    "statuses"
                ][
                    "file"
                ][
                    "scan_success"
                ]

                hits.extend(
                    file_hits
                )

                logger.info(
                    "%s | hits=%d | EOF=%s | sha256=%s.",
                    accession,
                    len(file_hits),
                    eof_reached,
                    digest,
                )

            except Exception as exc:

                manifest[
                    "scan_success"
                ] = False

                manifest[
                    "eof_reached"
                ] = False

                manifest[
                    "scan_error"
                ] = (
                    f"{type(exc).__name__}: {exc}"
                )

                manifest[
                    "file_status"
                ] = config[
                    "statuses"
                ][
                    "file"
                ][
                    "scan_failed"
                ]

                logger.warning(
                    "%s scan failed | %s.",
                    accession,
                    exc,
                )

            manifests.append(
                manifest
            )

        return (
            pd.DataFrame(
                manifests,
                columns=FILE_COLUMNS,
            ),
            pd.DataFrame(
                hits,
                columns=HIT_COLUMNS,
            ),
        )

    # ======================================================================
    # Persist
    # ======================================================================

    def _persist(
        self,
        *,
        study_inventory: pd.DataFrame,
        file_manifest: pd.DataFrame,
        exact_hits: pd.DataFrame,
        result: Any,
        config: dict[str, Any],
    ) -> dict[str, str]:

        outputs = config[
            "outputs"
        ]

        directory = (
            self.project_root
            /
            outputs[
                "processed_directory"
            ]
        )

        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        paths = {
            "study_inventory":
                directory
                /
                outputs[
                    "study_inventory_filename"
                ],

            "file_manifest":
                directory
                /
                outputs[
                    "file_manifest_filename"
                ],

            "exact_variant_hits":
                directory
                /
                outputs[
                    "exact_variant_hits_filename"
                ],

            "study_candidate_audit":
                directory
                /
                outputs[
                    "study_candidate_audit_filename"
                ],

            "candidate_resolution":
                directory
                /
                outputs[
                    "candidate_resolution_filename"
                ],

            "summary":
                directory
                /
                outputs[
                    "summary_filename"
                ],
        }

        write_parquet(
            study_inventory,
            paths["study_inventory"],
            index=False,
        )

        write_parquet(
            file_manifest,
            paths["file_manifest"],
            index=False,
        )

        write_parquet(
            exact_hits,
            paths["exact_variant_hits"],
            index=False,
        )

        write_parquet(
            result.study_candidate_audit,
            paths["study_candidate_audit"],
            index=False,
        )

        write_parquet(
            result.candidate_resolution,
            paths["candidate_resolution"],
            index=False,
        )

        write_parquet(
            result.summary,
            paths["summary"],
            index=False,
        )

        return {
            key:
                str(value)
            for key, value in paths.items()
        }

    # ======================================================================
    # Main
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:

        config = self._load_config()

        logger.info(
            "Starting M7.5B."
        )

        upstream = self._validate_upstream(
            config=config,
        )

        studies = self._discover_studies(
            config=config,
        )

        study_inventory = self._build_study_inventory(
            studies=studies,
            config=config,
        )

        file_manifest, exact_hits = self._audit_files(
            study_inventory=study_inventory,
            config=config,
        )

        result = assess_full_summary_statistics_audit(
            study_inventory=study_inventory,
            file_manifest=file_manifest,
            exact_variant_hits=exact_hits,
            config=config,
        )

        paths = self._persist(
            study_inventory=study_inventory,
            file_manifest=file_manifest,
            exact_hits=exact_hits,
            result=result,
            config=config,
        )

        report = {
            "milestone":
                "M7.5B",

            "stage":
                "full_summary_statistics_audit",

            "upstream_validation":
                upstream,

            "summary":
                _records(
                    result.summary
                )[0],

            "study_inventory":
                _records(
                    study_inventory
                ),

            "file_manifest":
                _records(
                    file_manifest
                ),

            "exact_variant_hits":
                _records(
                    exact_hits
                ),

            "study_candidate_audit":
                _records(
                    result.study_candidate_audit
                ),

            "candidate_resolution":
                _records(
                    result.candidate_resolution
                ),

            "scientific_policy":
                config[
                    "scientific_policy"
                ],

            "outputs":
                paths,
        }

        qc_path = (
            self.project_root
            /
            config[
                "outputs"
            ][
                "qc_relative_path"
            ]
        )

        qc_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        qc_path.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
                allow_nan=False,
            ),
            encoding="utf-8",
        )

        logger.info(
            "M7.5B complete | overall=%s.",
            report[
                "summary"
            ][
                "overall_status"
            ],
        )

        return report