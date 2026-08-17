"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/assess_variant_validation_feasibility.py

Description:
    Runner for M7.5 Variant-Level Validation Feasibility.

    The runner queries:

        GET /v2/associations

    using:

        rs_id=<candidate>
        page=<zero-based page>
        size=<configured page size>

    The runner follows page.number and page.total_pages deterministically and
    merges all returned content records into one immutable candidate snapshot.

    Service failures are recorded as acquisition limitations and are never
    interpreted as negative genetic evidence.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd
import yaml

from pcatrfqtl.analysis.m7.variant_validation_feasibility import (
    assess_variant_validation_feasibility,
    extract_v2_content,
    extract_v2_page,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


# ============================================================================
# HTTP result
# ============================================================================


@dataclass(frozen=True)
class HttpResult:
    """One REST API V2 HTTP request result."""

    success: bool
    status_code: int | None
    payload: dict[str, Any] | None
    error_type: str | None
    error_message: str | None


# ============================================================================
# Utility helpers
# ============================================================================


def _sha256_bytes(
    payload: bytes,
) -> str:
    """Return SHA-256 digest."""

    return hashlib.sha256(
        payload
    ).hexdigest()


def _json_safe_row(
    row: pd.Series,
) -> dict[str, Any]:
    """Convert a pandas row to JSON-safe values."""

    output: dict[str, Any] = {}

    for key, value in row.to_dict().items():

        if value is None:

            output[
                str(
                    key
                )
            ] = None

            continue

        try:

            if pd.isna(
                value
            ):

                output[
                    str(
                        key
                    )
                ] = None

                continue

        except (
            TypeError,
            ValueError,
        ):
            pass

        if hasattr(
            value,
            "item",
        ):

            value = value.item()

        output[
            str(
                key
            )
        ] = value

    return output


def _records(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Convert dataframe to JSON-safe records."""

    return [
        _json_safe_row(
            row
        )
        for _, row in dataframe.iterrows()
    ]


# ============================================================================
# Runner
# ============================================================================


class M75VariantValidationFeasibilityRunner:
    """Execute M7.5."""

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
        """Load and validate configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M7.5 config not found: {self.config_path}"
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

            raise RuntimeError(
                "M7.5 config root must be a mapping."
            )

        if config.get(
            "milestone"
        ) != "M7.5":

            raise RuntimeError(
                "Unexpected M7.5 milestone."
            )

        if config.get(
            "stage"
        ) != "variant_level_validation_feasibility":

            raise RuntimeError(
                "Unexpected M7.5 stage."
            )

        candidates = config.get(
            "candidates"
        )

        if not isinstance(
            candidates,
            list,
        ):

            raise RuntimeError(
                "M7.5 candidates must be a list."
            )

        expected_count = int(
            config[
                "validation"
            ][
                "expected_candidate_count"
            ]
        )

        if len(
            candidates
        ) != expected_count:

            raise RuntimeError(
                "M7.5 candidate count mismatch."
            )

        return config

    # ======================================================================
    # Locked upstream
    # ======================================================================

    def _load_qc(
        self,
        *,
        relative_path: str,
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """Load one upstream QC file."""

        path = (
            self.project_root
            /
            relative_path
        )

        if not path.exists():

            raise FileNotFoundError(
                f"Upstream QC not found: {path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            payload = json.load(
                handle
            )

        return (
            payload,
            path,
        )

    def _validate_required_state(
        self,
        *,
        label: str,
        summary: dict[str, Any],
        expected: dict[str, Any],
        path: Path,
    ) -> None:
        """Validate locked summary fields."""

        for field, expected_value in expected.items():

            if field not in summary:

                raise RuntimeError(
                    f"{label} missing locked field {field!r}.\n"
                    f"Path: {path}"
                )

            observed = summary[
                field
            ]

            if observed != expected_value:

                raise RuntimeError(
                    f"{label} locked state mismatch.\n"
                    f"Field: {field!r}\n"
                    f"Observed: {observed!r}\n"
                    f"Expected: {expected_value!r}"
                )

    def _validate_upstream_item(
        self,
        *,
        label: str,
        specification: dict[str, Any],
    ) -> dict[str, Any]:
        """Validate one locked upstream artifact."""

        payload, path = self._load_qc(
            relative_path=str(
                specification[
                    "relative_path"
                ]
            )
        )

        if payload.get(
            "milestone"
        ) != specification[
            "expected_milestone"
        ]:

            raise RuntimeError(
                f"{label} milestone mismatch."
            )

        if payload.get(
            "stage"
        ) != specification[
            "expected_stage"
        ]:

            raise RuntimeError(
                f"{label} stage mismatch."
            )

        summary = payload.get(
            "summary"
        )

        if not isinstance(
            summary,
            dict,
        ):

            raise RuntimeError(
                f"{label} summary must be a mapping."
            )

        status_field = str(
            specification[
                "status_field"
            ]
        )

        observed_status = summary.get(
            status_field
        )

        if observed_status != specification[
            "expected_status"
        ]:

            raise RuntimeError(
                f"{label} status mismatch.\n"
                f"Field: {status_field}\n"
                f"Observed: {observed_status!r}\n"
                f"Expected: {specification['expected_status']!r}"
            )

        self._validate_required_state(
            label=label,
            summary=summary,
            expected=specification[
                "required_summary_state"
            ],
            path=path,
        )

        logger.info(
            "Validated %s | stage=%s | %s=%s.",
            label,
            payload[
                "stage"
            ],
            status_field,
            observed_status,
        )

        return {
            "path":
                str(
                    path
                ),

            "milestone":
                payload[
                    "milestone"
                ],

            "stage":
                payload[
                    "stage"
                ],

            "status_field":
                status_field,

            "status":
                observed_status,
        }

    def _validate_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> dict[str, Any]:
        """Validate M7.3C and M7.4D."""

        return {
            "m7_3c":
                self._validate_upstream_item(
                    label="M7.3C",
                    specification=config[
                        "upstream_qc"
                    ][
                        "controlled_access"
                    ],
                ),

            "m7_4d":
                self._validate_upstream_item(
                    label="M7.4D",
                    specification=config[
                        "upstream_qc"
                    ][
                        "external_clinical"
                    ],
                ),
        }

    # ======================================================================
    # HTTP
    # ======================================================================

    def _request_json(
        self,
        *,
        url: str,
        config: dict[str, Any],
    ) -> HttpResult:
        """Execute bounded-retry REST API request."""

        api = config[
            "gwas_catalog"
        ][
            "rest_api"
        ]

        last_status: int | None = None
        last_type: str | None = None
        last_message: str | None = None

        max_retries = int(
            api[
                "max_retries"
            ]
        )

        for attempt in range(
            1,
            max_retries + 1,
        ):

            request = Request(
                url,
                headers={
                    "Accept":
                        "application/json",

                    "User-Agent":
                        str(
                            api[
                                "user_agent"
                            ]
                        ),
                },
            )

            try:

                with urlopen(
                    request,
                    timeout=int(
                        api[
                            "timeout_seconds"
                        ]
                    ),
                ) as response:

                    raw = response.read()

                    status_code = int(
                        response.status
                    )

                if not raw:

                    last_status = status_code
                    last_type = "EMPTY_RESPONSE"
                    last_message = "API returned an empty response."

                else:

                    payload = json.loads(
                        raw.decode(
                            "utf-8"
                        )
                    )

                    if not isinstance(
                        payload,
                        dict,
                    ):

                        raise RuntimeError(
                            "REST API V2 response root is not a mapping."
                        )

                    return HttpResult(
                        success=True,
                        status_code=status_code,
                        payload=payload,
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

                if (
                    500
                    <=
                    exc.code
                    <=
                    599
                    and
                    attempt
                    <
                    max_retries
                ):

                    logger.warning(
                        "GWAS V2 server error | attempt=%d/%d | status=%d.",
                        attempt,
                        max_retries,
                        exc.code,
                    )

                    time.sleep(
                        int(
                            api[
                                "retry_backoff_seconds"
                            ]
                        )
                        *
                        attempt
                    )

                    continue

            except URLError as exc:

                last_type = "CONNECTION_ERROR"

                last_message = str(
                    exc.reason
                )

                if attempt < max_retries:

                    time.sleep(
                        int(
                            api[
                                "retry_backoff_seconds"
                            ]
                        )
                        *
                        attempt
                    )

                    continue

            except (
                UnicodeDecodeError,
                json.JSONDecodeError,
            ) as exc:

                last_type = "INVALID_JSON"
                last_message = str(
                    exc
                )

            break

        return HttpResult(
            success=False,
            status_code=last_status,
            payload=None,
            error_type=last_type,
            error_message=last_message,
        )

    # ======================================================================
    # URL
    # ======================================================================

    def _association_url(
        self,
        *,
        rsid: str,
        page: int,
        config: dict[str, Any],
    ) -> str:
        """Build exact V2 association query URL."""

        api = config[
            "gwas_catalog"
        ][
            "rest_api"
        ]

        query = api[
            "query"
        ]

        base = (
            str(
                api[
                    "base_url"
                ]
            ).rstrip(
                "/"
            )
            +
            str(
                api[
                    "endpoint"
                ]
            )
        )

        parameters = {
            str(
                query[
                    "exact_variant_parameter"
                ]
            ):
                rsid,

            str(
                query[
                    "page_parameter"
                ]
            ):
                int(
                    page
                ),

            str(
                query[
                    "size_parameter"
                ]
            ):
                int(
                    query[
                        "page_size"
                    ]
                ),
        }

        return (
            base
            +
            "?"
            +
            urlencode(
                parameters
            )
        )

    # ======================================================================
    # Snapshot
    # ======================================================================

    def _snapshot_paths(
        self,
        *,
        rsid: str,
        config: dict[str, Any],
    ) -> tuple[
        Path,
        Path,
    ]:
        """Resolve immutable raw snapshot paths."""

        directory = (
            self.project_root
            /
            config[
                "gwas_catalog"
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
            f"{rsid}.associations.v2.json",

            directory
            /
            f"{rsid}.associations.v2.metadata.json",
        )

    def _load_existing_snapshot(
        self,
        *,
        rsid: str,
        config: dict[str, Any],
    ) -> dict[str, Any] | None:
        """Reuse a previously stored successful immutable snapshot."""

        raw_path, metadata_path = self._snapshot_paths(
            rsid=rsid,
            config=config,
        )

        if not raw_path.exists():

            return None

        raw_bytes = raw_path.read_bytes()

        payload = json.loads(
            raw_bytes.decode(
                "utf-8"
            )
        )

        if metadata_path.exists():

            metadata = json.loads(
                metadata_path.read_text(
                    encoding="utf-8"
                )
            )

            if metadata.get(
                "sha256"
            ) != _sha256_bytes(
                raw_bytes
            ):

                raise RuntimeError(
                    f"Snapshot checksum mismatch: {raw_path}"
                )

        return payload

    def _persist_snapshot(
        self,
        *,
        rsid: str,
        payload: dict[str, Any],
        config: dict[str, Any],
    ) -> None:
        """Persist one successful immutable combined candidate snapshot."""

        raw_path, metadata_path = self._snapshot_paths(
            rsid=rsid,
            config=config,
        )

        if raw_path.exists():

            raise RuntimeError(
                f"Refusing to overwrite immutable snapshot: {raw_path}"
            )

        raw_bytes = json.dumps(
            payload,
            indent=2,
            ensure_ascii=False,
            sort_keys=True,
        ).encode(
            "utf-8"
        )

        digest = _sha256_bytes(
            raw_bytes
        )

        raw_path.write_bytes(
            raw_bytes
        )

        metadata = {
            "milestone":
                "M7.5",

            "source":
                "NHGRI-EBI GWAS Catalog REST API V2",

            "endpoint":
                "/v2/associations",

            "query_rsid":
                rsid,

            "retrieved_at_utc":
                datetime.now(
                    timezone.utc
                ).isoformat(),

            "sha256":
                digest,

            "immutable":
                True,

            "zero_result_interpreted_as_negative_evidence":
                False,

            "zero_result_interpreted_as_variant_absence":
                False,
        }

        metadata_path.write_text(
            json.dumps(
                metadata,
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )

    # ======================================================================
    # Candidate query
    # ======================================================================

    def _query_candidate(
        self,
        *,
        rsid: str,
        config: dict[str, Any],
    ) -> dict[str, Any]:
        """Query all V2 pages for one candidate."""

        existing = self._load_existing_snapshot(
            rsid=rsid,
            config=config,
        )

        if existing is not None:

            page = extract_v2_page(
                existing
            )

            records = extract_v2_content(
                existing
            )

            logger.info(
                "Reusing immutable V2 snapshot | %s | records=%d.",
                rsid,
                len(
                    records
                ),
            )

            return {
                "request_success":
                    True,

                "http_status":
                    200,

                "payload":
                    existing,

                "pages_retrieved":
                    max(
                        int(
                            page[
                                "total_pages"
                            ]
                        ),
                        1,
                    ),

                "records_retrieved":
                    len(
                        records
                    ),

                "snapshot_reused":
                    True,

                "error_type":
                    None,

                "error_message":
                    None,
            }

        page_number = int(
            config[
                "gwas_catalog"
            ][
                "rest_api"
            ][
                "query"
            ][
                "start_page"
            ]
        )

        combined_content: list[
            dict[str, Any]
        ] = []

        pages_retrieved = 0

        total_elements: int | None = None
        total_pages: int | None = None
        page_size: int | None = None

        first_links: list[Any] = []

        while True:

            url = self._association_url(
                rsid=rsid,
                page=page_number,
                config=config,
            )

            logger.info(
                "GWAS Catalog V2 | %s | page=%d.",
                rsid,
                page_number,
            )

            response = self._request_json(
                url=url,
                config=config,
            )

            if not response.success:

                return {
                    "request_success":
                        False,

                    "http_status":
                        response.status_code,

                    "payload":
                        None,

                    "pages_retrieved":
                        pages_retrieved,

                    "records_retrieved":
                        0,

                    "snapshot_reused":
                        False,

                    "error_type":
                        response.error_type,

                    "error_message":
                        response.error_message,

                    "acquisition_status":
                        config[
                            "statuses"
                        ][
                            "acquisition"
                        ][
                            "service_unavailable"
                        ],
                }

            payload = response.payload or {}

            try:

                page = extract_v2_page(
                    payload
                )

                content = extract_v2_content(
                    payload
                )

            except RuntimeError as exc:

                return {
                    "request_success":
                        False,

                    "http_status":
                        response.status_code,

                    "payload":
                        None,

                    "pages_retrieved":
                        pages_retrieved,

                    "records_retrieved":
                        0,

                    "snapshot_reused":
                        False,

                    "error_type":
                        "INVALID_RESPONSE_SCHEMA",

                    "error_message":
                        str(
                            exc
                        ),

                    "acquisition_status":
                        config[
                            "statuses"
                        ][
                            "acquisition"
                        ][
                            "invalid_response"
                        ],
                }

            if pages_retrieved == 0:

                links = payload.get(
                    "links",
                    [],
                )

                if isinstance(
                    links,
                    list,
                ):

                    first_links = links

                total_elements = int(
                    page[
                        "total_elements"
                    ]
                )

                total_pages = int(
                    page[
                        "total_pages"
                    ]
                )

                page_size = int(
                    page[
                        "size"
                    ]
                )

            else:

                if int(
                    page[
                        "total_elements"
                    ]
                ) != total_elements:

                    raise RuntimeError(
                        "GWAS V2 total_elements changed during pagination."
                    )

                if int(
                    page[
                        "total_pages"
                    ]
                ) != total_pages:

                    raise RuntimeError(
                        "GWAS V2 total_pages changed during pagination."
                    )

            combined_content.extend(
                content
            )

            pages_retrieved += 1

            if total_pages is None:

                raise RuntimeError(
                    "GWAS V2 total_pages was not initialized."
                )

            # Successful zero result.
            if total_pages == 0:

                break

            if page_number + 1 >= total_pages:

                break

            page_number += 1

            # Conservative request spacing.
            queries_per_second = float(
                config[
                    "gwas_catalog"
                ][
                    "rest_api"
                ][
                    "client_queries_per_second"
                ]
            )

            if queries_per_second > 0:

                time.sleep(
                    1.0
                    /
                    queries_per_second
                )

        if total_elements is None:

            total_elements = 0

        if total_pages is None:

            total_pages = 0

        if page_size is None:

            page_size = int(
                config[
                    "gwas_catalog"
                ][
                    "rest_api"
                ][
                    "query"
                ][
                    "page_size"
                ]
            )

        if len(
            combined_content
        ) != total_elements:

            raise RuntimeError(
                "GWAS V2 retrieved record count does not match "
                f"page.total_elements for {rsid}: "
                f"{len(combined_content)} != {total_elements}"
            )

        combined_payload = {
            "links":
                first_links,

            "content":
                combined_content,

            "page": {
                "size":
                    page_size,

                "total_elements":
                    total_elements,

                "total_pages":
                    total_pages,

                "number":
                    0,
            },

            "_pcatrfqtl_metadata": {
                "query_rsid":
                    rsid,

                "pages_retrieved":
                    pages_retrieved,

                "records_retrieved":
                    len(
                        combined_content
                    ),

                "pagination_collapsed":
                    True,

                "source_api":
                    "GWAS_CATALOG_REST_API_V2",
            },
        }

        self._persist_snapshot(
            rsid=rsid,
            payload=combined_payload,
            config=config,
        )

        return {
            "request_success":
                True,

            "http_status":
                200,

            "payload":
                combined_payload,

            "pages_retrieved":
                pages_retrieved,

            "records_retrieved":
                len(
                    combined_content
                ),

            "snapshot_reused":
                False,

            "error_type":
                None,

            "error_message":
                None,
        }

    # ======================================================================
    # Acquisition
    # ======================================================================

    def _acquire_candidates(
        self,
        *,
        config: dict[str, Any],
    ) -> dict[str, dict[str, Any]]:
        """Query all configured candidate rsIDs."""

        results: dict[
            str,
            dict[str, Any]
        ] = {}

        candidates = config[
            "candidates"
        ]

        for index, candidate in enumerate(
            candidates,
            start=1,
        ):

            rsid = str(
                candidate[
                    "rsid"
                ]
            )

            logger.info(
                "M7.5 REST API V2 lookup %d/%d | %s.",
                index,
                len(
                    candidates
                ),
                rsid,
            )

            result = self._query_candidate(
                rsid=rsid,
                config=config,
            )

            results[
                rsid
            ] = result

            logger.info(
                "%s | success=%s | records=%d.",
                rsid,
                result[
                    "request_success"
                ],
                result[
                    "records_retrieved"
                ],
            )

        return results

    # ======================================================================
    # Persist processed outputs
    # ======================================================================

    def _persist_outputs(
        self,
        *,
        result: Any,
        config: dict[str, Any],
    ) -> dict[str, str]:
        """Write M7.5 processed parquet artifacts."""

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
            "acquisition_status":
                directory
                /
                outputs[
                    "acquisition_status_filename"
                ],

            "association_inventory":
                directory
                /
                outputs[
                    "association_inventory_filename"
                ],

            "candidate_feasibility":
                directory
                /
                outputs[
                    "candidate_feasibility_filename"
                ],

            "route_inventory":
                directory
                /
                outputs[
                    "route_inventory_filename"
                ],

            "summary":
                directory
                /
                outputs[
                    "summary_filename"
                ],
        }

        write_parquet(
            result.acquisition_status,
            paths[
                "acquisition_status"
            ],
            index=False,
        )

        write_parquet(
            result.association_inventory,
            paths[
                "association_inventory"
            ],
            index=False,
        )

        write_parquet(
            result.candidate_feasibility,
            paths[
                "candidate_feasibility"
            ],
            index=False,
        )

        write_parquet(
            result.route_inventory,
            paths[
                "route_inventory"
            ],
            index=False,
        )

        write_parquet(
            result.summary,
            paths[
                "summary"
            ],
            index=False,
        )

        return {
            key:
                str(
                    value
                )
            for key, value in paths.items()
        }

    # ======================================================================
    # QC
    # ======================================================================

    def _write_qc(
        self,
        *,
        result: Any,
        config: dict[str, Any],
        upstream_validation: dict[str, Any],
        output_paths: dict[str, str],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """Write M7.5 QC report."""

        report = {
            "milestone":
                "M7.5",

            "stage":
                "variant_level_validation_feasibility",

            "api":
                {
                    "provider":
                        config[
                            "gwas_catalog"
                        ][
                            "provider"
                        ],

                    "version":
                        config[
                            "gwas_catalog"
                        ][
                            "api_version"
                        ],

                    "endpoint":
                        "/v2/associations",

                    "exact_variant_parameter":
                        "rs_id",
                },

            "upstream_validation":
                upstream_validation,

            "summary":
                _json_safe_row(
                    result.summary.iloc[
                        0
                    ]
                ),

            "acquisition_status":
                _records(
                    result.acquisition_status
                ),

            "association_inventory":
                _records(
                    result.association_inventory
                ),

            "candidate_feasibility":
                _records(
                    result.candidate_feasibility
                ),

            "route_inventory":
                _records(
                    result.route_inventory
                ),

            "classification_policy":
                config[
                    "classification"
                ],

            "scientific_policy":
                config[
                    "scientific_policy"
                ],

            "outputs":
                output_paths,
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

        with qc_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
                allow_nan=False,
            )

        return (
            report,
            qc_path,
        )

    # ======================================================================
    # Main
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M7.5."""

        config = self._load_config()

        logger.info(
            "Starting M7.5 Variant-Level Validation Feasibility."
        )

        upstream_validation = self._validate_upstream(
            config=config,
        )

        logger.info(
            "Locked upstream validation passed."
        )

        api_results = self._acquire_candidates(
            config=config,
        )

        result = assess_variant_validation_feasibility(
            api_results=api_results,
            config=config,
        )

        output_paths = self._persist_outputs(
            result=result,
            config=config,
        )

        report, qc_path = self._write_qc(
            result=result,
            config=config,
            upstream_validation=upstream_validation,
            output_paths=output_paths,
        )

        summary = report[
            "summary"
        ]

        logger.info(
            "M7.5 complete."
        )

        logger.info(
            "Candidates=%d | successful queries=%d | zero-result=%d.",
            summary[
                "candidates_assessed"
            ],
            summary[
                "api_queries_successful"
            ],
            summary[
                "api_queries_zero_curated_results"
            ],
        )

        for row in report[
            "candidate_feasibility"
        ]:

            logger.info(
                "%s | curated=%d | exact=%d | prostate=%d | "
                "reused=%d | nonreused=%d | status=%s.",
                row[
                    "rsid"
                ],
                row[
                    "curated_records_returned"
                ],
                row[
                    "exact_candidate_associations"
                ],
                row[
                    "prostate_trait_candidate_associations"
                ],
                row[
                    "reused_prostate_studies"
                ],
                row[
                    "nonreused_prostate_studies"
                ],
                row[
                    "candidate_status"
                ],
            )

        logger.info(
            "Independent replication verified=%d | claimed=%s.",
            summary[
                "independent_replications_verified"
            ],
            summary[
                "independent_replication_claimed"
            ],
        )

        logger.info(
            "Overall=%s | next=%s.",
            summary[
                "overall_status"
            ],
            summary[
                "next_stage"
            ],
        )

        logger.info(
            "QC=%s.",
            qc_path,
        )

        return report