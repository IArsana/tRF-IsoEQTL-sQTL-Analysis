"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/runners/audit_external_trf_raw_read_acquisition.py

Description:
    Runner for M7.4C.1 External tRF Raw-Read Acquisition Audit.

    This runner retrieves SRA RunInfo metadata for GSE80400/SRP073456,
    stores an immutable raw metadata snapshot, validates the run inventory,
    audits SRA Toolkit availability, and creates the downstream acquisition
    manifest.

    No sequencing reads are downloaded in this stage.

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
import io
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd
import yaml

from pcatrfqtl.analysis.m7.external_trf_raw_read_acquisition import (
    audit_external_trf_raw_read_acquisition,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


# ============================================================================
# Helpers
# ============================================================================


def _sha256_bytes(
    payload: bytes,
) -> str:
    """Return SHA-256 digest for raw metadata bytes."""

    return hashlib.sha256(
        payload
    ).hexdigest()


def _records_with_json_nulls(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Convert DataFrame records into JSON-safe structures."""

    if dataframe.empty:
        return []

    cleaned = dataframe.astype(object).where(
        pd.notna(dataframe),
        None,
    )

    records = cleaned.to_dict(
        orient="records"
    )

    for record in records:

        for key, value in list(
            record.items()
        ):

            if hasattr(
                value,
                "item",
            ):
                record[key] = value.item()

    return records


def _json_safe_row(
    row: pd.Series,
) -> dict[str, Any]:
    """Convert one pandas row to JSON-safe scalars."""

    output: dict[str, Any] = {}

    for key, value in row.to_dict().items():

        if value is None:
            output[str(key)] = None
            continue

        try:
            if pd.isna(value):
                output[str(key)] = None
                continue
        except (
            TypeError,
            ValueError,
        ):
            pass

        if hasattr(value, "item"):
            value = value.item()

        output[str(key)] = value

    return output


# ============================================================================
# Runner
# ============================================================================


class M74C1ExternalTrfRawReadAcquisitionRunner:
    """Execute M7.4C.1."""

    def __init__(
        self,
        *,
        project_root: str | Path,
        config_path: str | Path,
    ) -> None:

        self.project_root = Path(
            project_root
        )

        self.config_path = Path(
            config_path
        )

    # ======================================================================
    # Config
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:

        if not self.config_path.exists():

            raise FileNotFoundError(
                f"M7.4C.1 config not found: {self.config_path}"
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
                "M7.4C.1 YAML root must be a mapping."
            )

        if config.get(
            "milestone"
        ) != "M7.4C.1":

            raise RuntimeError(
                "Unexpected M7.4C.1 milestone."
            )

        if config.get(
            "stage"
        ) != "external_trf_raw_read_acquisition_audit":

            raise RuntimeError(
                "Unexpected M7.4C.1 stage."
            )

        return config

    # ======================================================================
    # Upstream
    # ======================================================================

    def _load_upstream(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:

        upstream = (
            config[
                "upstream_qc"
            ][
                "technical_eligibility"
            ]
        )

        path = (
            self.project_root
            /
            str(
                upstream[
                    "relative_path"
                ]
            )
        )

        if not path.exists():

            raise FileNotFoundError(
                f"M7.4B QC not found: {path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            payload = json.load(
                handle
            )

        if payload.get(
            "milestone"
        ) != upstream[
            "expected_milestone"
        ]:

            raise RuntimeError(
                "Unexpected upstream milestone."
            )

        if payload.get(
            "stage"
        ) != upstream[
            "expected_stage"
        ]:

            raise RuntimeError(
                "Unexpected upstream stage."
            )

        observed_status = (
            payload
            .get(
                "summary",
                {},
            )
            .get(
                "overall_technical_status"
            )
        )

        if observed_status != upstream[
            "expected_overall_status"
        ]:

            raise RuntimeError(
                "Unexpected upstream M7.4B status: "
                f"{observed_status!r}"
            )

        pilot_key = config[
            "pilot_dataset"
        ][
            "dataset_key"
        ]

        rows = payload.get(
            "dataset_technical_eligibility",
            [],
        )

        pilot_rows = [
            row
            for row in rows
            if row.get(
                "dataset_key"
            )
            ==
            pilot_key
        ]

        if len(
            pilot_rows
        ) != 1:

            raise RuntimeError(
                "Pilot dataset not uniquely resolved in M7.4B."
            )

        pilot = pilot_rows[
            0
        ]

        if not bool(
            pilot.get(
                "technical_eligibility",
                False,
            )
        ):

            raise RuntimeError(
                "Pilot dataset is not technically eligible."
            )

        return (
            payload,
            path,
        )

    # ======================================================================
    # SRA metadata acquisition
    # ======================================================================

    def _fetch_runinfo(
        self,
        *,
        config: dict[str, Any],
    ) -> tuple[
        bytes,
        str,
    ]:

        metadata = config[
            "sra_metadata"
        ]

        query = urlencode(
            {
                metadata[
                    "accession_parameter"
                ]:
                    metadata[
                        "query_accession"
                    ],
            }
        )

        url = (
            f"{metadata['endpoint']}?"
            f"{query}"
        )

        request = Request(
            url,
            headers={
                "User-Agent":
                    metadata[
                        "user_agent"
                    ],
            },
        )

        with urlopen(
            request,
            timeout=int(
                metadata[
                    "timeout_seconds"
                ]
            ),
        ) as response:

            payload = response.read()

        if not payload:

            raise RuntimeError(
                "NCBI SRA RunInfo returned empty content."
            )

        return (
            payload,
            url,
        )

    def _persist_raw_snapshot(
        self,
        *,
        payload: bytes,
        source_url: str,
        config: dict[str, Any],
    ) -> tuple[
        Path,
        Path,
        str,
    ]:
        """
        Store immutable SRA metadata snapshot.

        Existing raw snapshot is never silently overwritten.
        """

        metadata = config[
            "sra_metadata"
        ]

        raw_path = (
            self.project_root
            /
            metadata[
                "raw_snapshot_relative_path"
            ]
        )

        metadata_path = (
            self.project_root
            /
            metadata[
                "raw_snapshot_metadata_relative_path"
            ]
        )

        raw_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        digest = _sha256_bytes(
            payload
        )

        if raw_path.exists():

            existing = raw_path.read_bytes()

            existing_digest = _sha256_bytes(
                existing
            )

            if existing_digest != digest:

                raise RuntimeError(
                    "Immutable raw SRA metadata snapshot already exists "
                    "with different content. Do not overwrite it. "
                    "Create a versioned snapshot instead."
                )

        else:

            raw_path.write_bytes(
                payload
            )

        if not metadata_path.exists():

            snapshot_metadata = {
                "source":
                    "NCBI SRA Database Backend",

                "source_url":
                    source_url,

                "query_accession":
                    metadata[
                        "query_accession"
                    ],

                "retrieved_at_utc":
                    datetime.now(
                        timezone.utc
                    ).isoformat(),

                "sha256":
                    digest,

                "immutable":
                    True,
            }

            metadata_path.write_text(
                json.dumps(
                    snapshot_metadata,
                    indent=2,
                    ensure_ascii=False,
                ),
                encoding="utf-8",
            )

        return (
            raw_path,
            metadata_path,
            digest,
        )

    # ======================================================================
    # Toolchain
    # ======================================================================

    def _audit_toolchain(
        self,
    ) -> dict[str, bool]:

        return {
            "prefetch":
                shutil.which(
                    "prefetch"
                )
                is not None,

            "fasterq-dump":
                shutil.which(
                    "fasterq-dump"
                )
                is not None,
        }

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:

        config = self._load_config()

        (
            m7_4b_qc,
            m7_4b_path,
        ) = self._load_upstream(
            config=config,
        )

        logger.info(
            "Loaded locked M7.4B technical eligibility artifact."
        )

        payload, source_url = self._fetch_runinfo(
            config=config,
        )

        (
            raw_snapshot_path,
            raw_metadata_path,
            raw_sha256,
        ) = self._persist_raw_snapshot(
            payload=payload,
            source_url=source_url,
            config=config,
        )

        logger.info(
            "SRA RunInfo snapshot stored/verified: %s.",
            raw_snapshot_path,
        )

        runinfo = pd.read_csv(
            io.BytesIO(
                payload
            ),
            low_memory=False,
        )

        toolchain = self._audit_toolchain()

        result = audit_external_trf_raw_read_acquisition(
            runinfo=runinfo,
            toolchain=toolchain,
            config=config,
        )

        outputs = config[
            "outputs"
        ]

        processed_dir = (
            self.project_root
            /
            outputs[
                "processed_directory"
            ]
        )

        processed_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        qc_path = (
            self.project_root
            /
            outputs[
                "qc_relative_path"
            ]
        )

        qc_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        run_inventory_path = (
            processed_dir
            /
            outputs[
                "run_inventory_filename"
            ]
        )

        manifest_path = (
            processed_dir
            /
            outputs[
                "acquisition_manifest_filename"
            ]
        )

        summary_path = (
            processed_dir
            /
            outputs[
                "summary_filename"
            ]
        )

        write_parquet(
            result.run_inventory,
            run_inventory_path,
            index=False,
        )

        write_parquet(
            result.acquisition_manifest,
            manifest_path,
            index=False,
        )

        write_parquet(
            result.summary,
            summary_path,
            index=False,
        )

        summary = _json_safe_row(
            result.summary.iloc[
                0
            ]
        )

        report = {
            "milestone":
                "M7.4C.1",

            "stage":
                "external_trf_raw_read_acquisition_audit",

            "candidate":
                config[
                    "candidate"
                ],

            "pilot_dataset":
                config[
                    "pilot_dataset"
                ],

            "upstream_validation": {
                "m7_4b_loaded":
                    True,

                "m7_4b_path":
                    str(
                        m7_4b_path
                    ),

                "m7_4b_status":
                    (
                        m7_4b_qc
                        .get(
                            "summary",
                            {},
                        )
                        .get(
                            "overall_technical_status"
                        )
                    ),
            },

            "raw_metadata_snapshot": {
                "path":
                    str(
                        raw_snapshot_path
                    ),

                "metadata_path":
                    str(
                        raw_metadata_path
                    ),

                "sha256":
                    raw_sha256,

                "source_url":
                    source_url,

                "immutable":
                    True,
            },

            "toolchain":
                toolchain,

            "summary":
                summary,

            "run_inventory":
                _records_with_json_nulls(
                    result.run_inventory
                ),

            "acquisition_manifest":
                _records_with_json_nulls(
                    result.acquisition_manifest
                ),

            "scientific_policy":
                config[
                    "scientific_policy"
                ],

            "outputs": {
                "run_inventory":
                    str(
                        run_inventory_path
                    ),

                "acquisition_manifest":
                    str(
                        manifest_path
                    ),

                "summary":
                    str(
                        summary_path
                    ),
            },
        }

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

        logger.info(
            "M7.4C.1 complete."
        )

        logger.info(
            "Runs=%d | total spots=%d | total bases=%d.",
            summary[
                "runs_identified"
            ],
            summary[
                "total_spots"
            ],
            summary[
                "total_bases"
            ],
        )

        logger.info(
            "prefetch=%s | fasterq-dump=%s | toolkit_ready=%s.",
            summary[
                "prefetch_available"
            ],
            summary[
                "fasterq_dump_available"
            ],
            summary[
                "sra_toolkit_ready"
            ],
        )

        logger.info(
            "Raw download=%s | FASTQ generation=%s | "
            "sequence search=%s | exact counting=%s.",
            summary[
                "raw_sequence_download_performed"
            ],
            summary[
                "fastq_generation_performed"
            ],
            summary[
                "candidate_sequence_search_performed"
            ],
            summary[
                "exact_match_counting_performed"
            ],
        )

        logger.info(
            "Acquisition status=%s.",
            summary[
                "acquisition_status"
            ],
        )

        logger.info(
            "Overall status=%s.",
            summary[
                "overall_status"
            ],
        )

        logger.info(
            "Next stage=%s.",
            summary[
                "next_stage"
            ],
        )

        logger.info(
            "QC output=%s.",
            qc_path,
        )

        return report