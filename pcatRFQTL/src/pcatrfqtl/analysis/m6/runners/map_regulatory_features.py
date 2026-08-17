"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/runners/map_regulatory_features.py

Description:
    Runner for M6.2A Candidate → Regulatory Feature Resolution.

    Inputs:
        - locked M6.1 candidate prioritization QC;
        - Moradi Supplementary Table S1;
        - Moradi Supplementary Table S2;
        - Moradi Supplementary Table S7.

    Outputs:
        - regulatory feature evidence;
        - candidate → regulatory feature mapping;
        - M6.2A summary;
        - standards-compliant QC JSON.

    Important safeguards:
        - Raw source files are read-only.
        - Candidate-feature assignments are inherited from M6.1.
        - Missing values are serialized as JSON null, never NaN.
        - Parent genes are not inferred.
        - Nearest-gene mapping is not performed.
        - Positional gene mapping is not performed.
        - Regulatory evidence is not interpreted as causal.
        - No formal colocalization is performed.
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
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from pcatrfqtl.analysis.m6.regulatory_feature_mapping import (
    resolve_regulatory_features,
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
# JSON helpers
# ============================================================================


def _records_with_json_nulls(
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """
    Convert a DataFrame to JSON-safe records.

    Pandas and NumPy missing values are explicitly converted to Python None
    so that the final JSON artifact contains standards-compliant ``null``
    rather than non-standard ``NaN`` values.

    This transformation affects only JSON serialization and does not alter
    the analytical DataFrame or Parquet outputs.
    """

    if dataframe.empty:

        return []

    cleaned = dataframe.astype(
        object
    ).where(
        pd.notna(
            dataframe
        ),
        None,
    )

    return cleaned.to_dict(
        orient="records"
    )


def _json_safe_value(
    value: Any,
) -> Any:
    """
    Convert one scalar value to a JSON-safe Python representation.

    Missing scalar values are returned as None.
    """

    if value is None:

        return None

    try:

        if pd.isna(
            value
        ):

            return None

    except (
        TypeError,
        ValueError,
    ):

        pass

    return value


# ============================================================================
# Runner
# ============================================================================


class M62ARegulatoryFeatureMappingRunner:
    """Execute M6.2A candidate-regulatory feature mapping."""

    def __init__(
        self,
        *,
        project_root: str | Path,
        config_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.project_root = Path(
            project_root
        )

        self.config_path = Path(
            config_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    # ======================================================================
    # Outputs
    # ======================================================================

    @property
    def feature_evidence_output(
        self,
    ) -> Path:
        """Return regulatory-feature evidence output path."""

        return (
            self.output_directory
            / "regulatory_feature_evidence.parquet"
        )

    @property
    def candidate_mapping_output(
        self,
    ) -> Path:
        """Return candidate-feature mapping output path."""

        return (
            self.output_directory
            / "candidate_regulatory_feature_mapping.parquet"
        )

    @property
    def summary_output(
        self,
    ) -> Path:
        """Return M6.2A summary output path."""

        return (
            self.output_directory
            / "regulatory_feature_mapping_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:
        """Return M6.2A QC JSON output path."""

        return (
            self.qc_directory
            / "m6_2a_regulatory_feature_mapping.json"
        )

    # ======================================================================
    # Configuration
    # ======================================================================

    def _load_config(
        self,
    ) -> dict[str, Any]:
        """Load and validate M6.2A configuration."""

        if not self.config_path.exists():

            raise FileNotFoundError(
                "M6.2A configuration not found: "
                f"{self.config_path}"
            )

        if not self.config_path.is_file():

            raise RuntimeError(
                "M6.2A configuration path is not a file: "
                f"{self.config_path}"
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
                "M6.2A configuration must be a YAML mapping."
            )

        required_sections = {
            "upstream_qc",
            "moradi_sources",
            "disease_context",
            "feature_policy",
            "evidence_policy",
            "resolution_status",
            "scientific_policy",
        }

        missing_sections = (
            required_sections
            -
            set(
                config
            )
        )

        if missing_sections:

            raise ValueError(
                "M6.2A configuration missing required sections: "
                f"{sorted(missing_sections)}"
            )

        return config

    # ======================================================================
    # Generic input helpers
    # ======================================================================

    def _resolve_path(
        self,
        relative_path: str,
    ) -> Path:
        """Resolve project-relative path."""

        return (
            self.project_root
            /
            str(
                relative_path
            )
        )

    @staticmethod
    def _require_columns(
        dataframe: pd.DataFrame,
        columns: list[str],
        *,
        source_name: str,
    ) -> None:
        """Validate required source columns."""

        missing = [
            column
            for column
            in columns
            if column not in dataframe.columns
        ]

        if missing:

            raise RuntimeError(
                f"{source_name} missing required columns: "
                f"{missing}"
            )

    # ======================================================================
    # Upstream M6.1 QC
    # ======================================================================

    def _load_m6_1_qc(
        self,
        config: dict[str, Any],
    ) -> tuple[
        dict[str, Any],
        Path,
    ]:
        """Load and validate locked M6.1 QC."""

        properties = config[
            "upstream_qc"
        ][
            "candidate_prioritization"
        ]

        path = self._resolve_path(
            str(
                properties[
                    "relative_path"
                ]
            )
        )

        if not path.exists():

            raise FileNotFoundError(
                "Locked M6.1 QC not found: "
                f"{path}"
            )

        if not path.is_file():

            raise RuntimeError(
                "M6.1 QC path is not a file: "
                f"{path}"
            )

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            payload = json.load(
                handle
            )

        if not isinstance(
            payload,
            dict,
        ):

            raise ValueError(
                "M6.1 QC must contain a JSON object."
            )

        expected_milestone = str(
            properties[
                "expected_milestone"
            ]
        )

        expected_stage = str(
            properties[
                "expected_stage"
            ]
        )

        observed_milestone = str(
            payload.get(
                "milestone"
            )
        )

        observed_stage = str(
            payload.get(
                "stage"
            )
        )

        if observed_milestone != expected_milestone:

            raise RuntimeError(
                "Unexpected upstream milestone: "
                f"{observed_milestone}; "
                f"expected {expected_milestone}."
            )

        if observed_stage != expected_stage:

            raise RuntimeError(
                "Unexpected upstream M6.1 stage: "
                f"{observed_stage}; "
                f"expected {expected_stage}."
            )

        candidate_priority = payload.get(
            "candidate_priority"
        )

        if not isinstance(
            candidate_priority,
            list,
        ):

            raise RuntimeError(
                "M6.1 QC does not contain candidate_priority records."
            )

        if not candidate_priority:

            raise RuntimeError(
                "M6.1 candidate_priority is empty."
            )

        return (
            payload,
            path,
        )

    # ======================================================================
    # Moradi S1
    # ======================================================================

    def _load_s1(
        self,
        config: dict[str, Any],
    ) -> tuple[
        pd.DataFrame,
        Path,
    ]:
        """Load Moradi cis intron-retention sQTL table."""

        source = config[
            "moradi_sources"
        ][
            "s1_cis_qtl"
        ]

        path = self._resolve_path(
            str(
                source[
                    "relative_path"
                ]
            )
        )

        if not path.exists():

            raise FileNotFoundError(
                "Moradi S1 not found: "
                f"{path}"
            )

        frame = pd.read_excel(
            path,
            sheet_name=str(
                source[
                    "sheet"
                ]
            ),
        )

        self._require_columns(
            frame,
            [
                str(
                    source[
                        "feature_column"
                    ]
                ),
                str(
                    source[
                        "variant_column"
                    ]
                ),
                str(
                    source[
                        "variant_position_column"
                    ]
                ),
                str(
                    source[
                        "p_value_column"
                    ]
                ),
                str(
                    source[
                        "fdr_column"
                    ]
                ),
                str(
                    source[
                        "beta_column"
                    ]
                ),
                str(
                    source[
                        "maf_column"
                    ]
                ),
            ],
            source_name="Moradi S1",
        )

        return (
            frame,
            path,
        )

    # ======================================================================
    # Moradi S2
    # ======================================================================

    def _load_s2(
        self,
        config: dict[str, Any],
    ) -> tuple[
        pd.DataFrame,
        Path,
    ]:
        """
        Load Moradi cis intron GWAS-LD table.

        Important:
            The actual S2 column values indicate that:

                sQTL_SNP      = genomic coordinate
                sQTL_SNP-pos  = rsID

                tag_SNP       = genomic coordinate
                tag_SNP_pos   = rsID

            Configuration follows actual source values rather than the
            apparent semantics of the original column names.
        """

        source = config[
            "moradi_sources"
        ][
            "s2_gwas_ld"
        ]

        path = self._resolve_path(
            str(
                source[
                    "relative_path"
                ]
            )
        )

        if not path.exists():

            raise FileNotFoundError(
                "Moradi S2 not found: "
                f"{path}"
            )

        frame = pd.read_excel(
            path,
            sheet_name=str(
                source[
                    "sheet"
                ]
            ),
        )

        self._require_columns(
            frame,
            [
                str(
                    source[
                        "feature_column"
                    ]
                ),
                str(
                    source[
                        "qtl_variant_column"
                    ]
                ),
                str(
                    source[
                        "qtl_position_column"
                    ]
                ),
                str(
                    source[
                        "tag_variant_column"
                    ]
                ),
                str(
                    source[
                        "tag_position_column"
                    ]
                ),
                str(
                    source[
                        "ld_column"
                    ]
                ),
                str(
                    source[
                        "disease_column"
                    ]
                ),
            ],
            source_name="Moradi S2",
        )

        return (
            frame,
            path,
        )

    # ======================================================================
    # Moradi S7
    # ======================================================================

    def _load_s7(
        self,
        config: dict[str, Any],
    ) -> tuple[
        pd.DataFrame,
        Path,
    ]:
        """Load Moradi differential cis intron-retention table."""

        source = config[
            "moradi_sources"
        ][
            "s7_differential_intron"
        ]

        path = self._resolve_path(
            str(
                source[
                    "relative_path"
                ]
            )
        )

        if not path.exists():

            raise FileNotFoundError(
                "Moradi S7 not found: "
                f"{path}"
            )

        frame = pd.read_csv(
            path,
            low_memory=False,
        )

        self._require_columns(
            frame,
            [
                str(
                    source[
                        "feature_column"
                    ]
                ),
                str(
                    source[
                        "base_mean_column"
                    ]
                ),
                str(
                    source[
                        "log2_fold_change_column"
                    ]
                ),
                str(
                    source[
                        "standard_error_column"
                    ]
                ),
                str(
                    source[
                        "statistic_column"
                    ]
                ),
                str(
                    source[
                        "p_value_column"
                    ]
                ),
                str(
                    source[
                        "adjusted_p_value_column"
                    ]
                ),
            ],
            source_name="Moradi S7",
        )

        return (
            frame,
            path,
        )

    # ======================================================================
    # QC report construction
    # ======================================================================

    def _build_qc_report(
        self,
        *,
        config: dict[str, Any],
        m6_1_qc: dict[str, Any],
        m6_1_path: Path,
        s1: pd.DataFrame,
        s1_path: Path,
        s2: pd.DataFrame,
        s2_path: Path,
        s7: pd.DataFrame,
        s7_path: Path,
        result: Any,
    ) -> dict[str, Any]:
        """Build standards-compliant M6.2A QC report."""

        if result.summary.empty:

            raise RuntimeError(
                "M6.2A summary output is unexpectedly empty."
            )

        summary = result.summary.iloc[
            0
        ]

        report = {
            "milestone":
                "M6.2A",

            "stage":
                "candidate_to_regulatory_feature_resolution",

            "policy":
                dict(
                    config[
                        "scientific_policy"
                    ]
                ),

            "upstream_validation": {
                "m6_1_loaded":
                    True,

                "m6_1_path":
                    str(
                        m6_1_path
                    ),

                "m6_1_highest_priority_candidate":
                    _json_safe_value(
                        m6_1_qc.get(
                            "summary",
                            {},
                        ).get(
                            "highest_priority_candidate"
                        )
                    ),
            },

            "inputs": {
                "moradi_s1":
                    str(
                        s1_path
                    ),

                "moradi_s2":
                    str(
                        s2_path
                    ),

                "moradi_s7":
                    str(
                        s7_path
                    ),

                "moradi_s1_rows":
                    int(
                        len(
                            s1
                        )
                    ),

                "moradi_s2_rows":
                    int(
                        len(
                            s2
                        )
                    ),

                "moradi_s7_rows":
                    int(
                        len(
                            s7
                        )
                    ),
            },

            "summary": {
                "candidates_assessed":
                    int(
                        summary[
                            "candidates_assessed"
                        ]
                    ),

                "candidates_with_regulatory_feature":
                    int(
                        summary[
                            "candidates_with_regulatory_feature"
                        ]
                    ),

                "candidates_without_regulatory_feature":
                    int(
                        summary[
                            "candidates_without_regulatory_feature"
                        ]
                    ),

                "unique_regulatory_features_assessed":
                    int(
                        summary[
                            "unique_regulatory_features_assessed"
                        ]
                    ),

                "fully_confirmed_regulatory_features":
                    int(
                        summary[
                            "fully_confirmed_regulatory_features"
                        ]
                    ),

                "features_with_parent_gene_resolved":
                    int(
                        summary[
                            "features_with_parent_gene_resolved"
                        ]
                    ),

                "features_with_parent_gene_pending":
                    int(
                        summary[
                            "features_with_parent_gene_pending"
                        ]
                    ),

                "parent_gene_inference_performed":
                    bool(
                        summary[
                            "parent_gene_inference_performed"
                        ]
                    ),

                "nearest_gene_mapping_performed":
                    bool(
                        summary[
                            "nearest_gene_mapping_performed"
                        ]
                    ),

                "formal_colocalization_performed":
                    bool(
                        summary[
                            "formal_colocalization_performed"
                        ]
                    ),

                "causal_inference_performed":
                    bool(
                        summary[
                            "causal_inference_performed"
                        ]
                    ),

                "next_stage":
                    str(
                        summary[
                            "next_stage"
                        ]
                    ),
            },

            "feature_evidence":
                _records_with_json_nulls(
                    result.feature_evidence
                ),

            "candidate_mapping":
                _records_with_json_nulls(
                    result.candidate_mapping
                ),

            "outputs": {
                "feature_evidence":
                    str(
                        self.feature_evidence_output
                    ),

                "candidate_mapping":
                    str(
                        self.candidate_mapping_output
                    ),

                "summary":
                    str(
                        self.summary_output
                    ),
            },
        }

        return report

    # ======================================================================
    # Run
    # ======================================================================

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute complete M6.2A workflow."""

        config = self._load_config()

        # ------------------------------------------------------------------
        # Locked upstream candidate prioritization
        # ------------------------------------------------------------------

        (
            m6_1_qc,
            m6_1_path,
        ) = self._load_m6_1_qc(
            config
        )

        logger.info(
            "Loaded M6.1 candidate prioritization QC."
        )

        # ------------------------------------------------------------------
        # Raw Moradi sources
        # ------------------------------------------------------------------

        (
            s1,
            s1_path,
        ) = self._load_s1(
            config
        )

        (
            s2,
            s2_path,
        ) = self._load_s2(
            config
        )

        (
            s7,
            s7_path,
        ) = self._load_s7(
            config
        )

        logger.info(
            "Loaded Moradi S1: %d rows.",
            len(
                s1
            ),
        )

        logger.info(
            "Loaded Moradi S2: %d rows.",
            len(
                s2
            ),
        )

        logger.info(
            "Loaded Moradi S7: %d rows.",
            len(
                s7
            ),
        )

        # ------------------------------------------------------------------
        # Core analysis
        # ------------------------------------------------------------------

        result = resolve_regulatory_features(
            m6_1_qc=m6_1_qc,
            s1=s1,
            s2=s2,
            s7=s7,
            config=config,
        )

        if result.candidate_mapping.empty:

            raise RuntimeError(
                "M6.2A candidate mapping is unexpectedly empty."
            )

        # ------------------------------------------------------------------
        # Output directories
        # ------------------------------------------------------------------

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        # ------------------------------------------------------------------
        # Analytical Parquet outputs
        # ------------------------------------------------------------------

        write_parquet(
            result.feature_evidence,
            self.feature_evidence_output,
            index=False,
        )

        write_parquet(
            result.candidate_mapping,
            self.candidate_mapping_output,
            index=False,
        )

        write_parquet(
            result.summary,
            self.summary_output,
            index=False,
        )

        # ------------------------------------------------------------------
        # QC
        # ------------------------------------------------------------------

        report = self._build_qc_report(
            config=config,
            m6_1_qc=m6_1_qc,
            m6_1_path=m6_1_path,
            s1=s1,
            s1_path=s1_path,
            s2=s2,
            s2_path=s2_path,
            s7=s7,
            s7_path=s7_path,
            result=result,
        )

        with self.qc_output.open(
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

        # ------------------------------------------------------------------
        # Logging
        # ------------------------------------------------------------------

        logger.info(
            "M6.2A complete."
        )

        logger.info(
            "Candidates=%d | with feature=%d | without feature=%d.",
            report[
                "summary"
            ][
                "candidates_assessed"
            ],
            report[
                "summary"
            ][
                "candidates_with_regulatory_feature"
            ],
            report[
                "summary"
            ][
                "candidates_without_regulatory_feature"
            ],
        )

        logger.info(
            "Confirmed regulatory features=%d | "
            "parent genes resolved=%d | pending=%d.",
            report[
                "summary"
            ][
                "fully_confirmed_regulatory_features"
            ],
            report[
                "summary"
            ][
                "features_with_parent_gene_resolved"
            ],
            report[
                "summary"
            ][
                "features_with_parent_gene_pending"
            ],
        )

        for row in report[
            "candidate_mapping"
        ]:

            logger.info(
                "%s | feature=%s | status=%s | "
                "continue_to_gene_annotation=%s.",
                row[
                    "lead_rsid"
                ],
                row.get(
                    "regulatory_feature_id"
                ),
                row[
                    "resolution_status"
                ],
                row[
                    "continue_to_gene_annotation"
                ],
            )

        logger.info(
            "Next stage: %s.",
            report[
                "summary"
            ][
                "next_stage"
            ],
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report