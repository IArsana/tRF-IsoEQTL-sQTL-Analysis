"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m4/runners/match_trfqtl_moradi.py

Description:
    Runner for M4.2 exact rsID matching between PRAD tRF-QTL
    associations and all parts of the M4.1 Unified Moradi QTL Index.

    Matching modes:

        - QTL_RSID_EXACT
        - SOURCE_TAG_RSID_EXACT

    Results from all twelve Moradi index parts are combined into one
    long-form evidence table.

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
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from pcatrfqtl.analysis.m4.moradi_index import (
    MORADI_INDEX_SPECS,
)
from pcatrfqtl.analysis.m4.variant_match import (
    TRFQTLMoradiVariantMatcher,
)
from pcatrfqtl.io.parquet import (
    read_parquet,
    write_parquet,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


logger = get_logger(
    __name__
)


@dataclass(frozen=True)
class M42TRFQTLMoradiInputs:
    """Input resources for M4.2."""

    prostate_trfqtl: Path

    moradi_index_directory: Path


class M42TRFQTLMoradiRunner:
    """Run exact tRF-QTL to Moradi variant integration."""

    OUTPUT_FILENAME = (
        "trfqtl_moradi_matches.parquet"
    )

    SUMMARY_FILENAME = (
        "m4_2_trfqtl_moradi_match_summary.json"
    )

    def __init__(
        self,
        *,
        inputs: M42TRFQTLMoradiInputs,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:
        """Initialize M4.2 runner."""

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _require_file(
        path: Path,
    ) -> None:
        """Require an existing file."""

        if not path.exists():
            raise FileNotFoundError(
                f"M4.2 required file not found: {path}"
            )

    @staticmethod
    def _safe_nunique(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count unique non-missing values."""

        if column not in dataframe.columns:
            return 0

        return int(
            dataframe[
                column
            ]
            .dropna()
            .nunique()
        )

    @staticmethod
    def _value_counts(
        dataframe: pd.DataFrame,
        column: str,
    ) -> dict[str, int]:
        """Return JSON-safe value counts."""

        if column not in dataframe.columns:
            return {}

        values = (
            dataframe[
                column
            ]
            .dropna()
            .astype(str)
            .value_counts()
        )

        return {
            str(
                key
            ): int(
                value
            )
            for (
                key,
                value,
            ) in values.items()
        }

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M4.2 exact rsID matching."""

        prostate_trfqtl_path = Path(
            self.inputs.prostate_trfqtl
        )

        self._require_file(
            prostate_trfqtl_path
        )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        trfqtl = read_parquet(
            prostate_trfqtl_path
        )

        prepared = (
            TRFQTLMoradiVariantMatcher
            .prepare_trfqtl(
                trfqtl
            )
        )

        logger.info(
            "Starting M4.2 exact tRF-QTL ↔ Moradi matching."
        )

        logger.info(
            "Eligible PRAD tRF-QTL rows: %d; unique rsIDs: %d.",
            len(
                prepared
            ),
            prepared[
                "harm_variant_rsid"
            ]
            .dropna()
            .nunique(),
        )

        matched_parts: list[
            pd.DataFrame
        ] = []

        table_reports: dict[
            str,
            Any,
        ] = {}

        for table_id in MORADI_INDEX_SPECS:

            part_path = (
                Path(
                    self.inputs.moradi_index_directory
                )
                / f"{table_id}.parquet"
            )

            self._require_file(
                part_path
            )

            moradi = read_parquet(
                part_path
            )

            matched = (
                TRFQTLMoradiVariantMatcher
                .match_table(
                    prepared,
                    moradi,
                )
            )

            if matched.empty:

                table_reports[
                    table_id
                ] = {
                    "moradi_rows":
                        len(
                            moradi
                        ),

                    "matched_evidence_rows":
                        0,

                    "matched_trfqtl_rsids":
                        0,

                    "match_modes":
                        {},
                }

                logger.info(
                    "M4.2 %s: no exact rsID evidence.",
                    table_id,
                )

                del moradi

                continue

            matched[
                "m4_index_part"
            ] = table_id

            matched_parts.append(
                matched
            )

            table_reports[
                table_id
            ] = {
                "moradi_rows":
                    len(
                        moradi
                    ),

                "matched_evidence_rows":
                    len(
                        matched
                    ),

                "matched_trfqtl_rsids":
                    self._safe_nunique(
                        matched,
                        "trfqtl_rsid",
                    ),

                "matched_regulatory_features":
                    self._safe_nunique(
                        matched,
                        "feature_identity_key",
                    ),

                "match_modes":
                    self._value_counts(
                        matched,
                        "m4_match_mode",
                    ),
            }

            logger.info(
                "M4.2 %s: %d evidence rows across %d tRF-QTL rsIDs.",
                table_id,
                len(
                    matched
                ),
                self._safe_nunique(
                    matched,
                    "trfqtl_rsid",
                ),
            )

            del moradi

        # --------------------------------------------------------------
        # Combined output
        # --------------------------------------------------------------

        if matched_parts:

            result = pd.concat(
                matched_parts,
                ignore_index=True,
                sort=False,
            )

            result = result.sort_values(
                by=[
                    "trfqtl_rsid",
                    "m4_match_mode",
                    "m4_source_table",
                    "m4_source_row",
                ],
                kind="stable",
            ).reset_index(
                drop=True
            )

            result[
                "m4_evidence_id"
            ] = [
                f"M4E{index:07d}"
                for index
                in range(
                    1,
                    len(
                        result
                    )
                    + 1,
                )
            ]

        else:

            result = pd.DataFrame(
                columns=[
                    "m4_evidence_id",
                    "trfqtl_rsid",
                    "trf_id",
                    "m4_match_mode",
                ]
            )

        output_path = (
            self.output_directory
            / self.OUTPUT_FILENAME
        )

        write_parquet(
            result,
            output_path,
            index=False,
        )

        # --------------------------------------------------------------
        # Candidate-level coverage
        # --------------------------------------------------------------

        eligible_rsids = set(
            prepared[
                "harm_variant_rsid"
            ]
            .dropna()
            .astype(str)
        )

        matched_rsids = set(
            result[
                "trfqtl_rsid"
            ]
            .dropna()
            .astype(str)
        ) if (
            "trfqtl_rsid"
            in result.columns
        ) else set()

        unmatched_rsids = (
            eligible_rsids
            - matched_rsids
        )

        direct_qtl_rsids = set()

        source_tag_rsids = set()

        if not result.empty:

            direct_mask = (
                result[
                    "m4_match_mode"
                ]
                == "QTL_RSID_EXACT"
            )

            tag_mask = (
                result[
                    "m4_match_mode"
                ]
                == "SOURCE_TAG_RSID_EXACT"
            )

            direct_qtl_rsids = set(
                result.loc[
                    direct_mask,
                    "trfqtl_rsid",
                ]
                .dropna()
                .astype(str)
            )

            source_tag_rsids = set(
                result.loc[
                    tag_mask,
                    "trfqtl_rsid",
                ]
                .dropna()
                .astype(str)
            )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M4.2",

            "stage":
                "exact_trfqtl_moradi_variant_matching",

            "inputs": {
                "prostate_trfqtl":
                    str(
                        prostate_trfqtl_path
                    ),

                "moradi_index_directory":
                    str(
                        self.inputs.moradi_index_directory
                    ),
            },

            "output":
                str(
                    output_path
                ),

            "policy": {
                "matching_key":
                    "canonical_rsid",

                "coordinate_matching_performed":
                    False,

                "liftover_performed":
                    False,

                "new_ld_inference_performed":
                    False,

                "source_tag_evidence_preserved":
                    True,

                "source_tag_evidence_interpreted_as_new_ld":
                    False,

                "colocalization_performed":
                    False,

                "statistical_filtering_performed":
                    False,

                "candidate_ranking_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "trfqtl_input_rows":
                    len(
                        trfqtl
                    ),

                "trfqtl_match_eligible_rows":
                    len(
                        prepared
                    ),

                "trfqtl_unique_eligible_rsids":
                    len(
                        eligible_rsids
                    ),

                "matched_evidence_rows":
                    len(
                        result
                    ),

                "matched_trfqtl_rsids":
                    len(
                        matched_rsids
                    ),

                "unmatched_trfqtl_rsids":
                    len(
                        unmatched_rsids
                    ),

                "direct_qtl_match_rsids":
                    len(
                        direct_qtl_rsids
                    ),

                "source_tag_match_rsids":
                    len(
                        source_tag_rsids
                    ),

                "unique_regulatory_features":
                    self._safe_nunique(
                        result,
                        "feature_identity_key",
                    ),

                "match_mode_counts":
                    self._value_counts(
                        result,
                        "m4_match_mode",
                    ),

                "regulatory_scope_counts":
                    self._value_counts(
                        result,
                        "m4_regulatory_scope",
                    ),

                "feature_class_counts":
                    self._value_counts(
                        result,
                        "m4_feature_class",
                    ),
            },

            "matched_rsids":
                sorted(
                    matched_rsids
                ),

            "unmatched_rsids":
                sorted(
                    unmatched_rsids
                ),

            "direct_qtl_match_rsids":
                sorted(
                    direct_qtl_rsids
                ),

            "source_tag_match_rsids":
                sorted(
                    source_tag_rsids
                ),

            "tables":
                table_reports,
        }

        summary_path = (
            self.qc_directory
            / self.SUMMARY_FILENAME
        )

        report[
            "report_path"
        ] = str(
            summary_path
        )

        with summary_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        logger.info(
            "M4.2 completed: %d exact evidence rows; "
            "%d/%d tRF-QTL rsIDs matched.",
            len(
                result
            ),
            len(
                matched_rsids
            ),
            len(
                eligible_rsids
            ),
        )

        logger.info(
            "M4.2 QC report: %s",
            summary_path,
        )

        return report