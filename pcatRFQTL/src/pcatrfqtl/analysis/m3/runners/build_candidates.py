"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/runners/build_candidates.py

Description:
    Dataset-level runner for M3.4 integrated SNP-tRF candidate table.

    Inputs:

        data/processed/m3/prostate_gwas.parquet
        data/processed/m3/prostate_trfqtl.parquet

    Output:

        data/processed/m3/candidate_snptrf.parquet

    QC:

        data/processed/m3/qc/m3_4_candidate_summary.json

    M3.4 preserves every prostate-cancer tRF-QTL association while
    adding direct GWAS evidence derived from exact canonical-rsID
    identity.

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

from pcatrfqtl.analysis.m3.candidates import (
    CandidateDirectEvidence,
    SNPTRFCandidateBuilder,
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
class M34CandidateInputs:
    """Input paths required for M3.4."""

    prostate_gwas: Path
    prostate_trfqtl: Path


class M34CandidateRunner:
    """Execute M3.4 integrated SNP-tRF candidate generation."""

    OUTPUT_FILENAME = (
        "candidate_snptrf.parquet"
    )

    SUMMARY_FILENAME = (
        "m3_4_candidate_summary.json"
    )

    def __init__(
        self,
        inputs: M34CandidateInputs,
        output_directory: str | Path,
    ) -> None:
        """Initialize M3.4 runner."""

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

    @staticmethod
    def _require_file(
        path: Path,
    ) -> None:
        """Require a valid M3 input file."""

        if not path.exists():
            raise FileNotFoundError(
                f"M3.4 input file not found: {path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                f"M3.4 expected file but found: {path}"
            )

    @staticmethod
    def _safe_nunique(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count unique non-missing values."""

        if (
            dataframe.empty
            or column
            not in dataframe.columns
        ):
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

        if (
            dataframe.empty
            or column
            not in dataframe.columns
        ):
            return {}

        counts = (
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
            ) in counts.items()
        }

    @staticmethod
    def _count_true(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count explicit True values."""

        if (
            column
            not in dataframe.columns
        ):
            return 0

        return int(
            dataframe[
                column
            ]
            .fillna(False)
            .eq(True)
            .sum()
        )

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M3.4."""

        self._require_file(
            self.inputs.prostate_gwas
        )

        self._require_file(
            self.inputs.prostate_trfqtl
        )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        qc_directory = (
            self.output_directory
            / "qc"
        )

        qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        logger.info(
            "Starting M3.4 integrated SNP-tRF candidate generation."
        )

        gwas = read_parquet(
            self.inputs.prostate_gwas
        )

        trfqtl = read_parquet(
            self.inputs.prostate_trfqtl
        )

        candidates = (
            SNPTRFCandidateBuilder
            .build(
                gwas=gwas,
                trfqtl=trfqtl,
            )
        )

        if (
            len(
                candidates
            )
            != len(
                trfqtl
            )
        ):
            raise RuntimeError(
                "M3.4 failed candidate cardinality invariant."
            )

        output_path = (
            self.output_directory
            / self.OUTPUT_FILENAME
        )

        write_parquet(
            candidates,
            output_path,
            index=False,
        )

        direct_match_rows = (
            self._count_true(
                candidates,
                "m3_direct_gwas_match",
            )
        )

        no_direct_match_rows = int(
            candidates[
                "m3_candidate_direct_evidence"
            ]
            .eq(
                CandidateDirectEvidence
                .NO_DIRECT_GWAS_MATCH
                .value
            )
            .sum()
        )

        ineligible_rows = int(
            candidates[
                "m3_candidate_direct_evidence"
            ]
            .eq(
                CandidateDirectEvidence
                .DIRECT_MATCH_INELIGIBLE
                .value
            )
            .sum()
        )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M3.4",

            "stage":
                "integrated_snptrf_candidate_table",

            "inputs": {
                "gwas": {
                    "path":
                        str(
                            self.inputs.prostate_gwas
                        ),
                    "rows":
                        int(
                            len(
                                gwas
                            )
                        ),
                },

                "trfqtl": {
                    "path":
                        str(
                            self.inputs.prostate_trfqtl
                        ),
                    "rows":
                        int(
                            len(
                                trfqtl
                            )
                        ),
                },
            },

            "output": {
                "path":
                    str(
                        output_path
                    ),

                "rows":
                    int(
                        len(
                            candidates
                        )
                    ),

                "unique_candidate_ids":
                    self._safe_nunique(
                        candidates,
                        "m3_candidate_id",
                    ),

                "unique_trfqtl_rsids":
                    self._safe_nunique(
                        candidates,
                        "harm_variant_rsid",
                    ),

                "unique_trfs":
                    self._safe_nunique(
                        candidates,
                        "harm_feature_id",
                    ),

                "unique_variant_trf_pairs":
                    int(
                        candidates[
                            [
                                "harm_variant_rsid",
                                "harm_feature_identity_key",
                            ]
                        ]
                        .dropna()
                        .drop_duplicates()
                        .shape[
                            0
                        ]
                    ),

                "direct_gwas_match_rows":
                    direct_match_rows,

                "no_direct_gwas_match_rows":
                    no_direct_match_rows,

                "direct_match_ineligible_rows":
                    ineligible_rows,
            },

            "method": {
                "candidate_unit":
                    "prostate_cancer_snptrf_association",

                "direct_gwas_evidence_key":
                    "canonical_rsid",

                "candidate_cardinality_preserved":
                    True,

                "trfqtl_candidates_removed":
                    False,

                "ld_analysis_performed":
                    False,

                "distance_analysis_performed":
                    False,

                "coordinate_matching_performed":
                    False,

                "liftover_performed":
                    False,

                "colocalization_performed":
                    False,

                "candidate_ranking_performed":
                    False,
            },

            "qc": {
                "candidate_rows_equal_trfqtl_rows":
                    bool(
                        len(
                            candidates
                        )
                        == len(
                            trfqtl
                        )
                    ),

                "candidate_ids_unique":
                    bool(
                        candidates[
                            "m3_candidate_id"
                        ]
                        .is_unique
                    ),

                "candidate_direct_evidence_counts":
                    self._value_counts(
                        candidates,
                        "m3_candidate_direct_evidence",
                    ),

                "gwas_association_count_distribution":
                    self._value_counts(
                        candidates,
                        "m3_gwas_association_count",
                    ),
            },
        }

        summary_path = (
            qc_directory
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
            "M3.4 complete: %d SNP-tRF candidates.",
            len(
                candidates
            ),
        )

        logger.info(
            "M3.4 direct GWAS matches: %d.",
            direct_match_rows,
        )

        logger.info(
            "M3.4 candidates without direct GWAS match: %d.",
            no_direct_match_rows,
        )

        logger.info(
            "M3.4 QC report: %s",
            summary_path,
        )

        return report