"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/build_candidate_loci.py

Description:
    Combined runner for:

        M5.1 - LD Reference Strategy
        M5.2 - Candidate Locus Construction

    The runner writes the LD analysis policy separately from the
    candidate locus table.

    No LD computation is performed.

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

from pcatrfqtl.analysis.m5.candidate_loci import (
    M52CandidateLocusBuilder,
)
from pcatrfqtl.analysis.m5.ld_strategy import (
    default_m5_ld_strategy,
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
class M51M52Inputs:
    """Input datasets required by M5.1-M5.2."""

    bridge_readiness: Path

    prostate_trfqtl: Path


class M51M52CandidateLocusRunner:
    """
    Write LD strategy and construct physical candidate loci.
    """

    LOCUS_FILENAME = (
        "candidate_loci.parquet"
    )

    STRATEGY_FILENAME = (
        "m5_1_ld_reference_strategy.json"
    )

    QC_FILENAME = (
        "m5_2_candidate_loci_summary.json"
    )

    def __init__(
        self,
        *,
        inputs: M51M52Inputs,
        output_directory: str | Path,
        qc_directory: str | Path,
        screening_window_bp: int = 500_000,
        primary_r2_threshold: float = 0.8,
        secondary_r2_threshold: float = 0.5,
    ) -> None:
        """Initialize M5.1 + M5.2."""

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

        self.screening_window_bp = (
            screening_window_bp
        )

        self.primary_r2_threshold = (
            primary_r2_threshold
        )

        self.secondary_r2_threshold = (
            secondary_r2_threshold
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _require_file(
        path: Path,
    ) -> None:
        """Require an existing input file."""

        if not path.exists():
            raise FileNotFoundError(
                f"M5 required input not found: {path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                f"M5 expected file but received: {path}"
            )

    @staticmethod
    def _count_true(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count true values."""

        if column not in dataframe.columns:
            return 0

        return int(
            dataframe[
                column
            ]
            .fillna(False)
            .eq(True)
            .sum()
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
            str(key): int(value)
            for key, value
            in values.items()
        }

    # ------------------------------------------------------------------
    # Runner
    # ------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M5.1 and M5.2."""

        self._require_file(
            self.inputs.bridge_readiness
        )

        self._require_file(
            self.inputs.prostate_trfqtl
        )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        # --------------------------------------------------------------
        # M5.1 strategy
        # --------------------------------------------------------------

        strategy = default_m5_ld_strategy(
            screening_window_bp=(
                self.screening_window_bp
            ),
            primary_r2_threshold=(
                self.primary_r2_threshold
            ),
            secondary_r2_threshold=(
                self.secondary_r2_threshold
            ),
        )

        strategy_path = (
            self.qc_directory
            / self.STRATEGY_FILENAME
        )

        strategy_report = {
            "milestone":
                "M5.1",

            "stage":
                "ld_reference_strategy",

            "strategy":
                strategy.to_dict(),

            "interpretation_safeguards": {
                "screening_window_is_ld_block":
                    False,

                "physical_proximity_implies_ld":
                    False,

                "source_reported_ld_equals_recomputed_ld":
                    False,

                "ld_implies_causality":
                    False,
            },
        }

        with strategy_path.open(
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                strategy_report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        logger.info(
            "M5.1 LD reference strategy written: %s",
            strategy_path,
        )

        # --------------------------------------------------------------
        # M5.2 loci
        # --------------------------------------------------------------

        bridge = read_parquet(
            self.inputs.bridge_readiness
        )

        trfqtl = read_parquet(
            self.inputs.prostate_trfqtl
        )

        logger.info(
            "Starting M5.2 candidate locus construction."
        )

        loci = (
            M52CandidateLocusBuilder
            .build(
                bridge,
                trfqtl,
                strategy=strategy,
            )
        )

        if (
            len(
                loci
            )
            != len(
                bridge
            )
        ):
            raise RuntimeError(
                "M5.2 changed candidate cardinality: "
                f"{len(bridge)} -> {len(loci)}"
            )

        output_path = (
            self.output_directory
            / self.LOCUS_FILENAME
        )

        write_parquet(
            loci,
            output_path,
            index=False,
        )

        qc_report: dict[str, Any] = {
            "milestone":
                "M5.2",

            "stage":
                "candidate_locus_construction",

            "inputs": {
                "bridge_readiness":
                    str(
                        self.inputs.bridge_readiness
                    ),

                "prostate_trfqtl":
                    str(
                        self.inputs.prostate_trfqtl
                    ),
            },

            "output":
                str(
                    output_path
                ),

            "strategy_reference":
                str(
                    strategy_path
                ),

            "policy": {
                "candidate_cardinality_preserved":
                    True,

                "screening_window_bp":
                    strategy.screening_window_bp,

                "screening_window_interpreted_as_ld_block":
                    False,

                "coordinate_matching_performed":
                    False,

                "external_ld_calculated":
                    False,

                "liftover_performed":
                    False,

                "colocalization_performed":
                    False,

                "candidate_ranking_performed":
                    False,
            },

            "summary": {
                "candidate_rows":
                    len(
                        loci
                    ),

                "unique_lead_rsids":
                    int(
                        loci[
                            "trfqtl_rsid"
                        ]
                        .dropna()
                        .nunique()
                    ),

                "lead_coordinate_usable":
                    self._count_true(
                        loci,
                        "lead_coordinate_usable",
                    ),

                "ld_query_ready":
                    self._count_true(
                        loci,
                        "m5_ld_query_ready",
                    ),

                "locus_status_counts":
                    self._value_counts(
                        loci,
                        "m5_locus_status",
                    ),

                "chromosome_counts":
                    self._value_counts(
                        loci,
                        "lead_chromosome",
                    ),
            },
        }

        qc_path = (
            self.qc_directory
            / self.QC_FILENAME
        )

        qc_report[
            "report_path"
        ] = str(
            qc_path
        )

        with qc_path.open(
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                qc_report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        logger.info(
            "M5.2 completed: %d candidate loci; "
            "%d LD-query ready.",
            len(
                loci
            ),
            self._count_true(
                loci,
                "m5_ld_query_ready",
            ),
        )

        logger.info(
            "M5.2 QC report: %s",
            qc_path,
        )

        return {
            "m5_1":
                strategy_report,

            "m5_2":
                qc_report,
        }