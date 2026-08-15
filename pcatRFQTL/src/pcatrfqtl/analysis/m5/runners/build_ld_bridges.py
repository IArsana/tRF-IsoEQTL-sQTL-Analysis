"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/build_ld_bridges.py

Description:
    Combined runner for:

        M5.3 - LD-aware GWAS bridge
        M5.4 - LD-aware Moradi regulatory bridge

    This runner consumes population-specific external LD evidence from:

        data/processed/m5/ld/<POPULATION>/

    Example:

        data/processed/m5/ld/EAS/
        data/processed/m5/ld/EUR/
        data/processed/m5/ld/SAS/

    The default primary exploratory population is EAS.

    Important behavior:

        - Only Parquet LD evidence from the selected population directory
          is loaded.

        - Lead variants that are unavailable in the selected LD reference
          population are not fabricated as empty LD evidence rows.

        - Candidate cardinality is preserved upstream in candidate_loci;
          M5.3 and M5.4 are long-form evidence tables and therefore contain
          only actual LD-mediated bridge evidence.

        - Imported LD evidence is normalized once and reused by both M5.3
          and M5.4.

        - Moradi direct QTL variants and source-reported tag variants remain
          distinct evidence roles.

        - GWAS matching is performed using canonical rsID identity and does
          not assume coordinate-build compatibility with GWAS Catalog.

        - No colocalization, causal inference, or candidate ranking is
          performed.

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
from pcatrfqtl.analysis.m5.gwas_ld_bridge import (
    M53GWASLDBridgeBuilder,
)
from pcatrfqtl.analysis.m5.ld_evidence import (
    M5LDEvidenceNormalizer,
)
from pcatrfqtl.analysis.m5.moradi_ld_bridge import (
    M54MoradiLDBridgeBuilder,
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


# ============================================================================
# Data classes
# ============================================================================


@dataclass(frozen=True)
class M53M54Inputs:
    """
    Inputs required for M5.3 and M5.4.

    ld_root_directory points to:

        data/processed/m5/ld/

    The selected population subdirectory is resolved automatically.
    """

    candidate_loci: Path

    ld_root_directory: Path

    prostate_gwas: Path

    moradi_index_directory: Path


# ============================================================================
# Runner
# ============================================================================


class M53M54LDBridgeRunner:
    """Build population-specific GWAS and Moradi LD-aware bridges."""

    NORMALIZED_LD_FILENAME = (
        "normalized_ld_evidence.parquet"
    )

    GWAS_OUTPUT_TEMPLATE = (
        "gwas_ld_bridge_{population}.parquet"
    )

    MORADI_OUTPUT_TEMPLATE = (
        "moradi_ld_bridge_{population}.parquet"
    )

    M53_QC_TEMPLATE = (
        "m5_3_gwas_ld_bridge_{population}_summary.json"
    )

    M54_QC_TEMPLATE = (
        "m5_4_moradi_ld_bridge_{population}_summary.json"
    )

    def __init__(
        self,
        *,
        inputs: M53M54Inputs,
        output_directory: str | Path,
        qc_directory: str | Path,
        population: str = "EAS",
        primary_r2_threshold: float = 0.8,
        secondary_r2_threshold: float = 0.5,
    ) -> None:
        """Initialize population-specific M5.3/M5.4 analysis."""

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

        self.population = (
            str(
                population
            )
            .strip()
            .upper()
        )

        self.primary_r2_threshold = (
            primary_r2_threshold
        )

        self.secondary_r2_threshold = (
            secondary_r2_threshold
        )

        if self.population not in {
            "EAS",
            "EUR",
            "SAS",
        }:

            raise ValueError(
                "Unsupported M5 LD population: "
                f"{self.population}. "
                "Expected one of EAS, EUR, SAS."
            )

        if not (
            0.0
            <= self.secondary_r2_threshold
            <= self.primary_r2_threshold
            <= 1.0
        ):

            raise ValueError(
                "LD thresholds must satisfy "
                "0 <= secondary <= primary <= 1."
            )

    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------

    @property
    def ld_population_directory(
        self,
    ) -> Path:
        """Return selected population-specific LD directory."""

        return (
            Path(
                self.inputs.ld_root_directory
            )
            / self.population
        )

    @property
    def normalized_ld_output_path(
        self,
    ) -> Path:
        """Return normalized LD evidence output path."""

        return (
            self.output_directory
            / "ld"
            / self.population
            / self.NORMALIZED_LD_FILENAME
        )

    @property
    def gwas_output_path(
        self,
    ) -> Path:
        """Return M5.3 output path."""

        return (
            self.output_directory
            / self.GWAS_OUTPUT_TEMPLATE.format(
                population=self.population
            )
        )

    @property
    def moradi_output_path(
        self,
    ) -> Path:
        """Return M5.4 output path."""

        return (
            self.output_directory
            / self.MORADI_OUTPUT_TEMPLATE.format(
                population=self.population
            )
        )

    @property
    def m53_qc_path(
        self,
    ) -> Path:
        """Return M5.3 QC path."""

        return (
            self.qc_directory
            / self.M53_QC_TEMPLATE.format(
                population=self.population
            )
        )

    @property
    def m54_qc_path(
        self,
    ) -> Path:
        """Return M5.4 QC path."""

        return (
            self.qc_directory
            / self.M54_QC_TEMPLATE.format(
                population=self.population
            )
        )

    # ------------------------------------------------------------------
    # Generic helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _require_file(
        path: Path,
    ) -> None:
        """Require an existing file."""

        if not path.exists():

            raise FileNotFoundError(
                f"Required M5 file not found: {path}"
            )

        if not path.is_file():

            raise FileNotFoundError(
                f"M5 expected a file but received: {path}"
            )

    @staticmethod
    def _require_directory(
        path: Path,
    ) -> None:
        """Require an existing directory."""

        if not path.exists():

            raise FileNotFoundError(
                f"Required M5 directory not found: {path}"
            )

        if not path.is_dir():

            raise NotADirectoryError(
                f"M5 expected a directory but received: {path}"
            )

    @staticmethod
    def _value_counts(
        dataframe: pd.DataFrame,
        column: str,
    ) -> dict[str, int]:
        """Return JSON-safe value counts."""

        if column not in dataframe.columns:
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
            str(key): int(value)
            for key, value
            in counts.items()
        }

    @staticmethod
    def _safe_nunique(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Return unique non-missing values."""

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
    def _count_true(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count true values safely."""

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

    # ------------------------------------------------------------------
    # LD loading
    # ------------------------------------------------------------------

    def _load_ld_evidence(
        self,
    ) -> tuple[
        pd.DataFrame,
        dict[str, Any],
    ]:
        """
        Load and normalize LD evidence for one population.

        Missing lead-specific Parquet files are allowed because some lead
        variants may be monoallelic or otherwise unavailable in the selected
        reference population.
        """

        directory = (
            self.ld_population_directory
        )

        self._require_directory(
            directory
        )

        paths = sorted(
            path
            for path
            in directory.glob(
                "*.parquet"
            )
            if (
                path.name
                != self.NORMALIZED_LD_FILENAME
            )
        )

        if not paths:

            raise RuntimeError(
                "No external/precomputed LD Parquet files found for "
                f"population {self.population}: {directory}. "
                "M5.3/M5.4 must not fabricate LD evidence."
            )

        tables: list[
            tuple[
                pd.DataFrame,
                str,
            ]
        ] = []

        source_reports: dict[
            str,
            Any,
        ] = {}

        for path in paths:

            logger.info(
                "Loading %s LD evidence: %s",
                self.population,
                path.name,
            )

            dataframe = read_parquet(
                path
            )

            tables.append(
                (
                    dataframe,
                    path.name,
                )
            )

            source_reports[
                path.name
            ] = {
                "rows":
                    len(
                        dataframe
                    ),

                "lead_rsids":
                    (
                        sorted(
                            dataframe[
                                "lead_rsid"
                            ]
                            .dropna()
                            .astype(str)
                            .unique()
                            .tolist()
                        )
                        if (
                            "lead_rsid"
                            in dataframe.columns
                        )
                        else []
                    ),

                "population_values":
                    (
                        sorted(
                            dataframe[
                                "population"
                            ]
                            .dropna()
                            .astype(str)
                            .unique()
                            .tolist()
                        )
                        if (
                            "population"
                            in dataframe.columns
                        )
                        else []
                    ),
            }

        normalized = (
            M5LDEvidenceNormalizer
            .combine(
                tables,
                primary_r2_threshold=(
                    self.primary_r2_threshold
                ),
                secondary_r2_threshold=(
                    self.secondary_r2_threshold
                ),
            )
        )

        if normalized.empty:

            raise RuntimeError(
                "Population-specific LD evidence normalized to zero rows "
                f"for {self.population}."
            )

        # --------------------------------------------------------------
        # Validate population provenance
        # --------------------------------------------------------------

        if (
            "ld_population"
            in normalized.columns
        ):

            observed_populations = set(
                normalized[
                    "ld_population"
                ]
                .dropna()
                .astype(str)
                .str.upper()
            )

            unexpected = (
                observed_populations
                - {
                    self.population
                }
            )

            if unexpected:

                raise RuntimeError(
                    "LD evidence contains population values outside "
                    f"the selected population {self.population}: "
                    f"{sorted(unexpected)}"
                )

        report = {
            "population":
                self.population,

            "directory":
                str(
                    directory
                ),

            "source_files":
                len(
                    paths
                ),

            "source_file_details":
                source_reports,

            "normalized_rows":
                len(
                    normalized
                ),

            "normalized_unique_lead_rsids":
                self._safe_nunique(
                    normalized,
                    "lead_rsid",
                ),

            "bridge_eligible_rows":
                self._count_true(
                    normalized,
                    "m5_bridge_eligible",
                ),

            "primary_ld_rows":
                self._count_true(
                    normalized,
                    "m5_primary_ld_support",
                ),

            "ld_class_counts":
                self._value_counts(
                    normalized,
                    "ld_evidence_class",
                ),
        }

        return (
            normalized,
            report,
        )

    # ------------------------------------------------------------------
    # Candidate availability
    # ------------------------------------------------------------------

    def _build_ld_availability_summary(
        self,
        *,
        loci: pd.DataFrame,
        ld_evidence: pd.DataFrame,
    ) -> dict[str, Any]:
        """
        Compare full M5 candidate set against available LD evidence.

        This prevents population-specific LD unavailability from being
        interpreted as candidate removal.
        """

        candidate_leads = set(
            loci[
                "trfqtl_rsid"
            ]
            .dropna()
            .astype(str)
            .str.lower()
        )

        evidence_leads = set(
            ld_evidence[
                "lead_rsid"
            ]
            .dropna()
            .astype(str)
            .str.lower()
        )

        unavailable = (
            candidate_leads
            - evidence_leads
        )

        unexpected = (
            evidence_leads
            - candidate_leads
        )

        if unexpected:

            raise RuntimeError(
                "LD evidence contains lead rsIDs outside the M5 candidate "
                f"set: {sorted(unexpected)}"
            )

        return {
            "candidate_lead_rsids":
                sorted(
                    candidate_leads
                ),

            "candidate_lead_count":
                len(
                    candidate_leads
                ),

            "ld_available_rsids":
                sorted(
                    evidence_leads
                ),

            "ld_available_count":
                len(
                    evidence_leads
                ),

            "ld_unavailable_rsids":
                sorted(
                    unavailable
                ),

            "ld_unavailable_count":
                len(
                    unavailable
                ),

            "candidate_set_reduced":
                False,
        }

    # ------------------------------------------------------------------
    # M5.3
    # ------------------------------------------------------------------

    def _run_gwas_bridge(
        self,
        *,
        loci: pd.DataFrame,
        ld_evidence: pd.DataFrame,
    ) -> tuple[
        pd.DataFrame,
        dict[str, Any],
    ]:
        """Execute M5.3 LD-aware GWAS bridge."""

        self._require_file(
            self.inputs.prostate_gwas
        )

        prostate_gwas = read_parquet(
            self.inputs.prostate_gwas
        )

        logger.info(
            "Starting M5.3 LD-aware GWAS bridge "
            "for population %s.",
            self.population,
        )

        result = (
            M53GWASLDBridgeBuilder
            .build(
                loci,
                ld_evidence,
                prostate_gwas,
            )
        )

        if not result.empty:

            result[
                "m5_ld_analysis_population"
            ] = self.population

            result = result.sort_values(
                by=[
                    "trfqtl_rsid",
                    "r2",
                    "neighbor_rsid",
                ],
                ascending=[
                    True,
                    False,
                    True,
                ],
                kind="stable",
            ).reset_index(
                drop=True
            )

            result[
                "m5_gwas_bridge_id"
            ] = [
                (
                    f"M5G-{self.population}-"
                    f"{index:07d}"
                )
                for index
                in range(
                    1,
                    len(
                        result
                    )
                    + 1,
                )
            ]

        write_parquet(
            result,
            self.gwas_output_path,
            index=False,
        )

        matched_leads = (
            set(
                result[
                    "trfqtl_rsid"
                ]
                .dropna()
                .astype(str)
            )
            if (
                not result.empty
                and "trfqtl_rsid"
                in result.columns
            )
            else set()
        )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M5.3",

            "stage":
                "ld_aware_gwas_bridge",

            "population":
                self.population,

            "output":
                str(
                    self.gwas_output_path
                ),

            "policy": {
                "matching_key":
                    "canonical_rsid",

                "gwas_coordinate_matching":
                    False,

                "external_ld_used":
                    True,

                "population_specific_ld":
                    True,

                "primary_r2_threshold":
                    self.primary_r2_threshold,

                "secondary_r2_threshold":
                    self.secondary_r2_threshold,

                "physical_proximity_alone_accepted":
                    False,

                "missing_population_ld_imputed":
                    False,

                "colocalization_performed":
                    False,

                "causal_inference_performed":
                    False,

                "candidate_ranking_performed":
                    False,
            },

            "summary": {
                "bridge_rows":
                    len(
                        result
                    ),

                "lead_rsids_with_gwas_bridge":
                    len(
                        matched_leads
                    ),

                "unique_gwas_proxy_rsids":
                    self._safe_nunique(
                        result,
                        "gwas_rsid",
                    ),

                "primary_high_ld_bridge_rows":
                    self._count_true(
                        result,
                        "m5_primary_gwas_ld_bridge",
                    ),

                "ld_class_counts":
                    self._value_counts(
                        result,
                        "ld_evidence_class",
                    ),
            },

            "matched_lead_rsids":
                sorted(
                    matched_leads
                ),
        }

        return (
            result,
            report,
        )

    # ------------------------------------------------------------------
    # M5.4
    # ------------------------------------------------------------------

    def _run_moradi_bridge(
        self,
        *,
        loci: pd.DataFrame,
        ld_evidence: pd.DataFrame,
    ) -> tuple[
        pd.DataFrame,
        dict[str, Any],
    ]:
        """Execute M5.4 LD-aware Moradi regulatory bridge."""

        logger.info(
            "Starting M5.4 LD-aware Moradi regulatory bridge "
            "for population %s.",
            self.population,
        )

        frames: list[
            pd.DataFrame
        ] = []

        table_reports: dict[
            str,
            Any,
        ] = {}

        for table_id in MORADI_INDEX_SPECS:

            path = (
                Path(
                    self.inputs.moradi_index_directory
                )
                / f"{table_id}.parquet"
            )

            self._require_file(
                path
            )

            moradi = read_parquet(
                path
            )

            matched = (
                M54MoradiLDBridgeBuilder
                .build_table(
                    loci,
                    ld_evidence,
                    moradi,
                )
            )

            if not matched.empty:

                matched[
                    "m5_index_part"
                ] = table_id

                matched[
                    "m5_ld_analysis_population"
                ] = self.population

                frames.append(
                    matched
                )

            table_reports[
                table_id
            ] = {
                "source_rows":
                    len(
                        moradi
                    ),

                "bridge_rows":
                    len(
                        matched
                    ),

                "matched_lead_rsids":
                    self._safe_nunique(
                        matched,
                        "trfqtl_rsid",
                    ),

                "matched_proxy_rsids":
                    self._safe_nunique(
                        matched,
                        "neighbor_rsid",
                    ),

                "matched_features":
                    self._safe_nunique(
                        matched,
                        "feature_identity_key",
                    ),

                "bridge_type_counts":
                    self._value_counts(
                        matched,
                        "m5_bridge_type",
                    ),
            }

            if not matched.empty:

                logger.info(
                    "M5.4 %s [%s]: %d regulatory bridge rows.",
                    table_id,
                    self.population,
                    len(
                        matched
                    ),
                )

            del moradi

        if frames:

            result = pd.concat(
                frames,
                ignore_index=True,
                sort=False,
            )

            sort_columns = [
                column
                for column
                in (
                    "trfqtl_rsid",
                    "r2",
                    "m5_bridge_type",
                    "m4_source_table",
                    "m4_source_row",
                )
                if column
                in result.columns
            ]

            ascending_map = {
                "trfqtl_rsid":
                    True,

                "r2":
                    False,

                "m5_bridge_type":
                    True,

                "m4_source_table":
                    True,

                "m4_source_row":
                    True,
            }

            result = result.sort_values(
                by=sort_columns,
                ascending=[
                    ascending_map[
                        column
                    ]
                    for column
                    in sort_columns
                ],
                kind="stable",
            ).reset_index(
                drop=True
            )

            result[
                "m5_regulatory_bridge_id"
            ] = [
                (
                    f"M5R-{self.population}-"
                    f"{index:07d}"
                )
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

            result = pd.DataFrame()

        write_parquet(
            result,
            self.moradi_output_path,
            index=False,
        )

        matched_leads = (
            set(
                result[
                    "trfqtl_rsid"
                ]
                .dropna()
                .astype(str)
            )
            if (
                not result.empty
                and "trfqtl_rsid"
                in result.columns
            )
            else set()
        )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M5.4",

            "stage":
                "ld_aware_moradi_regulatory_bridge",

            "population":
                self.population,

            "output":
                str(
                    self.moradi_output_path
                ),

            "policy": {
                "external_ld_used":
                    True,

                "population_specific_ld":
                    True,

                "moradi_source_reported_ld_recomputed":
                    False,

                "moradi_direct_qtl_and_tag_distinguished":
                    True,

                "cis_and_trans_distinguished":
                    True,

                "primary_r2_threshold":
                    self.primary_r2_threshold,

                "secondary_r2_threshold":
                    self.secondary_r2_threshold,

                "missing_population_ld_imputed":
                    False,

                "colocalization_performed":
                    False,

                "causal_inference_performed":
                    False,

                "candidate_ranking_performed":
                    False,
            },

            "summary": {
                "bridge_rows":
                    len(
                        result
                    ),

                "lead_rsids_with_regulatory_bridge":
                    len(
                        matched_leads
                    ),

                "unique_bridge_proxy_rsids":
                    self._safe_nunique(
                        result,
                        "neighbor_rsid",
                    ),

                "unique_regulatory_features":
                    self._safe_nunique(
                        result,
                        "feature_identity_key",
                    ),

                "bridge_type_counts":
                    self._value_counts(
                        result,
                        "m5_bridge_type",
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

                "ld_class_counts":
                    self._value_counts(
                        result,
                        "ld_evidence_class",
                    ),
            },

            "matched_lead_rsids":
                sorted(
                    matched_leads
                ),

            "tables":
                table_reports,
        }

        return (
            result,
            report,
        )

    # ------------------------------------------------------------------
    # Public runner
    # ------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute population-specific M5.3 and M5.4."""

        self._require_file(
            self.inputs.candidate_loci
        )

        self._require_directory(
            self.inputs.moradi_index_directory
        )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.normalized_ld_output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        loci = read_parquet(
            self.inputs.candidate_loci
        )

        logger.info(
            "Starting M5.3 + M5.4 for LD population %s.",
            self.population,
        )

        (
            ld_evidence,
            ld_report,
        ) = self._load_ld_evidence()

        availability = (
            self._build_ld_availability_summary(
                loci=loci,
                ld_evidence=ld_evidence,
            )
        )

        logger.info(
            "M5 LD availability [%s]: %d/%d lead variants available.",
            self.population,
            availability[
                "ld_available_count"
            ],
            availability[
                "candidate_lead_count"
            ],
        )

        if availability[
            "ld_unavailable_rsids"
        ]:

            logger.info(
                "LD unavailable [%s]: %s",
                self.population,
                ", ".join(
                    availability[
                        "ld_unavailable_rsids"
                    ]
                ),
            )

        # --------------------------------------------------------------
        # Save normalized LD evidence
        # --------------------------------------------------------------

        write_parquet(
            ld_evidence,
            self.normalized_ld_output_path,
            index=False,
        )

        # --------------------------------------------------------------
        # M5.3
        # --------------------------------------------------------------

        (
            _,
            m53_report,
        ) = self._run_gwas_bridge(
            loci=loci,
            ld_evidence=ld_evidence,
        )

        # --------------------------------------------------------------
        # M5.4
        # --------------------------------------------------------------

        (
            _,
            m54_report,
        ) = self._run_moradi_bridge(
            loci=loci,
            ld_evidence=ld_evidence,
        )

        # --------------------------------------------------------------
        # Common provenance
        # --------------------------------------------------------------

        common_metadata = {
            "population":
                self.population,

            "candidate_ld_availability":
                availability,

            "ld_evidence":
                ld_report,

            "normalized_ld_evidence":
                str(
                    self.normalized_ld_output_path
                ),
        }

        m53_report[
            "candidate_ld_availability"
        ] = availability

        m53_report[
            "ld_evidence"
        ] = ld_report

        m53_report[
            "normalized_ld_evidence"
        ] = str(
            self.normalized_ld_output_path
        )

        m54_report[
            "candidate_ld_availability"
        ] = availability

        m54_report[
            "ld_evidence"
        ] = ld_report

        m54_report[
            "normalized_ld_evidence"
        ] = str(
            self.normalized_ld_output_path
        )

        # --------------------------------------------------------------
        # Write QC
        # --------------------------------------------------------------

        with self.m53_qc_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                m53_report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        with self.m54_qc_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                m54_report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        logger.info(
            "M5.3 + M5.4 completed for population %s.",
            self.population,
        )

        logger.info(
            "M5.3 output: %s",
            self.gwas_output_path,
        )

        logger.info(
            "M5.4 output: %s",
            self.moradi_output_path,
        )

        return {
            "population":
                self.population,

            "common":
                common_metadata,

            "m5_3":
                m53_report,

            "m5_4":
                m54_report,
        }