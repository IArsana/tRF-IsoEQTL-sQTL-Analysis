"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/integrate_disease_ld_bridges.py

Description:
    Runner for M5.3C.3F disease-signal ↔ LD bridge integration.

    Inputs:
        - M5.3C.3E suggestive disease-variant recurrence table.
        - Population-specific normalized LD evidence.
        - Existing population-specific Moradi regulatory LD bridge tables.

    Populations:
        - EAS
        - EUR
        - SAS

    Outputs:
        - Disease ↔ candidate-lead LD reference matches by population.
        - Disease ↔ regulatory shared-proxy bridges by population.
        - Cross-population LD synthesis.
        - Candidate-level bridge summary.
        - QC JSON report.

    Important terminology:
        RAW LD REFERENCE MATCH
            A disease-associated rsID is present in the lead-specific LD
            reference table. This does NOT imply that the r² value passes
            any biological LD threshold.

        MODERATE-OR-HIGH LD BRIDGE
            r² >= 0.50.

        HIGH LD BRIDGE
            r² >= 0.80.

        DISEASE-REGULATORY SHARED PROXY
            A disease-associated variant is the same canonical rsID as a
            previously identified Moradi regulatory LD proxy for the same
            candidate lead and reference population.

    Build policy:
        - Disease-side harmonised GWAS summary statistics are GRCh38.
        - Existing LD and Moradi evidence are GRCh37/hg19.
        - No cross-build coordinate join is performed.
        - No silent liftover is performed.
        - Cross-build integration is canonical-rsID based only.

    Scientific safeguards:
        - Raw LD lookup matches are not interpreted as LD bridges.
        - Physical proximity is not interpreted as LD.
        - Moderate/high LD does not imply causality.
        - Shared proxy identity does not prove colocalization.
        - Cross-population consistency is sensitivity evidence, not
          independent biological replication.
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
import pyarrow.parquet as pq

from pcatrfqtl.analysis.m5.disease_ld_bridge import (
    HIGH_LD_R2,
    MODERATE_LD_R2,
    build_disease_ld_bridge,
    build_disease_regulatory_bridge,
    synthesize_disease_ld_bridges,
)
from pcatrfqtl.io.parquet import write_parquet
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


POPULATIONS = (
    "EAS",
    "EUR",
    "SAS",
)


# ============================================================================
# Parquet compatibility
# ============================================================================


def _read_parquet_compat(
    path: Path,
) -> pd.DataFrame:
    """Read M5 Parquet artifacts without pandas reconstruction metadata."""

    if not path.exists():

        raise FileNotFoundError(
            f"Input Parquet not found: {path}"
        )

    table = pq.read_table(
        path,
        use_pandas_metadata=False,
    )

    return table.to_pandas(
        ignore_metadata=True,
    )


# ============================================================================
# Helpers
# ============================================================================


def _unique_variant_count(
    dataframe: pd.DataFrame,
    *,
    column: str = "disease_variant_rsid",
) -> int:
    """Count unique non-missing disease variants."""

    if dataframe.empty:
        return 0

    if column not in dataframe.columns:
        return 0

    return int(
        dataframe[
            column
        ]
        .dropna()
        .astype("string")
        .nunique()
    )


def _unique_lead_variant_pair_count(
    dataframe: pd.DataFrame,
) -> int:
    """Count unique candidate-lead × disease-variant pairs."""

    if dataframe.empty:
        return 0

    required = {
        "lead_rsid",
        "disease_variant_rsid",
    }

    if not required.issubset(
        dataframe.columns
    ):
        return 0

    return int(
        dataframe[
            [
                "lead_rsid",
                "disease_variant_rsid",
            ]
        ]
        .dropna()
        .drop_duplicates()
        .shape[0]
    )


def _safe_bool_series(
    series: pd.Series,
) -> pd.Series:
    """Convert nullable boolean-like data to a strict bool Series."""

    return (
        series
        .fillna(False)
        .astype(bool)
    )


# ============================================================================
# Runner
# ============================================================================


class M53C3FDiseaseLDBridgeRunner:
    """Execute M5.3C.3F disease-signal ↔ LD bridge integration."""

    def __init__(
        self,
        *,
        recurrence_path: str | Path,
        m5_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.recurrence_path = Path(
            recurrence_path
        )

        self.m5_directory = Path(
            m5_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    # ------------------------------------------------------------------
    # Existing input paths
    # ------------------------------------------------------------------

    def ld_path(
        self,
        population: str,
    ) -> Path:
        """Return normalized LD evidence path for one population."""

        return (
            self.m5_directory
            / "ld"
            / population
            / "normalized_ld_evidence.parquet"
        )

    def moradi_bridge_path(
        self,
        population: str,
    ) -> Path:
        """Return existing Moradi LD bridge path for one population."""

        return (
            self.m5_directory
            / f"moradi_ld_bridge_{population}.parquet"
        )

    # ------------------------------------------------------------------
    # Output paths
    # ------------------------------------------------------------------

    def disease_ld_output(
        self,
        population: str,
    ) -> Path:
        """Return per-population disease-LD output path."""

        return (
            self.m5_directory
            / f"disease_ld_bridge_{population}.parquet"
        )

    def disease_regulatory_output(
        self,
        population: str,
    ) -> Path:
        """Return per-population disease-regulatory bridge output path."""

        return (
            self.m5_directory
            / f"disease_regulatory_bridge_{population}.parquet"
        )

    @property
    def synthesis_output(
        self,
    ) -> Path:
        """Return cross-population LD synthesis path."""

        return (
            self.m5_directory
            / "disease_ld_bridge_synthesis.parquet"
        )

    @property
    def candidate_summary_output(
        self,
    ) -> Path:
        """Return candidate-level disease-LD summary path."""

        return (
            self.m5_directory
            / "disease_ld_candidate_summary.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:
        """Return M5.3C.3F QC JSON path."""

        return (
            self.qc_directory
            / "m5_3c_3f_disease_ld_bridge.json"
        )

    # ------------------------------------------------------------------
    # Main
    # ------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Execute M5.3C.3F."""

        recurrence = _read_parquet_compat(
            self.recurrence_path
        )

        logger.info(
            "Loaded M5.3C.3E suggestive variant identities: %d.",
            len(recurrence),
        )

        if recurrence.empty:

            raise RuntimeError(
                "M5.3C.3F received an empty disease-variant recurrence table."
            )

        required_recurrence_columns = {
            "lead_rsid",
            "variant_identity_key",
            "rsid",
            "minimum_p_value",
            "study_count",
        }

        missing_recurrence_columns = (
            required_recurrence_columns
            - set(recurrence.columns)
        )

        if missing_recurrence_columns:

            raise ValueError(
                "M5.3C.3E recurrence table missing required columns: "
                f"{sorted(missing_recurrence_columns)}"
            )

        self.m5_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        disease_bridges: dict[
            str,
            pd.DataFrame,
        ] = {}

        regulatory_bridges: dict[
            str,
            pd.DataFrame,
        ] = {}

        population_records: list[
            dict[str, Any]
        ] = []

        # ==================================================================
        # Per-population integration
        # ==================================================================

        for population in POPULATIONS:

            logger.info(
                "Processing population %s.",
                population,
            )

            ld_path = self.ld_path(
                population
            )

            moradi_path = self.moradi_bridge_path(
                population
            )

            if not ld_path.exists():

                raise FileNotFoundError(
                    "Normalized LD evidence missing for "
                    f"{population}: {ld_path}"
                )

            if not moradi_path.exists():

                raise FileNotFoundError(
                    "Moradi LD bridge missing for "
                    f"{population}: {moradi_path}"
                )

            ld = _read_parquet_compat(
                ld_path
            )

            moradi = _read_parquet_compat(
                moradi_path
            )

            # --------------------------------------------------------------
            # Disease ↔ candidate-lead LD lookup
            #
            # IMPORTANT:
            # This output contains raw rsID matches across the lead-specific
            # LD reference tables. Rows are not automatically considered an
            # LD bridge until r² thresholds are applied.
            # --------------------------------------------------------------

            disease_result = build_disease_ld_bridge(
                recurrence=recurrence,
                ld_evidence=ld,
                population=population,
            )

            disease_bridge = (
                disease_result.bridge
            )

            disease_bridges[
                population
            ] = disease_bridge

            write_parquet(
                disease_bridge,
                self.disease_ld_output(
                    population
                ),
                index=False,
            )

            # --------------------------------------------------------------
            # Disease ↔ regulatory same-proxy identity
            # --------------------------------------------------------------

            regulatory_result = (
                build_disease_regulatory_bridge(
                    disease_ld_bridge=disease_bridge,
                    regulatory_bridge=moradi,
                    population=population,
                )
            )

            regulatory_bridge = (
                regulatory_result.bridge
            )

            regulatory_bridges[
                population
            ] = regulatory_bridge

            write_parquet(
                regulatory_bridge,
                self.disease_regulatory_output(
                    population
                ),
                index=False,
            )

            # --------------------------------------------------------------
            # Population-specific threshold counts
            # --------------------------------------------------------------

            if disease_bridge.empty:

                moderate_or_high_rows = 0
                moderate_or_high_variants = 0
                high_rows = 0
                high_variants = 0

            else:

                moderate_or_high_mask = _safe_bool_series(
                    disease_bridge[
                        "passes_moderate_ld"
                    ]
                )

                high_mask = _safe_bool_series(
                    disease_bridge[
                        "passes_high_ld"
                    ]
                )

                moderate_or_high = disease_bridge.loc[
                    moderate_or_high_mask
                ]

                high = disease_bridge.loc[
                    high_mask
                ]

                moderate_or_high_rows = int(
                    len(
                        moderate_or_high
                    )
                )

                moderate_or_high_variants = (
                    _unique_variant_count(
                        moderate_or_high
                    )
                )

                high_rows = int(
                    len(
                        high
                    )
                )

                high_variants = (
                    _unique_variant_count(
                        high
                    )
                )

            population_record = {
                "population":
                    population,

                "input_disease_variant_identities":
                    disease_result.input_disease_variants,

                "disease_variant_identities_with_rsid":
                    disease_result.disease_variants_with_rsid,

                "disease_variant_identities_without_rsid":
                    disease_result.disease_variants_without_rsid,

                "raw_ld_reference_match_rows":
                    disease_result.matched_ld_rows,

                "disease_variants_with_raw_ld_reference_match":
                    (
                        disease_result
                        .matched_disease_variant_identities
                    ),

                "moderate_or_high_ld_bridge_rows":
                    moderate_or_high_rows,

                "disease_variants_with_moderate_or_high_ld_bridge":
                    moderate_or_high_variants,

                "high_ld_bridge_rows":
                    high_rows,

                "disease_variants_with_high_ld_bridge":
                    high_variants,

                "regulatory_bridge_rows":
                    regulatory_result.regulatory_rows,

                "disease_regulatory_bridge_rows":
                    regulatory_result.integrated_rows,

                "disease_variants_with_regulatory_bridge":
                    (
                        regulatory_result
                        .integrated_disease_variants
                    ),
            }

            population_records.append(
                population_record
            )

            logger.info(
                "%s: raw LD matches=%d variants=%d | "
                "r2>=%.2f variants=%d | "
                "r2>=%.2f variants=%d | "
                "disease-regulatory rows=%d.",
                population,
                disease_result.matched_ld_rows,
                disease_result.matched_disease_variant_identities,
                MODERATE_LD_R2,
                moderate_or_high_variants,
                HIGH_LD_R2,
                high_variants,
                regulatory_result.integrated_rows,
            )

        # ==================================================================
        # Cross-population synthesis
        # ==================================================================

        synthesis = synthesize_disease_ld_bridges(
            disease_bridges
        )

        write_parquet(
            synthesis,
            self.synthesis_output,
            index=False,
        )

        population_qc = pd.DataFrame(
            population_records
        )

        # ==================================================================
        # Candidate-level synthesis
        # ==================================================================

        candidate_records: list[
            dict[str, Any]
        ] = []

        candidate_leads = sorted(
            recurrence[
                "lead_rsid"
            ]
            .dropna()
            .astype("string")
            .unique()
            .tolist()
        )

        for lead in candidate_leads:

            disease_signal_rows = recurrence.loc[
                recurrence[
                    "lead_rsid"
                ].astype("string")
                .eq(
                    lead
                )
            ]

            disease_signal_count = int(
                len(
                    disease_signal_rows
                )
            )

            if synthesis.empty:

                candidate_synthesis = synthesis

            else:

                candidate_synthesis = synthesis.loc[
                    synthesis[
                        "lead_rsid"
                    ]
                    .astype("string")
                    .eq(
                        lead
                    )
                ]

            # --------------------------------------------------------------
            # Raw reference matches
            # --------------------------------------------------------------

            raw_ld_match_variant_count = (
                _unique_variant_count(
                    candidate_synthesis
                )
            )

            # --------------------------------------------------------------
            # Moderate-or-high LD
            #
            # Cross-population synthesis stores one row per lead × disease
            # variant and the number of populations satisfying r² >= 0.50.
            # --------------------------------------------------------------

            if candidate_synthesis.empty:

                moderate_or_high_ld_variant_count = 0
                high_ld_variant_count = 0

                multipop_moderate_or_high_ld_count = 0
                multipop_high_ld_count = 0

            else:

                moderate_or_high_mask = (
                    pd.to_numeric(
                        candidate_synthesis[
                            "moderate_or_high_ld_population_count"
                        ],
                        errors="coerce",
                    )
                    .fillna(0)
                    .gt(0)
                )

                high_mask = _safe_bool_series(
                    candidate_synthesis[
                        "high_ld_in_any_population"
                    ]
                )

                multipop_moderate_mask = (
                    pd.to_numeric(
                        candidate_synthesis[
                            "moderate_or_high_ld_population_count"
                        ],
                        errors="coerce",
                    )
                    .fillna(0)
                    .ge(2)
                )

                multipop_high_mask = _safe_bool_series(
                    candidate_synthesis[
                        "high_ld_in_multiple_populations"
                    ]
                )

                moderate_or_high_ld_variant_count = (
                    _unique_variant_count(
                        candidate_synthesis.loc[
                            moderate_or_high_mask
                        ]
                    )
                )

                high_ld_variant_count = (
                    _unique_variant_count(
                        candidate_synthesis.loc[
                            high_mask
                        ]
                    )
                )

                multipop_moderate_or_high_ld_count = (
                    _unique_variant_count(
                        candidate_synthesis.loc[
                            multipop_moderate_mask
                        ]
                    )
                )

                multipop_high_ld_count = (
                    _unique_variant_count(
                        candidate_synthesis.loc[
                            multipop_high_mask
                        ]
                    )
                )

            # --------------------------------------------------------------
            # Regulatory shared-proxy integration
            # --------------------------------------------------------------

            regulatory_variant_set: set[
                str
            ] = set()

            regulatory_row_count = 0

            regulatory_populations: list[
                str
            ] = []

            for population in POPULATIONS:

                frame = regulatory_bridges[
                    population
                ]

                if frame.empty:
                    continue

                if "lead_rsid" not in frame.columns:
                    continue

                current = frame.loc[
                    frame[
                        "lead_rsid"
                    ]
                    .astype("string")
                    .eq(
                        lead
                    )
                ]

                if current.empty:
                    continue

                regulatory_row_count += int(
                    len(
                        current
                    )
                )

                if (
                    "disease_variant_rsid"
                    in current.columns
                ):

                    regulatory_variant_set.update(
                        current[
                            "disease_variant_rsid"
                        ]
                        .dropna()
                        .astype("string")
                        .tolist()
                    )

                regulatory_populations.append(
                    population
                )

            # --------------------------------------------------------------
            # Candidate integration classification
            #
            # Precedence:
            #   1. disease-regulatory shared proxy
            #   2. high LD disease bridge
            #   3. moderate LD disease bridge
            #   4. raw lookup match below threshold
            #   5. regional signal without LD lookup bridge
            # --------------------------------------------------------------

            if regulatory_variant_set:

                integration_class = (
                    "DISEASE_REGULATORY_SHARED_LD_PROXY"
                )

            elif high_ld_variant_count > 0:

                integration_class = (
                    "HIGH_LD_DISEASE_BRIDGE"
                )

            elif moderate_or_high_ld_variant_count > 0:

                integration_class = (
                    "MODERATE_LD_DISEASE_BRIDGE"
                )

            elif raw_ld_match_variant_count > 0:

                integration_class = (
                    "LD_REFERENCE_MATCH_BELOW_THRESHOLD"
                )

            elif disease_signal_count > 0:

                integration_class = (
                    "REGIONAL_SIGNAL_NO_LD_BRIDGE"
                )

            else:

                integration_class = (
                    "NO_REGIONAL_DISEASE_SIGNAL"
                )

            candidate_records.append(
                {
                    "lead_rsid":
                        lead,

                    "suggestive_disease_variant_identities":
                        disease_signal_count,

                    "disease_variants_with_raw_ld_reference_match":
                        raw_ld_match_variant_count,

                    "disease_variants_with_moderate_or_high_ld_bridge":
                        moderate_or_high_ld_variant_count,

                    "disease_variants_with_high_ld_bridge":
                        high_ld_variant_count,

                    "disease_variants_with_moderate_or_high_ld_in_multiple_populations":
                        multipop_moderate_or_high_ld_count,

                    "disease_variants_with_high_ld_in_multiple_populations":
                        multipop_high_ld_count,

                    "disease_variants_with_regulatory_shared_proxy":
                        int(
                            len(
                                regulatory_variant_set
                            )
                        ),

                    "disease_regulatory_bridge_rows":
                        regulatory_row_count,

                    "regulatory_bridge_populations":
                        "|".join(
                            sorted(
                                set(
                                    regulatory_populations
                                )
                            )
                        ),

                    "integration_class":
                        integration_class,

                    "colocalization_assessed":
                        False,

                    "fine_mapping_performed":
                        False,

                    "causal_inference_performed":
                        False,
                }
            )

        candidate_summary = pd.DataFrame(
            candidate_records
        )

        write_parquet(
            candidate_summary,
            self.candidate_summary_output,
            index=False,
        )

        # ==================================================================
        # Global threshold summaries
        # ==================================================================

        raw_ld_reference_match_rows = int(
            population_qc[
                "raw_ld_reference_match_rows"
            ].sum()
        )

        raw_ld_reference_match_variants = (
            _unique_lead_variant_pair_count(
                synthesis
            )
        )

        if synthesis.empty:

            moderate_or_high_ld_disease_variants = 0

            high_ld_disease_variants = 0

            multi_population_moderate_or_high_ld_disease_variants = 0

            multi_population_high_ld_disease_variants = 0

        else:

            moderate_or_high_mask = (
                pd.to_numeric(
                    synthesis[
                        "moderate_or_high_ld_population_count"
                    ],
                    errors="coerce",
                )
                .fillna(0)
                .gt(0)
            )

            high_mask = _safe_bool_series(
                synthesis[
                    "high_ld_in_any_population"
                ]
            )

            multipop_moderate_mask = (
                pd.to_numeric(
                    synthesis[
                        "moderate_or_high_ld_population_count"
                    ],
                    errors="coerce",
                )
                .fillna(0)
                .ge(2)
            )

            multipop_high_mask = _safe_bool_series(
                synthesis[
                    "high_ld_in_multiple_populations"
                ]
            )

            moderate_or_high_ld_disease_variants = (
                _unique_lead_variant_pair_count(
                    synthesis.loc[
                        moderate_or_high_mask
                    ]
                )
            )

            high_ld_disease_variants = (
                _unique_lead_variant_pair_count(
                    synthesis.loc[
                        high_mask
                    ]
                )
            )

            multi_population_moderate_or_high_ld_disease_variants = (
                _unique_lead_variant_pair_count(
                    synthesis.loc[
                        multipop_moderate_mask
                    ]
                )
            )

            multi_population_high_ld_disease_variants = (
                _unique_lead_variant_pair_count(
                    synthesis.loc[
                        multipop_high_mask
                    ]
                )
            )

        # ==================================================================
        # Candidate-level global counts
        # ==================================================================

        candidate_leads_with_any_raw_ld_reference_match = int(
            candidate_summary[
                "disease_variants_with_raw_ld_reference_match"
            ]
            .gt(0)
            .sum()
        )

        candidate_leads_with_moderate_or_high_ld_bridge = int(
            candidate_summary[
                "disease_variants_with_moderate_or_high_ld_bridge"
            ]
            .gt(0)
            .sum()
        )

        candidate_leads_with_high_ld_bridge = int(
            candidate_summary[
                "disease_variants_with_high_ld_bridge"
            ]
            .gt(0)
            .sum()
        )

        candidate_leads_with_disease_regulatory_shared_proxy = int(
            candidate_summary[
                "disease_variants_with_regulatory_shared_proxy"
            ]
            .gt(0)
            .sum()
        )

        # ==================================================================
        # QC report
        # ==================================================================

        report = {
            "milestone":
                "M5.3C.3F",

            "stage":
                "disease_signal_ld_bridge_integration",

            "policy": {
                "disease_source":
                    "M5.3C.3E suggestive variant recurrence",

                "ld_populations":
                    list(
                        POPULATIONS
                    ),

                "ld_source_build":
                    "GRCh37/hg19",

                "disease_gwas_build":
                    "GRCh38",

                "cross_build_coordinate_join_performed":
                    False,

                "cross_build_matching_method":
                    "canonical_rsid_only",

                "silent_liftover_performed":
                    False,

                "raw_ld_reference_match_interpreted_as_ld_bridge":
                    False,

                "moderate_ld_threshold":
                    MODERATE_LD_R2,

                "high_ld_threshold":
                    HIGH_LD_R2,

                "moderate_or_high_ld_definition":
                    "r2 >= 0.50",

                "high_ld_definition":
                    "r2 >= 0.80",

                "physical_distance_interpreted_as_ld":
                    False,

                "population_recurrence_interpreted_as_replication":
                    False,

                "same_proxy_identity_interpreted_as_colocalization":
                    False,

                "colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "input_suggestive_variant_identities":
                    int(
                        len(
                            recurrence
                        )
                    ),

                "candidate_leads":
                    int(
                        len(
                            candidate_leads
                        )
                    ),

                "populations":
                    int(
                        len(
                            POPULATIONS
                        )
                    ),

                "raw_ld_reference_match_rows":
                    raw_ld_reference_match_rows,

                "unique_cross_population_raw_ld_matches":
                    raw_ld_reference_match_variants,

                "moderate_or_high_ld_disease_variants":
                    moderate_or_high_ld_disease_variants,

                "high_ld_disease_variants":
                    high_ld_disease_variants,

                "multi_population_moderate_or_high_ld_disease_variants":
                    (
                        multi_population_moderate_or_high_ld_disease_variants
                    ),

                "multi_population_high_ld_disease_variants":
                    multi_population_high_ld_disease_variants,

                "disease_regulatory_bridge_rows":
                    int(
                        population_qc[
                            "disease_regulatory_bridge_rows"
                        ].sum()
                    ),

                "candidate_leads_with_any_raw_ld_reference_match":
                    candidate_leads_with_any_raw_ld_reference_match,

                "candidate_leads_with_moderate_or_high_ld_bridge":
                    candidate_leads_with_moderate_or_high_ld_bridge,

                "candidate_leads_with_high_ld_bridge":
                    candidate_leads_with_high_ld_bridge,

                "candidate_leads_with_disease_regulatory_shared_proxy":
                    (
                        candidate_leads_with_disease_regulatory_shared_proxy
                    ),
            },

            "population_qc":
                population_records,

            "candidate_summary":
                candidate_summary.to_dict(
                    orient="records"
                ),

            "outputs": {
                "cross_population_synthesis":
                    str(
                        self.synthesis_output
                    ),

                "candidate_summary":
                    str(
                        self.candidate_summary_output
                    ),

                "per_population_disease_ld":
                    {
                        population:
                            str(
                                self.disease_ld_output(
                                    population
                                )
                            )
                        for population
                        in POPULATIONS
                    },

                "per_population_disease_regulatory":
                    {
                        population:
                            str(
                                self.disease_regulatory_output(
                                    population
                                )
                            )
                        for population
                        in POPULATIONS
                    },
            },
        }

        with self.qc_output.open(
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
            "M5.3C.3F complete."
        )

        logger.info(
            "Raw LD reference matches: %d rows / %d unique "
            "lead-variant pairs.",
            raw_ld_reference_match_rows,
            raw_ld_reference_match_variants,
        )

        logger.info(
            "Moderate-or-high LD disease variants (r2>=%.2f): %d.",
            MODERATE_LD_R2,
            moderate_or_high_ld_disease_variants,
        )

        logger.info(
            "High-LD disease variants (r2>=%.2f): %d.",
            HIGH_LD_R2,
            high_ld_disease_variants,
        )

        logger.info(
            "Multi-population moderate-or-high LD variants: %d | "
            "multi-population high-LD variants: %d.",
            multi_population_moderate_or_high_ld_disease_variants,
            multi_population_high_ld_disease_variants,
        )

        logger.info(
            "Candidate leads with raw LD reference match: %d/%d.",
            candidate_leads_with_any_raw_ld_reference_match,
            len(
                candidate_leads
            ),
        )

        logger.info(
            "Candidate leads with moderate-or-high LD bridge: %d/%d.",
            candidate_leads_with_moderate_or_high_ld_bridge,
            len(
                candidate_leads
            ),
        )

        logger.info(
            "Candidate leads with high-LD bridge: %d/%d.",
            candidate_leads_with_high_ld_bridge,
            len(
                candidate_leads
            ),
        )

        logger.info(
            "Candidate leads with disease-regulatory shared proxy: %d/%d.",
            candidate_leads_with_disease_regulatory_shared_proxy,
            len(
                candidate_leads
            ),
        )

        logger.info(
            "Cross-population synthesis: %s",
            self.synthesis_output,
        )

        logger.info(
            "Candidate summary: %s",
            self.candidate_summary_output,
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report