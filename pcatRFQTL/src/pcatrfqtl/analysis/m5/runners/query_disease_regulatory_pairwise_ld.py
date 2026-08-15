"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/query_disease_regulatory_pairwise_ld.py

Description:
    Full runner for M5.3C.3G disease-proxy ↔ regulatory-proxy pairwise LD.

    The runner:
        1. loads the validated M5.3C.3G pairwise candidate scaffold;
        2. queries pairwise LD for EAS, EUR, and SAS;
        3. uses a persistent SQLite cache through CachedPairwiseLDProvider;
        4. preserves unavailable states explicitly;
        5. writes one pairwise LD artifact per population;
        6. synthesizes pairwise evidence across reference populations;
        7. writes detailed QC including cache diagnostics and LD thresholds.

    Reference framework:
        Provider:
            NIH LDlink LDpair

        Reference panel:
            1000 Genomes Project

        Genome build:
            GRCh37

        Populations:
            EAS
            EUR
            SAS

    LD interpretation:
        HIGH:
            r² >= 0.80

        MODERATE:
            0.50 <= r² < 0.80

        LOW:
            r² < 0.50

        UNAVAILABLE:
            no usable pairwise r² returned.

    Scientific safeguards:
        - Missing LD is not represented as r² = 0.
        - Physical proximity is not interpreted as LD.
        - Pairwise LD does not prove colocalization.
        - Cross-population recurrence is sensitivity evidence, not
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
import os
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq
from dotenv import load_dotenv

from pcatrfqtl.analysis.m5.cached_pairwise_ld_provider import (
    CachedPairwiseLDProvider,
)
from pcatrfqtl.analysis.m5.disease_regulatory_pairwise_ld import (
    HIGH_LD_R2,
    MODERATE_LD_R2,
    classify_pairwise_ld,
    synthesize_pairwise_ld,
)
from pcatrfqtl.analysis.m5.ldlink_pairwise_provider import (
    LDlinkPairwiseProvider,
    PROVIDER_NAME,
)
from pcatrfqtl.analysis.m5.pairwise_ld_cache import (
    PairwiseLDCache,
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
    """Read M5 Parquet artifact without pandas metadata reconstruction."""

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


def _safe_bool(
    value: Any,
) -> bool:
    """Convert nullable value to bool."""

    if pd.isna(
        value
    ):
        return False

    return bool(
        value
    )


def _status_counts(
    dataframe: pd.DataFrame,
) -> dict[str, int]:
    """Return pairwise query status counts."""

    if dataframe.empty:

        return {}

    if "query_status" not in dataframe.columns:

        return {}

    counts = (
        dataframe[
            "query_status"
        ]
        .astype("string")
        .fillna("MISSING")
        .value_counts(
            dropna=False
        )
    )

    return {
        str(
            key
        ):
            int(
                value
            )
        for key, value in counts.items()
    }


def _ld_class_counts(
    dataframe: pd.DataFrame,
) -> dict[str, int]:
    """Return LD-class counts."""

    if dataframe.empty:

        return {}

    if "ld_class" not in dataframe.columns:

        return {}

    counts = (
        dataframe[
            "ld_class"
        ]
        .astype("string")
        .fillna("MISSING")
        .value_counts(
            dropna=False
        )
    )

    return {
        str(
            key
        ):
            int(
                value
            )
        for key, value in counts.items()
    }


# ============================================================================
# Runner
# ============================================================================


class M53C3GPairwiseLDRunner:
    """Execute full M5.3C.3G pairwise LD integration."""

    def __init__(
        self,
        *,
        scaffold_path: str | Path,
        cache_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
        env_path: str | Path,
    ) -> None:

        self.scaffold_path = Path(
            scaffold_path
        )

        self.cache_path = Path(
            cache_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

        self.env_path = Path(
            env_path
        )

    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------

    def population_output(
        self,
        population: str,
    ) -> Path:

        return (
            self.output_directory
            / (
                "disease_regulatory_pairwise_ld_"
                f"{population}.parquet"
            )
        )

    @property
    def synthesis_output(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "disease_regulatory_pairwise_ld_synthesis.parquet"
        )

    @property
    def qc_output(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_3c_3g_disease_regulatory_pairwise_ld.json"
        )

    # ------------------------------------------------------------------
    # Environment
    # ------------------------------------------------------------------

    def _load_token(
        self,
    ) -> str:
        """Load LDLINK_TOKEN from project .env."""

        if not self.env_path.exists():

            raise FileNotFoundError(
                f".env file not found: {self.env_path}"
            )

        load_dotenv(
            dotenv_path=self.env_path,
            override=False,
        )

        token = os.environ.get(
            "LDLINK_TOKEN"
        )

        if token is None:

            raise RuntimeError(
                "LDLINK_TOKEN was not found."
            )

        token = token.strip()

        if not token:

            raise RuntimeError(
                "LDLINK_TOKEN is empty."
            )

        return token

    # ------------------------------------------------------------------
    # Main
    # ------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Run M5.3C.3G."""

        token = self._load_token()

        scaffold = _read_parquet_compat(
            self.scaffold_path
        )

        if scaffold.empty:

            raise RuntimeError(
                "M5.3C.3G candidate scaffold is empty."
            )

        required_columns = {
            "lead_rsid",
            "pair_id",
            "disease_proxy_rsid",
            "regulatory_proxy_rsid",
            "minimum_p_value",
            "study_count",
            "genome_wide_in_any_study",
            "regulatory_feature_ids",
            "regulatory_source_partitions",
        }

        missing_columns = (
            required_columns
            - set(
                scaffold.columns
            )
        )

        if missing_columns:

            raise ValueError(
                "Pairwise candidate scaffold missing required columns: "
                f"{sorted(missing_columns)}"
            )

        # ------------------------------------------------------------------
        # Validate scaffold cardinality.
        # ------------------------------------------------------------------

        candidate_pairs = int(
            scaffold[
                [
                    "disease_proxy_rsid",
                    "regulatory_proxy_rsid",
                ]
            ]
            .drop_duplicates()
            .shape[
                0
            ]
        )

        disease_variants = int(
            scaffold[
                "disease_proxy_rsid"
            ].nunique()
        )

        regulatory_variants = int(
            scaffold[
                "regulatory_proxy_rsid"
            ].nunique()
        )

        if candidate_pairs != len(
            scaffold
        ):

            raise RuntimeError(
                "Pairwise scaffold contains duplicate disease-regulatory "
                "variant combinations: "
                f"rows={len(scaffold)} "
                f"unique_pairs={candidate_pairs}"
            )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        logger.info(
            "Loaded pairwise scaffold: disease=%d regulatory=%d pairs=%d.",
            disease_variants,
            regulatory_variants,
            candidate_pairs,
        )

        expected_pair_population_evaluations = int(
            candidate_pairs
            * len(
                POPULATIONS
            )
        )

        population_results: dict[
            str,
            pd.DataFrame,
        ] = {}

        population_qc: list[
            dict[str, Any]
        ] = []

        # ==================================================================
        # Persistent cache lifecycle
        # ==================================================================

        with PairwiseLDCache(
            self.cache_path
        ) as cache:

            cache_before = int(
                cache.count()
            )

            logger.info(
                "Pairwise LD cache before run: %d records.",
                cache_before,
            )

            raw_provider = LDlinkPairwiseProvider(
                token=token,
            )

            cached_provider = (
                CachedPairwiseLDProvider(
                    provider=raw_provider,
                    cache=cache,
                    provider_name=PROVIDER_NAME,
                )
            )

            # ==============================================================
            # Population loop
            # ==============================================================

            for population in POPULATIONS:

                logger.info(
                    "Processing M5.3C.3G population %s.",
                    population,
                )

                diagnostic_before = (
                    cached_provider
                    .diagnostics()
                    .copy()
                )

                pairs = [
                    (
                        str(
                            row[
                                "disease_proxy_rsid"
                            ]
                        ),
                        str(
                            row[
                                "regulatory_proxy_rsid"
                            ]
                        ),
                    )
                    for _, row in scaffold.iterrows()
                ]

                pairwise_records = (
                    cached_provider.query_pairs(
                        pairs=pairs,
                        population=population,
                    )
                )

                if len(
                    pairwise_records
                ) != candidate_pairs:

                    raise RuntimeError(
                        "Pairwise provider cardinality mismatch for "
                        f"{population}: expected={candidate_pairs} "
                        f"received={len(pairwise_records)}"
                    )

                # ----------------------------------------------------------
                # Attach provider results to scaffold rows in exact order.
                # ----------------------------------------------------------

                records: list[
                    dict[str, Any]
                ] = []

                for (
                    (_, scaffold_row),
                    result,
                ) in zip(
                    scaffold.iterrows(),
                    pairwise_records,
                    strict=True,
                ):

                    r2 = (
                        None
                        if result.r2 is None
                        else float(
                            result.r2
                        )
                    )

                    d_prime = (
                        None
                        if result.d_prime is None
                        else float(
                            result.d_prime
                        )
                    )

                    ld_class = (
                        classify_pairwise_ld(
                            r2
                        )
                    )

                    records.append(
                        {
                            "lead_rsid":
                                str(
                                    scaffold_row[
                                        "lead_rsid"
                                    ]
                                ),

                            "pair_id":
                                str(
                                    scaffold_row[
                                        "pair_id"
                                    ]
                                ),

                            "disease_proxy_rsid":
                                str(
                                    scaffold_row[
                                        "disease_proxy_rsid"
                                    ]
                                ),

                            "regulatory_proxy_rsid":
                                str(
                                    scaffold_row[
                                        "regulatory_proxy_rsid"
                                    ]
                                ),

                            "population":
                                population,

                            "r2":
                                r2,

                            "d_prime":
                                d_prime,

                            "query_status":
                                str(
                                    result.status
                                ),

                            "provider":
                                str(
                                    result.provider
                                ),

                            "reason":
                                result.reason,

                            "ld_class":
                                ld_class,

                            "passes_moderate_or_high_ld":
                                bool(
                                    r2 is not None
                                    and
                                    r2 >= MODERATE_LD_R2
                                ),

                            "passes_high_ld":
                                bool(
                                    r2 is not None
                                    and
                                    r2 >= HIGH_LD_R2
                                ),

                            "minimum_disease_p_value":
                                (
                                    None
                                    if pd.isna(
                                        scaffold_row[
                                            "minimum_p_value"
                                        ]
                                    )
                                    else float(
                                        scaffold_row[
                                            "minimum_p_value"
                                        ]
                                    )
                                ),

                            "disease_study_count":
                                (
                                    None
                                    if pd.isna(
                                        scaffold_row[
                                            "study_count"
                                        ]
                                    )
                                    else int(
                                        scaffold_row[
                                            "study_count"
                                        ]
                                    )
                                ),

                            "genome_wide_in_any_disease_study":
                                _safe_bool(
                                    scaffold_row[
                                        "genome_wide_in_any_study"
                                    ]
                                ),

                            "regulatory_feature_ids":
                                scaffold_row[
                                    "regulatory_feature_ids"
                                ],

                            "regulatory_source_partitions":
                                scaffold_row[
                                    "regulatory_source_partitions"
                                ],
                        }
                    )

                frame = pd.DataFrame(
                    records
                )

                # ----------------------------------------------------------
                # Cardinality safeguard.
                # ----------------------------------------------------------

                if len(
                    frame
                ) != candidate_pairs:

                    raise RuntimeError(
                        "Population output cardinality mismatch for "
                        f"{population}: "
                        f"expected={candidate_pairs} "
                        f"actual={len(frame)}"
                    )

                population_results[
                    population
                ] = frame

                write_parquet(
                    frame,
                    self.population_output(
                        population
                    ),
                    index=False,
                )

                # ----------------------------------------------------------
                # Population-specific QC.
                # ----------------------------------------------------------

                numeric_r2 = pd.to_numeric(
                    frame[
                        "r2"
                    ],
                    errors="coerce",
                )

                numeric_r2_count = int(
                    numeric_r2.notna().sum()
                )

                unavailable_count = int(
                    len(
                        frame
                    )
                    - numeric_r2_count
                )

                moderate_or_high_count = int(
                    numeric_r2.ge(
                        MODERATE_LD_R2
                    ).sum()
                )

                high_count = int(
                    numeric_r2.ge(
                        HIGH_LD_R2
                    ).sum()
                )

                diagnostic_after = (
                    cached_provider
                    .diagnostics()
                    .copy()
                )

                diagnostic_delta = {
                    key:
                        int(
                            diagnostic_after[
                                key
                            ]
                            -
                            diagnostic_before[
                                key
                            ]
                        )
                    for key
                    in diagnostic_after
                }

                population_record = {
                    "population":
                        population,

                    "candidate_pairs":
                        candidate_pairs,

                    "output_rows":
                        int(
                            len(
                                frame
                            )
                        ),

                    "cardinality_preserved":
                        bool(
                            len(
                                frame
                            )
                            == candidate_pairs
                        ),

                    "pairs_with_numeric_r2":
                        numeric_r2_count,

                    "unavailable_or_failed_pairs":
                        unavailable_count,

                    "moderate_or_high_ld_pairs":
                        moderate_or_high_count,

                    "high_ld_pairs":
                        high_count,

                    "status_counts":
                        _status_counts(
                            frame
                        ),

                    "ld_class_counts":
                        _ld_class_counts(
                            frame
                        ),

                    "cache_diagnostics":
                        diagnostic_delta,
                }

                population_qc.append(
                    population_record
                )

                logger.info(
                    "%s: numeric r2=%d/%d | "
                    "r2>=%.2f=%d | r2>=%.2f=%d | "
                    "cache hits=%d misses=%d remote_pairs=%d.",
                    population,
                    numeric_r2_count,
                    candidate_pairs,
                    MODERATE_LD_R2,
                    moderate_or_high_count,
                    HIGH_LD_R2,
                    high_count,
                    diagnostic_delta[
                        "cache_hits"
                    ],
                    diagnostic_delta[
                        "cache_misses"
                    ],
                    diagnostic_delta[
                        "remote_pairs"
                    ],
                )

            # ==============================================================
            # Cross-population synthesis
            # ==============================================================

            synthesis_result = synthesize_pairwise_ld(
                population_results
            )

            synthesis = (
                synthesis_result.dataframe
            )

            if len(
                synthesis
            ) != candidate_pairs:

                raise RuntimeError(
                    "Cross-population synthesis cardinality mismatch: "
                    f"expected={candidate_pairs} "
                    f"actual={len(synthesis)}"
                )

            write_parquet(
                synthesis,
                self.synthesis_output,
                index=False,
            )

            cache_after = int(
                cache.count()
            )

            final_cache_diagnostics = (
                cached_provider
                .diagnostics()
            )

            cache_by_status = (
                cache.count_by_status()
            )

            cache_by_population = (
                cache.count_by_population()
            )

            # ==============================================================
            # Cross-population result summaries
            # ==============================================================

            moderate_any = int(
                synthesis[
                    "moderate_or_high_ld_in_any_population"
                ]
                .fillna(
                    False
                )
                .sum()
                if not synthesis.empty
                else 0
            )

            high_any = int(
                synthesis[
                    "high_ld_in_any_population"
                ]
                .fillna(
                    False
                )
                .sum()
                if not synthesis.empty
                else 0
            )

            moderate_multi = int(
                synthesis[
                    "moderate_or_high_ld_in_multiple_populations"
                ]
                .fillna(
                    False
                )
                .sum()
                if not synthesis.empty
                else 0
            )

            high_multi = int(
                synthesis[
                    "high_ld_in_multiple_populations"
                ]
                .fillna(
                    False
                )
                .sum()
                if not synthesis.empty
                else 0
            )

            complete_numeric_all_populations = int(
                synthesis[
                    "unavailable_population_count"
                ]
                .eq(
                    0
                )
                .sum()
                if not synthesis.empty
                else 0
            )

            any_unavailable = int(
                synthesis[
                    "unavailable_population_count"
                ]
                .gt(
                    0
                )
                .sum()
                if not synthesis.empty
                else 0
            )

            # ==============================================================
            # Report
            # ==============================================================

            report = {
                "milestone":
                    "M5.3C.3G",

                "stage":
                    "disease_proxy_regulatory_proxy_pairwise_ld",

                "policy": {
                    "candidate_lead":
                        "rs10216902",

                    "provider":
                        PROVIDER_NAME,

                    "reference_panel":
                        "1000 Genomes Project",

                    "genome_build":
                        "GRCh37",

                    "populations":
                        list(
                            POPULATIONS
                        ),

                    "moderate_ld_threshold":
                        MODERATE_LD_R2,

                    "high_ld_threshold":
                        HIGH_LD_R2,

                    "same_proxy_pairs_included":
                        False,

                    "pairwise_ld_directly_queried":
                        True,

                    "lead_centered_ld_reused_as_pairwise_ld":
                        False,

                    "missing_ld_interpreted_as_zero":
                        False,

                    "persistent_cache_used":
                        True,

                    "physical_distance_interpreted_as_ld":
                        False,

                    "cross_population_consistency_interpreted_as_replication":
                        False,

                    "colocalization_performed":
                        False,

                    "fine_mapping_performed":
                        False,

                    "causal_inference_performed":
                        False,
                },

                "summary": {
                    "disease_variants":
                        disease_variants,

                    "regulatory_variants":
                        regulatory_variants,

                    "candidate_variant_pairs":
                        candidate_pairs,

                    "populations":
                        int(
                            len(
                                POPULATIONS
                            )
                        ),

                    "expected_pair_population_evaluations":
                        expected_pair_population_evaluations,

                    "population_output_rows":
                        int(
                            sum(
                                record[
                                    "output_rows"
                                ]
                                for record
                                in population_qc
                            )
                        ),

                    "synthesized_variant_pairs":
                        int(
                            len(
                                synthesis
                            )
                        ),

                    "pairs_numeric_in_all_populations":
                        complete_numeric_all_populations,

                    "pairs_unavailable_in_at_least_one_population":
                        any_unavailable,

                    "moderate_or_high_ld_variant_pairs":
                        moderate_any,

                    "high_ld_variant_pairs":
                        high_any,

                    "multi_population_moderate_or_high_ld_pairs":
                        moderate_multi,

                    "multi_population_high_ld_pairs":
                        high_multi,

                    "cache_records_before":
                        cache_before,

                    "cache_records_after":
                        cache_after,

                    "cache_records_added":
                        int(
                            cache_after
                            - cache_before
                        ),
                },

                "cache_diagnostics":
                    final_cache_diagnostics,

                "cache_by_status":
                    cache_by_status,

                "cache_by_population":
                    cache_by_population,

                "population_qc":
                    population_qc,

                "outputs": {
                    "pairwise_synthesis":
                        str(
                            self.synthesis_output
                        ),

                    "per_population":
                        {
                            population:
                                str(
                                    self.population_output(
                                        population
                                    )
                                )
                            for population
                            in POPULATIONS
                        },

                    "cache":
                        str(
                            self.cache_path
                        ),
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

        logger.info(
            "M5.3C.3G complete."
        )

        logger.info(
            "Pairs synthesized: %d/%d.",
            report[
                "summary"
            ][
                "synthesized_variant_pairs"
            ],
            candidate_pairs,
        )

        logger.info(
            "Moderate/high LD pairs: %d | high-LD pairs: %d.",
            report[
                "summary"
            ][
                "moderate_or_high_ld_variant_pairs"
            ],
            report[
                "summary"
            ][
                "high_ld_variant_pairs"
            ],
        )

        logger.info(
            "Cache: before=%d after=%d added=%d.",
            report[
                "summary"
            ][
                "cache_records_before"
            ],
            report[
                "summary"
            ][
                "cache_records_after"
            ],
            report[
                "summary"
            ][
                "cache_records_added"
            ],
        )

        logger.info(
            "QC output: %s",
            self.qc_output,
        )

        return report