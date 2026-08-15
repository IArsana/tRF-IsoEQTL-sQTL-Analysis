"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/standardize_gwas_loci.py

Description:
    Runner for M5.3C.3D GWAS locus standardization and descriptive
    disease-signal quality control.

    Inputs:
        - M5.3C.3C locus retrieval manifest.
        - Retrieved harmonised GWAS locus Parquet files.

    Outputs:
        - Unified standardized GWAS locus dataset.
        - Study × candidate disease-signal summary.
        - M5.3C.3D QC report.

    Source-level repeated variant representations are preserved and classified
    as exact duplicates or conflicting representations. No rows are removed.

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

from pcatrfqtl.analysis.m5.gwas_locus_standardization import (
    GENOME_WIDE_SIGNIFICANCE,
    SUGGESTIVE_SIGNIFICANCE,
    standardize_gwas_locus,
    summarize_disease_signal,
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
# Parquet reader
# ============================================================================


def _read_parquet_compat(
    path: Path,
) -> pd.DataFrame:
    """Read M5 Parquet artifacts while ignoring pandas metadata."""

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
# Runner
# ============================================================================


class M53C3DGWASLocusStandardizationRunner:
    """Execute M5.3C.3D."""

    def __init__(
        self,
        *,
        retrieval_manifest_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.retrieval_manifest_path = Path(
            retrieval_manifest_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    # ------------------------------------------------------------------
    # Output paths
    # ------------------------------------------------------------------

    @property
    def standardized_output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "standardized"
            / "gwas_locus_standardized.parquet"
        )

    @property
    def disease_signal_output_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "gwas_locus_disease_signal_summary.parquet"
        )

    @property
    def qc_output_path(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_3c_3d_gwas_locus_standardization.json"
        )

    # ------------------------------------------------------------------
    # Main
    # ------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Run GWAS locus standardization and disease-signal QC."""

        manifest = _read_parquet_compat(
            self.retrieval_manifest_path
        )

        required = {
            "study_accession",
            "lead_rsid",
            "chromosome",
            "lead_position",
            "status",
            "variant_rows",
            "output_path",
        }

        missing = (
            required
            - set(
                manifest.columns
            )
        )

        if missing:

            raise ValueError(
                "Retrieval manifest missing required columns: "
                f"{sorted(missing)}"
            )

        retrieved = manifest.loc[
            manifest[
                "status"
            ].eq(
                "LOCUS_RETRIEVED"
            )
        ].copy()

        if retrieved.empty:

            raise RuntimeError(
                "M5.3C.3D found zero retrieved loci."
            )

        self.standardized_output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        standardized_frames: list[
            pd.DataFrame
        ] = []

        signal_records: list[
            dict[str, Any]
        ] = []

        locus_qc_records: list[
            dict[str, Any]
        ] = []

        # ==================================================================
        # Per-locus processing
        # ==================================================================

        for _, row in retrieved.iterrows():

            accession = str(
                row[
                    "study_accession"
                ]
            )

            lead_rsid = str(
                row[
                    "lead_rsid"
                ]
            )

            chromosome = str(
                row[
                    "chromosome"
                ]
            )

            lead_position = int(
                row[
                    "lead_position"
                ]
            )

            source_path = Path(
                str(
                    row[
                        "output_path"
                    ]
                )
            )

            logger.info(
                "Standardizing %s %s.",
                accession,
                lead_rsid,
            )

            source = _read_parquet_compat(
                source_path
            )

            result = standardize_gwas_locus(
                source,
                study_accession=accession,
                lead_rsid=lead_rsid,
                lead_position=lead_position,
                expected_chromosome=chromosome,
            )

            standardized = (
                result.dataframe
            )

            standardized_frames.append(
                standardized
            )

            signal = summarize_disease_signal(
                standardized,
                study_accession=accession,
                lead_rsid=lead_rsid,
            )

            signal_records.append(
                signal
            )

            manifest_rows = int(
                row[
                    "variant_rows"
                ]
            )

            locus_qc_records.append(
                {
                    "study_accession":
                        accession,

                    "lead_rsid":
                        lead_rsid,

                    "manifest_rows":
                        manifest_rows,

                    "input_rows":
                        result.input_rows,

                    "output_rows":
                        result.output_rows,

                    "cardinality_matches_manifest":
                        (
                            manifest_rows
                            == result.input_rows
                            == result.output_rows
                        ),

                    "malformed_coordinate_rows":
                        result.malformed_coordinate_rows,

                    "missing_p_value_rows":
                        result.missing_p_value_rows,

                    "missing_beta_rows":
                        result.missing_beta_rows,

                    "missing_standard_error_rows":
                        result.missing_standard_error_rows,

                    "missing_effect_allele_frequency_rows":
                        result.missing_effect_allele_frequency_rows,

                    "repeated_variant_identity_rows":
                        result.repeated_variant_identity_rows,

                    "repeated_variant_identity_keys":
                        result.repeated_variant_identity_keys,

                    "exact_duplicate_rows":
                        result.exact_duplicate_rows,

                    "exact_duplicate_keys":
                        result.exact_duplicate_keys,

                    "conflicting_variant_representation_rows":
                        result.conflicting_variant_representation_rows,

                    "conflicting_variant_representation_keys":
                        result.conflicting_variant_representation_keys,

                    "lead_present":
                        signal[
                            "lead_present"
                        ],

                    "lead_p_value":
                        signal[
                            "lead_p_value"
                        ],

                    "minimum_p_value":
                        signal[
                            "minimum_p_value"
                        ],

                    "genome_wide_significant_variants":
                        signal[
                            "genome_wide_significant_variants"
                        ],

                    "suggestive_significant_variants":
                        signal[
                            "suggestive_significant_variants"
                        ],

                    "suggestive_only_variants":
                        signal[
                            "suggestive_only_variants"
                        ],
                }
            )

        # ==================================================================
        # Combine
        # ==================================================================

        combined = pd.concat(
            standardized_frames,
            ignore_index=True,
        )

        signals = pd.DataFrame(
            signal_records
        )

        locus_qc = pd.DataFrame(
            locus_qc_records
        )

        # ==================================================================
        # Cardinality guards
        # ==================================================================

        manifest_total_rows = int(
            pd.to_numeric(
                retrieved[
                    "variant_rows"
                ],
                errors="coerce",
            )
            .fillna(
                0
            )
            .sum()
        )

        if len(
            combined
        ) != manifest_total_rows:

            raise RuntimeError(
                "M5.3C.3D total cardinality mismatch: "
                f"manifest={manifest_total_rows}, "
                f"standardized={len(combined)}"
            )

        cardinality_failures = int(
            (
                ~locus_qc[
                    "cardinality_matches_manifest"
                ]
            ).sum()
        )

        if cardinality_failures:

            raise RuntimeError(
                "M5.3C.3D detected locus cardinality failures: "
                f"{cardinality_failures}"
            )

        # ==================================================================
        # Persist
        # ==================================================================

        write_parquet(
            combined,
            self.standardized_output_path,
            index=False,
        )

        write_parquet(
            signals,
            self.disease_signal_output_path,
            index=False,
        )

        # ==================================================================
        # Global QC
        # ==================================================================

        study_lead_pairs = int(
            len(
                signals
            )
        )

        lead_present_pairs = int(
            signals[
                "lead_present"
            ]
            .fillna(
                False
            )
            .sum()
        )

        genome_wide_signal_pairs = int(
            signals[
                "has_genome_wide_signal"
            ]
            .fillna(
                False
            )
            .sum()
        )

        suggestive_signal_pairs = int(
            signals[
                "has_suggestive_signal"
            ]
            .fillna(
                False
            )
            .sum()
        )

        suggestive_only_signal_pairs = int(
            signals[
                "has_suggestive_only_signal"
            ]
            .fillna(
                False
            )
            .sum()
        )

        # ==================================================================
        # Per-lead synthesis
        # ==================================================================

        per_lead_summary = []

        for lead, group in signals.groupby(
            "lead_rsid",
            sort=True,
        ):

            valid_minimums = group[
                "minimum_p_value"
            ].dropna()

            valid_lead_p = group[
                "lead_p_value"
            ].dropna()

            per_lead_summary.append(
                {
                    "lead_rsid":
                        str(
                            lead
                        ),

                    "studies_with_locus":
                        int(
                            len(
                                group
                            )
                        ),

                    "studies_with_direct_lead":
                        int(
                            group[
                                "lead_present"
                            ]
                            .fillna(
                                False
                            )
                            .sum()
                        ),

                    "studies_with_genome_wide_signal":
                        int(
                            group[
                                "has_genome_wide_signal"
                            ]
                            .fillna(
                                False
                            )
                            .sum()
                        ),

                    "studies_with_suggestive_signal":
                        int(
                            group[
                                "has_suggestive_signal"
                            ]
                            .fillna(
                                False
                            )
                            .sum()
                        ),

                    "studies_with_suggestive_only_signal":
                        int(
                            group[
                                "has_suggestive_only_signal"
                            ]
                            .fillna(
                                False
                            )
                            .sum()
                        ),

                    "minimum_p_value_across_studies":
                        (
                            float(
                                valid_minimums.min()
                            )
                            if not valid_minimums.empty
                            else None
                        ),

                    "minimum_direct_lead_p_value_across_studies":
                        (
                            float(
                                valid_lead_p.min()
                            )
                            if not valid_lead_p.empty
                            else None
                        ),
                }
            )

        # ==================================================================
        # QC report
        # ==================================================================

        report = {
            "milestone":
                "M5.3C.3D",

            "stage":
                "gwas_locus_standardization_and_disease_signal_qc",

            "policy": {
                "all_variants_preserved":
                    True,

                "significance_filtering_performed":
                    False,

                "source_representations_removed":
                    False,

                "deduplication_performed":
                    False,

                "ld_pruning_performed":
                    False,

                "genome_wide_significance_threshold":
                    GENOME_WIDE_SIGNIFICANCE,

                "suggestive_significance_threshold":
                    SUGGESTIVE_SIGNIFICANCE,

                "suggestive_only_definition":
                    "5e-8 <= p < 1e-5",

                "significance_thresholds_descriptive_only":
                    True,

                "repeated_identity_scope":
                    "within_study_lead_locus",

                "repeated_identity_classification":
                    [
                        "exact_duplicate_representation",
                        "conflicting_variant_representation",
                    ],

                "physical_distance_interpreted_as_ld":
                    False,

                "colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,
            },

            "summary": {
                "retrieved_loci":
                    int(
                        len(
                            retrieved
                        )
                    ),

                "studies":
                    int(
                        combined[
                            "study_accession"
                        ].nunique()
                    ),

                "candidate_leads":
                    int(
                        combined[
                            "lead_rsid"
                        ].nunique()
                    ),

                "study_lead_pairs":
                    study_lead_pairs,

                "manifest_variant_rows":
                    manifest_total_rows,

                "standardized_variant_rows":
                    int(
                        len(
                            combined
                        )
                    ),

                "cardinality_failures":
                    cardinality_failures,

                "malformed_coordinate_rows":
                    int(
                        locus_qc[
                            "malformed_coordinate_rows"
                        ].sum()
                    ),

                "missing_p_value_rows":
                    int(
                        combined[
                            "p_value"
                        ].isna().sum()
                    ),

                "missing_beta_rows":
                    int(
                        combined[
                            "beta"
                        ].isna().sum()
                    ),

                "missing_standard_error_rows":
                    int(
                        combined[
                            "standard_error"
                        ].isna().sum()
                    ),

                "missing_effect_allele_frequency_rows":
                    int(
                        combined[
                            "effect_allele_frequency"
                        ].isna().sum()
                    ),

                "repeated_variant_identity_rows":
                    int(
                        locus_qc[
                            "repeated_variant_identity_rows"
                        ].sum()
                    ),

                "repeated_variant_identity_keys":
                    int(
                        locus_qc[
                            "repeated_variant_identity_keys"
                        ].sum()
                    ),

                "exact_duplicate_rows":
                    int(
                        locus_qc[
                            "exact_duplicate_rows"
                        ].sum()
                    ),

                "exact_duplicate_keys":
                    int(
                        locus_qc[
                            "exact_duplicate_keys"
                        ].sum()
                    ),

                "conflicting_variant_representation_rows":
                    int(
                        locus_qc[
                            "conflicting_variant_representation_rows"
                        ].sum()
                    ),

                "conflicting_variant_representation_keys":
                    int(
                        locus_qc[
                            "conflicting_variant_representation_keys"
                        ].sum()
                    ),

                "lead_present_pairs":
                    lead_present_pairs,

                "lead_absent_pairs":
                    (
                        study_lead_pairs
                        - lead_present_pairs
                    ),

                "pairs_with_genome_wide_signal":
                    genome_wide_signal_pairs,

                "pairs_with_suggestive_signal":
                    suggestive_signal_pairs,

                "pairs_with_suggestive_only_signal":
                    suggestive_only_signal_pairs,

                "total_genome_wide_significant_variants":
                    int(
                        signals[
                            "genome_wide_significant_variants"
                        ].sum()
                    ),

                "total_suggestive_significant_variants":
                    int(
                        signals[
                            "suggestive_significant_variants"
                        ].sum()
                    ),

                "total_suggestive_only_variants":
                    int(
                        signals[
                            "suggestive_only_variants"
                        ].sum()
                    ),
            },

            "per_lead_summary":
                per_lead_summary,

            "locus_qc":
                locus_qc.to_dict(
                    orient="records"
                ),

            "standardized_output":
                str(
                    self.standardized_output_path
                ),

            "disease_signal_output":
                str(
                    self.disease_signal_output_path
                ),
        }

        with self.qc_output_path.open(
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
            "M5.3C.3D complete."
        )

        logger.info(
            "Standardized rows: %d/%d.",
            len(
                combined
            ),
            manifest_total_rows,
        )

        logger.info(
            "Repeated identities: %d rows / %d keys | "
            "exact duplicate: %d rows / %d keys | "
            "conflicting: %d rows / %d keys.",
            report[
                "summary"
            ][
                "repeated_variant_identity_rows"
            ],
            report[
                "summary"
            ][
                "repeated_variant_identity_keys"
            ],
            report[
                "summary"
            ][
                "exact_duplicate_rows"
            ],
            report[
                "summary"
            ][
                "exact_duplicate_keys"
            ],
            report[
                "summary"
            ][
                "conflicting_variant_representation_rows"
            ],
            report[
                "summary"
            ][
                "conflicting_variant_representation_keys"
            ],
        )

        logger.info(
            "Disease signals — genome-wide pairs: %d | "
            "p<1e-5 pairs: %d | suggestive-only pairs: %d.",
            genome_wide_signal_pairs,
            suggestive_signal_pairs,
            suggestive_only_signal_pairs,
        )

        logger.info(
            "Direct lead present: %d/%d.",
            lead_present_pairs,
            study_lead_pairs,
        )

        logger.info(
            "Standardized output: %s",
            self.standardized_output_path,
        )

        logger.info(
            "Disease-signal summary: %s",
            self.disease_signal_output_path,
        )

        logger.info(
            "QC output: %s",
            self.qc_output_path,
        )

        return report