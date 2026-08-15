"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/runners/synthesize_gwas_disease_signal.py

Description:
    Runner for M5.3C.3E cross-study prostate cancer GWAS disease-signal
    synthesis.

    Inputs:
        - M5.3C.3D standardized GWAS locus dataset.
        - M5.3C.3D disease-signal summary.

    Outputs:
        - Study × candidate-locus evidence table.
        - Candidate-level cross-study synthesis.
        - Suggestive-variant recurrence table.
        - QC JSON report.

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

from pcatrfqtl.analysis.m5.gwas_cross_study_signal import (
    synthesize_cross_study_signal,
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
# Parquet compatibility
# ============================================================================


def _read_parquet_compat(
    path: Path,
) -> pd.DataFrame:
    """Read M5 Parquet while ignoring pandas reconstruction metadata."""

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


class M53C3ECrossStudyDiseaseSignalRunner:
    """Execute M5.3C.3E."""

    def __init__(
        self,
        *,
        standardized_path: str | Path,
        disease_signal_path: str | Path,
        output_directory: str | Path,
        qc_directory: str | Path,
    ) -> None:

        self.standardized_path = Path(
            standardized_path
        )

        self.disease_signal_path = Path(
            disease_signal_path
        )

        self.output_directory = Path(
            output_directory
        )

        self.qc_directory = Path(
            qc_directory
        )

    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------

    @property
    def locus_evidence_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "gwas_cross_study_locus_evidence.parquet"
        )

    @property
    def lead_synthesis_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "gwas_cross_study_lead_synthesis.parquet"
        )

    @property
    def recurrence_path(
        self,
    ) -> Path:

        return (
            self.output_directory
            / "gwas_cross_study_variant_recurrence.parquet"
        )

    @property
    def qc_path(
        self,
    ) -> Path:

        return (
            self.qc_directory
            / "m5_3c_3e_cross_study_disease_signal.json"
        )

    # ------------------------------------------------------------------
    # Main
    # ------------------------------------------------------------------

    def run(
        self,
    ) -> dict[str, Any]:
        """Run M5.3C.3E."""

        standardized = _read_parquet_compat(
            self.standardized_path
        )

        disease_signal = _read_parquet_compat(
            self.disease_signal_path
        )

        logger.info(
            "Loaded standardized GWAS rows: %d.",
            len(
                standardized
            ),
        )

        logger.info(
            "Loaded disease-signal loci: %d.",
            len(
                disease_signal
            ),
        )

        result = synthesize_cross_study_signal(
            standardized=standardized,
            disease_signal=disease_signal,
        )

        locus_evidence = (
            result.locus_evidence
        )

        lead_synthesis = (
            result.lead_synthesis
        )

        recurrence = (
            result.variant_recurrence
        )

        self.output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        self.qc_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        write_parquet(
            locus_evidence,
            self.locus_evidence_path,
            index=False,
        )

        write_parquet(
            lead_synthesis,
            self.lead_synthesis_path,
            index=False,
        )

        write_parquet(
            recurrence,
            self.recurrence_path,
            index=False,
        )

        # ==================================================================
        # QC
        # ==================================================================

        expected_pairs = int(
            disease_signal[
                [
                    "study_accession",
                    "lead_rsid",
                ]
            ]
            .drop_duplicates()
            .shape[
                0
            ]
        )

        actual_pairs = int(
            locus_evidence[
                [
                    "study_accession",
                    "lead_rsid",
                ]
            ]
            .drop_duplicates()
            .shape[
                0
            ]
        )

        candidate_leads = int(
            lead_synthesis[
                "lead_rsid"
            ].nunique()
        )

        genome_wide_leads = int(
            lead_synthesis[
                "studies_with_genome_wide_regional_signal"
            ]
            .gt(
                0
            )
            .sum()
        )

        recurrent_suggestive_leads = int(
            lead_synthesis[
                "studies_with_suggestive_regional_signal"
            ]
            .ge(
                2
            )
            .sum()
        )

        direct_suggestive_leads = int(
            lead_synthesis[
                "studies_with_suggestive_direct_lead"
            ]
            .gt(
                0
            )
            .sum()
        )

        ld_required_leads = int(
            lead_synthesis[
                "ld_analysis_required"
            ]
            .fillna(
                False
            )
            .sum()
        )

        recurrent_variants = (
            recurrence.loc[
                recurrence[
                    "study_count"
                ].ge(
                    2
                )
            ]
            if not recurrence.empty
            else recurrence
        )

        report = {
            "milestone":
                "M5.3C.3E",

            "stage":
                "cross_study_disease_signal_synthesis",

            "policy": {
                "source_standardized_stage":
                    "M5.3C.3D",

                "significance_filtering_of_master_dataset":
                    False,

                "suggestive_recurrence_threshold":
                    1e-5,

                "variant_recurrence_interpreted_as_independent_replication":
                    False,

                "physical_distance_interpreted_as_ld":
                    False,

                "ld_analysis_performed":
                    False,

                "colocalization_performed":
                    False,

                "fine_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,

                "candidate_ranking":
                    "descriptive_only",
            },

            "summary": {
                "standardized_variant_rows":
                    int(
                        len(
                            standardized
                        )
                    ),

                "source_study_lead_pairs":
                    expected_pairs,

                "synthesized_study_lead_pairs":
                    actual_pairs,

                "study_lead_pair_cardinality_preserved":
                    bool(
                        expected_pairs
                        == actual_pairs
                    ),

                "candidate_leads":
                    candidate_leads,

                "candidate_leads_with_genome_wide_regional_signal":
                    genome_wide_leads,

                "candidate_leads_with_recurrent_suggestive_regional_signal":
                    recurrent_suggestive_leads,

                "candidate_leads_with_suggestive_direct_lead_signal":
                    direct_suggestive_leads,

                "candidate_leads_requiring_ld_analysis":
                    ld_required_leads,

                "suggestive_variant_identities":
                    int(
                        len(
                            recurrence
                        )
                    ),

                "recurrent_suggestive_variant_identities":
                    int(
                        len(
                            recurrent_variants
                        )
                    ),

                "maximum_cross_study_variant_recurrence":
                    (
                        int(
                            recurrence[
                                "study_count"
                            ].max()
                        )
                        if not recurrence.empty
                        else 0
                    ),
            },

            "lead_synthesis":
                lead_synthesis.to_dict(
                    orient="records"
                ),

            "outputs": {
                "locus_evidence":
                    str(
                        self.locus_evidence_path
                    ),

                "lead_synthesis":
                    str(
                        self.lead_synthesis_path
                    ),

                "variant_recurrence":
                    str(
                        self.recurrence_path
                    ),
            },
        }

        with self.qc_path.open(
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
            "M5.3C.3E complete."
        )

        logger.info(
            "Study-lead pairs synthesized: %d/%d.",
            actual_pairs,
            expected_pairs,
        )

        logger.info(
            "Candidate leads: %d | genome-wide regional: %d | "
            "recurrent suggestive: %d.",
            candidate_leads,
            genome_wide_leads,
            recurrent_suggestive_leads,
        )

        logger.info(
            "Direct suggestive candidate leads: %d | "
            "LD analysis required: %d.",
            direct_suggestive_leads,
            ld_required_leads,
        )

        logger.info(
            "Suggestive variant identities: %d | recurrent: %d.",
            len(
                recurrence
            ),
            len(
                recurrent_variants
            ),
        )

        logger.info(
            "Lead synthesis: %s",
            self.lead_synthesis_path,
        )

        logger.info(
            "Variant recurrence: %s",
            self.recurrence_path,
        )

        logger.info(
            "QC output: %s",
            self.qc_path,
        )

        return report