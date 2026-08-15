"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m3/candidates.py

Description:
    Core scientific logic for M3.4 integrated SNP-tRF candidate table.

    M3.4 constructs one analysis row per prostate-cancer tRF-QTL
    association identified during M3.2 and annotates whether the tRF-QTL
    variant has a direct canonical-rsID match among prostate-cancer GWAS
    associations from M3.1.

    The resulting table is a candidate scaffold rather than a ranked
    candidate list.

    M3.4 intentionally does not:

        - perform LD proxy analysis
        - calculate genomic distance to GWAS loci
        - perform genome liftover
        - perform colocalization
        - infer causal relationships
        - rank candidates
        - remove candidates lacking direct GWAS overlap
        - collapse distinct SNP-tRF pairs

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from enum import Enum

import pandas as pd


class CandidateDirectEvidence(str, Enum):
    """Direct GWAS evidence status for one SNP-tRF candidate."""

    DIRECT_GWAS_MATCH = "DIRECT_GWAS_MATCH"
    NO_DIRECT_GWAS_MATCH = "NO_DIRECT_GWAS_MATCH"
    DIRECT_MATCH_INELIGIBLE = "DIRECT_MATCH_INELIGIBLE"


class SNPTRFCandidateBuilder:
    """
    Build integrated prostate-cancer SNP-tRF candidate records.

    One output row corresponds to one input tRF-QTL association row.
    """

    GWAS_REQUIRED_COLUMNS = {
        "harm_variant_rsid",
        "m3_direct_overlap_eligible",
    }

    TRFQTL_REQUIRED_COLUMNS = {
        "harm_variant_rsid",
        "harm_variant_identity_key",
        "harm_feature_id",
        "harm_feature_identity_key",
        "harm_feature_usable",
        "m3_direct_overlap_eligible",
    }

    @classmethod
    def validate_gwas_schema(
        cls,
        dataframe: pd.DataFrame,
    ) -> None:
        """Require M3.1 GWAS fields used by M3.4."""

        missing = (
            cls.GWAS_REQUIRED_COLUMNS
            - set(
                dataframe.columns
            )
        )

        if missing:
            raise ValueError(
                "M3.1 GWAS input is missing required M3.4 columns: "
                f"{sorted(missing)}"
            )

    @classmethod
    def validate_trfqtl_schema(
        cls,
        dataframe: pd.DataFrame,
    ) -> None:
        """Require M3.2 tRF-QTL fields used by M3.4."""

        missing = (
            cls.TRFQTL_REQUIRED_COLUMNS
            - set(
                dataframe.columns
            )
        )

        if missing:
            raise ValueError(
                "M3.2 tRF-QTL input is missing required M3.4 columns: "
                f"{sorted(missing)}"
            )

    @classmethod
    def build_gwas_evidence_table(
        cls,
        gwas: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Collapse eligible GWAS association rows to one row per rsID.

        Only association multiplicity is summarized at this stage.
        Individual GWAS association rows remain preserved in the M3.1
        dataset and are not modified.
        """

        cls.validate_gwas_schema(
            gwas
        )

        eligible = (
            gwas.loc[
                gwas[
                    "m3_direct_overlap_eligible"
                ]
                .fillna(False)
                .eq(True)
                & gwas[
                    "harm_variant_rsid"
                ]
                .notna()
            ]
            .copy()
        )

        if eligible.empty:
            return pd.DataFrame(
                {
                    "harm_variant_rsid":
                        pd.Series(
                            dtype="string"
                        ),
                    "m3_gwas_association_count":
                        pd.Series(
                            dtype="int64"
                        ),
                }
            )

        evidence = (
            eligible.groupby(
                "harm_variant_rsid",
                dropna=False,
                sort=False,
            )
            .size()
            .rename(
                "m3_gwas_association_count"
            )
            .reset_index()
        )

        evidence[
            "m3_gwas_association_count"
        ] = (
            evidence[
                "m3_gwas_association_count"
            ]
            .astype(
                "int64"
            )
        )

        return evidence

    @classmethod
    def build(
        cls,
        *,
        gwas: pd.DataFrame,
        trfqtl: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Build the M3.4 integrated SNP-tRF candidate table.

        Candidate cardinality equals the number of M3.2 tRF-QTL rows.
        """

        cls.validate_gwas_schema(
            gwas
        )

        cls.validate_trfqtl_schema(
            trfqtl
        )

        candidates = (
            trfqtl.copy()
        )

        input_rows = len(
            candidates
        )

        gwas_evidence = (
            cls.build_gwas_evidence_table(
                gwas
            )
        )

        candidates = (
            candidates.merge(
                gwas_evidence,
                on="harm_variant_rsid",
                how="left",
                validate="many_to_one",
                sort=False,
            )
        )

        if (
            len(
                candidates
            )
            != input_rows
        ):
            raise RuntimeError(
                "M3.4 candidate integration changed tRF-QTL "
                f"cardinality: {input_rows} -> {len(candidates)}"
            )

        candidates[
            "m3_gwas_association_count"
        ] = (
            candidates[
                "m3_gwas_association_count"
            ]
            .fillna(
                0
            )
            .astype(
                "int64"
            )
        )

        candidates[
            "m3_direct_gwas_match"
        ] = (
            candidates[
                "m3_gwas_association_count"
            ]
            .gt(
                0
            )
        )

        candidates[
            "m3_candidate_direct_evidence"
        ] = (
            CandidateDirectEvidence
            .NO_DIRECT_GWAS_MATCH
            .value
        )

        eligible_mask = (
            candidates[
                "m3_direct_overlap_eligible"
            ]
            .fillna(False)
            .eq(True)
        )

        direct_mask = (
            eligible_mask
            & candidates[
                "m3_direct_gwas_match"
            ]
        )

        ineligible_mask = (
            ~eligible_mask
        )

        candidates.loc[
            direct_mask,
            "m3_candidate_direct_evidence",
        ] = (
            CandidateDirectEvidence
            .DIRECT_GWAS_MATCH
            .value
        )

        candidates.loc[
            ineligible_mask,
            "m3_candidate_direct_evidence",
        ] = (
            CandidateDirectEvidence
            .DIRECT_MATCH_INELIGIBLE
            .value
        )

        candidates[
            "m3_candidate_id"
        ] = [
            f"M3C{index:06d}"
            for index in range(
                1,
                len(
                    candidates
                )
                + 1,
            )
        ]

        candidates[
            "m3_candidate_stage"
        ] = "M3.4"

        candidates[
            "m3_ld_evidence_evaluated"
        ] = False

        candidates[
            "m3_colocalization_evaluated"
        ] = False

        candidates[
            "m3_candidate_ranked"
        ] = False

        return candidates