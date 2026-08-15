"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_ld_bridge.py

Description:
    Builds M5.3 LD-aware bridges between tRF-QTL lead variants and
    prostate cancer GWAS variants.

    The bridge is defined as:

        tRF-QTL lead SNP
            -> externally supplied LD proxy
            -> prostate cancer GWAS rsID

    Matching to GWAS is performed by canonical rsID identity.

    GWAS coordinates are not used because their reference-build
    provenance is not assumed to match the hg19/GRCh37 locus system.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import pandas as pd


class M53GWASLDBridgeBuilder:
    """Construct LD-aware GWAS evidence."""

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _first_existing_column(
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> str | None:
        """Return first compatible column."""

        for column in candidates:
            if column in dataframe.columns:
                return column

        return None

    @classmethod
    def _resolve_gwas_rsid_column(
        cls,
        dataframe: pd.DataFrame,
    ) -> str:
        """Resolve canonical prostate GWAS rsID field."""

        column = cls._first_existing_column(
            dataframe,
            (
                "canonical_rsid",
                "harm_variant_rsid",
                "variant_rsid",
                "rsid",
                "SNPS",
            ),
        )

        if column is None:
            raise ValueError(
                "M5.3 could not resolve the GWAS canonical rsID column."
            )

        return column

    @staticmethod
    def _project_gwas(
        dataframe: pd.DataFrame,
        *,
        rsid_column: str,
    ) -> pd.DataFrame:
        """Project useful GWAS evidence while retaining provenance."""

        preferred = (
            rsid_column,
            "association_eligible",
            "harm_variant_status",
            "variant_status",
            "DISEASE/TRAIT",
            "disease_trait",
            "P-VALUE",
            "p_value",
            "PVALUE_MLOG",
            "pvalue_mlog",
            "OR or BETA",
            "or_or_beta",
            "MAPPED_GENE",
            "mapped_gene",
            "STUDY ACCESSION",
            "study_accession",
            "PUBMEDID",
            "pubmedid",
        )

        columns = []

        for column in preferred:
            if (
                column
                in dataframe.columns
                and column
                not in columns
            ):
                columns.append(
                    column
                )

        result = dataframe[
            columns
        ].copy()

        result = result.rename(
            columns={
                rsid_column:
                    "gwas_rsid",
            }
        )

        result[
            "gwas_rsid"
        ] = (
            result[
                "gwas_rsid"
            ]
            .astype(
                "string"
            )
            .str.strip()
            .str.lower()
        )

        return result

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @classmethod
    def build(
        cls,
        loci: pd.DataFrame,
        ld_evidence: pd.DataFrame,
        prostate_gwas: pd.DataFrame,
    ) -> pd.DataFrame:
        """Build M5.3 long-form GWAS LD bridge."""

        required_loci = {
            "m5_locus_id",
            "m4_candidate_id",
            "trfqtl_rsid",
            "trf_id",
        }

        missing_loci = (
            required_loci
            - set(
                loci.columns
            )
        )

        if missing_loci:
            raise ValueError(
                "M5.3 loci missing columns: "
                f"{sorted(missing_loci)}"
            )

        required_ld = {
            "lead_rsid",
            "neighbor_rsid",
            "r2",
            "ld_evidence_class",
            "m5_bridge_eligible",
        }

        missing_ld = (
            required_ld
            - set(
                ld_evidence.columns
            )
        )

        if missing_ld:
            raise ValueError(
                "M5.3 LD evidence missing columns: "
                f"{sorted(missing_ld)}"
            )

        gwas_rsid_column = (
            cls._resolve_gwas_rsid_column(
                prostate_gwas
            )
        )

        gwas = cls._project_gwas(
            prostate_gwas,
            rsid_column=gwas_rsid_column,
        )

        candidate_columns = [
            column
            for column
            in (
                "m5_locus_id",
                "m4_candidate_id",
                "trfqtl_rsid",
                "trf_id",
                "trf_identity_key",
                "lead_chromosome",
                "lead_position",
                "locus_start",
                "locus_end",
            )
            if column
            in loci.columns
        ]

        candidates = loci[
            candidate_columns
        ].copy()

        candidates[
            "trfqtl_rsid"
        ] = (
            candidates[
                "trfqtl_rsid"
            ]
            .astype(
                "string"
            )
            .str.lower()
        )

        usable_ld = (
            ld_evidence.loc[
                ld_evidence[
                    "m5_bridge_eligible"
                ]
                .fillna(False)
                .eq(True)
            ]
            .copy()
        )

        if usable_ld.empty:
            return pd.DataFrame()

        candidate_ld = candidates.merge(
            usable_ld,
            how="inner",
            left_on="trfqtl_rsid",
            right_on="lead_rsid",
            validate="one_to_many",
        )

        if candidate_ld.empty:
            return pd.DataFrame()

        result = candidate_ld.merge(
            gwas,
            how="inner",
            left_on="neighbor_rsid",
            right_on="gwas_rsid",
            validate="many_to_many",
        )

        if result.empty:
            return result

        result[
            "m5_bridge_type"
        ] = "LD_TO_GWAS"

        result[
            "m5_bridge_variant_rsid"
        ] = result[
            "neighbor_rsid"
        ]

        result[
            "m5_exact_lead_gwas_match"
        ] = False

        result[
            "m5_bridge_supported_by_ld"
        ] = True

        result[
            "m5_primary_gwas_ld_bridge"
        ] = result[
            "m5_primary_ld_support"
        ].astype(
            "boolean"
        )

        result[
            "m5_colocalization_performed"
        ] = False

        result[
            "m5_causal_inference_performed"
        ] = False

        return result