"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_sumstats_discovery.py

Description:
    M5.3C prostate cancer GWAS summary-statistics study discovery.

    This module identifies prostate-related GWAS Catalog studies from the
    standardized local GWAS association dataset and prioritizes studies for
    full summary-statistics retrieval.

    The discovery stage does not download genome-wide summary statistics.

    Important safeguards:
        - Existing M3 disease harmonization is preserved.
        - Broad prostate-related text matching is used only as an audit and
          discovery mechanism.
        - Study accessions are deduplicated only at study-discovery level.
        - Association rows are not rewritten.
        - No locus association inference is performed.
        - No colocalization is performed.
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

import re
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


PROSTATE_PATTERN = re.compile(
    r"\bprostate\b",
    flags=re.IGNORECASE,
)

MALIGNANCY_PATTERN = re.compile(
    r"\b("
    r"cancer"
    r"|carcinoma"
    r"|neoplasm"
    r"|malignan\w*"
    r"|tumou?r"
    r")\b",
    flags=re.IGNORECASE,
)


class M53CProstateGWASStudyDiscovery:
    """Discover prostate-related GWAS studies from standardized Catalog data."""

    STUDY_COLUMNS = (
        "STUDY ACCESSION",
        "study_accession",
        "study_id",
    )

    TRAIT_COLUMNS = (
        "DISEASE/TRAIT",
        "disease_trait",
        "reported_trait",
        "trait",
    )

    MAPPED_TRAIT_COLUMNS = (
        "MAPPED_TRAIT",
        "mapped_trait",
    )

    TITLE_COLUMNS = (
        "STUDY",
        "study",
        "title",
    )

    PMID_COLUMNS = (
        "PUBMEDID",
        "pubmedid",
        "pmid",
    )

    AUTHOR_COLUMNS = (
        "FIRST AUTHOR",
        "first_author",
        "author",
    )

    @staticmethod
    def _first_existing_column(
        dataframe: pd.DataFrame,
        candidates: tuple[str, ...],
    ) -> str | None:
        """Resolve the first case-insensitive compatible column."""

        lookup = {
            str(column).strip().lower():
                column
            for column
            in dataframe.columns
        }

        for candidate in candidates:

            resolved = lookup.get(
                candidate.strip().lower()
            )

            if resolved is not None:
                return resolved

        return None

    @classmethod
    def _extract_metadata(
        cls,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Extract standardized study-level discovery fields."""

        study_column = cls._first_existing_column(
            dataframe,
            cls.STUDY_COLUMNS,
        )

        if study_column is None:
            raise ValueError(
                "GWAS data does not contain a study accession column."
            )

        trait_column = cls._first_existing_column(
            dataframe,
            cls.TRAIT_COLUMNS,
        )

        mapped_trait_column = cls._first_existing_column(
            dataframe,
            cls.MAPPED_TRAIT_COLUMNS,
        )

        title_column = cls._first_existing_column(
            dataframe,
            cls.TITLE_COLUMNS,
        )

        pmid_column = cls._first_existing_column(
            dataframe,
            cls.PMID_COLUMNS,
        )

        author_column = cls._first_existing_column(
            dataframe,
            cls.AUTHOR_COLUMNS,
        )

        result = pd.DataFrame(
            index=dataframe.index
        )

        result[
            "study_accession"
        ] = (
            dataframe[
                study_column
            ]
            .astype("string")
            .str.strip()
        )

        result[
            "reported_trait"
        ] = (
            dataframe[
                trait_column
            ].astype("string")
            if trait_column
            else pd.Series(
                pd.NA,
                index=dataframe.index,
                dtype="string",
            )
        )

        result[
            "mapped_trait"
        ] = (
            dataframe[
                mapped_trait_column
            ].astype("string")
            if mapped_trait_column
            else pd.Series(
                pd.NA,
                index=dataframe.index,
                dtype="string",
            )
        )

        result[
            "study_title"
        ] = (
            dataframe[
                title_column
            ].astype("string")
            if title_column
            else pd.Series(
                pd.NA,
                index=dataframe.index,
                dtype="string",
            )
        )

        result[
            "pmid"
        ] = (
            dataframe[
                pmid_column
            ].astype("string")
            if pmid_column
            else pd.Series(
                pd.NA,
                index=dataframe.index,
                dtype="string",
            )
        )

        result[
            "first_author"
        ] = (
            dataframe[
                author_column
            ].astype("string")
            if author_column
            else pd.Series(
                pd.NA,
                index=dataframe.index,
                dtype="string",
            )
        )

        return result

    @classmethod
    def discover(
        cls,
        gwas_paths: Iterable[Path],
    ) -> pd.DataFrame:
        """Discover candidate prostate cancer GWAS studies."""

        frames: list[pd.DataFrame] = []

        for path in gwas_paths:

            path = Path(
                path
            )

            source = pd.read_parquet(
                path
            )

            metadata = cls._extract_metadata(
                source
            )

            frames.append(
                metadata
            )

            del source

        combined = pd.concat(
            frames,
            ignore_index=True,
        )

        audit_text = (
            combined[
                [
                    "reported_trait",
                    "mapped_trait",
                    "study_title",
                ]
            ]
            .fillna("")
            .astype("string")
            .agg(
                " | ".join,
                axis=1,
            )
        )

        combined[
            "prostate_text_hit"
        ] = (
            audit_text
            .str.contains(
                PROSTATE_PATTERN,
                regex=True,
                na=False,
            )
            .astype(
                "boolean"
            )
        )

        combined[
            "malignancy_text_hit"
        ] = (
            audit_text
            .str.contains(
                MALIGNANCY_PATTERN,
                regex=True,
                na=False,
            )
            .astype(
                "boolean"
            )
        )

        combined[
            "prostate_cancer_candidate"
        ] = (
            combined[
                "prostate_text_hit"
            ]
            &
            combined[
                "malignancy_text_hit"
            ]
        ).astype(
            "boolean"
        )

        candidates = (
            combined.loc[
                combined[
                    "prostate_cancer_candidate"
                ]
                .fillna(False)
            ]
            .copy()
        )

        if candidates.empty:
            return candidates

        # One representative row per GCST, while retaining trait aggregation.
        grouped = (
            candidates
            .groupby(
                "study_accession",
                dropna=False,
            )
            .agg(
                reported_traits=(
                    "reported_trait",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
                mapped_traits=(
                    "mapped_trait",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
                study_titles=(
                    "study_title",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
                pmids=(
                    "pmid",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
                first_authors=(
                    "first_author",
                    lambda values:
                        sorted(
                            {
                                str(value)
                                for value
                                in values.dropna()
                            }
                        ),
                ),
                catalog_association_rows=(
                    "study_accession",
                    "size",
                ),
            )
            .reset_index()
        )

        grouped[
            "m53c_discovery_status"
        ] = "CANDIDATE_FOR_SUMSTATS_AVAILABILITY_CHECK"

        return grouped