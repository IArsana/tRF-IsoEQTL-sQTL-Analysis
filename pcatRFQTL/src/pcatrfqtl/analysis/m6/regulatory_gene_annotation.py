"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/regulatory_gene_annotation.py

Description:
    Core logic for M6.2B source-native regulatory feature → gene annotation.

    This stage attempts to resolve retained Moradi regulatory features using
    source-compatible annotation artifacts.

    The preferred model follows FASE, where event annotation can provide:

        event
        chromosome
        start
        stop
        strand
        gene

    M6.2B does NOT:
        - decode INT identifiers as genomic coordinates;
        - use sQTL SNP coordinates as feature coordinates;
        - use GWAS/tag SNP coordinates as feature coordinates;
        - assign nearest genes;
        - infer a gene merely because an association is cis;
        - perform external GTF reconstruction.

    If source annotation cannot be identified, the feature remains unresolved
    and is forwarded to M6.2C.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import pyarrow.parquet as pq


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class RegulatoryGeneAnnotationResult:
    """Container for M6.2B outputs."""

    annotation_candidates: pd.DataFrame
    feature_resolution: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Helpers
# ============================================================================


def _normalize_text(
    value: Any,
) -> str | None:
    """Normalize a scalar text value."""

    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass

    value = str(value).strip()

    return value or None


def _normalize_column(
    value: Any,
) -> str:
    """Normalize a column name for semantic matching."""

    return (
        str(value)
        .strip()
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
    )


def _first_matching_column(
    columns: Iterable[Any],
    aliases: Iterable[str],
) -> str | None:
    """Return first source column matching configured aliases."""

    normalized_lookup = {
        _normalize_column(column): str(column)
        for column in columns
    }

    for alias in aliases:

        normalized_alias = _normalize_column(alias)

        if normalized_alias in normalized_lookup:
            return normalized_lookup[normalized_alias]

    return None


# ============================================================================
# Input candidate extraction
# ============================================================================


def extract_features_requiring_annotation(
    m6_2a_qc: dict[str, Any],
) -> pd.DataFrame:
    """Extract M6.2A regulatory features requiring gene annotation."""

    rows: list[dict[str, Any]] = []

    for candidate in m6_2a_qc.get(
        "candidate_mapping",
        [],
    ):

        if not bool(
            candidate.get(
                "continue_to_gene_annotation",
                False,
            )
        ):
            continue

        feature_id = _normalize_text(
            candidate.get(
                "regulatory_feature_id"
            )
        )

        if feature_id is None:
            continue

        rows.append(
            {
                "lead_rsid":
                    _normalize_text(
                        candidate.get(
                            "lead_rsid"
                        )
                    ),

                "priority_rank":
                    candidate.get(
                        "priority_rank"
                    ),

                "priority_class":
                    _normalize_text(
                        candidate.get(
                            "priority_class"
                        )
                    ),

                "regulatory_feature_id":
                    feature_id,

                "feature_class":
                    _normalize_text(
                        candidate.get(
                            "feature_class"
                        )
                    ),

                "regulatory_scope":
                    _normalize_text(
                        candidate.get(
                            "regulatory_scope"
                        )
                    ),

                "upstream_parent_gene_status":
                    _normalize_text(
                        candidate.get(
                            "parent_gene_status"
                        )
                    ),
            }
        )

    return pd.DataFrame(rows)


# ============================================================================
# Source reading
# ============================================================================


def _read_parquet_safely(
    path: Path,
) -> pd.DataFrame:
    """
    Read Parquet without pandas metadata reconstruction.

    This follows the M5/M6-safe reader policy for potentially nested columns.
    """

    table = pq.read_table(
        path,
        use_pandas_metadata=False,
    )

    return table.to_pandas(
        ignore_metadata=True,
    )


def _iter_tabular_objects(
    path: Path,
) -> Iterable[tuple[str | None, pd.DataFrame]]:
    """Yield tabular objects from supported local source files."""

    suffix = path.suffix.lower()

    if suffix == ".csv":

        yield (
            None,
            pd.read_csv(
                path,
                low_memory=False,
            ),
        )

        return

    if suffix in {
        ".tsv",
        ".txt",
    }:

        yield (
            None,
            pd.read_csv(
                path,
                sep="\t",
                low_memory=False,
            ),
        )

        return

    if suffix in {
        ".xlsx",
        ".xls",
    }:

        book = pd.ExcelFile(path)

        for sheet in book.sheet_names:

            yield (
                sheet,
                pd.read_excel(
                    book,
                    sheet_name=sheet,
                ),
            )

        return

    if suffix == ".parquet":

        yield (
            None,
            _read_parquet_safely(
                path
            ),
        )


# ============================================================================
# Annotation inspection
# ============================================================================


def inspect_annotation_table(
    *,
    dataframe: pd.DataFrame,
    path: Path,
    sheet: str | None,
    target_features: set[str],
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    """
    Search one table for feature → gene annotation.

    A valid annotation candidate requires:
        1. recognizable feature/event column;
        2. recognizable gene column;
        3. exact feature ID match.

    Location fields are retained when available but are not required for
    source gene resolution.
    """

    search_config = config[
        "source_annotation_search"
    ]

    feature_column = _first_matching_column(
        dataframe.columns,
        search_config[
            "feature_identifier_columns"
        ],
    )

    gene_column = _first_matching_column(
        dataframe.columns,
        search_config[
            "gene_identifier_columns"
        ],
    )

    if (
        feature_column is None
        or gene_column is None
    ):
        return []

    chromosome_column = _first_matching_column(
        dataframe.columns,
        search_config[
            "chromosome_columns"
        ],
    )

    start_column = _first_matching_column(
        dataframe.columns,
        search_config[
            "start_columns"
        ],
    )

    stop_column = _first_matching_column(
        dataframe.columns,
        search_config[
            "stop_columns"
        ],
    )

    strand_column = _first_matching_column(
        dataframe.columns,
        search_config[
            "strand_columns"
        ],
    )

    feature_values = (
        dataframe[
            feature_column
        ]
        .astype("string")
        .str.strip()
    )

    hits = dataframe.loc[
        feature_values.isin(
            target_features
        )
    ].copy()

    if hits.empty:
        return []

    results: list[dict[str, Any]] = []

    for _, row in hits.iterrows():

        feature_id = _normalize_text(
            row.get(
                feature_column
            )
        )

        gene = _normalize_text(
            row.get(
                gene_column
            )
        )

        if feature_id is None:
            continue

        results.append(
            {
                "regulatory_feature_id":
                    feature_id,

                "source_path":
                    str(path),

                "source_sheet":
                    sheet,

                "feature_column":
                    feature_column,

                "gene_column":
                    gene_column,

                "source_gene":
                    gene,

                "source_chromosome":
                    (
                        _normalize_text(
                            row.get(
                                chromosome_column
                            )
                        )
                        if chromosome_column
                        else None
                    ),

                "source_start":
                    (
                        row.get(
                            start_column
                        )
                        if start_column
                        else None
                    ),

                "source_stop":
                    (
                        row.get(
                            stop_column
                        )
                        if stop_column
                        else None
                    ),

                "source_strand":
                    (
                        _normalize_text(
                            row.get(
                                strand_column
                            )
                        )
                        if strand_column
                        else None
                    ),

                "annotation_class":
                    "SOURCE_COMPATIBLE_EVENT_ANNOTATION",

                "exact_feature_match":
                    True,

                "nearest_gene_inference":
                    False,

                "positional_gene_inference":
                    False,
            }
        )

    return results


# ============================================================================
# Project-wide source annotation search
# ============================================================================


def search_source_annotations(
    *,
    project_root: Path,
    target_features: set[str],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Search configured local roots for exact source annotation."""

    source_config = config[
        "source_annotation_search"
    ]

    allowed_extensions = {
        str(extension).lower()
        for extension
        in source_config[
            "allowed_extensions"
        ]
    }

    rows: list[dict[str, Any]] = []

    seen_paths: set[Path] = set()

    for relative_root in source_config[
        "search_roots"
    ]:

        root = (
            project_root
            / str(relative_root)
        )

        if not root.exists():
            continue

        for path in sorted(
            root.rglob("*")
        ):

            if not path.is_file():
                continue

            if path.suffix.lower() not in allowed_extensions:
                continue

            resolved = path.resolve()

            if resolved in seen_paths:
                continue

            seen_paths.add(
                resolved
            )

            try:

                for sheet, frame in _iter_tabular_objects(
                    path
                ):

                    rows.extend(
                        inspect_annotation_table(
                            dataframe=frame,
                            path=path,
                            sheet=sheet,
                            target_features=target_features,
                            config=config,
                        )
                    )

            except Exception:
                # M6.2B is an annotation audit across heterogeneous source
                # files. Unreadable/non-tabular artifacts do not constitute
                # biological negative evidence.
                continue

    if not rows:

        return pd.DataFrame(
            columns=[
                "regulatory_feature_id",
                "source_path",
                "source_sheet",
                "feature_column",
                "gene_column",
                "source_gene",
                "source_chromosome",
                "source_start",
                "source_stop",
                "source_strand",
                "annotation_class",
                "exact_feature_match",
                "nearest_gene_inference",
                "positional_gene_inference",
            ]
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Final resolution
# ============================================================================


def resolve_feature_genes(
    *,
    features: pd.DataFrame,
    annotation_candidates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Resolve each regulatory feature conservatively."""

    statuses = config[
        "resolution_status"
    ]

    output: list[
        dict[str, Any]
    ] = []

    for feature in features.to_dict(
        orient="records"
    ):

        feature_id = str(
            feature[
                "regulatory_feature_id"
            ]
        )

        if annotation_candidates.empty:

            hits = pd.DataFrame()

        else:

            hits = annotation_candidates.loc[
                annotation_candidates[
                    "regulatory_feature_id"
                ].astype(
                    str
                )
                ==
                feature_id
            ].copy()

        genes = sorted(
            {
                str(gene)
                for gene
                in hits.get(
                    "source_gene",
                    pd.Series(
                        dtype="object"
                    ),
                )
                .dropna()
                .tolist()
                if str(gene).strip()
            }
        )

        # ------------------------------------------------------------------
        # Exactly one source gene
        # ------------------------------------------------------------------

        if len(genes) == 1:

            gene = genes[
                0
            ]

            status = str(
                statuses[
                    "resolved"
                ]
            )

            resolved = True
            continue_to_reconstruction = False

        # ------------------------------------------------------------------
        # Multiple incompatible source genes
        # ------------------------------------------------------------------

        elif len(genes) > 1:

            gene = None

            status = str(
                statuses[
                    "ambiguous"
                ]
            )

            resolved = False
            continue_to_reconstruction = True

        # ------------------------------------------------------------------
        # No source-native mapping
        # ------------------------------------------------------------------

        else:

            gene = None

            status = str(
                statuses[
                    "source_artifact_missing"
                ]
            )

            resolved = False
            continue_to_reconstruction = True

        output.append(
            {
                "lead_rsid":
                    feature[
                        "lead_rsid"
                    ],

                "priority_rank":
                    feature[
                        "priority_rank"
                    ],

                "priority_class":
                    feature[
                        "priority_class"
                    ],

                "regulatory_feature_id":
                    feature_id,

                "feature_class":
                    feature[
                        "feature_class"
                    ],

                "regulatory_scope":
                    feature[
                        "regulatory_scope"
                    ],

                "source_annotation_matches":
                    int(
                        len(
                            hits
                        )
                    ),

                "unique_source_genes":
                    int(
                        len(
                            genes
                        )
                    ),

                "parent_gene":
                    gene,

                "parent_gene_resolved":
                    resolved,

                "gene_annotation_status":
                    status,

                "gene_assignment_method":
                    (
                        "SOURCE_REPORTED_EVENT_ANNOTATION"
                        if resolved
                        else None
                    ),

                "nearest_gene_assignment_performed":
                    False,

                "snp_proximity_gene_assignment_performed":
                    False,

                "continue_to_m6_2c":
                    continue_to_reconstruction,

                "causal_claim_allowed":
                    False,
            }
        )

    return pd.DataFrame(
        output
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    feature_resolution: pd.DataFrame,
) -> pd.DataFrame:
    """Build M6.2B summary."""

    assessed = int(
        len(
            feature_resolution
        )
    )

    resolved = int(
        feature_resolution[
            "parent_gene_resolved"
        ]
        .astype(
            bool
        )
        .sum()
        if not feature_resolution.empty
        else 0
    )

    unresolved = (
        assessed
        -
        resolved
    )

    next_stage = (
        "M6.3_CANDIDATE_GENE_INTEGRATION"
        if unresolved == 0
        else "M6.2C_SOURCE_COMPATIBLE_EVENT_RECONSTRUCTION"
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M6.2B",

                "features_assessed":
                    assessed,

                "features_source_gene_resolved":
                    resolved,

                "features_source_gene_unresolved":
                    unresolved,

                "nearest_gene_mapping_performed":
                    False,

                "snp_proximity_mapping_performed":
                    False,

                "external_gtf_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,

                "next_stage":
                    next_stage,
            }
        ]
    )


# ============================================================================
# Public API
# ============================================================================


def resolve_source_gene_annotations(
    *,
    project_root: Path,
    m6_2a_qc: dict[str, Any],
    config: dict[str, Any],
) -> RegulatoryGeneAnnotationResult:
    """Execute M6.2B."""

    features = extract_features_requiring_annotation(
        m6_2a_qc
    )

    targets = set(
        features[
            "regulatory_feature_id"
        ].astype(
            str
        )
    )

    annotation_candidates = search_source_annotations(
        project_root=project_root,
        target_features=targets,
        config=config,
    )

    resolution = resolve_feature_genes(
        features=features,
        annotation_candidates=annotation_candidates,
        config=config,
    )

    summary = build_summary(
        resolution
    )

    return RegulatoryGeneAnnotationResult(
        annotation_candidates=annotation_candidates,
        feature_resolution=resolution,
        summary=summary,
    )