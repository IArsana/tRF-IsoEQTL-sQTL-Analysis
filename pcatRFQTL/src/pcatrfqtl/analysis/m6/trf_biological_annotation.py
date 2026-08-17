"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/trf_biological_annotation.py

Description:
    Core logic for M6.4A Candidate/tRF Biological Annotation.

    This stage integrates:
        - retained candidate variants;
        - candidate-associated tRF identities;
        - available chromosome/locus information;
        - source-reported tRF annotations;
        - retained regulatory-feature context;
        - RNA-processing context.

    M6.4A does NOT:
        - infer genes from variant coordinates;
        - decode tRF identifiers;
        - infer parent tRNAs;
        - infer amino acids or anticodons;
        - perform nearest-gene mapping;
        - perform colocalization or fine-mapping;
        - make causal claims.

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
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class TrfBiologicalAnnotationResult:
    """Container for M6.4A outputs."""

    candidate_trf_evidence: pd.DataFrame
    candidate_biological_annotation: pd.DataFrame
    summary: pd.DataFrame


def _normalize_text(
    value: Any,
) -> str | None:
    """Normalize optional text."""

    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass

    value = str(value).strip()

    return value or None


def _find_column(
    dataframe: pd.DataFrame,
    aliases: list[str],
) -> str | None:
    """Find first available column matching configured aliases."""

    lookup = {
        str(column).strip().lower(): str(column)
        for column in dataframe.columns
    }

    for alias in aliases:

        normalized = str(alias).strip().lower()

        if normalized in lookup:
            return lookup[normalized]

    return None


def _optional_series(
    dataframe: pd.DataFrame,
    column: str | None,
) -> pd.Series:
    """Return source column or a null series."""

    if column is None:
        return pd.Series(
            [None] * len(dataframe),
            index=dataframe.index,
            dtype="object",
        )

    return dataframe[column]


def standardize_candidate_trf_source(
    *,
    dataframe: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Standardize candidate/tRF source using configured aliases."""

    aliases = config[
        "candidate_trf_source"
    ][
        "column_aliases"
    ]

    resolved = {
        key: _find_column(
            dataframe,
            list(values),
        )
        for key, values
        in aliases.items()
    }

    if resolved["lead_rsid"] is None:

        raise RuntimeError(
            "M6.4A candidate/tRF source does not contain a recognizable rsID "
            "column."
        )

    if resolved["trf_id"] is None:

        raise RuntimeError(
            "M6.4A candidate/tRF source does not contain a recognizable tRF "
            "identity column."
        )

    output = pd.DataFrame(
        {
            "lead_rsid":
                _optional_series(
                    dataframe,
                    resolved["lead_rsid"],
                ),

            "trf_id":
                _optional_series(
                    dataframe,
                    resolved["trf_id"],
                ),

            "chromosome":
                _optional_series(
                    dataframe,
                    resolved["chromosome"],
                ),

            "position":
                _optional_series(
                    dataframe,
                    resolved["position"],
                ),

            "cancer_type":
                _optional_series(
                    dataframe,
                    resolved["cancer_type"],
                ),

            "trf_class":
                _optional_series(
                    dataframe,
                    resolved["trf_class"],
                ),

            "parent_trna":
                _optional_series(
                    dataframe,
                    resolved["parent_trna"],
                ),

            "amino_acid":
                _optional_series(
                    dataframe,
                    resolved["amino_acid"],
                ),

            "anticodon":
                _optional_series(
                    dataframe,
                    resolved["anticodon"],
                ),

            "trf_qtl_p_value":
                _optional_series(
                    dataframe,
                    resolved["qtl_p_value"],
                ),

            "trf_qtl_beta":
                _optional_series(
                    dataframe,
                    resolved["qtl_beta"],
                ),
        }
    )

    for column in [
        "lead_rsid",
        "trf_id",
        "chromosome",
        "cancer_type",
        "trf_class",
        "parent_trna",
        "amino_acid",
        "anticodon",
    ]:

        output[column] = output[column].map(
            _normalize_text
        )

    output = output.loc[
        output["lead_rsid"].notna()
    ].copy()

    return output.reset_index(
        drop=True
    )


def extract_m6_3_candidates(
    *,
    m6_3_qc: dict[str, Any],
) -> pd.DataFrame:
    """Extract locked M6.3 candidate state."""

    rows: list[dict[str, Any]] = []

    for record in m6_3_qc.get(
        "final_candidate_gene_integration",
        [],
    ):

        rows.append(
            {
                "lead_rsid":
                    _normalize_text(
                        record.get(
                            "lead_rsid"
                        )
                    ),

                "priority_rank":
                    record.get(
                        "priority_rank"
                    ),

                "priority_class":
                    _normalize_text(
                        record.get(
                            "priority_class"
                        )
                    ),

                "regulatory_feature_id":
                    _normalize_text(
                        record.get(
                            "regulatory_feature_id"
                        )
                    ),

                "regulatory_feature_class":
                    _normalize_text(
                        record.get(
                            "regulatory_feature_class"
                        )
                    ),

                "regulatory_feature_supported":
                    bool(
                        record.get(
                            "regulatory_feature_supported",
                            False,
                        )
                    ),

                "integrated_gene":
                    _normalize_text(
                        record.get(
                            "integrated_gene"
                        )
                    ),

                "integrated_gene_resolved":
                    bool(
                        record.get(
                            "integrated_gene_resolved",
                            False,
                        )
                    ),

                "gene_resolution_status":
                    _normalize_text(
                        record.get(
                            "gene_resolution_status"
                        )
                    ),

                "candidate_retained":
                    bool(
                        record.get(
                            "candidate_retained",
                            True,
                        )
                    ),
            }
        )

    result = pd.DataFrame(rows)

    if result.empty:
        raise RuntimeError(
            "M6.3 contains no retained candidates."
        )

    return result


def build_candidate_trf_evidence(
    *,
    candidates: pd.DataFrame,
    trf_source: pd.DataFrame,
) -> pd.DataFrame:
    """Join locked candidates to source tRF annotations."""

    retained = set(
        candidates[
            "lead_rsid"
        ]
        .dropna()
        .astype(str)
    )

    source = trf_source.loc[
        trf_source[
            "lead_rsid"
        ]
        .astype(str)
        .isin(retained)
    ].copy()

    return source.sort_values(
        [
            "lead_rsid",
            "trf_id",
        ],
        kind="stable",
    ).reset_index(
        drop=True
    )


def _build_rna_context(
    *,
    regulatory_feature_supported: bool,
    regulatory_feature_class: str | None,
    config: dict[str, Any],
) -> str:
    """Create conservative RNA-processing context label."""

    labels = config[
        "rna_processing_context"
    ]

    if regulatory_feature_supported:

        if regulatory_feature_class == "INTRON_RETENTION":

            return str(
                labels[
                    "combined_context_label"
                ]
            )

        return str(
            labels[
                "regulatory_feature_context_label"
            ]
        )

    return str(
        labels[
            "no_regulatory_feature_label"
        ]
    )


def build_candidate_biological_annotation(
    *,
    candidates: pd.DataFrame,
    candidate_trf_evidence: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Create one-row-per-candidate biological annotation."""

    statuses = config[
        "annotation_status"
    ]

    rows: list[dict[str, Any]] = []

    for candidate in candidates.to_dict(
        orient="records"
    ):

        lead_rsid = str(
            candidate[
                "lead_rsid"
            ]
        )

        trf_rows = candidate_trf_evidence.loc[
            candidate_trf_evidence[
                "lead_rsid"
            ]
            .astype(str)
            ==
            lead_rsid
        ].copy()

        trf_ids = sorted(
            {
                str(value)
                for value
                in trf_rows[
                    "trf_id"
                ]
                .dropna()
                .tolist()
                if str(value).strip()
            }
        )

        chromosomes = sorted(
            {
                str(value)
                for value
                in trf_rows[
                    "chromosome"
                ]
                .dropna()
                .tolist()
                if str(value).strip()
            }
        )

        trf_classes = sorted(
            {
                str(value)
                for value
                in trf_rows[
                    "trf_class"
                ]
                .dropna()
                .tolist()
                if str(value).strip()
            }
        )

        parent_trnas = sorted(
            {
                str(value)
                for value
                in trf_rows[
                    "parent_trna"
                ]
                .dropna()
                .tolist()
                if str(value).strip()
            }
        )

        amino_acids = sorted(
            {
                str(value)
                for value
                in trf_rows[
                    "amino_acid"
                ]
                .dropna()
                .tolist()
                if str(value).strip()
            }
        )

        anticodons = sorted(
            {
                str(value)
                for value
                in trf_rows[
                    "anticodon"
                ]
                .dropna()
                .tolist()
                if str(value).strip()
            }
        )

        feature_supported = bool(
            candidate[
                "regulatory_feature_supported"
            ]
        )

        rna_context = _build_rna_context(
            regulatory_feature_supported=feature_supported,
            regulatory_feature_class=_normalize_text(
                candidate.get(
                    "regulatory_feature_class"
                )
            ),
            config=config,
        )

        if not trf_ids:

            status = str(
                statuses[
                    "missing_trf"
                ]
            )

        elif (
            trf_classes
            or parent_trnas
            or amino_acids
            or anticodons
        ):

            status = str(
                statuses[
                    "complete"
                ]
            )

        else:

            status = str(
                statuses[
                    "partial"
                ]
            )

        rows.append(
            {
                "lead_rsid":
                    lead_rsid,

                "priority_rank":
                    int(
                        candidate[
                            "priority_rank"
                        ]
                    ),

                "priority_class":
                    str(
                        candidate[
                            "priority_class"
                        ]
                    ),

                "trf_count":
                    int(
                        len(
                            trf_ids
                        )
                    ),

                "trf_ids":
                    trf_ids,

                "chromosomes":
                    chromosomes,

                "trf_classes":
                    trf_classes,

                "parent_trnas":
                    parent_trnas,

                "amino_acids":
                    amino_acids,

                "anticodons":
                    anticodons,

                "regulatory_feature_id":
                    _normalize_text(
                        candidate.get(
                            "regulatory_feature_id"
                        )
                    ),

                "regulatory_feature_class":
                    _normalize_text(
                        candidate.get(
                            "regulatory_feature_class"
                        )
                    ),

                "regulatory_feature_supported":
                    feature_supported,

                "integrated_gene":
                    _normalize_text(
                        candidate.get(
                            "integrated_gene"
                        )
                    ),

                "integrated_gene_resolved":
                    bool(
                        candidate.get(
                            "integrated_gene_resolved",
                            False,
                        )
                    ),

                "gene_resolution_status":
                    _normalize_text(
                        candidate.get(
                            "gene_resolution_status"
                        )
                    ),

                "rna_processing_context":
                    rna_context,

                "annotation_status":
                    status,

                "candidate_retained":
                    True,

                "trf_identifier_decoded":
                    False,

                "parent_trna_inferred":
                    False,

                "nearest_gene_assignment_performed":
                    False,

                "causal_claim_allowed":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "priority_rank",
            "lead_rsid",
        ],
        kind="stable",
    ).reset_index(
        drop=True
    )


def build_summary(
    *,
    annotation: pd.DataFrame,
) -> pd.DataFrame:
    """Build M6.4A summary."""

    with_trf = int(
        (
            annotation[
                "trf_count"
            ]
            >
            0
        ).sum()
    )

    with_regulatory_context = int(
        annotation[
            "regulatory_feature_supported"
        ]
        .astype(bool)
        .sum()
    )

    resolved_gene = int(
        annotation[
            "integrated_gene_resolved"
        ]
        .astype(bool)
        .sum()
    )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M6.4A",

                "candidates_assessed":
                    int(
                        len(
                            annotation
                        )
                    ),

                "candidates_with_trf_identity":
                    with_trf,

                "candidates_with_regulatory_feature_context":
                    with_regulatory_context,

                "candidates_with_resolved_gene":
                    resolved_gene,

                "trf_identifier_inference_performed":
                    False,

                "nearest_gene_mapping_performed":
                    False,

                "causal_inference_performed":
                    False,

                "next_stage":
                    "M6.4B_PATHWAY_READINESS_ASSESSMENT",
            }
        ]
    )


def annotate_candidate_trf_biology(
    *,
    m6_3_qc: dict[str, Any],
    candidate_trf_source: pd.DataFrame,
    config: dict[str, Any],
) -> TrfBiologicalAnnotationResult:
    """Execute M6.4A."""

    candidates = extract_m6_3_candidates(
        m6_3_qc=m6_3_qc
    )

    standardized_source = standardize_candidate_trf_source(
        dataframe=candidate_trf_source,
        config=config,
    )

    evidence = build_candidate_trf_evidence(
        candidates=candidates,
        trf_source=standardized_source,
    )

    annotation = build_candidate_biological_annotation(
        candidates=candidates,
        candidate_trf_evidence=evidence,
        config=config,
    )

    summary = build_summary(
        annotation=annotation
    )

    return TrfBiologicalAnnotationResult(
        candidate_trf_evidence=evidence,
        candidate_biological_annotation=annotation,
        summary=summary,
    )