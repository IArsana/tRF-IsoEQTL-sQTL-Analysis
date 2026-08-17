"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m6/event_reconstruction_readiness.py

Description:
    Core logic for M6.2C.1 FASE/Event Reconstruction Feasibility Audit.

    This stage evaluates whether a retained regulatory event such as
    INT98200 can be mapped back to its source parent gene using an exact or
    source-compatible reconstruction of the FASE event namespace.

    A critical distinction is enforced:

        EVENT IDENTIFIER OBSERVED
            !=
        EVENT NAMESPACE REPRODUCIBLE

    The appearance of INT98200 in Moradi supplementary tables proves that
    the event identifier exists in the published source outputs.

    It does NOT establish that:
        - the source FASE annotation matrix is available;
        - the original event-generation namespace is available;
        - the event-to-gene mapping can be reproduced;
        - INT98200 can be decoded from its numeric suffix.

    M6.2C.1 is an audit only.

    It does NOT:
        - assign a parent gene;
        - decode INT identifiers as genomic coordinates;
        - interpret the numeric suffix as a GTF row number;
        - interpret the numeric suffix as event ordering;
        - use QTL SNP positions as event coordinates;
        - use tag SNP positions as event coordinates;
        - perform nearest-gene mapping;
        - perform approximate GTF-based gene assignment.

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
from typing import Any

import pandas as pd


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class EventReconstructionReadinessResult:
    """Container for M6.2C.1 outputs."""

    components: pd.DataFrame
    feature_readiness: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Generic helpers
# ============================================================================


def _normalize_text(
    value: Any,
) -> str | None:
    """Normalize optional text value."""

    if value is None:
        return None

    try:

        if pd.isna(
            value
        ):
            return None

    except (
        TypeError,
        ValueError,
    ):
        pass

    normalized = str(
        value
    ).strip()

    return normalized or None


def _contains_term(
    text: str,
    term: str,
    *,
    case_sensitive: bool,
) -> bool:
    """Perform configured substring matching."""

    if not case_sensitive:

        text = text.lower()
        term = term.lower()

    return term in text


# ============================================================================
# Upstream feature extraction
# ============================================================================


def extract_reconstruction_targets(
    *,
    m6_2b_qc: dict[str, Any],
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Extract unresolved regulatory features forwarded by M6.2B.
    """

    require_continue = bool(
        config[
            "target_policy"
        ][
            "require_continue_to_m6_2c"
        ]
    )

    supported_classes = {
        str(value)
        for value
        in config[
            "target_policy"
        ][
            "supported_feature_classes"
        ]
    }

    rows: list[
        dict[str, Any]
    ] = []

    for record in m6_2b_qc.get(
        "feature_resolution",
        [],
    ):

        if (
            require_continue
            and
            not bool(
                record.get(
                    "continue_to_m6_2c",
                    False,
                )
            )
        ):
            continue

        feature_id = _normalize_text(
            record.get(
                "regulatory_feature_id"
            )
        )

        feature_class = _normalize_text(
            record.get(
                "feature_class"
            )
        )

        if feature_id is None:
            continue

        if (
            feature_class is not None
            and
            feature_class not in supported_classes
        ):
            continue

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
                    feature_id,

                "feature_class":
                    feature_class,

                "regulatory_scope":
                    _normalize_text(
                        record.get(
                            "regulatory_scope"
                        )
                    ),

                "upstream_gene_annotation_status":
                    _normalize_text(
                        record.get(
                            "gene_annotation_status"
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Local file discovery
# ============================================================================


def discover_local_files(
    *,
    project_root: Path,
    config: dict[str, Any],
) -> list[Path]:
    """
    Collect unique candidate files under configured local search roots.

    These files are used only for generic reconstruction-component discovery.
    Event namespace is audited separately.
    """

    search_config = config[
        "local_search"
    ]

    extensions = {
        str(extension).lower()
        for extension
        in search_config[
            "inspect_extensions"
        ]
    }

    paths: dict[
        str,
        Path,
    ] = {}

    for relative_root in search_config[
        "roots"
    ]:

        root = (
            project_root
            /
            str(
                relative_root
            )
        )

        if not root.exists():
            continue

        for path in root.rglob("*"):

            if not path.is_file():
                continue

            simple_suffix = path.suffix.lower()

            compound_suffix = "".join(
                path.suffixes
            ).lower()

            accepted = (
                simple_suffix in extensions
                or
                compound_suffix in extensions
                or
                (
                    simple_suffix == ".gz"
                    and
                    ".gz" in extensions
                )
            )

            if not accepted:
                continue

            paths[
                str(
                    path.resolve()
                )
            ] = path

    return sorted(
        paths.values(),
        key=lambda item: str(
            item
        ),
    )


def _find_filename_hits(
    *,
    files: list[Path],
    search_terms: list[str],
    case_sensitive: bool,
) -> list[Path]:
    """Find local artifacts using explicit filename/path terms."""

    hits: list[
        Path
    ] = []

    for path in files:

        text = str(
            path
        )

        if any(
            _contains_term(
                text,
                term,
                case_sensitive=case_sensitive,
            )
            for term
            in search_terms
        ):

            hits.append(
                path
            )

    return hits


# ============================================================================
# Event namespace audit
# ============================================================================


def audit_event_namespace(
    *,
    project_root: Path,
    component: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    """
    Audit exact event-namespace evidence for INT98200.

    Presence of the target identifier in Moradi supplementary tables confirms
    only that the event ID was reported.

    It does NOT make the event namespace reproducible.
    """

    target_event = str(
        component[
            "exact_target_event"
        ]
    ).strip()

    source_tables = component.get(
        "source_event_tables",
        [],
    )

    observed_locations: list[
        str
    ] = []

    for source in source_tables:

        path = (
            project_root
            /
            str(
                source[
                    "relative_path"
                ]
            )
        )

        if not path.exists():
            continue

        suffix = path.suffix.lower()

        try:

            if suffix in {
                ".xlsx",
                ".xls",
            }:

                sheet = source.get(
                    "sheet"
                )

                if sheet is None:
                    continue

                frame = pd.read_excel(
                    path,
                    sheet_name=str(
                        sheet
                    ),
                )

            elif suffix == ".csv":

                frame = pd.read_csv(
                    path,
                    low_memory=False,
                )

            elif suffix in {
                ".tsv",
                ".txt",
            }:

                frame = pd.read_csv(
                    path,
                    sep="\t",
                    low_memory=False,
                )

            else:

                continue

        except Exception:

            continue

        feature_column = str(
            source[
                "feature_column"
            ]
        )

        if feature_column not in frame.columns:
            continue

        values = (
            frame[
                feature_column
            ]
            .astype("string")
            .str.strip()
        )

        if not values.eq(
            target_event
        ).any():

            continue

        location = str(
            path
        )

        if source.get(
            "sheet"
        ) is not None:

            location += (
                "::"
                +
                str(
                    source[
                        "sheet"
                    ]
                )
            )

        observed_locations.append(
            location
        )

    states = config[
        "component_states"
    ]

    target_observed = bool(
        observed_locations
    )

    # ----------------------------------------------------------------------
    # Crucial distinction:
    #
    # We can observe INT98200 in source tables while still lacking the
    # original namespace / annotation artifact needed to reproduce how
    # INT98200 maps to its parent gene.
    # ----------------------------------------------------------------------

    if target_observed:

        state = str(
            states[
                "event_observed"
            ]
        )

    else:

        state = str(
            states[
                "not_identified"
            ]
        )

    return {
        "component_id":
            "event_namespace",

        "component_label":
            str(
                component[
                    "label"
                ]
            ),

        "required":
            bool(
                component[
                    "required"
                ]
            ),

        "expected_role":
            str(
                component[
                    "expected_role"
                ]
            ),

        "reproducibility_requirement":
            str(
                component[
                    "reproducibility_requirement"
                ]
            ),

        "local_artifact_available":
            False,

        "local_artifact_count":
            0,

        "local_artifact_paths":
            [],

        "target_event_identifier":
            target_event,

        "target_event_observed":
            target_observed,

        "target_event_observed_locations":
            observed_locations,

        "namespace_reproducible":
            False,

        "component_state":
            state,
    }


# ============================================================================
# Generic component audit
# ============================================================================


def audit_reconstruction_components(
    *,
    project_root: Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Audit reconstruction components.

    Generic filename discovery is used for ordinary reconstruction artifacts.

    event_namespace is handled separately using exact source-table inspection
    to prevent false-positive namespace detection.
    """

    files = discover_local_files(
        project_root=project_root,
        config=config,
    )

    case_sensitive = bool(
        config[
            "local_search"
        ][
            "filename_search_case_sensitive"
        ]
    )

    states = config[
        "component_states"
    ]

    components = config[
        "reconstruction_components"
    ]

    rows: list[
        dict[str, Any]
    ] = []

    # ----------------------------------------------------------------------
    # Methodology is identifiable for these components from the known FASE
    # workflow even when source-specific project artifacts are unavailable.
    # ----------------------------------------------------------------------

    method_identified_components = {
        "fase_implementation",
        "reference_gtf",
        "junction_matrix",
        "read_membership_matrix",
        "intron_membership_matrix",
        "annotation_matrix",
    }

    for component_id, component in components.items():

        # ==================================================================
        # Event namespace requires dedicated exact audit.
        # ==================================================================

        if component_id == "event_namespace":

            rows.append(
                audit_event_namespace(
                    project_root=project_root,
                    component=component,
                    config=config,
                )
            )

            continue

        # ==================================================================
        # Generic reconstruction artifact audit
        # ==================================================================

        search_terms = [
            str(value)
            for value
            in component.get(
                "local_search_terms",
                [],
            )
        ]

        local_extensions = {
            str(value).lower()
            for value
            in component.get(
                "local_extensions",
                [],
            )
        }

        if search_terms:

            hits = _find_filename_hits(
                files=files,
                search_terms=search_terms,
                case_sensitive=case_sensitive,
            )

        elif local_extensions:

            hits = [
                path
                for path
                in files
                if (
                    path.suffix.lower()
                    in local_extensions
                    or
                    "".join(
                        path.suffixes
                    ).lower()
                    in local_extensions
                )
            ]

        else:

            hits = []

        # ------------------------------------------------------------------
        # Exclude generated analysis outputs that contain words such as
        # annotation merely because this pipeline produced them.
        # ------------------------------------------------------------------

        filtered_hits: list[
            Path
        ] = []

        for path in hits:

            path_text = str(
                path
            ).lower()

            if component_id == "annotation_matrix":

                generated_tokens = (
                    "regulatory_gene_annotation",
                    "source_annotation_candidates",
                    "gene_resolution",
                    "event_reconstruction",
                )

                if any(
                    token in path_text
                    for token
                    in generated_tokens
                ):
                    continue

            if component_id == "fase_implementation":

                generated_tokens = (
                    "event_reconstruction",
                    "regulatory_gene_annotation",
                )

                if any(
                    token in path_text
                    for token
                    in generated_tokens
                ):
                    continue

            filtered_hits.append(
                path
            )

        hits = filtered_hits

        # ------------------------------------------------------------------
        # Determine component state.
        # ------------------------------------------------------------------

        if hits:

            state = str(
                states[
                    "available"
                ]
            )

            artifact_available = True

        elif component_id in method_identified_components:

            state = str(
                states[
                    "method_identified"
                ]
            )

            artifact_available = False

        else:

            state = str(
                states[
                    "not_identified"
                ]
            )

            artifact_available = False

        rows.append(
            {
                "component_id":
                    component_id,

                "component_label":
                    str(
                        component[
                            "label"
                        ]
                    ),

                "required":
                    bool(
                        component[
                            "required"
                        ]
                    ),

                "expected_role":
                    str(
                        component[
                            "expected_role"
                        ]
                    ),

                "reproducibility_requirement":
                    str(
                        component[
                            "reproducibility_requirement"
                        ]
                    ),

                "local_artifact_available":
                    artifact_available,

                "local_artifact_count":
                    int(
                        len(
                            hits
                        )
                    ),

                "local_artifact_paths":
                    [
                        str(
                            path
                        )
                        for path
                        in hits
                    ],

                "component_state":
                    state,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Reconstruction classification
# ============================================================================


def classify_reconstruction(
    *,
    components: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[
    str,
    bool,
    bool,
]:
    """
    Classify exact/source-compatible reconstruction feasibility.

    Event-identifier observation alone does not count as an available event
    namespace.
    """

    classification = config[
        "classification"
    ]

    available_state = str(
        config[
            "component_states"
        ][
            "available"
        ]
    )

    required = components.loc[
        components[
            "required"
        ].astype(
            bool
        )
    ]

    all_required_available = bool(
        not required.empty
        and
        (
            required[
                "component_state"
            ]
            ==
            available_state
        ).all()
    )

    exact_allowed = bool(
        all_required_available
    )

    state_lookup = {
        str(
            row[
                "component_id"
            ]
        ):
            str(
                row[
                    "component_state"
                ]
            )
        for row
        in components.to_dict(
            orient="records"
        )
    }

    event_namespace_available = (
        state_lookup.get(
            "event_namespace"
        )
        ==
        available_state
    )

    annotation_matrix_available = (
        state_lookup.get(
            "annotation_matrix"
        )
        ==
        available_state
    )

    reference_gtf_available = (
        state_lookup.get(
            "reference_gtf"
        )
        ==
        available_state
    )

    # ----------------------------------------------------------------------
    # Source-compatible reconstruction requires all three critical identity
    # anchors:
    #
    #   reference GTF
    #   source annotation matrix
    #   reproducible event namespace
    # ----------------------------------------------------------------------

    source_compatible_allowed = bool(
        event_namespace_available
        and
        annotation_matrix_available
        and
        reference_gtf_available
    )

    if exact_allowed:

        return (
            str(
                classification[
                    "exact"
                ]
            ),
            True,
            True,
        )

    if source_compatible_allowed:

        return (
            str(
                classification[
                    "source_compatible"
                ]
            ),
            False,
            True,
        )

    # ----------------------------------------------------------------------
    # Approximate reconstruction is used only when at least one genuine
    # reconstruction anchor exists locally.
    #
    # Merely seeing INT98200 in published source tables does NOT satisfy this
    # condition.
    # ----------------------------------------------------------------------

    if (
        bool(
            config[
                "decision_policy"
            ][
                "approximate_reconstruction_may_use_external_gtf"
            ]
        )
        and
        (
            reference_gtf_available
            or
            annotation_matrix_available
            or
            event_namespace_available
        )
    ):

        return (
            str(
                classification[
                    "approximate_only"
                ]
            ),
            False,
            False,
        )

    return (
        str(
            classification[
                "not_reproducible"
            ]
        ),
        False,
        False,
    )


# ============================================================================
# Feature readiness
# ============================================================================


def build_feature_readiness(
    *,
    targets: pd.DataFrame,
    reconstruction_classification: str,
    exact_allowed: bool,
    source_compatible_allowed: bool,
) -> pd.DataFrame:
    """Apply reconstruction decision to unresolved regulatory features."""

    rows: list[
        dict[str, Any]
    ] = []

    for target in targets.to_dict(
        orient="records"
    ):

        assignment_allowed = bool(
            exact_allowed
            or
            source_compatible_allowed
        )

        if assignment_allowed:

            next_action = (
                "RUN_M6_2C_2_EVENT_GENE_RECONSTRUCTION"
            )

        else:

            next_action = (
                "DO_NOT_ASSIGN_PARENT_GENE_FROM_CURRENT_RECONSTRUCTION"
            )

        rows.append(
            {
                "lead_rsid":
                    target[
                        "lead_rsid"
                    ],

                "priority_rank":
                    target[
                        "priority_rank"
                    ],

                "priority_class":
                    target[
                        "priority_class"
                    ],

                "regulatory_feature_id":
                    target[
                        "regulatory_feature_id"
                    ],

                "feature_class":
                    target[
                        "feature_class"
                    ],

                "upstream_gene_annotation_status":
                    target[
                        "upstream_gene_annotation_status"
                    ],

                "reconstruction_classification":
                    reconstruction_classification,

                "exact_reconstruction_allowed":
                    exact_allowed,

                "source_compatible_reconstruction_allowed":
                    source_compatible_allowed,

                "parent_gene_assignment_from_reconstruction_allowed":
                    assignment_allowed,

                "parent_gene":
                    None,

                "parent_gene_assigned":
                    False,

                "nearest_gene_assignment_performed":
                    False,

                "approximate_mapping_used":
                    False,

                "next_action":
                    next_action,

                "causal_claim_allowed":
                    False,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    components: pd.DataFrame,
    targets: pd.DataFrame,
    feature_readiness: pd.DataFrame,
    reconstruction_classification: str,
    exact_allowed: bool,
    source_compatible_allowed: bool,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M6.2C.1 summary."""

    required = components.loc[
        components[
            "required"
        ].astype(
            bool
        )
    ]

    available_state = str(
        config[
            "component_states"
        ][
            "available"
        ]
    )

    event_observed_state = str(
        config[
            "component_states"
        ][
            "event_observed"
        ]
    )

    required_available = int(
        (
            required[
                "component_state"
            ]
            ==
            available_state
        ).sum()
    )

    event_namespace_row = components.loc[
        components[
            "component_id"
        ]
        ==
        "event_namespace"
    ]

    event_identifier_observed = False
    event_namespace_reproducible = False

    if not event_namespace_row.empty:

        event_record = event_namespace_row.iloc[
            0
        ]

        if "target_event_observed" in event_namespace_row.columns:

            event_identifier_observed = bool(
                event_record.get(
                    "target_event_observed",
                    False,
                )
            )

        event_namespace_reproducible = (
            str(
                event_record[
                    "component_state"
                ]
            )
            ==
            available_state
        )

        # Defensive consistency check.
        if (
            str(
                event_record[
                    "component_state"
                ]
            )
            ==
            event_observed_state
        ):

            event_namespace_reproducible = False

    features_ready = int(
        feature_readiness[
            "parent_gene_assignment_from_reconstruction_allowed"
        ]
        .astype(
            bool
        )
        .sum()
        if not feature_readiness.empty
        else 0
    )

    if (
        exact_allowed
        or
        source_compatible_allowed
    ):

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_exact_or_source_compatible"
            ]
        )

    elif reconstruction_classification == str(
        config[
            "classification"
        ][
            "approximate_only"
        ]
    ):

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_approximate_only"
            ]
        )

    else:

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_not_reproducible"
            ]
        )

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M6.2C.1",

                "features_assessed":
                    int(
                        len(
                            targets
                        )
                    ),

                "reconstruction_components_assessed":
                    int(
                        len(
                            components
                        )
                    ),

                "required_components":
                    int(
                        len(
                            required
                        )
                    ),

                "required_components_available":
                    required_available,

                "required_components_unavailable":
                    int(
                        len(
                            required
                        )
                        -
                        required_available
                    ),

                "target_event_identifier_observed":
                    event_identifier_observed,

                "event_namespace_reproducible":
                    event_namespace_reproducible,

                "reconstruction_classification":
                    reconstruction_classification,

                "exact_reconstruction_allowed":
                    exact_allowed,

                "source_compatible_reconstruction_allowed":
                    source_compatible_allowed,

                "features_ready_for_parent_gene_assignment":
                    features_ready,

                "parent_gene_assignment_performed":
                    False,

                "nearest_gene_mapping_performed":
                    False,

                "approximate_mapping_used":
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


def assess_event_reconstruction_readiness(
    *,
    project_root: Path,
    m6_2b_qc: dict[str, Any],
    config: dict[str, Any],
) -> EventReconstructionReadinessResult:
    """Execute M6.2C.1 event-reconstruction feasibility audit."""

    targets = extract_reconstruction_targets(
        m6_2b_qc=m6_2b_qc,
        config=config,
    )

    if targets.empty:

        raise RuntimeError(
            "M6.2C.1 received no unresolved regulatory features "
            "from M6.2B."
        )

    components = audit_reconstruction_components(
        project_root=project_root,
        config=config,
    )

    (
        reconstruction_classification,
        exact_allowed,
        source_compatible_allowed,
    ) = classify_reconstruction(
        components=components,
        config=config,
    )

    feature_readiness = build_feature_readiness(
        targets=targets,
        reconstruction_classification=reconstruction_classification,
        exact_allowed=exact_allowed,
        source_compatible_allowed=source_compatible_allowed,
    )

    summary = build_summary(
        components=components,
        targets=targets,
        feature_readiness=feature_readiness,
        reconstruction_classification=reconstruction_classification,
        exact_allowed=exact_allowed,
        source_compatible_allowed=source_compatible_allowed,
        config=config,
    )

    return EventReconstructionReadinessResult(
        components=components,
        feature_readiness=feature_readiness,
        summary=summary,
    )