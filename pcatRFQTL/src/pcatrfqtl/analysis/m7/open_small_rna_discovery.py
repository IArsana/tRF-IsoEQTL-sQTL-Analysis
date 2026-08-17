"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/open_small_rna_discovery.py

Description:
    Core logic for M7.4A Open/Public Small-RNA Dataset Discovery.

    This stage identifies public human prostate cancer small-RNA /
    non-coding-RNA sequencing datasets that may support independent
    external validation of the prioritized tRF candidate.

    Dataset metadata are read from the canonical nested schema in
    configs/datasets.yaml.

    M7.4A performs discovery-level eligibility and hierarchical
    prioritization only.

    It does NOT:
        - download FASTQ or SRA sequencing files;
        - inspect raw sequencing reads;
        - search for the candidate tRF sequence;
        - quantify candidate tRF abundance;
        - establish technical quantification eligibility;
        - perform differential expression;
        - perform clinical association;
        - interpret missing explicit tRF annotation as biological absence;
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


# ============================================================================
# Result model
# ============================================================================


@dataclass(frozen=True)
class OpenSmallRnaDiscoveryResult:
    """Container for M7.4A output tables."""

    dataset_discovery: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Generic helpers
# ============================================================================


def _normalize_text(
    value: Any,
) -> str | None:
    """Normalize optional scalar text."""

    if value is None:
        return None

    try:
        if pd.isna(value):
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


def _as_bool(
    value: Any,
) -> bool:
    """
    Convert a value conservatively to bool.

    Missing values evaluate to False.
    """

    if value is None:
        return False

    try:
        if pd.isna(value):
            return False

    except (
        TypeError,
        ValueError,
    ):
        pass

    if isinstance(
        value,
        bool,
    ):
        return value

    if isinstance(
        value,
        str,
    ):

        normalized = value.strip().lower()

        if normalized in {
            "true",
            "yes",
            "1",
        }:
            return True

        if normalized in {
            "false",
            "no",
            "0",
            "",
        }:
            return False

    return bool(
        value
    )


def _as_int(
    value: Any,
    *,
    default: int = 0,
) -> int:
    """Convert optional numeric value to integer safely."""

    if value is None:
        return default

    try:
        if pd.isna(value):
            return default

    except (
        TypeError,
        ValueError,
    ):
        pass

    try:
        return int(
            value
        )

    except (
        TypeError,
        ValueError,
    ):
        return default


def _nested_mapping(
    mapping: dict[str, Any],
    key: str,
) -> dict[str, Any]:
    """
    Return a nested mapping or an empty dict.

    This keeps the core tolerant to optional nested sections while
    preserving strict validation for fields needed for classification.
    """

    value = mapping.get(
        key
    )

    if isinstance(
        value,
        dict,
    ):
        return value

    return {}


def _first_truthy(
    values: list[Any],
) -> bool:
    """Return True when at least one supplied value is explicitly truthy."""

    return any(
        _as_bool(
            value
        )
        for value in values
    )


# ============================================================================
# Registry extraction
# ============================================================================


def extract_registry_datasets(
    *,
    dataset_registry: dict[str, Any],
    dataset_keys: list[str],
) -> list[dict[str, Any]]:
    """
    Extract requested datasets from the canonical dataset registry.

    The supplied dataset_registry must already correspond to the top-level
    "datasets" mapping from configs/datasets.yaml.
    """

    if not isinstance(
        dataset_registry,
        dict,
    ):
        raise RuntimeError(
            "M7.4A dataset registry must be a mapping."
        )

    if not isinstance(
        dataset_keys,
        list,
    ) or not dataset_keys:

        raise RuntimeError(
            "M7.4A requires at least one dataset key."
        )

    records: list[
        dict[str, Any]
    ] = []

    for dataset_key in dataset_keys:

        if dataset_key not in dataset_registry:

            raise RuntimeError(
                "M7.4A dataset key missing from canonical registry: "
                f"{dataset_key}"
            )

        dataset = dataset_registry[
            dataset_key
        ]

        if not isinstance(
            dataset,
            dict,
        ):

            raise RuntimeError(
                f"Dataset {dataset_key!r} must be a mapping."
            )

        records.append(
            {
                "dataset_key":
                    dataset_key,

                "dataset":
                    dataset,
            }
        )

    return records


# ============================================================================
# Nested dataset parsing
# ============================================================================


def _parse_dataset(
    *,
    dataset_key: str,
    dataset: dict[str, Any],
) -> dict[str, Any]:
    """
    Parse one nested datasets.yaml entry into analysis-ready metadata.

    This function performs no priority classification.
    """

    source = _nested_mapping(
        dataset,
        "source",
    )

    organism = _nested_mapping(
        dataset,
        "organism",
    )

    disease = _nested_mapping(
        dataset,
        "disease",
    )

    experiment = _nested_mapping(
        dataset,
        "experiment",
    )

    cohort = _nested_mapping(
        dataset,
        "cohort",
    )

    sequencing = _nested_mapping(
        dataset,
        "sequencing",
    )

    m7 = _nested_mapping(
        dataset,
        "m7",
    )

    discovery_stage = _nested_mapping(
        m7,
        "discovery_stage",
    )

    technical_stage = _nested_mapping(
        m7,
        "technical_audit_stage",
    )

    interpretation = _nested_mapping(
        dataset,
        "interpretation",
    )

    provenance = _nested_mapping(
        dataset,
        "provenance",
    )

    intended_use = dataset.get(
        "intended_use",
        [],
    )

    if not isinstance(
        intended_use,
        list,
    ):
        intended_use = []

    dataset_id = (
        _normalize_text(
            source.get(
                "accession"
            )
        )
        or
        _normalize_text(
            dataset.get(
                "id"
            )
        )
        or
        dataset_key
    )

    sample_count = _as_int(
        cohort.get(
            "sample_count_series"
        )
    )

    patient_count_raw = cohort.get(
        "patient_count_reported"
    )

    patient_count = (
        _as_int(
            patient_count_raw
        )
        if patient_count_raw is not None
        else None
    )

    return {
        "dataset_key":
            dataset_key,

        "registry_id":
            _normalize_text(
                dataset.get(
                    "id"
                )
            ),

        "dataset_id":
            dataset_id,

        "name":
            _normalize_text(
                dataset.get(
                    "name"
                )
            ),

        "description":
            _normalize_text(
                dataset.get(
                    "description"
                )
            ),

        # ------------------------------------------------------------------
        # Source
        # ------------------------------------------------------------------

        "source_organization":
            _normalize_text(
                source.get(
                    "organization"
                )
            ),

        "repository":
            _normalize_text(
                source.get(
                    "database"
                )
            ),

        "bioproject":
            _normalize_text(
                source.get(
                    "bioproject"
                )
            ),

        "sra_study":
            _normalize_text(
                source.get(
                    "sra_study"
                )
            ),

        "source_data_type":
            _normalize_text(
                source.get(
                    "data_type"
                )
            ),

        "source_access_type":
            _normalize_text(
                source.get(
                    "access_type"
                )
            ),

        "source_url":
            _normalize_text(
                source.get(
                    "source_url"
                )
            ),

        # ------------------------------------------------------------------
        # Biological context
        # ------------------------------------------------------------------

        "organism":
            _normalize_text(
                organism.get(
                    "species"
                )
            ),

        "taxon_id":
            organism.get(
                "taxon_id"
            ),

        "disease_context":
            _normalize_text(
                disease.get(
                    "primary"
                )
            ),

        # ------------------------------------------------------------------
        # Experiment
        # ------------------------------------------------------------------

        "experiment_type":
            _normalize_text(
                experiment.get(
                    "experiment_type"
                )
            ),

        "library_scope":
            _normalize_text(
                experiment.get(
                    "library_scope"
                )
            ),

        "small_rna_library":
            _as_bool(
                experiment.get(
                    "small_rna_library"
                )
            ),

        "trna_derived_rna_reported":
            _as_bool(
                experiment.get(
                    "trna_derived_rna_reported"
                )
            ),

        "raw_sequence_available":
            _as_bool(
                experiment.get(
                    "raw_sequence_available"
                )
            ),

        "raw_sequence_repository":
            _normalize_text(
                experiment.get(
                    "raw_sequence_repository"
                )
            ),

        "raw_sequence_access_type":
            _normalize_text(
                experiment.get(
                    "raw_sequence_access_type"
                )
            ),

        # ------------------------------------------------------------------
        # Cohort
        # ------------------------------------------------------------------

        "sample_count_series":
            sample_count,

        "patient_count_reported":
            patient_count,

        "prostate_tissue":
            _as_bool(
                cohort.get(
                    "prostate_tissue"
                )
            ),

        "tumor_samples_present":
            _as_bool(
                cohort.get(
                    "tumor_samples_present"
                )
            ),

        "normal_or_benign_samples_present":
            _as_bool(
                cohort.get(
                    "normal_or_benign_samples_present"
                )
            ),

        "adjacent_normal_samples_present":
            _as_bool(
                cohort.get(
                    "adjacent_normal_samples_present"
                )
            ),

        "benign_samples_present":
            _as_bool(
                cohort.get(
                    "benign_samples_present"
                )
            ),

        "clinical_metadata_reported":
            _as_bool(
                cohort.get(
                    "clinical_metadata_reported"
                )
            ),

        "extensive_clinicopathologic_metadata_reported":
            _as_bool(
                cohort.get(
                    "extensive_clinicopathologic_metadata_reported"
                )
            ),

        "biochemical_recurrence_metadata_reported":
            _as_bool(
                cohort.get(
                    "biochemical_recurrence_metadata_reported"
                )
            ),

        "multiple_disease_stages_reported":
            _as_bool(
                cohort.get(
                    "multiple_disease_stages_reported"
                )
            ),

        "ffpe":
            _as_bool(
                cohort.get(
                    "ffpe"
                )
            ),

        # ------------------------------------------------------------------
        # Sequencing
        # ------------------------------------------------------------------

        "platform":
            _normalize_text(
                sequencing.get(
                    "platform"
                )
            ),

        "instrument":
            _normalize_text(
                sequencing.get(
                    "instrument"
                )
            ),

        "size_fractionated_library_reported":
            _as_bool(
                sequencing.get(
                    "size_fractionated_library_reported"
                )
            ),

        "read_length_compatibility_verified":
            _as_bool(
                sequencing.get(
                    "read_length_compatibility_verified"
                )
            ),

        "adapter_information_verified":
            _as_bool(
                sequencing.get(
                    "adapter_information_verified"
                )
            ),

        "read_structure_verified":
            _as_bool(
                sequencing.get(
                    "read_structure_verified"
                )
            ),

        # ------------------------------------------------------------------
        # Registry intent / M7 state
        # ------------------------------------------------------------------

        "intended_use":
            [
                str(
                    value
                )
                for value in intended_use
            ],

        "registry_discovery_eligible":
            _as_bool(
                discovery_stage.get(
                    "discovery_eligible"
                )
            ),

        "registry_discovery_milestone":
            _normalize_text(
                discovery_stage.get(
                    "milestone"
                )
            ),

        "registry_technical_audit_status":
            _normalize_text(
                technical_stage.get(
                    "status"
                )
            ),

        # ------------------------------------------------------------------
        # Interpretation / provenance
        # ------------------------------------------------------------------

        "explicit_trna_derived_rna_context":
            _as_bool(
                interpretation.get(
                    "explicit_trna_derived_rna_context"
                )
            ),

        "candidate_trf_detected_registry":
            _as_bool(
                interpretation.get(
                    "candidate_trf_detected"
                )
            ),

        "candidate_trf_absent_registry":
            _as_bool(
                interpretation.get(
                    "candidate_trf_absent"
                )
            ),

        "direct_quantification_claimed_registry":
            _as_bool(
                interpretation.get(
                    "direct_quantification_claimed"
                )
            ),

        "clinical_validation_claimed_registry":
            _as_bool(
                interpretation.get(
                    "clinical_validation_claimed"
                )
            ),

        "evidence_status":
            _normalize_text(
                provenance.get(
                    "evidence_status"
                )
            ),
    }


# ============================================================================
# Discovery eligibility
# ============================================================================


def _evaluate_discovery_eligibility(
    *,
    parsed: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, bool]:
    """
    Evaluate M7.4A discovery inclusion criteria.

    The criteria determine only whether a dataset proceeds to M7.4B.
    """

    criteria = config[
        "inclusion_criteria"
    ]

    organism = (
        parsed.get(
            "organism"
        )
        or ""
    ).strip().lower()

    disease = (
        parsed.get(
            "disease_context"
        )
        or ""
    ).strip().lower()

    source_access_type = (
        parsed.get(
            "source_access_type"
        )
        or ""
    ).strip().lower()

    raw_access_type = (
        parsed.get(
            "raw_sequence_access_type"
        )
        or ""
    ).strip().lower()

    human_pass = bool(
        organism
        ==
        "homo sapiens"
    )

    prostate_pass = bool(
        parsed[
            "prostate_tissue"
        ]
    )

    prostate_cancer_pass = bool(
        "prostate cancer"
        in disease
    )

    small_rna_pass = bool(
        parsed[
            "small_rna_library"
        ]
    )

    public_record_pass = bool(
        parsed[
            "repository"
        ]
        and
        source_access_type
        ==
        "public"
    )

    raw_sequence_route_pass = bool(
        parsed[
            "raw_sequence_available"
        ]
        and
        parsed[
            "raw_sequence_repository"
        ]
        and
        raw_access_type
        ==
        "public"
    )

    registry_discovery_pass = bool(
        parsed[
            "registry_discovery_eligible"
        ]
    )

    checks = {
        "human_criterion_pass":
            human_pass,

        "prostate_context_criterion_pass":
            prostate_pass,

        "prostate_cancer_criterion_pass":
            prostate_cancer_pass,

        "small_rna_sequencing_criterion_pass":
            small_rna_pass,

        "public_dataset_record_criterion_pass":
            public_record_pass,

        "raw_sequence_route_criterion_pass":
            raw_sequence_route_pass,

        "registry_discovery_eligible":
            registry_discovery_pass,
    }

    required_passes: list[
        bool
    ] = []

    if _as_bool(
        criteria.get(
            "require_human"
        )
    ):
        required_passes.append(
            human_pass
        )

    if _as_bool(
        criteria.get(
            "require_prostate_context"
        )
    ):
        required_passes.append(
            prostate_pass
        )

    if _as_bool(
        criteria.get(
            "require_prostate_cancer_context"
        )
    ):
        required_passes.append(
            prostate_cancer_pass
        )

    if _as_bool(
        criteria.get(
            "require_high_throughput_small_or_non_coding_rna_sequencing"
        )
    ):
        required_passes.append(
            small_rna_pass
        )

    if _as_bool(
        criteria.get(
            "require_public_dataset_record"
        )
    ):
        required_passes.append(
            public_record_pass
        )

    if _as_bool(
        criteria.get(
            "require_raw_sequence_route"
        )
    ):
        required_passes.append(
            raw_sequence_route_pass
        )

    # Registry eligibility is additionally required because datasets.yaml
    # is the canonical source of discovery state.
    required_passes.append(
        registry_discovery_pass
    )

    checks[
        "discovery_eligible"
    ] = bool(
        all(
            required_passes
        )
    )

    return checks


# ============================================================================
# Priority classification
# ============================================================================


def _classify_priority(
    *,
    parsed: dict[str, Any],
    discovery_eligible: bool,
    config: dict[str, Any],
) -> tuple[
    int | None,
    str,
    str,
    str,
]:
    """
    Apply hierarchical M7.4A prioritization.

    Returns:
        priority_rank,
        priority_class,
        prioritization_basis,
        discovery_status

    No arbitrary additive score is used.
    """

    statuses = config[
        "statuses"
    ][
        "dataset"
    ]

    if not discovery_eligible:

        return (
            None,
            "DISCOVERY_EXCLUDED",
            "FAILED_DISCOVERY_INCLUSION_CRITERIA",
            str(
                statuses[
                    "discovery_excluded"
                ]
            ),
        )

    explicit_trna_context = bool(
        parsed[
            "trna_derived_rna_reported"
        ]
        or
        parsed[
            "explicit_trna_derived_rna_context"
        ]
    )

    comparator_available = _first_truthy(
        [
            parsed[
                "normal_or_benign_samples_present"
            ],
            parsed[
                "adjacent_normal_samples_present"
            ],
            parsed[
                "benign_samples_present"
            ],
        ]
    )

    clinicopathologic_metadata = bool(
        parsed[
            "clinical_metadata_reported"
        ]
        or
        parsed[
            "extensive_clinicopathologic_metadata_reported"
        ]
    )

    recurrence_metadata = bool(
        parsed[
            "biochemical_recurrence_metadata_reported"
        ]
    )

    large_cohort = bool(
        parsed[
            "sample_count_series"
        ]
        >=
        50
    )

    raw_route = bool(
        parsed[
            "raw_sequence_available"
        ]
    )

    # ----------------------------------------------------------------------
    # Priority 1:
    # explicit tRNA-derived RNA context.
    # ----------------------------------------------------------------------

    if explicit_trna_context:

        return (
            1,
            "HIGH_PRIORITY",
            "EXPLICIT_TRNA_DERIVED_RNA_CONTEXT",
            str(
                statuses[
                    "high_priority"
                ]
            ),
        )

    # ----------------------------------------------------------------------
    # Priority 2:
    # large clinical cohort with recurrence data.
    # ----------------------------------------------------------------------

    if (
        raw_route
        and
        clinicopathologic_metadata
        and
        recurrence_metadata
        and
        large_cohort
    ):

        return (
            2,
            "HIGH_PRIORITY",
            "LARGE_CLINICAL_COHORT_WITH_RECURRENCE_METADATA",
            str(
                statuses[
                    "high_priority"
                ]
            ),
        )

    # ----------------------------------------------------------------------
    # Priority 3:
    # large cohort with tumor comparator.
    # ----------------------------------------------------------------------

    if (
        raw_route
        and
        comparator_available
        and
        large_cohort
    ):

        return (
            3,
            "INTERMEDIATE_PRIORITY",
            "LARGE_COHORT_WITH_TUMOR_COMPARATOR",
            str(
                statuses[
                    "intermediate_priority"
                ]
            ),
        )

    # ----------------------------------------------------------------------
    # Priority 4:
    # discovery-eligible but primarily exploratory.
    # ----------------------------------------------------------------------

    return (
        4,
        "EXPLORATORY",
        "PUBLIC_RAW_SMALL_RNA_ROUTE_WITH_LIMITED_VALIDATION_DEPTH",
        str(
            statuses[
                "exploratory"
            ]
        ),
    )


# ============================================================================
# Per-dataset classification
# ============================================================================


def classify_dataset(
    *,
    dataset_key: str,
    dataset: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    """
    Classify one canonical nested dataset entry.

    This produces discovery-level evidence only.
    """

    parsed = _parse_dataset(
        dataset_key=dataset_key,
        dataset=dataset,
    )

    eligibility = _evaluate_discovery_eligibility(
        parsed=parsed,
        config=config,
    )

    (
        priority_rank,
        priority_class,
        prioritization_basis,
        discovery_status,
    ) = _classify_priority(
        parsed=parsed,
        discovery_eligible=eligibility[
            "discovery_eligible"
        ],
        config=config,
    )

    comparator_available = _first_truthy(
        [
            parsed[
                "normal_or_benign_samples_present"
            ],
            parsed[
                "adjacent_normal_samples_present"
            ],
            parsed[
                "benign_samples_present"
            ],
        ]
    )

    clinicopathologic_metadata = bool(
        parsed[
            "clinical_metadata_reported"
        ]
        or
        parsed[
            "extensive_clinicopathologic_metadata_reported"
        ]
    )

    explicit_trna_context = bool(
        parsed[
            "trna_derived_rna_reported"
        ]
        or
        parsed[
            "explicit_trna_derived_rna_context"
        ]
    )

    large_cohort = bool(
        parsed[
            "sample_count_series"
        ]
        >=
        50
    )

    # ----------------------------------------------------------------------
    # Registry safeguards
    #
    # At M7.4A, these must remain false.
    # ----------------------------------------------------------------------

    if parsed[
        "candidate_trf_detected_registry"
    ]:

        raise RuntimeError(
            f"{parsed['dataset_id']} registry already claims candidate "
            "tRF detection before M7.4B/M7.4C."
        )

    if parsed[
        "candidate_trf_absent_registry"
    ]:

        raise RuntimeError(
            f"{parsed['dataset_id']} registry claims candidate tRF absence "
            "before raw-read analysis."
        )

    if parsed[
        "direct_quantification_claimed_registry"
    ]:

        raise RuntimeError(
            f"{parsed['dataset_id']} registry claims direct quantification "
            "before M7.4C."
        )

    if parsed[
        "clinical_validation_claimed_registry"
    ]:

        raise RuntimeError(
            f"{parsed['dataset_id']} registry claims clinical validation "
            "before M7.4D."
        )

    return {
        # ------------------------------------------------------------------
        # Identity
        # ------------------------------------------------------------------

        "dataset_key":
            parsed[
                "dataset_key"
            ],

        "registry_id":
            parsed[
                "registry_id"
            ],

        "dataset_id":
            parsed[
                "dataset_id"
            ],

        "name":
            parsed[
                "name"
            ],

        "repository":
            parsed[
                "repository"
            ],

        "source_organization":
            parsed[
                "source_organization"
            ],

        "source_url":
            parsed[
                "source_url"
            ],

        "bioproject":
            parsed[
                "bioproject"
            ],

        "sra_study":
            parsed[
                "sra_study"
            ],

        "evidence_status":
            parsed[
                "evidence_status"
            ],

        # ------------------------------------------------------------------
        # Biological / cohort context
        # ------------------------------------------------------------------

        "organism":
            parsed[
                "organism"
            ],

        "disease_context":
            parsed[
                "disease_context"
            ],

        "sample_count_series":
            parsed[
                "sample_count_series"
            ],

        "patient_count_reported":
            parsed[
                "patient_count_reported"
            ],

        "prostate_tissue":
            parsed[
                "prostate_tissue"
            ],

        "tumor_samples_present":
            parsed[
                "tumor_samples_present"
            ],

        "tumor_comparator_available":
            comparator_available,

        "adjacent_normal_samples_present":
            parsed[
                "adjacent_normal_samples_present"
            ],

        "benign_samples_present":
            parsed[
                "benign_samples_present"
            ],

        "clinicopathologic_metadata_reported":
            clinicopathologic_metadata,

        "recurrence_metadata_reported":
            parsed[
                "biochemical_recurrence_metadata_reported"
            ],

        "multiple_disease_stages_reported":
            parsed[
                "multiple_disease_stages_reported"
            ],

        "ffpe":
            parsed[
                "ffpe"
            ],

        "large_cohort":
            large_cohort,

        # ------------------------------------------------------------------
        # Sequencing
        # ------------------------------------------------------------------

        "experiment_type":
            parsed[
                "experiment_type"
            ],

        "library_scope":
            parsed[
                "library_scope"
            ],

        "small_rna_library":
            parsed[
                "small_rna_library"
            ],

        "platform":
            parsed[
                "platform"
            ],

        "instrument":
            parsed[
                "instrument"
            ],

        "size_fractionated_library_reported":
            parsed[
                "size_fractionated_library_reported"
            ],

        "raw_sequence_available":
            parsed[
                "raw_sequence_available"
            ],

        "raw_sequence_repository":
            parsed[
                "raw_sequence_repository"
            ],

        "raw_sequence_access_type":
            parsed[
                "raw_sequence_access_type"
            ],

        "raw_sequence_route":
            eligibility[
                "raw_sequence_route_criterion_pass"
            ],

        # ------------------------------------------------------------------
        # tRNA context
        # ------------------------------------------------------------------

        "trna_derived_rna_context_reported":
            explicit_trna_context,

        "trna_context_not_reported_does_not_mean_absent":
            True,

        # ------------------------------------------------------------------
        # Inclusion criteria
        # ------------------------------------------------------------------

        "human_criterion_pass":
            eligibility[
                "human_criterion_pass"
            ],

        "prostate_context_criterion_pass":
            eligibility[
                "prostate_context_criterion_pass"
            ],

        "prostate_cancer_criterion_pass":
            eligibility[
                "prostate_cancer_criterion_pass"
            ],

        "small_rna_sequencing_criterion_pass":
            eligibility[
                "small_rna_sequencing_criterion_pass"
            ],

        "public_dataset_record_criterion_pass":
            eligibility[
                "public_dataset_record_criterion_pass"
            ],

        "raw_sequence_route_criterion_pass":
            eligibility[
                "raw_sequence_route_criterion_pass"
            ],

        "registry_discovery_eligible":
            eligibility[
                "registry_discovery_eligible"
            ],

        "discovery_eligible":
            eligibility[
                "discovery_eligible"
            ],

        # ------------------------------------------------------------------
        # Priority
        # ------------------------------------------------------------------

        "priority_rank":
            priority_rank,

        "priority_class":
            priority_class,

        "prioritization_basis":
            prioritization_basis,

        "discovery_status":
            discovery_status,

        # ------------------------------------------------------------------
        # Technical audit state
        #
        # These values deliberately remain unresolved until M7.4B.
        # ------------------------------------------------------------------

        "read_length_compatibility_verified":
            parsed[
                "read_length_compatibility_verified"
            ],

        "adapter_information_verified":
            parsed[
                "adapter_information_verified"
            ],

        "read_structure_verified":
            parsed[
                "read_structure_verified"
            ],

        "technical_eligibility_established":
            False,

        # ------------------------------------------------------------------
        # M7.4A safeguards
        # ------------------------------------------------------------------

        "fastq_download_performed":
            False,

        "candidate_sequence_search_performed":
            False,

        "candidate_detected":
            False,

        "candidate_absent":
            False,

        "trf_quantification_performed":
            False,

        "differential_expression_performed":
            False,

        "clinical_association_performed":
            False,

        "causal_inference_performed":
            False,
    }


# ============================================================================
# Discovery table
# ============================================================================


def build_dataset_discovery(
    *,
    dataset_registry: dict[str, Any],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build the M7.4A dataset discovery table."""

    dataset_keys = list(
        config[
            "dataset_registry"
        ][
            "dataset_keys"
        ]
    )

    registry_records = extract_registry_datasets(
        dataset_registry=dataset_registry,
        dataset_keys=dataset_keys,
    )

    rows = [
        classify_dataset(
            dataset_key=record[
                "dataset_key"
            ],
            dataset=record[
                "dataset"
            ],
            config=config,
        )
        for record in registry_records
    ]

    result = pd.DataFrame(
        rows
    )

    if result.empty:

        raise RuntimeError(
            "M7.4A produced an empty dataset-discovery table."
        )

    # ----------------------------------------------------------------------
    # Identity validation
    # ----------------------------------------------------------------------

    if result[
        "dataset_key"
    ].isna().any():

        raise RuntimeError(
            "M7.4A discovery output contains missing dataset keys."
        )

    if result[
        "dataset_key"
    ].duplicated().any():

        duplicates = (
            result.loc[
                result[
                    "dataset_key"
                ].duplicated(
                    keep=False
                ),
                "dataset_key",
            ]
            .astype(
                str
            )
            .unique()
            .tolist()
        )

        raise RuntimeError(
            "M7.4A duplicate dataset keys detected: "
            f"{duplicates}"
        )

    if result[
        "dataset_id"
    ].isna().any():

        raise RuntimeError(
            "M7.4A discovery output contains missing dataset IDs."
        )

    # ----------------------------------------------------------------------
    # Discovery-only safeguards
    # ----------------------------------------------------------------------

    forbidden_true_columns = [
        "technical_eligibility_established",
        "fastq_download_performed",
        "candidate_sequence_search_performed",
        "candidate_detected",
        "candidate_absent",
        "trf_quantification_performed",
        "differential_expression_performed",
        "clinical_association_performed",
        "causal_inference_performed",
    ]

    for column in forbidden_true_columns:

        if result[
            column
        ].fillna(
            False
        ).astype(
            bool
        ).any():

            raise RuntimeError(
                "M7.4A discovery-only safeguard violated: "
                f"{column}=True."
            )

    # ----------------------------------------------------------------------
    # Stable ordering
    # ----------------------------------------------------------------------

    return result.sort_values(
        [
            "discovery_eligible",
            "priority_rank",
            "sample_count_series",
            "dataset_id",
        ],
        ascending=[
            False,
            True,
            False,
            True,
        ],
        kind="stable",
        na_position="last",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    discovery: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build the one-row M7.4A summary artifact."""

    eligible = discovery.loc[
        discovery[
            "discovery_eligible"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
    ].copy()

    high_priority = int(
        (
            eligible[
                "priority_class"
            ]
            ==
            "HIGH_PRIORITY"
        ).sum()
    )

    intermediate_priority = int(
        (
            eligible[
                "priority_class"
            ]
            ==
            "INTERMEDIATE_PRIORITY"
        ).sum()
    )

    exploratory = int(
        (
            eligible[
                "priority_class"
            ]
            ==
            "EXPLORATORY"
        ).sum()
    )

    raw_routes = int(
        discovery[
            "raw_sequence_route"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    explicit_trna_contexts = int(
        discovery[
            "trna_derived_rna_context_reported"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    recurrence_datasets = int(
        discovery[
            "recurrence_metadata_reported"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    comparator_datasets = int(
        discovery[
            "tumor_comparator_available"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    ffpe_datasets = int(
        discovery[
            "ffpe"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    if not eligible.empty:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall"
            ][
                "candidates_identified"
            ]
        )

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_candidates_identified"
            ]
        )

        leading_dataset = _normalize_text(
            eligible.iloc[
                0
            ][
                "dataset_id"
            ]
        )

        leading_dataset_key = _normalize_text(
            eligible.iloc[
                0
            ][
                "dataset_key"
            ]
        )

    else:

        overall_status = str(
            config[
                "statuses"
            ][
                "overall"
            ][
                "none_identified"
            ]
        )

        next_stage = str(
            config[
                "next_stage"
            ][
                "if_none_identified"
            ]
        )

        leading_dataset = None
        leading_dataset_key = None

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.4A",

                "candidate_trf":
                    _normalize_text(
                        config[
                            "discovery_scope"
                        ].get(
                            "target_candidate"
                        )
                    ),

                "candidate_sequence":
                    _normalize_text(
                        config[
                            "discovery_scope"
                        ].get(
                            "exact_candidate_sequence"
                        )
                    ),

                "candidate_sequence_length_nt":
                    _as_int(
                        config[
                            "discovery_scope"
                        ].get(
                            "candidate_sequence_length_nt"
                        )
                    ),

                # ----------------------------------------------------------
                # Dataset counts
                # ----------------------------------------------------------

                "datasets_assessed":
                    int(
                        len(
                            discovery
                        )
                    ),

                "datasets_discovery_eligible":
                    int(
                        len(
                            eligible
                        )
                    ),

                "datasets_discovery_excluded":
                    int(
                        len(
                            discovery
                        )
                        -
                        len(
                            eligible
                        )
                    ),

                "high_priority_datasets":
                    high_priority,

                "intermediate_priority_datasets":
                    intermediate_priority,

                "exploratory_datasets":
                    exploratory,

                # ----------------------------------------------------------
                # Evidence characteristics
                # ----------------------------------------------------------

                "datasets_with_raw_sequence_route":
                    raw_routes,

                "datasets_with_trna_context":
                    explicit_trna_contexts,

                "datasets_with_tumor_comparator":
                    comparator_datasets,

                "datasets_with_recurrence_metadata":
                    recurrence_datasets,

                "ffpe_datasets":
                    ffpe_datasets,

                # ----------------------------------------------------------
                # Leading dataset
                # ----------------------------------------------------------

                "leading_dataset_key_for_technical_audit":
                    leading_dataset_key,

                "leading_dataset_for_technical_audit":
                    leading_dataset,

                # ----------------------------------------------------------
                # M7.4A safeguards
                # ----------------------------------------------------------

                "technical_eligibility_established":
                    False,

                "raw_fastq_download_performed":
                    False,

                "candidate_sequence_search_performed":
                    False,

                "candidate_detection_claimed":
                    False,

                "candidate_absence_claimed":
                    False,

                "trf_quantification_performed":
                    False,

                "differential_expression_performed":
                    False,

                "clinical_association_performed":
                    False,

                "causal_inference_performed":
                    False,

                # ----------------------------------------------------------
                # Resolution
                # ----------------------------------------------------------

                "overall_discovery_status":
                    overall_status,

                "next_stage":
                    next_stage,
            }
        ]
    )


# ============================================================================
# Final consistency checks
# ============================================================================


def _validate_summary_consistency(
    *,
    discovery: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    """Validate consistency between dataset-level and summary artifacts."""

    if summary.empty or len(
        summary
    ) != 1:

        raise RuntimeError(
            "M7.4A summary must contain exactly one row."
        )

    row = summary.iloc[
        0
    ]

    if int(
        row[
            "datasets_assessed"
        ]
    ) != len(
        discovery
    ):

        raise RuntimeError(
            "M7.4A summary dataset count is inconsistent "
            "with discovery table."
        )

    observed_eligible = int(
        discovery[
            "discovery_eligible"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    if int(
        row[
            "datasets_discovery_eligible"
        ]
    ) != observed_eligible:

        raise RuntimeError(
            "M7.4A eligible-dataset count is inconsistent."
        )

    forbidden_true_fields = [
        "technical_eligibility_established",
        "raw_fastq_download_performed",
        "candidate_sequence_search_performed",
        "candidate_detection_claimed",
        "candidate_absence_claimed",
        "trf_quantification_performed",
        "differential_expression_performed",
        "clinical_association_performed",
        "causal_inference_performed",
    ]

    for field in forbidden_true_fields:

        if bool(
            row[
                field
            ]
        ):

            raise RuntimeError(
                "M7.4A summary violates discovery-only safeguard: "
                f"{field}=True."
            )


# ============================================================================
# Public API
# ============================================================================


def discover_open_small_rna_datasets(
    *,
    dataset_registry: dict[str, Any],
    config: dict[str, Any],
) -> OpenSmallRnaDiscoveryResult:
    """
    Execute M7.4A Open/Public Small-RNA Dataset Discovery.

    Parameters
    ----------
    dataset_registry:
        The top-level "datasets" mapping extracted from configs/datasets.yaml.

    config:
        Parsed configs/m7_open_small_rna_discovery.yaml.

    Returns
    -------
    OpenSmallRnaDiscoveryResult
        Dataset-level discovery table and one-row summary.
    """

    discovery = build_dataset_discovery(
        dataset_registry=dataset_registry,
        config=config,
    )

    summary = build_summary(
        discovery=discovery,
        config=config,
    )

    _validate_summary_consistency(
        discovery=discovery,
        summary=summary,
    )

    return OpenSmallRnaDiscoveryResult(
        dataset_discovery=discovery,
        summary=summary,
    )