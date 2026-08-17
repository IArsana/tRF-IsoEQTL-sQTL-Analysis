"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m7/variant_validation_feasibility.py

Description:
    Core logic for M7.5 Variant-Level Validation Feasibility.

    M7.5 audits literature-curated exact-variant associations returned by the
    NHGRI-EBI GWAS Catalog REST API V2.

    The live API and the published API reference may expose slightly different
    response naming conventions. This module therefore accepts both:

        Reference-style:
            content
            links
            page.total_elements
            page.total_pages

        Live HAL-style:
            _embedded.associations
            _links
            page.totalElements
            page.totalPages

    Both representations are normalized internally to:

        content
        page.size
        page.total_elements
        page.total_pages
        page.number

    Scientific interpretation remains conservative:

        - successful zero-result queries are not negative genetic evidence;
        - API service failure is not negative genetic evidence;
        - absence of a curated association is not variant absence;
        - exact rsID is confirmed from snp_allele[].rs_id when association
          records are returned;
        - previously used GWAS studies are not counted as new validation;
        - non-reused GWAS studies do not establish independent replication;
        - full summary statistics are not analyzed in M7.5;
        - TCGA somatic calls cannot substitute germline genotype.

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
# Result container
# ============================================================================


@dataclass(frozen=True)
class VariantValidationFeasibilityResult:
    """Container for processed M7.5 results."""

    acquisition_status: pd.DataFrame
    association_inventory: pd.DataFrame
    candidate_feasibility: pd.DataFrame
    route_inventory: pd.DataFrame
    summary: pd.DataFrame


# ============================================================================
# Generic scalar helpers
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

    text = str(
        value
    ).strip()

    return text or None


def _as_float(
    value: Any,
) -> float | None:
    """Convert an optional scalar to float."""

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

    try:
        return float(
            value
        )
    except (
        TypeError,
        ValueError,
    ):
        return None


def _as_int(
    value: Any,
) -> int | None:
    """Convert an optional scalar to integer."""

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

    try:
        return int(
            value
        )
    except (
        TypeError,
        ValueError,
    ):
        return None


def _normalize_text_list(
    value: Any,
) -> list[str]:
    """Normalize a list-like field into unique strings."""

    if not isinstance(
        value,
        list,
    ):
        return []

    output: list[str] = []

    for item in value:

        normalized = _normalize_text(
            item
        )

        if normalized is not None:
            output.append(
                normalized
            )

    return list(
        dict.fromkeys(
            output
        )
    )


# ============================================================================
# GWAS Catalog REST API V2 schema compatibility
# ============================================================================


def _first_existing_key(
    mapping: dict[str, Any],
    *,
    keys: tuple[str, ...],
) -> tuple[str | None, Any]:
    """
    Return the first existing key/value pair.

    Used only for explicit known V2 schema aliases.
    """

    for key in keys:

        if key in mapping:
            return (
                key,
                mapping[
                    key
                ],
            )

    return (
        None,
        None,
    )


def extract_v2_page(
    payload: dict[str, Any],
) -> dict[str, int]:
    """
    Extract GWAS Catalog REST API V2 pagination metadata.

    Supported page representations:

        Reference/OpenAPI:
            {
                "size": 20,
                "total_elements": 10,
                "total_pages": 1,
                "number": 0
            }

        Live endpoint:
            {
                "size": 20,
                "totalElements": 10,
                "totalPages": 1,
                "number": 0
            }

    Returned representation is always canonical snake_case.
    """

    if not isinstance(
        payload,
        dict,
    ):
        raise RuntimeError(
            "GWAS Catalog V2 payload must be a mapping."
        )

    page = payload.get(
        "page"
    )

    if not isinstance(
        page,
        dict,
    ):
        raise RuntimeError(
            "GWAS Catalog V2 payload is missing a valid 'page' object."
        )

    # ----------------------------------------------------------------------
    # size
    # ----------------------------------------------------------------------

    if "size" not in page:
        raise RuntimeError(
            "GWAS Catalog V2 page object is missing 'size'."
        )

    # ----------------------------------------------------------------------
    # number
    # ----------------------------------------------------------------------

    if "number" not in page:
        raise RuntimeError(
            "GWAS Catalog V2 page object is missing 'number'."
        )

    # ----------------------------------------------------------------------
    # total elements
    #
    # API reference:
    #   total_elements
    #
    # Live endpoint:
    #   totalElements
    # ----------------------------------------------------------------------

    (
        total_elements_key,
        total_elements_raw,
    ) = _first_existing_key(
        page,
        keys=(
            "total_elements",
            "totalElements",
        ),
    )

    if total_elements_key is None:
        raise RuntimeError(
            "GWAS Catalog V2 page object is missing "
            "'total_elements' or 'totalElements'. "
            f"Available keys: {sorted(page.keys())}"
        )

    # ----------------------------------------------------------------------
    # total pages
    # ----------------------------------------------------------------------

    (
        total_pages_key,
        total_pages_raw,
    ) = _first_existing_key(
        page,
        keys=(
            "total_pages",
            "totalPages",
        ),
    )

    if total_pages_key is None:
        raise RuntimeError(
            "GWAS Catalog V2 page object is missing "
            "'total_pages' or 'totalPages'. "
            f"Available keys: {sorted(page.keys())}"
        )

    # ----------------------------------------------------------------------
    # Integer conversion
    # ----------------------------------------------------------------------

    try:

        result = {
            "size":
                int(
                    page[
                        "size"
                    ]
                ),

            "total_elements":
                int(
                    total_elements_raw
                ),

            "total_pages":
                int(
                    total_pages_raw
                ),

            "number":
                int(
                    page[
                        "number"
                    ]
                ),
        }

    except (
        TypeError,
        ValueError,
    ) as exc:

        raise RuntimeError(
            "GWAS Catalog V2 pagination fields could not be "
            "converted to integers."
        ) from exc

    # ----------------------------------------------------------------------
    # Basic sanity checks
    # ----------------------------------------------------------------------

    for field in (
        "size",
        "total_elements",
        "total_pages",
        "number",
    ):

        if result[
            field
        ] < 0:
            raise RuntimeError(
                "GWAS Catalog V2 pagination value cannot be negative: "
                f"{field}={result[field]}"
            )

    # ----------------------------------------------------------------------
    # Additional consistency
    # ----------------------------------------------------------------------

    if (
        result[
            "total_elements"
        ]
        ==
        0
        and
        result[
            "total_pages"
        ]
        not in (
            0,
            1,
        )
    ):
        raise RuntimeError(
            "GWAS Catalog V2 returned an inconsistent zero-result "
            "pagination state."
        )

    return result


def extract_v2_content(
    payload: dict[str, Any],
) -> list[dict[str, Any]]:
    """
    Extract association records from GWAS Catalog REST API V2.

    Supported response representations:

    1. Reference/OpenAPI style:

        {
            "content": [...]
        }

    2. HAL-style live response:

        {
            "_embedded": {
                "associations": [...]
            }
        }

    3. Valid zero-result response:

        {
            "_links": {...},
            "page": {
                "totalElements": 0,
                ...
            }
        }

       In that case no content collection is required and [] is returned.

    This function does not infer absence from an empty collection.
    """

    if not isinstance(
        payload,
        dict,
    ):
        raise RuntimeError(
            "GWAS Catalog V2 payload must be a mapping."
        )

    # ----------------------------------------------------------------------
    # OpenAPI/reference representation
    # ----------------------------------------------------------------------

    if "content" in payload:

        content = payload[
            "content"
        ]

        if not isinstance(
            content,
            list,
        ):
            raise RuntimeError(
                "GWAS Catalog V2 'content' field must be a list."
            )

        invalid_items = [
            index
            for index, item in enumerate(
                content
            )
            if not isinstance(
                item,
                dict,
            )
        ]

        if invalid_items:
            raise RuntimeError(
                "GWAS Catalog V2 'content' contains non-object "
                f"records at indices: {invalid_items[:10]}"
            )

        return content

    # ----------------------------------------------------------------------
    # Live HAL-style representation
    # ----------------------------------------------------------------------

    if "_embedded" in payload:

        embedded = payload[
            "_embedded"
        ]

        if not isinstance(
            embedded,
            dict,
        ):
            raise RuntimeError(
                "GWAS Catalog V2 '_embedded' field must be a mapping."
            )

        # Association endpoint is expected to expose an association
        # collection. The exact collection name is accepted from the two
        # common spellings below.
        (
            association_key,
            associations,
        ) = _first_existing_key(
            embedded,
            keys=(
                "associations",
                "association",
            ),
        )

        if association_key is None:

            page = extract_v2_page(
                payload
            )

            if page[
                "total_elements"
            ] == 0:
                return []

            raise RuntimeError(
                "GWAS Catalog V2 reports non-zero results but "
                "'_embedded' does not contain an association collection. "
                f"Available keys: {sorted(embedded.keys())}"
            )

        if isinstance(
            associations,
            list,
        ):

            invalid_items = [
                index
                for index, item in enumerate(
                    associations
                )
                if not isinstance(
                    item,
                    dict,
                )
            ]

            if invalid_items:
                raise RuntimeError(
                    "GWAS Catalog V2 '_embedded.associations' contains "
                    f"non-object records at indices: {invalid_items[:10]}"
                )

            return associations

        if isinstance(
            associations,
            dict,
        ):

            values = list(
                associations.values()
            )

            if not all(
                isinstance(
                    item,
                    dict,
                )
                for item in values
            ):
                raise RuntimeError(
                    "GWAS Catalog V2 '_embedded.associations' mapping "
                    "contains non-object values."
                )

            return values

        raise RuntimeError(
            "GWAS Catalog V2 '_embedded.associations' must be a "
            "list or mapping."
        )

    # ----------------------------------------------------------------------
    # No content collection present.
    #
    # This is valid when total_elements == 0.
    # ----------------------------------------------------------------------

    page = extract_v2_page(
        payload
    )

    if page[
        "total_elements"
    ] == 0:
        return []

    raise RuntimeError(
        "GWAS Catalog V2 reports "
        f"{page['total_elements']} association(s), but no supported "
        "association collection was found. Expected 'content' or "
        "'_embedded.associations'. "
        f"Top-level keys: {sorted(payload.keys())}"
    )


def extract_v2_links(
    payload: dict[str, Any],
) -> Any:
    """
    Extract response links without relying on one representation.

    Supported:
        links
        _links

    Links are currently retained for provenance only. Pagination in the
    M7.5 runner remains deterministic using page and size query parameters.
    """

    if not isinstance(
        payload,
        dict,
    ):
        raise RuntimeError(
            "GWAS Catalog V2 payload must be a mapping."
        )

    if "links" in payload:
        return payload[
            "links"
        ]

    if "_links" in payload:
        return payload[
            "_links"
        ]

    return None


# ============================================================================
# Association parsers
# ============================================================================


def extract_efo_traits(
    association: dict[str, Any],
) -> tuple[list[str], list[str]]:
    """Extract EFO IDs and EFO labels."""

    values = association.get(
        "efo_traits",
        [],
    )

    if not isinstance(
        values,
        list,
    ):
        return (
            [],
            [],
        )

    efo_ids: list[str] = []
    efo_labels: list[str] = []

    for value in values:

        if not isinstance(
            value,
            dict,
        ):
            continue

        efo_id = _normalize_text(
            value.get(
                "efo_id"
            )
        )

        efo_trait = _normalize_text(
            value.get(
                "efo_trait"
            )
        )

        if efo_id is not None:
            efo_ids.append(
                efo_id
            )

        if efo_trait is not None:
            efo_labels.append(
                efo_trait
            )

    return (
        list(
            dict.fromkeys(
                efo_ids
            )
        ),
        list(
            dict.fromkeys(
                efo_labels
            )
        ),
    )


def extract_bg_efo_traits(
    association: dict[str, Any],
) -> tuple[list[str], list[str]]:
    """Extract background EFO IDs and labels."""

    values = association.get(
        "bg_efo_traits",
        [],
    )

    if not isinstance(
        values,
        list,
    ):
        return (
            [],
            [],
        )

    efo_ids: list[str] = []
    efo_labels: list[str] = []

    for value in values:

        if not isinstance(
            value,
            dict,
        ):
            continue

        efo_id = _normalize_text(
            value.get(
                "efo_id"
            )
        )

        efo_trait = _normalize_text(
            value.get(
                "efo_trait"
            )
        )

        if efo_id is not None:
            efo_ids.append(
                efo_id
            )

        if efo_trait is not None:
            efo_labels.append(
                efo_trait
            )

    return (
        list(
            dict.fromkeys(
                efo_ids
            )
        ),
        list(
            dict.fromkeys(
                efo_labels
            )
        ),
    )


def extract_candidate_allele(
    association: dict[str, Any],
    *,
    rsid: str,
) -> tuple[bool, str | None]:
    """
    Confirm exact candidate rsID and obtain its effect allele.

    Exact confirmation uses:

        snp_allele[].rs_id
    """

    values = association.get(
        "snp_allele",
        [],
    )

    if not isinstance(
        values,
        list,
    ):
        return (
            False,
            None,
        )

    for value in values:

        if not isinstance(
            value,
            dict,
        ):
            continue

        observed_rsid = _normalize_text(
            value.get(
                "rs_id"
            )
        )

        if observed_rsid != rsid:
            continue

        return (
            True,
            _normalize_text(
                value.get(
                    "effect_allele"
                )
            ),
        )

    return (
        False,
        None,
    )


def trait_is_prostate_relevant(
    *,
    reported_traits: list[str],
    efo_traits: list[str],
    config: dict[str, Any],
) -> bool:
    """
    Perform conservative prostate-trait text screening.

    This is screening only and does not establish phenotype compatibility.
    """

    trait_values = [
        *reported_traits,
        *efo_traits,
    ]

    keywords = [
        str(
            value
        )
        for value in config[
            "trait_relevance"
        ][
            "accepted_keywords"
        ]
    ]

    case_sensitive = bool(
        config[
            "trait_relevance"
        ][
            "case_sensitive"
        ]
    )

    if not case_sensitive:

        trait_values = [
            value.lower()
            for value in trait_values
        ]

        keywords = [
            value.lower()
            for value in keywords
        ]

    return any(
        keyword in trait_value
        for trait_value in trait_values
        for keyword in keywords
    )


# ============================================================================
# Acquisition status
# ============================================================================


def build_acquisition_status(
    *,
    api_results: dict[str, dict[str, Any]],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one API acquisition-status row per candidate."""

    rows: list[dict[str, Any]] = []

    acquisition_statuses = config[
        "statuses"
    ][
        "acquisition"
    ]

    for candidate in config[
        "candidates"
    ]:

        rsid = str(
            candidate[
                "rsid"
            ]
        )

        result = api_results.get(
            rsid,
            {},
        )

        request_success = bool(
            result.get(
                "request_success",
                False,
            )
        )

        payload = result.get(
            "payload"
        )

        total_elements = 0
        parsed_schema = False

        if request_success:

            if not isinstance(
                payload,
                dict,
            ):
                raise RuntimeError(
                    "Successful GWAS Catalog V2 request does not "
                    f"contain a valid payload for {rsid}."
                )

            page = extract_v2_page(
                payload
            )

            content = extract_v2_content(
                payload
            )

            total_elements = int(
                page[
                    "total_elements"
                ]
            )

            parsed_schema = True

            if len(
                content
            ) != total_elements:

                raise RuntimeError(
                    "GWAS Catalog V2 canonical payload record count "
                    f"mismatch for {rsid}: "
                    f"content={len(content)}, "
                    f"total_elements={total_elements}."
                )

            if total_elements > 0:

                acquisition_status = acquisition_statuses[
                    "success_with_results"
                ]

            else:

                acquisition_status = acquisition_statuses[
                    "success_zero_results"
                ]

        else:

            acquisition_status = (
                result.get(
                    "acquisition_status"
                )
                or
                acquisition_statuses[
                    "service_unavailable"
                ]
            )

        rows.append(
            {
                "rsid":
                    rsid,

                "priority_rank":
                    int(
                        candidate[
                            "priority_rank"
                        ]
                    ),

                "request_success":
                    request_success,

                "http_status":
                    result.get(
                        "http_status"
                    ),

                "response_schema_parsed":
                    parsed_schema,

                "total_elements":
                    total_elements,

                "pages_retrieved":
                    int(
                        result.get(
                            "pages_retrieved",
                            0,
                        )
                    ),

                "records_retrieved":
                    int(
                        result.get(
                            "records_retrieved",
                            0,
                        )
                    ),

                "acquisition_status":
                    acquisition_status,

                "error_type":
                    result.get(
                        "error_type"
                    ),

                "error_message":
                    result.get(
                        "error_message"
                    ),

                "snapshot_reused":
                    bool(
                        result.get(
                            "snapshot_reused",
                            False,
                        )
                    ),

                "negative_evidence":
                    False,

                "variant_absence_inferred":
                    False,
            }
        )

    result = pd.DataFrame(
        rows
    )

    return result.sort_values(
        "priority_rank",
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Association normalization
# ============================================================================


def normalize_association(
    *,
    query_rsid: str,
    association: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    """Normalize one GWAS Catalog REST API V2 association."""

    (
        exact_candidate_rsid,
        candidate_effect_allele,
    ) = extract_candidate_allele(
        association,
        rsid=query_rsid,
    )

    (
        efo_ids,
        efo_traits,
    ) = extract_efo_traits(
        association
    )

    (
        bg_efo_ids,
        bg_efo_traits,
    ) = extract_bg_efo_traits(
        association
    )

    reported_traits = _normalize_text_list(
        association.get(
            "reported_trait"
        )
    )

    prostate_trait_candidate = trait_is_prostate_relevant(
        reported_traits=reported_traits,
        efo_traits=efo_traits,
        config=config,
    )

    accession_id = _normalize_text(
        association.get(
            "accession_id"
        )
    )

    previously_used_accessions = {
        str(
            value
        )
        for value in config[
            "previously_used_gwas"
        ][
            "study_accessions"
        ]
    }

    previously_used_study = bool(
        accession_id is not None
        and
        accession_id in previously_used_accessions
    )

    statuses = config[
        "statuses"
    ][
        "association"
    ]

    if not exact_candidate_rsid:

        association_status = statuses[
            "exact_rsid_not_confirmed"
        ]

    elif (
        prostate_trait_candidate
        and
        previously_used_study
    ):

        association_status = statuses[
            "prostate_reused"
        ]

    elif prostate_trait_candidate:

        association_status = statuses[
            "prostate_nonreused"
        ]

    else:

        association_status = statuses[
            "nonprostate"
        ]

    return {
        "query_rsid":
            query_rsid,

        "association_id":
            _as_int(
                association.get(
                    "association_id"
                )
            ),

        "exact_candidate_rsid":
            exact_candidate_rsid,

        "candidate_effect_allele":
            candidate_effect_allele,

        "accession_id":
            accession_id,

        "pubmed_id":
            _normalize_text(
                association.get(
                    "pubmed_id"
                )
            ),

        "first_author":
            _normalize_text(
                association.get(
                    "first_author"
                )
            ),

        "reported_traits":
            reported_traits,

        "efo_ids":
            efo_ids,

        "efo_traits":
            efo_traits,

        "bg_efo_ids":
            bg_efo_ids,

        "bg_efo_traits":
            bg_efo_traits,

        "prostate_trait_candidate":
            prostate_trait_candidate,

        "phenotype_compatibility_verified":
            False,

        "previously_used_study":
            previously_used_study,

        "study_independence_verified":
            False,

        "risk_frequency":
            _normalize_text(
                association.get(
                    "risk_frequency"
                )
            ),

        "pvalue_description":
            _normalize_text(
                association.get(
                    "pvalue_description"
                )
            ),

        "pvalue_mantissa":
            _as_float(
                association.get(
                    "pvalue_mantissa"
                )
            ),

        "pvalue_exponent":
            _as_int(
                association.get(
                    "pvalue_exponent"
                )
            ),

        "p_value":
            _as_float(
                association.get(
                    "p_value"
                )
            ),

        "or_per_copy_num":
            _as_float(
                association.get(
                    "or_per_copy_num"
                )
            ),

        "or_value":
            _as_float(
                association.get(
                    "or_value"
                )
            ),

        "beta_num":
            _as_float(
                association.get(
                    "beta_num"
                )
            ),

        "beta":
            _as_float(
                association.get(
                    "beta"
                )
            ),

        "beta_unit":
            _normalize_text(
                association.get(
                    "beta_unit"
                )
            ),

        "beta_direction":
            _normalize_text(
                association.get(
                    "beta_direction"
                )
            ),

        "ci_lower":
            _as_float(
                association.get(
                    "ci_lower"
                )
            ),

        "ci_upper":
            _as_float(
                association.get(
                    "ci_upper"
                )
            ),

        "multi_snp_haplotype":
            association.get(
                "multi_snp_haplotype"
            ),

        "snp_interaction":
            association.get(
                "snp_interaction"
            ),

        "snp_type":
            _normalize_text(
                association.get(
                    "snp_type"
                )
            ),

        "description":
            _normalize_text(
                association.get(
                    "description"
                )
            ),

        "range":
            _normalize_text_list(
                association.get(
                    "range"
                )
            ),

        "locations":
            _normalize_text_list(
                association.get(
                    "locations"
                )
            ),

        "mapped_genes":
            _normalize_text_list(
                association.get(
                    "mapped_genes"
                )
            ),

        "snp_effect_alleles":
            _normalize_text_list(
                association.get(
                    "snp_effect_allele"
                )
            ),

        "association_status":
            association_status,

        "negative_evidence":
            False,

        "independent_replication_claimed":
            False,
    }


def build_association_inventory(
    *,
    api_results: dict[str, dict[str, Any]],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build normalized curated-association inventory."""

    rows: list[dict[str, Any]] = []

    for candidate in config[
        "candidates"
    ]:

        rsid = str(
            candidate[
                "rsid"
            ]
        )

        result = api_results.get(
            rsid,
            {},
        )

        if not bool(
            result.get(
                "request_success",
                False,
            )
        ):
            continue

        payload = result.get(
            "payload"
        )

        if not isinstance(
            payload,
            dict,
        ):
            continue

        content = extract_v2_content(
            payload
        )

        for association in content:

            rows.append(
                normalize_association(
                    query_rsid=rsid,
                    association=association,
                    config=config,
                )
            )

    columns = [
        "query_rsid",
        "association_id",
        "exact_candidate_rsid",
        "candidate_effect_allele",
        "accession_id",
        "pubmed_id",
        "first_author",
        "reported_traits",
        "efo_ids",
        "efo_traits",
        "bg_efo_ids",
        "bg_efo_traits",
        "prostate_trait_candidate",
        "phenotype_compatibility_verified",
        "previously_used_study",
        "study_independence_verified",
        "risk_frequency",
        "pvalue_description",
        "pvalue_mantissa",
        "pvalue_exponent",
        "p_value",
        "or_per_copy_num",
        "or_value",
        "beta_num",
        "beta",
        "beta_unit",
        "beta_direction",
        "ci_lower",
        "ci_upper",
        "multi_snp_haplotype",
        "snp_interaction",
        "snp_type",
        "description",
        "range",
        "locations",
        "mapped_genes",
        "snp_effect_alleles",
        "association_status",
        "negative_evidence",
        "independent_replication_claimed",
    ]

    return pd.DataFrame(
        rows,
        columns=columns,
    )


# ============================================================================
# Candidate feasibility
# ============================================================================


def build_candidate_feasibility(
    *,
    acquisition_status: pd.DataFrame,
    association_inventory: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Resolve candidate-level curated-association feasibility."""

    statuses = config[
        "statuses"
    ][
        "candidate"
    ]

    rows: list[dict[str, Any]] = []

    for candidate in config[
        "candidates"
    ]:

        rsid = str(
            candidate[
                "rsid"
            ]
        )

        acquisition_subset = acquisition_status.loc[
            acquisition_status[
                "rsid"
            ]
            ==
            rsid
        ]

        if len(
            acquisition_subset
        ) != 1:

            raise RuntimeError(
                f"Expected one M7.5 acquisition row for {rsid}."
            )

        acquisition = acquisition_subset.iloc[
            0
        ]

        associations = association_inventory.loc[
            association_inventory[
                "query_rsid"
            ]
            ==
            rsid
        ].copy()

        exact = associations.loc[
            associations[
                "exact_candidate_rsid"
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
        ].copy()

        prostate = exact.loc[
            exact[
                "prostate_trait_candidate"
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
        ].copy()

        reused_studies = sorted(
            {
                str(
                    value
                )
                for value in prostate.loc[
                    prostate[
                        "previously_used_study"
                    ]
                    .fillna(
                        False
                    )
                    .astype(
                        bool
                    ),
                    "accession_id",
                ].dropna()
            }
        )

        nonreused_studies = sorted(
            {
                str(
                    value
                )
                for value in prostate.loc[
                    ~prostate[
                        "previously_used_study"
                    ]
                    .fillna(
                        False
                    )
                    .astype(
                        bool
                    ),
                    "accession_id",
                ].dropna()
            }
        )

        request_success = bool(
            acquisition[
                "request_success"
            ]
        )

        if not request_success:

            candidate_status = statuses[
                "api_unavailable"
            ]

        elif nonreused_studies:

            candidate_status = statuses[
                "prostate_nonreused"
            ]

        elif reused_studies:

            candidate_status = statuses[
                "prostate_reused"
            ]

        elif len(
            exact
        ) > 0:

            candidate_status = statuses[
                "nonprostate_only"
            ]

        else:

            candidate_status = statuses[
                "no_curated_association"
            ]

        rows.append(
            {
                "rsid":
                    rsid,

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

                "chromosome":
                    str(
                        candidate[
                            "chromosome"
                        ]
                    ),

                "position":
                    int(
                        candidate[
                            "position"
                        ]
                    ),

                "api_query_successful":
                    request_success,

                "curated_records_returned":
                    int(
                        acquisition[
                            "records_retrieved"
                        ]
                    ),

                "exact_candidate_associations":
                    int(
                        len(
                            exact
                        )
                    ),

                "prostate_trait_candidate_associations":
                    int(
                        len(
                            prostate
                        )
                    ),

                "reused_prostate_studies":
                    int(
                        len(
                            reused_studies
                        )
                    ),

                "nonreused_prostate_studies":
                    int(
                        len(
                            nonreused_studies
                        )
                    ),

                "reused_study_accessions":
                    reused_studies,

                "nonreused_study_accessions":
                    nonreused_studies,

                "study_independence_verified":
                    False,

                "independent_replication_established":
                    False,

                "phenotype_compatibility_verified":
                    False,

                "ancestry_compatibility_verified":
                    False,

                "effect_allele_harmonization_performed":
                    False,

                # All candidates continue to M7.5B. Curated Catalog presence
                # alone is not sufficient for full variant-level validation.
                "full_summary_statistics_audit_required":
                    True,

                "tcga_germline_route_available_now":
                    False,

                "tcga_germline_deferred_controlled_access":
                    True,

                "tcga_somatic_surrogate_allowed":
                    False,

                "api_zero_result_interpreted_as_negative":
                    False,

                "api_failure_interpreted_as_negative":
                    False,

                "variant_absence_inferred":
                    False,

                "candidate_status":
                    candidate_status,
            }
        )

    result = pd.DataFrame(
        rows
    )

    return result.sort_values(
        "priority_rank",
        kind="stable",
    ).reset_index(
        drop=True
    )


# ============================================================================
# Route inventory
# ============================================================================


def build_route_inventory(
    *,
    acquisition_status: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build variant-validation route inventory."""

    all_queries_successful = bool(
        acquisition_status[
            "request_success"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .all()
    )

    any_query_successful = bool(
        acquisition_status[
            "request_success"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .any()
    )

    if all_queries_successful:

        curated_route_status = (
            "CURATED_ASSOCIATION_ROUTE_AUDITED_ALL_CANDIDATES"
        )

    elif any_query_successful:

        curated_route_status = (
            "CURATED_ASSOCIATION_ROUTE_PARTIALLY_AUDITED"
        )

    else:

        curated_route_status = (
            "CURATED_ASSOCIATION_ROUTE_SERVICE_OR_SCHEMA_UNAVAILABLE"
        )

    return pd.DataFrame(
        [
            {
                "route":
                    "GWAS_CATALOG_REST_API_V2_CURATED_ASSOCIATIONS",

                "technically_relevant":
                    True,

                "access_type":
                    "PUBLIC",

                "route_available_now":
                    True,

                "analysis_performed":
                    any_query_successful,

                "all_candidates_audited":
                    all_queries_successful,

                "status":
                    curated_route_status,
            },
            {
                "route":
                    "FULL_GWAS_SUMMARY_STATISTICS",

                "technically_relevant":
                    True,

                "access_type":
                    "PUBLIC_WHEN_AVAILABLE",

                "route_available_now":
                    True,

                "analysis_performed":
                    False,

                "all_candidates_audited":
                    False,

                "status":
                    "M7.5B_FULL_SUMMARY_STATISTICS_AUDIT_REQUIRED",
            },
            {
                "route":
                    "TCGA_GERMLINE",

                "technically_relevant":
                    True,

                "access_type":
                    "CONTROLLED",

                "route_available_now":
                    False,

                "analysis_performed":
                    False,

                "all_candidates_audited":
                    False,

                "status":
                    config[
                        "routes"
                    ][
                        "tcga_germline"
                    ][
                        "status"
                    ],
            },
            {
                "route":
                    "TCGA_SOMATIC_AS_GERMLINE",

                "technically_relevant":
                    False,

                "access_type":
                    "NOT_APPLICABLE",

                "route_available_now":
                    False,

                "analysis_performed":
                    False,

                "all_candidates_audited":
                    False,

                "status":
                    config[
                        "routes"
                    ][
                        "tcga_somatic"
                    ][
                        "status"
                    ],
            },
        ]
    )


# ============================================================================
# Summary
# ============================================================================


def build_summary(
    *,
    acquisition_status: pd.DataFrame,
    candidate_feasibility: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Build one-row M7.5 summary."""

    candidates_assessed = int(
        len(
            candidate_feasibility
        )
    )

    successful_queries = int(
        acquisition_status[
            "request_success"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .sum()
    )

    zero_curated_queries = int(
        (
            acquisition_status[
                "request_success"
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
            &
            (
                acquisition_status[
                    "total_elements"
                ]
                ==
                0
            )
        ).sum()
    )

    queries_with_curated_results = int(
        (
            acquisition_status[
                "request_success"
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
            &
            (
                acquisition_status[
                    "total_elements"
                ]
                >
                0
            )
        ).sum()
    )

    nonreused_candidates = int(
        (
            candidate_feasibility[
                "nonreused_prostate_studies"
            ]
            >
            0
        ).sum()
    )

    reused_only_candidates = int(
        (
            (
                candidate_feasibility[
                    "reused_prostate_studies"
                ]
                >
                0
            )
            &
            (
                candidate_feasibility[
                    "nonreused_prostate_studies"
                ]
                ==
                0
            )
        ).sum()
    )

    statuses = config[
        "statuses"
    ][
        "overall"
    ]

    if successful_queries == 0:

        overall_status = statuses[
            "api_unavailable"
        ]

    elif nonreused_candidates > 0:

        overall_status = statuses[
            "candidate_nonreused_evidence"
        ]

    elif reused_only_candidates > 0:

        overall_status = statuses[
            "reused_evidence_only"
        ]

    else:

        overall_status = statuses[
            "no_curated_evidence"
        ]

    return pd.DataFrame(
        [
            {
                "milestone":
                    "M7.5",

                "candidates_assessed":
                    candidates_assessed,

                "api_queries_successful":
                    successful_queries,

                "api_queries_with_curated_results":
                    queries_with_curated_results,

                "api_queries_zero_curated_results":
                    zero_curated_queries,

                "all_candidate_api_queries_successful":
                    bool(
                        successful_queries
                        ==
                        candidates_assessed
                    ),

                "candidates_with_nonreused_prostate_curated_evidence":
                    nonreused_candidates,

                "candidates_with_reused_only_prostate_curated_evidence":
                    reused_only_candidates,

                "full_summary_statistics_audit_required_candidates":
                    candidates_assessed,

                "independent_replications_verified":
                    0,

                "independent_replication_claimed":
                    False,

                "study_independence_verified":
                    False,

                "phenotype_compatibility_verified":
                    False,

                "ancestry_compatibility_verified":
                    False,

                "effect_allele_harmonization_performed":
                    False,

                "effect_direction_concordance_claimed":
                    False,

                "full_summary_statistics_audit_performed":
                    False,

                "tcga_germline_route_available_now":
                    False,

                "tcga_germline_deferred_controlled_access":
                    True,

                "tcga_somatic_used_as_germline":
                    False,

                "api_zero_result_means_negative_evidence":
                    False,

                "api_failure_means_negative_evidence":
                    False,

                "api_zero_result_means_variant_absent":
                    False,

                "overall_status":
                    overall_status,

                "next_stage":
                    config[
                        "next_stage"
                    ],
            }
        ]
    )


# ============================================================================
# Scientific validation
# ============================================================================


def _validate_result(
    *,
    acquisition_status: pd.DataFrame,
    association_inventory: pd.DataFrame,
    candidate_feasibility: pd.DataFrame,
    summary: pd.DataFrame,
    config: dict[str, Any],
) -> None:
    """Validate structural and scientific safeguards."""

    expected_count = int(
        config[
            "validation"
        ][
            "expected_candidate_count"
        ]
    )

    if len(
        candidate_feasibility
    ) != expected_count:

        raise RuntimeError(
            "M7.5 candidate count mismatch."
        )

    expected_rsids = {
        str(
            value
        )
        for value in config[
            "validation"
        ][
            "expected_rsids"
        ]
    }

    observed_rsids = {
        str(
            value
        )
        for value in candidate_feasibility[
            "rsid"
        ]
    }

    if observed_rsids != expected_rsids:

        raise RuntimeError(
            "M7.5 candidate rsID set mismatch.\n"
            f"Observed: {sorted(observed_rsids)}\n"
            f"Expected: {sorted(expected_rsids)}"
        )

    if len(
        acquisition_status
    ) != expected_count:

        raise RuntimeError(
            "M7.5 acquisition-status row count mismatch."
        )

    if len(
        summary
    ) != 1:

        raise RuntimeError(
            "M7.5 summary must contain exactly one row."
        )

    summary_row = summary.iloc[
        0
    ]

    forbidden_true = [
        "independent_replication_claimed",
        "study_independence_verified",
        "phenotype_compatibility_verified",
        "ancestry_compatibility_verified",
        "effect_allele_harmonization_performed",
        "effect_direction_concordance_claimed",
        "full_summary_statistics_audit_performed",
        "tcga_germline_route_available_now",
        "tcga_somatic_used_as_germline",
        "api_zero_result_means_negative_evidence",
        "api_failure_means_negative_evidence",
        "api_zero_result_means_variant_absent",
    ]

    for field in forbidden_true:

        if bool(
            summary_row[
                field
            ]
        ):

            raise RuntimeError(
                f"M7.5 safeguard violation: {field}=True."
            )

    if int(
        summary_row[
            "independent_replications_verified"
        ]
    ) != 0:

        raise RuntimeError(
            "M7.5 cannot verify independent replication."
        )

    if (
        acquisition_status[
            "negative_evidence"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .any()
    ):

        raise RuntimeError(
            "M7.5 API acquisition was interpreted as negative evidence."
        )

    if (
        acquisition_status[
            "variant_absence_inferred"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .any()
    ):

        raise RuntimeError(
            "M7.5 inferred variant absence from API state."
        )

    if not association_inventory.empty:

        if (
            association_inventory[
                "study_independence_verified"
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
            .any()
        ):

            raise RuntimeError(
                "M7.5 unexpectedly verified study independence."
            )

        if (
            association_inventory[
                "independent_replication_claimed"
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
            .any()
        ):

            raise RuntimeError(
                "M7.5 unexpectedly claimed independent replication."
            )

        if (
            association_inventory[
                "negative_evidence"
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
            .any()
        ):

            raise RuntimeError(
                "M7.5 association inventory contains negative-evidence "
                "interpretation."
            )

    if (
        candidate_feasibility[
            "independent_replication_established"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .any()
    ):

        raise RuntimeError(
            "M7.5 unexpectedly established independent replication."
        )

    if (
        candidate_feasibility[
            "tcga_somatic_surrogate_allowed"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .any()
    ):

        raise RuntimeError(
            "M7.5 allowed TCGA somatic calls as germline surrogate."
        )

    if not (
        candidate_feasibility[
            "full_summary_statistics_audit_required"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
        .all()
    ):

        raise RuntimeError(
            "M7.5 requires all candidates to proceed to full "
            "summary-statistics audit."
        )


# ============================================================================
# Public API
# ============================================================================


def assess_variant_validation_feasibility(
    *,
    api_results: dict[str, dict[str, Any]],
    config: dict[str, Any],
) -> VariantValidationFeasibilityResult:
    """Execute complete M7.5 core analysis."""

    acquisition_status = build_acquisition_status(
        api_results=api_results,
        config=config,
    )

    association_inventory = build_association_inventory(
        api_results=api_results,
        config=config,
    )

    candidate_feasibility = build_candidate_feasibility(
        acquisition_status=acquisition_status,
        association_inventory=association_inventory,
        config=config,
    )

    route_inventory = build_route_inventory(
        acquisition_status=acquisition_status,
        config=config,
    )

    summary = build_summary(
        acquisition_status=acquisition_status,
        candidate_feasibility=candidate_feasibility,
        config=config,
    )

    _validate_result(
        acquisition_status=acquisition_status,
        association_inventory=association_inventory,
        candidate_feasibility=candidate_feasibility,
        summary=summary,
        config=config,
    )

    return VariantValidationFeasibilityResult(
        acquisition_status=acquisition_status,
        association_inventory=association_inventory,
        candidate_feasibility=candidate_feasibility,
        route_inventory=route_inventory,
        summary=summary,
    )