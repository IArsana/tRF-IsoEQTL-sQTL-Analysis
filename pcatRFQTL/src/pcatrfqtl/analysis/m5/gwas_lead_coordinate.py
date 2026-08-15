"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/analysis/m5/gwas_lead_coordinate.py

Description:
    M5.3C.3B canonical GRCh38 coordinate resolver for candidate lead
    variants.

    Resolution policy:
        1. Query Ensembl REST variation endpoint.
        2. Retry transient Ensembl failures.
        3. If Ensembl remains unavailable or does not provide a usable
           chromosome-level GRCh38 mapping, query the NCBI dbSNP RefSNP API.
        4. Select a canonical GRCh38 chromosome placement.
        5. Preserve provider and fallback provenance.

    Coordinate conventions:
        - Ensembl chromosome mapping positions are 1-based genomic positions.
        - NCBI SPDI positions are 0-based interbase coordinates.
        - NCBI fallback therefore converts:

              position_1based = spdi_position + 1

    Important safeguards:
        - Genome assembly is explicitly GRCh38.
        - Only canonical chromosome mappings are accepted.
        - Scaffold, patch, alternate-locus, and non-chromosomal placements
          are excluded.
        - Conflicting canonical mappings are rejected.
        - GRCh37 coordinates are never silently reused.
        - No liftover is performed by this module.
        - Coordinate presence does not imply disease association.
        - Physical proximity does not imply LD.
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

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Any


# ============================================================================
# Constants
# ============================================================================


TARGET_ASSEMBLY = "GRCh38"

TARGET_SPECIES = "human"


ENSEMBL_REST_BASE = (
    "https://rest.ensembl.org"
)


NCBI_REFSNP_BASE = (
    "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp"
)


ENSEMBL_MAX_ATTEMPTS = 3

ENSEMBL_RETRY_SECONDS = 2.0


RETRYABLE_HTTP_CODES = {
    429,
    500,
    502,
    503,
    504,
}


# Canonical RefSeq chromosome accessions for GRCh38.
#
# Version suffixes are deliberately retained because NCBI placements
# explicitly reference sequence versions.
GRCH38_REFSEQ_CHROMOSOMES = {
    "NC_000001.11": "1",
    "NC_000002.12": "2",
    "NC_000003.12": "3",
    "NC_000004.12": "4",
    "NC_000005.10": "5",
    "NC_000006.12": "6",
    "NC_000007.14": "7",
    "NC_000008.11": "8",
    "NC_000009.12": "9",
    "NC_000010.11": "10",
    "NC_000011.10": "11",
    "NC_000012.12": "12",
    "NC_000013.11": "13",
    "NC_000014.9": "14",
    "NC_000015.10": "15",
    "NC_000016.10": "16",
    "NC_000017.11": "17",
    "NC_000018.10": "18",
    "NC_000019.10": "19",
    "NC_000020.11": "20",
    "NC_000021.9": "21",
    "NC_000022.11": "22",
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "MT",
}


# ============================================================================
# Data model
# ============================================================================


@dataclass(frozen=True)
class CanonicalLeadCoordinate:
    """Canonical GRCh38 coordinate for one candidate lead."""

    lead_rsid: str

    found: bool

    chromosome: str | None

    base_pair_location: int | None

    genome_assembly: str | None

    lookup_method: str

    source_url: str | None

    reason: str | None

    provider: str | None = None

    fallback_used: bool = False

    attempt_count: int = 0


# ============================================================================
# Generic helpers
# ============================================================================


def _normalize_rsid(
    value: str,
) -> str:
    """Normalize and validate one rsID."""

    rsid = (
        str(
            value
        )
        .strip()
        .lower()
    )

    if not re.fullmatch(
        r"rs\d+",
        rsid,
    ):

        raise ValueError(
            f"Invalid rsID: {value}"
        )

    return rsid


def _normalize_chromosome(
    value: Any,
) -> str:
    """Normalize chromosome notation."""

    chromosome = (
        str(
            value
        )
        .strip()
    )

    chromosome = re.sub(
        r"^chr",
        "",
        chromosome,
        flags=re.IGNORECASE,
    )

    chromosome = chromosome.upper()

    if chromosome in {
        "M",
        "MT",
    }:
        return "MT"

    return chromosome


def _is_primary_chromosome(
    value: Any,
) -> bool:
    """Return True for canonical human chromosomes."""

    chromosome = _normalize_chromosome(
        value
    )

    return chromosome in {
        *{
            str(
                number
            )
            for number
            in range(
                1,
                23,
            )
        },
        "X",
        "Y",
        "MT",
    }


def _fetch_json(
    url: str,
    *,
    timeout: int = 60,
    user_agent: str,
) -> dict[str, Any]:
    """Fetch one JSON object."""

    request = urllib.request.Request(
        url,
        headers={
            "Accept":
                "application/json",

            "Content-Type":
                "application/json",

            "User-Agent":
                user_agent,
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=timeout,
    ) as response:

        payload = json.load(
            response
        )

    if not isinstance(
        payload,
        dict,
    ):

        raise RuntimeError(
            "Remote variant endpoint returned a non-object JSON payload."
        )

    return payload


# ============================================================================
# Ensembl provider
# ============================================================================


def _ensembl_url(
    rsid: str,
) -> str:
    """Construct Ensembl variation endpoint URL."""

    return (
        f"{ENSEMBL_REST_BASE}"
        f"/variation/{TARGET_SPECIES}/"
        f"{urllib.parse.quote(rsid)}"
    )


def _select_ensembl_grch38_mapping(
    payload: dict[str, Any],
) -> tuple[
    str,
    int,
] | None:
    """Select canonical chromosome-level GRCh38 Ensembl mapping."""

    mappings = payload.get(
        "mappings",
        []
    )

    if not isinstance(
        mappings,
        list,
    ):

        return None

    candidates: set[
        tuple[
            str,
            int,
        ]
    ] = set()

    for mapping in mappings:

        if not isinstance(
            mapping,
            dict,
        ):

            continue

        assembly = (
            str(
                mapping.get(
                    "assembly_name",
                    ""
                )
            )
            .strip()
        )

        if assembly != TARGET_ASSEMBLY:
            continue

        coord_system = (
            str(
                mapping.get(
                    "coord_system",
                    ""
                )
            )
            .strip()
            .lower()
        )

        if coord_system != "chromosome":
            continue

        seq_region = mapping.get(
            "seq_region_name"
        )

        if seq_region is None:
            continue

        if not _is_primary_chromosome(
            seq_region
        ):

            continue

        try:

            start = int(
                mapping.get(
                    "start"
                )
            )

            end = int(
                mapping.get(
                    "end"
                )
            )

        except (
            TypeError,
            ValueError,
        ):

            continue

        if (
            start <= 0
            or end <= 0
        ):

            continue

        if start != end:
            continue

        chromosome = _normalize_chromosome(
            seq_region
        )

        candidates.add(
            (
                chromosome,
                start,
            )
        )

    if not candidates:

        return None

    if len(
        candidates
    ) > 1:

        raise RuntimeError(
            "Ensembl returned conflicting canonical GRCh38 mappings: "
            f"{sorted(candidates)}"
        )

    return next(
        iter(
            candidates
        )
    )


def _resolve_from_ensembl(
    rsid: str,
) -> CanonicalLeadCoordinate:
    """Resolve one rsID through Ensembl with bounded retries."""

    url = _ensembl_url(
        rsid
    )

    last_reason: str | None = None

    attempts_used = 0

    for attempt in range(
        1,
        ENSEMBL_MAX_ATTEMPTS
        + 1,
    ):

        attempts_used = attempt

        try:

            payload = _fetch_json(
                url,
                timeout=60,
                user_agent=(
                    "pcatRFQTL/1.0 "
                    "M5.3C.3-Ensembl-coordinate-resolution"
                ),
            )

        except urllib.error.HTTPError as exc:

            last_reason = (
                f"HTTPError: HTTP Error "
                f"{exc.code}: {exc.reason}"
            )

            if (
                exc.code
                in RETRYABLE_HTTP_CODES
                and attempt
                < ENSEMBL_MAX_ATTEMPTS
            ):

                time.sleep(
                    ENSEMBL_RETRY_SECONDS
                    * attempt
                )

                continue

            break

        except (
            urllib.error.URLError,
            TimeoutError,
        ) as exc:

            last_reason = (
                f"{type(exc).__name__}: {exc}"
            )

            if (
                attempt
                < ENSEMBL_MAX_ATTEMPTS
            ):

                time.sleep(
                    ENSEMBL_RETRY_SECONDS
                    * attempt
                )

                continue

            break

        except Exception as exc:

            last_reason = (
                f"{type(exc).__name__}: {exc}"
            )

            break

        # ------------------------------------------------------------------
        # Identity safeguard
        # ------------------------------------------------------------------

        returned_name = payload.get(
            "name"
        )

        if returned_name is not None:

            normalized_name = (
                str(
                    returned_name
                )
                .strip()
                .lower()
            )

            if normalized_name != rsid:

                return CanonicalLeadCoordinate(
                    lead_rsid=rsid,
                    found=False,
                    chromosome=None,
                    base_pair_location=None,
                    genome_assembly=None,
                    lookup_method="ENSEMBL_VARIATION_REST",
                    source_url=url,
                    reason=(
                        "Ensembl response identity does not match "
                        f"requested rsID: {normalized_name}"
                    ),
                    provider="ENSEMBL",
                    fallback_used=False,
                    attempt_count=attempts_used,
                )

        # ------------------------------------------------------------------
        # Coordinate selection
        # ------------------------------------------------------------------

        try:

            mapping = (
                _select_ensembl_grch38_mapping(
                    payload
                )
            )

        except Exception as exc:

            return CanonicalLeadCoordinate(
                lead_rsid=rsid,
                found=False,
                chromosome=None,
                base_pair_location=None,
                genome_assembly=None,
                lookup_method="ENSEMBL_VARIATION_REST",
                source_url=url,
                reason=(
                    f"{type(exc).__name__}: {exc}"
                ),
                provider="ENSEMBL",
                fallback_used=False,
                attempt_count=attempts_used,
            )

        if mapping is None:

            return CanonicalLeadCoordinate(
                lead_rsid=rsid,
                found=False,
                chromosome=None,
                base_pair_location=None,
                genome_assembly=None,
                lookup_method="ENSEMBL_VARIATION_REST",
                source_url=url,
                reason=(
                    "No canonical chromosome-level GRCh38 mapping "
                    "was returned by Ensembl."
                ),
                provider="ENSEMBL",
                fallback_used=False,
                attempt_count=attempts_used,
            )

        chromosome, position = (
            mapping
        )

        return CanonicalLeadCoordinate(
            lead_rsid=rsid,
            found=True,
            chromosome=chromosome,
            base_pair_location=position,
            genome_assembly=TARGET_ASSEMBLY,
            lookup_method="ENSEMBL_VARIATION_REST",
            source_url=url,
            reason=None,
            provider="ENSEMBL",
            fallback_used=False,
            attempt_count=attempts_used,
        )

    return CanonicalLeadCoordinate(
        lead_rsid=rsid,
        found=False,
        chromosome=None,
        base_pair_location=None,
        genome_assembly=None,
        lookup_method="ENSEMBL_VARIATION_REST",
        source_url=url,
        reason=last_reason,
        provider="ENSEMBL",
        fallback_used=False,
        attempt_count=attempts_used,
    )


# ============================================================================
# NCBI RefSNP provider
# ============================================================================


def _ncbi_refsnp_url(
    rsid: str,
) -> str:
    """Construct NCBI RefSNP API URL."""

    numeric = rsid.removeprefix(
        "rs"
    )

    return (
        f"{NCBI_REFSNP_BASE}/"
        f"{numeric}"
    )


def _extract_ncbi_sequence_id(
    placement: dict[str, Any],
) -> str | None:
    """Extract RefSeq sequence accession from one placement."""

    seq_id = placement.get(
        "seq_id"
    )

    if seq_id is not None:

        return str(
            seq_id
        )

    placement_annot = placement.get(
        "placement_annot",
        {}
    )

    if isinstance(
        placement_annot,
        dict,
    ):

        seq_id = placement_annot.get(
            "seq_id"
        )

        if seq_id is not None:

            return str(
                seq_id
            )

    return None


def _placement_is_grch38(
    placement: dict[str, Any],
) -> bool:
    """
    Identify a GRCh38 assembly placement.

    NCBI placement metadata can vary across RefSNP records. We therefore
    inspect placement annotations rather than relying on one single key.
    """

    placement_annot = placement.get(
        "placement_annot",
        {}
    )

    if not isinstance(
        placement_annot,
        dict,
    ):

        return False

    traits = placement_annot.get(
        "seq_id_traits_by_assembly",
        []
    )

    if isinstance(
        traits,
        list,
    ):

        for trait in traits:

            if not isinstance(
                trait,
                dict,
            ):

                continue

            assembly_name = (
                str(
                    trait.get(
                        "assembly_name",
                        ""
                    )
                )
                .strip()
            )

            if assembly_name.startswith(
                TARGET_ASSEMBLY
            ):

                return True

    assembly_name = (
        str(
            placement_annot.get(
                "assembly_name",
                ""
            )
        )
        .strip()
    )

    return assembly_name.startswith(
        TARGET_ASSEMBLY
    )


def _placement_is_top_level(
    placement: dict[str, Any],
) -> bool:
    """Require a chromosome/top-level RefSNP placement."""

    placement_annot = placement.get(
        "placement_annot",
        {}
    )

    if not isinstance(
        placement_annot,
        dict,
    ):

        return False

    if placement_annot.get(
        "is_top_level"
    ) is True:

        return True

    # Some RefSNP responses may omit the flag but still use canonical
    # RefSeq chromosome accessions. That is accepted only when the sequence
    # accession itself is one of the known GRCh38 chromosomes.
    seq_id = _extract_ncbi_sequence_id(
        placement
    )

    return (
        seq_id
        in GRCH38_REFSEQ_CHROMOSOMES
    )


def _extract_spdi_positions(
    placement: dict[str, Any],
) -> set[int]:
    """Extract candidate 0-based SPDI positions from one placement."""

    alleles = placement.get(
        "alleles",
        []
    )

    if not isinstance(
        alleles,
        list,
    ):

        return set()

    positions: set[
        int
    ] = set()

    for allele in alleles:

        if not isinstance(
            allele,
            dict,
        ):

            continue

        allele_data = allele.get(
            "allele",
            {}
        )

        if not isinstance(
            allele_data,
            dict,
        ):

            continue

        spdi = allele_data.get(
            "spdi",
            {}
        )

        if not isinstance(
            spdi,
            dict,
        ):

            continue

        try:

            position = int(
                spdi.get(
                    "position"
                )
            )

        except (
            TypeError,
            ValueError,
        ):

            continue

        if position < 0:
            continue

        positions.add(
            position
        )

    return positions


def _select_ncbi_grch38_mapping(
    payload: dict[str, Any],
) -> tuple[
    str,
    int,
] | None:
    """
    Select canonical GRCh38 chromosome placement from an NCBI RefSNP record.

    NCBI SPDI position is 0-based interbase. Returned position is converted
    to a 1-based genomic coordinate.
    """

    primary_snapshot = payload.get(
        "primary_snapshot_data",
        {}
    )

    if not isinstance(
        primary_snapshot,
        dict,
    ):

        return None

    placements = primary_snapshot.get(
        "placements_with_allele",
        []
    )

    if not isinstance(
        placements,
        list,
    ):

        return None

    candidates: set[
        tuple[
            str,
            int,
        ]
    ] = set()

    for placement in placements:

        if not isinstance(
            placement,
            dict,
        ):

            continue

        if not _placement_is_grch38(
            placement
        ):

            continue

        if not _placement_is_top_level(
            placement
        ):

            continue

        seq_id = _extract_ncbi_sequence_id(
            placement
        )

        if (
            seq_id is None
            or seq_id
            not in GRCH38_REFSEQ_CHROMOSOMES
        ):

            continue

        chromosome = (
            GRCH38_REFSEQ_CHROMOSOMES[
                seq_id
            ]
        )

        positions_0based = (
            _extract_spdi_positions(
                placement
            )
        )

        if not positions_0based:
            continue

        # All alleles of a simple RefSNP placement should share the same
        # starting locus. Multiple conflicting start positions are rejected.
        if len(
            positions_0based
        ) != 1:

            continue

        position_0based = next(
            iter(
                positions_0based
            )
        )

        position_1based = (
            position_0based
            + 1
        )

        candidates.add(
            (
                chromosome,
                position_1based,
            )
        )

    if not candidates:

        return None

    if len(
        candidates
    ) > 1:

        raise RuntimeError(
            "NCBI RefSNP returned conflicting canonical GRCh38 placements: "
            f"{sorted(candidates)}"
        )

    return next(
        iter(
            candidates
        )
    )


def _resolve_from_ncbi(
    rsid: str,
    *,
    previous_attempt_count: int,
) -> CanonicalLeadCoordinate:
    """Resolve one rsID through NCBI RefSNP."""

    url = _ncbi_refsnp_url(
        rsid
    )

    try:

        payload = _fetch_json(
            url,
            timeout=60,
            user_agent=(
                "pcatRFQTL/1.0 "
                "M5.3C.3-NCBI-RefSNP-coordinate-resolution"
            ),
        )

    except Exception as exc:

        return CanonicalLeadCoordinate(
            lead_rsid=rsid,
            found=False,
            chromosome=None,
            base_pair_location=None,
            genome_assembly=None,
            lookup_method="NCBI_REFSNP_API",
            source_url=url,
            reason=(
                f"{type(exc).__name__}: {exc}"
            ),
            provider="NCBI_DBSNP",
            fallback_used=True,
            attempt_count=(
                previous_attempt_count
                + 1
            ),
        )

    # ----------------------------------------------------------------------
    # rsID safeguard
    # ----------------------------------------------------------------------

    refsnp_id = payload.get(
        "refsnp_id"
    )

    if refsnp_id is not None:

        returned_rsid = (
            "rs"
            + str(
                refsnp_id
            )
        ).lower()

        if returned_rsid != rsid:

            return CanonicalLeadCoordinate(
                lead_rsid=rsid,
                found=False,
                chromosome=None,
                base_pair_location=None,
                genome_assembly=None,
                lookup_method="NCBI_REFSNP_API",
                source_url=url,
                reason=(
                    "NCBI RefSNP response identity does not match "
                    f"requested rsID: returned={returned_rsid}"
                ),
                provider="NCBI_DBSNP",
                fallback_used=True,
                attempt_count=(
                    previous_attempt_count
                    + 1
                ),
            )

    try:

        mapping = (
            _select_ncbi_grch38_mapping(
                payload
            )
        )

    except Exception as exc:

        return CanonicalLeadCoordinate(
            lead_rsid=rsid,
            found=False,
            chromosome=None,
            base_pair_location=None,
            genome_assembly=None,
            lookup_method="NCBI_REFSNP_API",
            source_url=url,
            reason=(
                f"{type(exc).__name__}: {exc}"
            ),
            provider="NCBI_DBSNP",
            fallback_used=True,
            attempt_count=(
                previous_attempt_count
                + 1
            ),
        )

    if mapping is None:

        return CanonicalLeadCoordinate(
            lead_rsid=rsid,
            found=False,
            chromosome=None,
            base_pair_location=None,
            genome_assembly=None,
            lookup_method="NCBI_REFSNP_API",
            source_url=url,
            reason=(
                "No canonical GRCh38 chromosome placement "
                "was found in the NCBI RefSNP record."
            ),
            provider="NCBI_DBSNP",
            fallback_used=True,
            attempt_count=(
                previous_attempt_count
                + 1
            ),
        )

    chromosome, position = (
        mapping
    )

    return CanonicalLeadCoordinate(
        lead_rsid=rsid,
        found=True,
        chromosome=chromosome,
        base_pair_location=position,
        genome_assembly=TARGET_ASSEMBLY,
        lookup_method="NCBI_REFSNP_API",
        source_url=url,
        reason=None,
        provider="NCBI_DBSNP",
        fallback_used=True,
        attempt_count=(
            previous_attempt_count
            + 1
        ),
    )


# ============================================================================
# Public resolver
# ============================================================================


def resolve_lead_coordinate(
    lead_rsid: str,
) -> CanonicalLeadCoordinate:
    """
    Resolve one candidate lead rsID to a canonical GRCh38 coordinate.

    Provider order:
        1. Ensembl REST
        2. NCBI dbSNP RefSNP API fallback
    """

    rsid = _normalize_rsid(
        lead_rsid
    )

    # ----------------------------------------------------------------------
    # Primary provider: Ensembl
    # ----------------------------------------------------------------------

    ensembl_result = (
        _resolve_from_ensembl(
            rsid
        )
    )

    if ensembl_result.found:

        return ensembl_result

    # ----------------------------------------------------------------------
    # Fallback provider: NCBI RefSNP
    # ----------------------------------------------------------------------

    ncbi_result = (
        _resolve_from_ncbi(
            rsid,
            previous_attempt_count=(
                ensembl_result.attempt_count
            ),
        )
    )

    if ncbi_result.found:

        return ncbi_result

    # ----------------------------------------------------------------------
    # Neither provider resolved the locus.
    #
    # Preserve both failure reasons in the final result.
    # ----------------------------------------------------------------------

    reason = (
        "Coordinate resolution failed across all providers. "
        f"Ensembl: {ensembl_result.reason}; "
        f"NCBI: {ncbi_result.reason}"
    )

    return CanonicalLeadCoordinate(
        lead_rsid=rsid,
        found=False,
        chromosome=None,
        base_pair_location=None,
        genome_assembly=None,
        lookup_method="ENSEMBL_THEN_NCBI_REFSNP",
        source_url=ncbi_result.source_url,
        reason=reason,
        provider=None,
        fallback_used=True,
        attempt_count=(
            ncbi_result.attempt_count
        ),
    )