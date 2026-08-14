"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/harmonization/variants.py

Description:
    Variant identity harmonization for M3.3.

    This module constructs a canonical variant identity from:

        - canonical dbSNP rsID
        - build-aware genomic coordinate
        - source variant identifier

    Variant matching follows two independent evidence paths:

        1. rsID identity
        2. verified build-aware coordinate identity

    The module intentionally does not:

        - resolve historical rsID aliases
        - query dbSNP
        - infer missing rsIDs from coordinates
        - perform genome liftover
        - correct malformed source identifiers
        - perform cross-dataset joins

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
from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any

from pcatrfqtl.harmonization.coordinates import (
    CanonicalCoordinate,
    CoordinateHarmonizer,
)


class VariantIdentityStatus(str, Enum):
    """Overall usability status of a harmonized variant identity."""

    RESOLVED = "RESOLVED"
    PARTIAL = "PARTIAL"
    UNRESOLVED = "UNRESOLVED"


class VariantMatchMethod(str, Enum):
    """Evidence used to establish variant identity."""

    RSID = "RSID"
    COORDINATE = "COORDINATE"
    RSID_AND_COORDINATE = "RSID_AND_COORDINATE"
    NONE = "NONE"


@dataclass(frozen=True)
class HarmonizedVariant:
    """
    Canonical representation of one source variant.

    Attributes:
        dataset:
            Source dataset identifier.

        source_variant_id:
            Original standardized/source variant identifier.

        canonical_rsid:
            Canonical rsID when valid.

        rsid_usable:
            Whether the rsID can participate in rsID-based matching.

        coordinate:
            Canonical build-aware coordinate metadata.

        coordinate_usable:
            Whether a valid chromosome-position representation exists.

        coordinate_join_allowed:
            Whether the dataset build permits coordinate matching.

        identity_key:
            Preferred stable identity for the current record.

        identity_method:
            Evidence used to construct the preferred identity.

        status:
            Whether the record is fully or partially identifiable.

        note:
            Optional provenance information.
    """

    dataset: str
    source_variant_id: str | None

    canonical_rsid: str | None
    rsid_usable: bool

    coordinate: CanonicalCoordinate

    coordinate_usable: bool
    coordinate_join_allowed: bool

    identity_key: str | None
    identity_method: VariantMatchMethod

    status: VariantIdentityStatus

    note: str | None = None

    def to_dict(
        self,
    ) -> dict[str, Any]:
        """Return JSON-serializable variant metadata."""

        result = asdict(
            self
        )

        result[
            "identity_method"
        ] = self.identity_method.value

        result[
            "status"
        ] = self.status.value

        return result


class VariantHarmonizer:
    """Construct and compare cross-source variant identities."""

    RSID_PATTERN = re.compile(
        r"^rs[0-9]+$",
        flags=re.IGNORECASE,
    )

    @classmethod
    def normalize_rsid(
        cls,
        value: Any,
    ) -> str | None:
        """
        Normalize a single dbSNP rsID.

        Examples:
            RS123 -> rs123
            rs123 -> rs123

        Unsupported values remain unresolved.
        """

        if value is None:
            return None

        text = str(
            value
        ).strip()

        if not text:
            return None

        if (
            cls.RSID_PATTERN
            .fullmatch(
                text
            )
            is None
        ):
            return None

        return (
            "rs"
            + text[
                2:
            ]
        )

    @staticmethod
    def make_rsid_identity_key(
        rsid: str | None,
    ) -> str | None:
        """Construct canonical rsID identity key."""

        if rsid is None:
            return None

        return (
            f"rsid:{rsid}"
        )

    @staticmethod
    def make_coordinate_identity_key(
        coordinate: CanonicalCoordinate,
    ) -> str | None:
        """Construct canonical build-aware coordinate identity key."""

        if not (
            coordinate.coordinate_usable
            and coordinate.coordinate_join_allowed
            and coordinate.build_aware_key
        ):
            return None

        return (
            "coord:"
            f"{coordinate.build_aware_key}"
        )

    @classmethod
    def harmonize(
        cls,
        *,
        dataset: str,
        source_variant_id: Any = None,
        rsid: Any = None,
        chromosome: Any = None,
        position: Any = None,
        coordinate: Any = None,
    ) -> HarmonizedVariant:
        """
        Harmonize one source variant.

        The caller may provide either:

            rsid
            chromosome + position
            coordinate

        or any combination of these.

        The function never infers an rsID from genomic coordinates.
        """

        canonical_rsid = (
            cls.normalize_rsid(
                rsid
            )
        )

        rsid_usable = (
            canonical_rsid
            is not None
        )

        coordinate_model = (
            CoordinateHarmonizer
            .harmonize(
                dataset=dataset,
                chromosome=chromosome,
                position=position,
                coordinate=coordinate,
            )
        )

        coordinate_identity_key = (
            cls.make_coordinate_identity_key(
                coordinate_model
            )
        )

        rsid_identity_key = (
            cls.make_rsid_identity_key(
                canonical_rsid
            )
        )

        if (
            rsid_identity_key
            is not None
            and coordinate_identity_key
            is not None
        ):
            identity_key = (
                rsid_identity_key
            )

            identity_method = (
                VariantMatchMethod
                .RSID_AND_COORDINATE
            )

            status = (
                VariantIdentityStatus
                .RESOLVED
            )

        elif (
            rsid_identity_key
            is not None
        ):
            identity_key = (
                rsid_identity_key
            )

            identity_method = (
                VariantMatchMethod
                .RSID
            )

            status = (
                VariantIdentityStatus
                .RESOLVED
            )

        elif (
            coordinate_identity_key
            is not None
        ):
            identity_key = (
                coordinate_identity_key
            )

            identity_method = (
                VariantMatchMethod
                .COORDINATE
            )

            status = (
                VariantIdentityStatus
                .PARTIAL
            )

        else:
            identity_key = None

            identity_method = (
                VariantMatchMethod
                .NONE
            )

            status = (
                VariantIdentityStatus
                .UNRESOLVED
            )

        normalized_source_variant_id = (
            None
            if source_variant_id is None
            else str(
                source_variant_id
            ).strip()
        )

        return HarmonizedVariant(
            dataset=dataset,
            source_variant_id=(
                normalized_source_variant_id
            ),
            canonical_rsid=(
                canonical_rsid
            ),
            rsid_usable=(
                rsid_usable
            ),
            coordinate=(
                coordinate_model
            ),
            coordinate_usable=(
                coordinate_model
                .coordinate_usable
            ),
            coordinate_join_allowed=(
                coordinate_model
                .coordinate_join_allowed
            ),
            identity_key=(
                identity_key
            ),
            identity_method=(
                identity_method
            ),
            status=(
                status
            ),
        )

    @classmethod
    def rsids_equal(
        cls,
        left: HarmonizedVariant,
        right: HarmonizedVariant,
    ) -> bool:
        """Compare variants using canonical rsID identity only."""

        if not (
            left.rsid_usable
            and right.rsid_usable
        ):
            return False

        return (
            left.canonical_rsid
            == right.canonical_rsid
        )

    @classmethod
    def coordinates_equal(
        cls,
        left: HarmonizedVariant,
        right: HarmonizedVariant,
    ) -> bool:
        """Compare variants using verified coordinate identity only."""

        if not (
            left.coordinate_join_allowed
            and right.coordinate_join_allowed
        ):
            return False

        return (
            CoordinateHarmonizer
            .coordinates_equal(
                left_dataset=(
                    left.dataset
                ),
                left_chromosome=(
                    left.coordinate
                    .chromosome
                ),
                left_position=(
                    left.coordinate
                    .position
                ),
                right_dataset=(
                    right.dataset
                ),
                right_chromosome=(
                    right.coordinate
                    .chromosome
                ),
                right_position=(
                    right.coordinate
                    .position
                ),
            )
        )

    @classmethod
    def variants_equal(
        cls,
        left: HarmonizedVariant,
        right: HarmonizedVariant,
    ) -> bool:
        """
        Determine whether two variants represent the same locus.

        Matching priority:

            1. canonical rsID
            2. verified build-aware coordinate

        Coordinate fallback is allowed only when the genome-build
        registry permits direct coordinate comparison.
        """

        if cls.rsids_equal(
            left,
            right,
        ):
            return True

        if cls.coordinates_equal(
            left,
            right,
        ):
            return True

        return False

    @classmethod
    def match_method(
        cls,
        left: HarmonizedVariant,
        right: HarmonizedVariant,
    ) -> VariantMatchMethod:
        """Return the evidence supporting a cross-source match."""

        rsid_match = (
            cls.rsids_equal(
                left,
                right,
            )
        )

        coordinate_match = (
            cls.coordinates_equal(
                left,
                right,
            )
        )

        if (
            rsid_match
            and coordinate_match
        ):
            return (
                VariantMatchMethod
                .RSID_AND_COORDINATE
            )

        if rsid_match:
            return (
                VariantMatchMethod
                .RSID
            )

        if coordinate_match:
            return (
                VariantMatchMethod
                .COORDINATE
            )

        return (
            VariantMatchMethod
            .NONE
        )