"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/harmonization/coordinates.py

Description:
    Canonical genomic coordinate representation for M3.2.

    This module converts standardized chromosome-position information
    into a common coordinate representation while preserving explicit
    genome-build provenance.

    M3.2 responsibilities:

        - canonicalize chromosome labels
        - validate genomic positions
        - construct chromosome-position keys
        - preserve source genome-build metadata
        - determine whether coordinate-level joins are permitted
        - distinguish source coordinates from future target coordinates

    This module intentionally does not:

        - perform genome liftover
        - infer unknown genome builds
        - resolve rsIDs
        - perform cross-source joins
        - silently repair invalid coordinates

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
from typing import Any

from pcatrfqtl.harmonization.genome import (
    GenomeBuildRegistry,
)
from pcatrfqtl.harmonization.models import (
    GenomeBuildStatus,
)


@dataclass(frozen=True)
class CanonicalCoordinate:
    """
    Canonical representation of one genomic coordinate.

    Attributes:
        chromosome:
            Canonical chromosome label without ``chr`` prefix.

        position:
            Positive 1-based genomic position.

        coordinate_key:
            Build-independent display key such as ``6:28958399``.

        build_aware_key:
            Build-aware coordinate key such as
            ``hg19:6:28958399``.

        source_genome_build:
            Source build label.

        source_assembly:
            Source reference assembly.

        source_coordinate:
            Original coordinate representation when available.

        coordinate_usable:
            Whether chromosome and position are syntactically valid.

        coordinate_join_allowed:
            Whether this source may currently participate in
            coordinate-based cross-source joins.

        liftover_required:
            Whether conversion is required to reach the canonical
            harmonization target.

        liftover_performed:
            Always False during M3.2.

        target_genome_build:
            Canonical future harmonization target.

        target_coordinate_key:
            None during M3.2 because liftover is not yet performed.
    """

    chromosome: str | None
    position: int | None
    coordinate_key: str | None
    build_aware_key: str | None

    source_genome_build: str
    source_assembly: str

    source_coordinate: str | None

    coordinate_usable: bool
    coordinate_join_allowed: bool

    liftover_required: bool
    liftover_performed: bool

    target_genome_build: str
    target_coordinate_key: str | None

    def to_dict(
        self,
    ) -> dict[str, Any]:
        """Return JSON-serializable coordinate metadata."""

        return asdict(
            self
        )


class CoordinateHarmonizer:
    """
    Build-aware coordinate harmonization utilities.

    M3.2 operates only on standardized source coordinates.
    """

    VALID_CHROMOSOMES = {
        *{
            str(value)
            for value in range(
                1,
                23,
            )
        },
        "X",
        "Y",
        "MT",
    }

    COORDINATE_PATTERN = re.compile(
        r"^(?:chr)?"
        r"(?P<chromosome>"
        r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)"
        r")"
        r":"
        r"(?P<position>[0-9]+)$",
        flags=re.IGNORECASE,
    )

    @classmethod
    def normalize_chromosome(
        cls,
        value: Any,
    ) -> str | None:
        """
        Normalize chromosome labels.

        Examples:
            chr1 -> 1
            01   -> 1
            chrX -> X
            M    -> MT
            chrM -> MT
        """

        if value is None:
            return None

        text = str(
            value
        ).strip()

        if not text:
            return None

        if text.lower().startswith(
            "chr"
        ):
            text = text[
                3:
            ]

        text = text.upper()

        if text == "M":
            text = "MT"

        if text.isdigit():
            text = str(
                int(
                    text
                )
            )

        if (
            text
            not in cls.VALID_CHROMOSOMES
        ):
            return None

        return text

    @staticmethod
    def normalize_position(
        value: Any,
    ) -> int | None:
        """
        Normalize genomic position.

        Only positive integral coordinates are accepted.
        """

        if value is None:
            return None

        if isinstance(
            value,
            bool,
        ):
            return None

        try:
            numeric = int(
                value
            )
        except (
            TypeError,
            ValueError,
        ):
            return None

        if numeric <= 0:
            return None

        return numeric

    @classmethod
    def parse_coordinate(
        cls,
        value: Any,
    ) -> tuple[
        str | None,
        int | None,
    ]:
        """
        Parse chromosome-position string.

        Returns:
            (chromosome, position)

        Invalid coordinates return:
            (None, None)
        """

        if value is None:
            return (
                None,
                None,
            )

        text = str(
            value
        ).strip()

        match = (
            cls.COORDINATE_PATTERN
            .fullmatch(
                text
            )
        )

        if match is None:
            return (
                None,
                None,
            )

        chromosome = (
            cls.normalize_chromosome(
                match.group(
                    "chromosome"
                )
            )
        )

        position = (
            cls.normalize_position(
                match.group(
                    "position"
                )
            )
        )

        if (
            chromosome is None
            or position is None
        ):
            return (
                None,
                None,
            )

        return (
            chromosome,
            position,
        )

    @staticmethod
    def make_coordinate_key(
        chromosome: str | None,
        position: int | None,
    ) -> str | None:
        """
        Construct chromosome-position key.

        Example:
            6 + 28958399 -> 6:28958399
        """

        if (
            chromosome is None
            or position is None
        ):
            return None

        return (
            f"{chromosome}:{position}"
        )

    @staticmethod
    def make_build_aware_key(
        *,
        genome_build: str,
        chromosome: str | None,
        position: int | None,
    ) -> str | None:
        """
        Construct build-aware genomic key.

        Example:
            hg19 + 6 + 28958399
                -> hg19:6:28958399
        """

        if (
            chromosome is None
            or position is None
        ):
            return None

        if not genome_build:
            return None

        return (
            f"{genome_build}:"
            f"{chromosome}:"
            f"{position}"
        )

    @classmethod
    def harmonize(
        cls,
        *,
        dataset: str,
        chromosome: Any = None,
        position: Any = None,
        coordinate: Any = None,
    ) -> CanonicalCoordinate:
        """
        Harmonize one standardized source coordinate.

        Either:

            chromosome + position

        or:

            coordinate="chr:position"

        may be supplied.

        Dataset genome-build provenance is taken exclusively from the
        M3.1 GenomeBuildRegistry.
        """

        build_metadata = (
            GenomeBuildRegistry.get(
                dataset
            )
        )

        source_coordinate = (
            None
            if coordinate is None
            else str(
                coordinate
            ).strip()
        )

        parsed_chromosome: str | None
        parsed_position: int | None

        if coordinate is not None:
            (
                parsed_chromosome,
                parsed_position,
            ) = cls.parse_coordinate(
                coordinate
            )

        else:
            parsed_chromosome = (
                cls.normalize_chromosome(
                    chromosome
                )
            )

            parsed_position = (
                cls.normalize_position(
                    position
                )
            )

        coordinate_usable = (
            parsed_chromosome
            is not None
            and parsed_position
            is not None
        )

        coordinate_join_allowed = (
            coordinate_usable
            and build_metadata.status
            is GenomeBuildStatus.VERIFIED
            and build_metadata.coordinate_join_allowed
        )

        coordinate_key = (
            cls.make_coordinate_key(
                parsed_chromosome,
                parsed_position,
            )
            if coordinate_usable
            else None
        )

        build_aware_key = (
            cls.make_build_aware_key(
                genome_build=(
                    build_metadata
                    .genome_build
                ),
                chromosome=(
                    parsed_chromosome
                ),
                position=(
                    parsed_position
                ),
            )
            if coordinate_usable
            else None
        )

        return CanonicalCoordinate(
            chromosome=(
                parsed_chromosome
            ),
            position=(
                parsed_position
            ),
            coordinate_key=(
                coordinate_key
            ),
            build_aware_key=(
                build_aware_key
            ),
            source_genome_build=(
                build_metadata
                .genome_build
            ),
            source_assembly=(
                build_metadata
                .assembly
                .value
            ),
            source_coordinate=(
                source_coordinate
            ),
            coordinate_usable=(
                coordinate_usable
            ),
            coordinate_join_allowed=(
                coordinate_join_allowed
            ),
            liftover_required=(
                GenomeBuildRegistry
                .requires_liftover_to_target(
                    dataset
                )
            ),
            liftover_performed=False,
            target_genome_build=(
                GenomeBuildRegistry
                .TARGET_BUILD
            ),
            target_coordinate_key=None,
        )

    @classmethod
    def are_coordinates_matchable(
        cls,
        *,
        left_dataset: str,
        left_chromosome: Any,
        left_position: Any,
        right_dataset: str,
        right_chromosome: Any,
        right_position: Any,
    ) -> bool:
        """
        Determine whether two source coordinates may be directly
        compared.

        Matching is allowed only when:

            - both coordinates are valid
            - both datasets permit coordinate joins
            - both datasets use the same verified assembly
        """

        if not (
            GenomeBuildRegistry
            .are_coordinate_compatible(
                left_dataset,
                right_dataset,
            )
        ):
            return False

        left = cls.harmonize(
            dataset=left_dataset,
            chromosome=left_chromosome,
            position=left_position,
        )

        right = cls.harmonize(
            dataset=right_dataset,
            chromosome=right_chromosome,
            position=right_position,
        )

        return (
            left.coordinate_join_allowed
            and right.coordinate_join_allowed
        )

    @classmethod
    def coordinates_equal(
        cls,
        *,
        left_dataset: str,
        left_chromosome: Any,
        left_position: Any,
        right_dataset: str,
        right_chromosome: Any,
        right_position: Any,
    ) -> bool:
        """
        Compare two genomic coordinates only when their source builds
        are directly compatible.
        """

        if not cls.are_coordinates_matchable(
            left_dataset=left_dataset,
            left_chromosome=left_chromosome,
            left_position=left_position,
            right_dataset=right_dataset,
            right_chromosome=right_chromosome,
            right_position=right_position,
        ):
            return False

        left = cls.harmonize(
            dataset=left_dataset,
            chromosome=left_chromosome,
            position=left_position,
        )

        right = cls.harmonize(
            dataset=right_dataset,
            chromosome=right_chromosome,
            position=right_position,
        )

        return (
            left.build_aware_key
            == right.build_aware_key
        )