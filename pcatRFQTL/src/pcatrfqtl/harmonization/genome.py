"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/harmonization/genome.py

Description:
    Genome-build verification registry for M3.1.

    This module records the reference assembly used by each standardized
    source dataset before any coordinate-level cross-source
    harmonization is attempted.

    Current verified sources:

        Cancer-tRFQTL
            hg19 / GRCh37

        Moradi QTL
            hg19 / GRCh37

    GWAS Catalog is retained as source-dependent at this stage and is
    therefore not assumed to be directly coordinate-compatible with
    hg19 datasets.

    This module intentionally does not perform liftover.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from typing import ClassVar

from pcatrfqtl.harmonization.models import (
    GenomeAssembly,
    GenomeBuildStatus,
    GenomeBuildVerification,
)


class GenomeBuildRegistry:
    """Registry of genome-build provenance for standardized datasets."""

    TARGET_BUILD: ClassVar[str] = "hg38"

    TARGET_ASSEMBLY: ClassVar[
        GenomeAssembly
    ] = GenomeAssembly.GRCH38

    _REGISTRY: ClassVar[
        dict[
            str,
            GenomeBuildVerification,
        ]
    ] = {
        "cancer_trfqtl": GenomeBuildVerification(
            dataset="cancer_trfqtl",
            genome_build="hg19",
            assembly=GenomeAssembly.GRCH37,
            status=GenomeBuildStatus.VERIFIED,
            coordinate_system="chromosome:position",
            evidence=(
                "Cancer-tRFQTL supplementary tables explicitly "
                "report SNP positions using hg19 coordinates."
            ),
            liftover_required=True,
            coordinate_join_allowed=True,
            note=(
                "Coordinate matching with other verified hg19 sources "
                "is permitted before any future conversion to hg38."
            ),
        ),

        "moradi_qtl": GenomeBuildVerification(
            dataset="moradi_qtl",
            genome_build="hg19",
            assembly=GenomeAssembly.GRCH37,
            status=GenomeBuildStatus.VERIFIED,
            coordinate_system="chromosome:position",
            evidence=(
                "The Moradi publication states that every sQTL and "
                "iso-eQTL table reports SNP genomic position in hg19."
            ),
            liftover_required=True,
            coordinate_join_allowed=True,
            note=(
                "Applies to the genomic positions represented in "
                "Moradi supplementary QTL tables S1-S4."
            ),
        ),

        "moradi_de": GenomeBuildVerification(
            dataset="moradi_de",
            genome_build="not_applicable",
            assembly=GenomeAssembly.UNKNOWN,
            status=GenomeBuildStatus.UNVERIFIED,
            coordinate_system="feature_identifier_only",
            evidence=(
                "Moradi differential-expression tables S5-S10 are "
                "feature-level outputs and do not provide genomic SNP "
                "coordinates used for cross-source variant matching."
            ),
            liftover_required=False,
            coordinate_join_allowed=False,
            note=(
                "Genome build is not required for direct DE feature "
                "statistics, but transcript identifiers may require "
                "separate feature harmonization."
            ),
        ),

        "gwas_catalog": GenomeBuildVerification(
            dataset="gwas_catalog",
            genome_build="source_dependent",
            assembly=GenomeAssembly.SOURCE_DEPENDENT,
            status=GenomeBuildStatus.SOURCE_DEPENDENT,
            coordinate_system="catalog_mapped_coordinates",
            evidence=(
                "GWAS Catalog coordinates are retained as "
                "source-dependent during M3.1 rather than assigning "
                "a genome build without record-level verification."
            ),
            liftover_required=False,
            coordinate_join_allowed=False,
            note=(
                "Coordinate-level matching against hg19 sources is "
                "disabled until the GWAS coordinate convention is "
                "explicitly resolved."
            ),
        ),
    }

    @classmethod
    def get(
        cls,
        dataset: str,
    ) -> GenomeBuildVerification:
        """Return genome-build metadata for one dataset."""

        key = str(
            dataset
        ).strip()

        if key not in cls._REGISTRY:
            raise KeyError(
                f"Genome-build metadata not registered for dataset: "
                f"{dataset}"
            )

        return cls._REGISTRY[
            key
        ]

    @classmethod
    def all(
        cls,
    ) -> dict[
        str,
        GenomeBuildVerification,
    ]:
        """Return a copy of all genome-build records."""

        return dict(
            cls._REGISTRY
        )

    @classmethod
    def are_coordinate_compatible(
        cls,
        left: str,
        right: str,
    ) -> bool:
        """
        Determine whether two datasets may currently be compared by
        genomic coordinate without liftover.

        Both datasets must:

            - have verified builds
            - explicitly allow coordinate joins
            - use the same genome assembly
        """

        left_build = cls.get(
            left
        )

        right_build = cls.get(
            right
        )

        if (
            left_build.status
            is not GenomeBuildStatus.VERIFIED
        ):
            return False

        if (
            right_build.status
            is not GenomeBuildStatus.VERIFIED
        ):
            return False

        if not (
            left_build.coordinate_join_allowed
            and right_build.coordinate_join_allowed
        ):
            return False

        return (
            left_build.assembly
            == right_build.assembly
        )

    @classmethod
    def requires_liftover_to_target(
        cls,
        dataset: str,
    ) -> bool:
        """
        Return whether the dataset requires conversion to the canonical
        harmonization target.
        """

        metadata = cls.get(
            dataset
        )

        if (
            metadata.status
            is not GenomeBuildStatus.VERIFIED
        ):
            return False

        return (
            metadata.assembly
            != cls.TARGET_ASSEMBLY
        )

    @classmethod
    def summary(
        cls,
    ) -> dict[
        str,
        dict,
    ]:
        """Return a JSON-serializable registry summary."""

        return {
            dataset: metadata.to_dict()
            for (
                dataset,
                metadata,
            ) in cls._REGISTRY.items()
        }