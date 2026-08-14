"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/harmonization/models.py

Description:
    Data models used by the cross-source harmonization stage.

    M3.1 introduces explicit genome-build metadata so genomic
    coordinates cannot be compared across datasets without first
    knowing their reference assembly.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any


class GenomeAssembly(str, Enum):
    """Canonical genome assemblies recognized by the pipeline."""

    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"
    SOURCE_DEPENDENT = "source_dependent"
    UNKNOWN = "unknown"


class GenomeBuildStatus(str, Enum):
    """Evidence status for a source genome build."""

    VERIFIED = "VERIFIED"
    SOURCE_DEPENDENT = "SOURCE_DEPENDENT"
    UNVERIFIED = "UNVERIFIED"


@dataclass(frozen=True)
class GenomeBuildVerification:
    """
    Genome-build provenance for one dataset.

    Attributes:
        dataset:
            Internal dataset identifier.

        genome_build:
            Common build label such as hg19 or hg38.

        assembly:
            Canonical reference assembly.

        status:
            Verification status.

        coordinate_system:
            Description of the coordinate convention represented by the
            source.

        evidence:
            Human-readable provenance supporting the assignment.

        liftover_required:
            Whether the dataset must eventually be converted before
            comparison against the selected harmonization target.

        coordinate_join_allowed:
            Whether coordinate-level cross-source matching is currently
            permitted.

        note:
            Optional additional source-specific information.
    """

    dataset: str
    genome_build: str
    assembly: GenomeAssembly
    status: GenomeBuildStatus
    coordinate_system: str
    evidence: str
    liftover_required: bool
    coordinate_join_allowed: bool
    note: str | None = None

    def to_dict(
        self,
    ) -> dict[str, Any]:
        """Return a JSON-serializable representation."""

        result = asdict(
            self
        )

        result[
            "assembly"
        ] = self.assembly.value

        result[
            "status"
        ] = self.status.value

        return result