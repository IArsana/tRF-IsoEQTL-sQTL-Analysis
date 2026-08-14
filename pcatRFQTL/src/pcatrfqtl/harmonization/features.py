"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/harmonization/features.py

Description:
    General biological feature identity harmonization for M3.4.

    This module provides canonical representations for feature
    identifiers occurring across the pcatRFQTL datasets.

    Supported feature classes include:

        - transcript isoforms
        - exon splicing events
        - intron-retention events
        - tRNA-derived fragments
        - genes
        - unknown/source-specific features

    Feature harmonization preserves source semantics. Identifiers are
    never converted into a different biological identifier class
    without explicit evidence.

    Examples:

        ENST00000264639.9
            -> transcript isoform

        EX12345
            -> exon event

        INT12345
            -> intron event

        INT1e+05
            -> preserved source identifier, but non-canonical

        source-specific tRF identifier
            -> preserved as tRF identity

    This module intentionally does not:

        - query external annotation databases
        - convert transcript IDs to gene IDs
        - infer exon genomic coordinates
        - infer intron genomic coordinates
        - repair malformed source identifiers
        - remove Ensembl transcript versions
        - map gene symbols to Ensembl gene IDs
        - normalize tRF nomenclature across external databases

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


class FeatureType(str, Enum):
    """Biological feature classes recognized by the pipeline."""

    TRANSCRIPT_ISOFORM = "TRANSCRIPT_ISOFORM"
    EXON_EVENT = "EXON_EVENT"
    INTRON_EVENT = "INTRON_EVENT"
    TRF = "TRF"
    GENE = "GENE"
    UNKNOWN = "UNKNOWN"


class FeatureIdentityStatus(str, Enum):
    """Canonicalization status of a feature identifier."""

    RESOLVED = "RESOLVED"
    PRESERVED = "PRESERVED"
    PARTIAL = "PARTIAL"
    UNRESOLVED = "UNRESOLVED"


class FeatureIdentifierSystem(str, Enum):
    """Identifier systems recognized during M3.4."""

    ENSEMBL_TRANSCRIPT = "ENSEMBL_TRANSCRIPT"
    MORADI_EXON_EVENT = "MORADI_EXON_EVENT"
    MORADI_INTRON_EVENT = "MORADI_INTRON_EVENT"
    SOURCE_TRF = "SOURCE_TRF"
    GENE_SYMBOL = "GENE_SYMBOL"
    ENSEMBL_GENE = "ENSEMBL_GENE"
    SOURCE_SPECIFIC = "SOURCE_SPECIFIC"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class HarmonizedFeature:
    """
    Canonical representation of one biological feature.

    Attributes:
        source_feature_id:
            Identifier exactly as represented by the standardized
            source dataset.

        canonical_feature_id:
            Canonical identifier usable for feature-level matching.

        feature_type:
            Biological feature class.

        identifier_system:
            Identifier namespace or source system.

        identity_key:
            Namespace-aware identity key.

        usable:
            Whether the feature may participate in canonical matching.

        status:
            Feature identity resolution status.

        version:
            Optional version component, primarily for Ensembl IDs.

        base_identifier:
            Identifier without version when applicable.

        source_anomaly:
            Whether the identifier is retained as a known source
            anomaly.

        note:
            Optional provenance note.
    """

    source_feature_id: str | None
    canonical_feature_id: str | None

    feature_type: FeatureType
    identifier_system: FeatureIdentifierSystem

    identity_key: str | None

    usable: bool
    status: FeatureIdentityStatus

    version: str | None = None
    base_identifier: str | None = None

    source_anomaly: bool = False
    note: str | None = None

    def to_dict(
        self,
    ) -> dict[str, Any]:
        """Return a JSON-serializable representation."""

        result = asdict(
            self
        )

        result["feature_type"] = (
            self.feature_type.value
        )

        result["identifier_system"] = (
            self.identifier_system.value
        )

        result["status"] = (
            self.status.value
        )

        return result


class FeatureHarmonizer:
    """Harmonize biological feature identifiers across source datasets."""

    ENSEMBL_TRANSCRIPT_PATTERN = re.compile(
        r"^(ENST[0-9]+)(?:\.([0-9]+))?$",
        flags=re.IGNORECASE,
    )

    ENSEMBL_GENE_PATTERN = re.compile(
        r"^(ENSG[0-9]+)(?:\.([0-9]+))?$",
        flags=re.IGNORECASE,
    )

    EXON_EVENT_PATTERN = re.compile(
        r"^EX[0-9]+$",
        flags=re.IGNORECASE,
    )

    INTRON_EVENT_PATTERN = re.compile(
        r"^INT[0-9]+$",
        flags=re.IGNORECASE,
    )

    SCIENTIFIC_INTRON_PATTERN = re.compile(
        r"^INT[0-9]+(?:\.[0-9]+)?[eE][+-]?[0-9]+$",
        flags=re.IGNORECASE,
    )

    @staticmethod
    def normalize_source_identifier(
        value: Any,
    ) -> str | None:
        """Normalize representation without changing identifier semantics."""

        if value is None:
            return None

        text = str(
            value
        ).strip()

        if not text:
            return None

        return text

    @classmethod
    def classify(
        cls,
        value: Any,
        *,
        feature_type_hint: FeatureType | None = None,
    ) -> FeatureType:
        """
        Classify a source identifier.

        Explicit feature type hints take precedence for identifiers whose
        nomenclature cannot safely be inferred, particularly tRFs and
        gene symbols.
        """

        text = cls.normalize_source_identifier(
            value
        )

        if text is None:
            return FeatureType.UNKNOWN

        if (
            cls.ENSEMBL_TRANSCRIPT_PATTERN
            .fullmatch(
                text
            )
        ):
            return (
                FeatureType
                .TRANSCRIPT_ISOFORM
            )

        if (
            cls.EXON_EVENT_PATTERN
            .fullmatch(
                text
            )
        ):
            return (
                FeatureType
                .EXON_EVENT
            )

        if (
            cls.INTRON_EVENT_PATTERN
            .fullmatch(
                text
            )
        ):
            return (
                FeatureType
                .INTRON_EVENT
            )

        if (
            cls.SCIENTIFIC_INTRON_PATTERN
            .fullmatch(
                text
            )
        ):
            return (
                FeatureType
                .INTRON_EVENT
            )

        if (
            cls.ENSEMBL_GENE_PATTERN
            .fullmatch(
                text
            )
        ):
            return (
                FeatureType
                .GENE
            )

        if feature_type_hint is not None:
            return feature_type_hint

        return FeatureType.UNKNOWN

    @classmethod
    def harmonize_transcript(
        cls,
        value: Any,
    ) -> HarmonizedFeature:
        """Harmonize an Ensembl transcript identifier."""

        source = cls.normalize_source_identifier(
            value
        )

        if source is None:
            return cls._unresolved(
                FeatureType.TRANSCRIPT_ISOFORM
            )

        match = (
            cls.ENSEMBL_TRANSCRIPT_PATTERN
            .fullmatch(
                source
            )
        )

        if match is None:
            return HarmonizedFeature(
                source_feature_id=source,
                canonical_feature_id=None,
                feature_type=(
                    FeatureType
                    .TRANSCRIPT_ISOFORM
                ),
                identifier_system=(
                    FeatureIdentifierSystem
                    .SOURCE_SPECIFIC
                ),
                identity_key=None,
                usable=False,
                status=(
                    FeatureIdentityStatus
                    .PARTIAL
                ),
                note=(
                    "Value was expected to represent an Ensembl "
                    "transcript but did not match canonical ENST syntax."
                ),
            )

        base_identifier = (
            match.group(
                1
            ).upper()
        )

        version = (
            match.group(
                2
            )
        )

        canonical = (
            base_identifier
            if version is None
            else (
                f"{base_identifier}."
                f"{version}"
            )
        )

        return HarmonizedFeature(
            source_feature_id=source,
            canonical_feature_id=canonical,
            feature_type=(
                FeatureType
                .TRANSCRIPT_ISOFORM
            ),
            identifier_system=(
                FeatureIdentifierSystem
                .ENSEMBL_TRANSCRIPT
            ),
            identity_key=(
                f"transcript:{canonical}"
            ),
            usable=True,
            status=(
                FeatureIdentityStatus
                .RESOLVED
            ),
            version=version,
            base_identifier=(
                base_identifier
            ),
        )

    @classmethod
    def harmonize_exon_event(
        cls,
        value: Any,
    ) -> HarmonizedFeature:
        """Harmonize Moradi exon-event identifiers."""

        source = cls.normalize_source_identifier(
            value
        )

        if source is None:
            return cls._unresolved(
                FeatureType.EXON_EVENT
            )

        if (
            cls.EXON_EVENT_PATTERN
            .fullmatch(
                source
            )
            is None
        ):
            return HarmonizedFeature(
                source_feature_id=source,
                canonical_feature_id=None,
                feature_type=(
                    FeatureType.EXON_EVENT
                ),
                identifier_system=(
                    FeatureIdentifierSystem
                    .SOURCE_SPECIFIC
                ),
                identity_key=None,
                usable=False,
                status=(
                    FeatureIdentityStatus
                    .PARTIAL
                ),
            )

        canonical = (
            source.upper()
        )

        return HarmonizedFeature(
            source_feature_id=source,
            canonical_feature_id=canonical,
            feature_type=(
                FeatureType.EXON_EVENT
            ),
            identifier_system=(
                FeatureIdentifierSystem
                .MORADI_EXON_EVENT
            ),
            identity_key=(
                f"exon:{canonical}"
            ),
            usable=True,
            status=(
                FeatureIdentityStatus
                .RESOLVED
            ),
            base_identifier=canonical,
        )

    @classmethod
    def harmonize_intron_event(
        cls,
        value: Any,
    ) -> HarmonizedFeature:
        """
        Harmonize Moradi intron-retention event identifiers.

        Known scientific-notation-like source anomalies remain preserved
        but are not promoted to canonical identifiers.
        """

        source = cls.normalize_source_identifier(
            value
        )

        if source is None:
            return cls._unresolved(
                FeatureType.INTRON_EVENT
            )

        if (
            cls.INTRON_EVENT_PATTERN
            .fullmatch(
                source
            )
        ):
            canonical = (
                source.upper()
            )

            return HarmonizedFeature(
                source_feature_id=source,
                canonical_feature_id=canonical,
                feature_type=(
                    FeatureType.INTRON_EVENT
                ),
                identifier_system=(
                    FeatureIdentifierSystem
                    .MORADI_INTRON_EVENT
                ),
                identity_key=(
                    f"intron:{canonical}"
                ),
                usable=True,
                status=(
                    FeatureIdentityStatus
                    .RESOLVED
                ),
                base_identifier=canonical,
            )

        if (
            cls.SCIENTIFIC_INTRON_PATTERN
            .fullmatch(
                source
            )
        ):
            return HarmonizedFeature(
                source_feature_id=source,
                canonical_feature_id=None,
                feature_type=(
                    FeatureType.INTRON_EVENT
                ),
                identifier_system=(
                    FeatureIdentifierSystem
                    .SOURCE_SPECIFIC
                ),
                identity_key=None,
                usable=False,
                status=(
                    FeatureIdentityStatus
                    .PARTIAL
                ),
                source_anomaly=True,
                note=(
                    "Scientific-notation-like intron identifier retained "
                    "without automatic correction."
                ),
            )

        return HarmonizedFeature(
            source_feature_id=source,
            canonical_feature_id=None,
            feature_type=(
                FeatureType.INTRON_EVENT
            ),
            identifier_system=(
                FeatureIdentifierSystem
                .SOURCE_SPECIFIC
            ),
            identity_key=None,
            usable=False,
            status=(
                FeatureIdentityStatus
                .PARTIAL
            ),
        )

    @classmethod
    def harmonize_trf(
        cls,
        value: Any,
    ) -> HarmonizedFeature:
        """
        Preserve a tRF source identifier.

        M3.4 does not assume that Cancer-tRFQTL identifiers are directly
        equivalent to external tRF nomenclature systems.
        """

        source = cls.normalize_source_identifier(
            value
        )

        if source is None:
            return cls._unresolved(
                FeatureType.TRF
            )

        return HarmonizedFeature(
            source_feature_id=source,
            canonical_feature_id=source,
            feature_type=(
                FeatureType.TRF
            ),
            identifier_system=(
                FeatureIdentifierSystem
                .SOURCE_TRF
            ),
            identity_key=(
                f"trf:source:{source}"
            ),
            usable=True,
            status=(
                FeatureIdentityStatus
                .PRESERVED
            ),
            base_identifier=source,
            note=(
                "tRF identifier preserved in source namespace; no "
                "cross-database nomenclature conversion performed."
            ),
        )

    @classmethod
    def harmonize_gene(
        cls,
        value: Any,
    ) -> HarmonizedFeature:
        """
        Harmonize gene identifiers conservatively.

        Ensembl gene identifiers can be canonicalized directly.
        Other non-empty values are preserved as gene symbols without
        asserting a symbol-to-Ensembl mapping.
        """

        source = cls.normalize_source_identifier(
            value
        )

        if source is None:
            return cls._unresolved(
                FeatureType.GENE
            )

        ensembl_match = (
            cls.ENSEMBL_GENE_PATTERN
            .fullmatch(
                source
            )
        )

        if ensembl_match is not None:
            base_identifier = (
                ensembl_match.group(
                    1
                ).upper()
            )

            version = (
                ensembl_match.group(
                    2
                )
            )

            canonical = (
                base_identifier
                if version is None
                else (
                    f"{base_identifier}."
                    f"{version}"
                )
            )

            return HarmonizedFeature(
                source_feature_id=source,
                canonical_feature_id=canonical,
                feature_type=(
                    FeatureType.GENE
                ),
                identifier_system=(
                    FeatureIdentifierSystem
                    .ENSEMBL_GENE
                ),
                identity_key=(
                    f"gene:ensembl:{canonical}"
                ),
                usable=True,
                status=(
                    FeatureIdentityStatus
                    .RESOLVED
                ),
                version=version,
                base_identifier=(
                    base_identifier
                ),
            )

        return HarmonizedFeature(
            source_feature_id=source,
            canonical_feature_id=source,
            feature_type=(
                FeatureType.GENE
            ),
            identifier_system=(
                FeatureIdentifierSystem
                .GENE_SYMBOL
            ),
            identity_key=(
                f"gene:symbol:{source}"
            ),
            usable=True,
            status=(
                FeatureIdentityStatus
                .PRESERVED
            ),
            base_identifier=source,
            note=(
                "Gene symbol preserved without external "
                "symbol-to-Ensembl resolution."
            ),
        )

    @classmethod
    def harmonize(
        cls,
        value: Any,
        *,
        feature_type: FeatureType | str | None = None,
    ) -> HarmonizedFeature:
        """
        Harmonize any supported biological feature.

        Explicit feature type should normally be supplied from dataset
        metadata. Automatic classification is used only when the
        identifier syntax provides reliable evidence.
        """

        if isinstance(
            feature_type,
            str,
        ):
            feature_type = (
                FeatureType(
                    feature_type
                )
            )

        resolved_type = (
            feature_type
            if feature_type is not None
            else cls.classify(
                value
            )
        )

        if (
            resolved_type
            is FeatureType.TRANSCRIPT_ISOFORM
        ):
            return (
                cls.harmonize_transcript(
                    value
                )
            )

        if (
            resolved_type
            is FeatureType.EXON_EVENT
        ):
            return (
                cls.harmonize_exon_event(
                    value
                )
            )

        if (
            resolved_type
            is FeatureType.INTRON_EVENT
        ):
            return (
                cls.harmonize_intron_event(
                    value
                )
            )

        if (
            resolved_type
            is FeatureType.TRF
        ):
            return (
                cls.harmonize_trf(
                    value
                )
            )

        if (
            resolved_type
            is FeatureType.GENE
        ):
            return (
                cls.harmonize_gene(
                    value
                )
            )

        source = (
            cls.normalize_source_identifier(
                value
            )
        )

        if source is None:
            return cls._unresolved(
                FeatureType.UNKNOWN
            )

        return HarmonizedFeature(
            source_feature_id=source,
            canonical_feature_id=source,
            feature_type=(
                FeatureType.UNKNOWN
            ),
            identifier_system=(
                FeatureIdentifierSystem
                .SOURCE_SPECIFIC
            ),
            identity_key=(
                f"source:{source}"
            ),
            usable=False,
            status=(
                FeatureIdentityStatus
                .PRESERVED
            ),
            base_identifier=source,
            note=(
                "Identifier preserved because its biological feature "
                "class could not be determined safely."
            ),
        )

    @staticmethod
    def _unresolved(
        feature_type: FeatureType,
    ) -> HarmonizedFeature:
        """Construct an unresolved feature result."""

        return HarmonizedFeature(
            source_feature_id=None,
            canonical_feature_id=None,
            feature_type=feature_type,
            identifier_system=(
                FeatureIdentifierSystem
                .UNKNOWN
            ),
            identity_key=None,
            usable=False,
            status=(
                FeatureIdentityStatus
                .UNRESOLVED
            ),
        )

    @staticmethod
    def features_equal(
        left: HarmonizedFeature,
        right: HarmonizedFeature,
    ) -> bool:
        """
        Compare two canonical feature identities.

        Feature classes and namespaces must agree. A transcript cannot
        match an exon event merely because source strings happen to
        resemble each other.
        """

        if not (
            left.usable
            and right.usable
        ):
            return False

        if (
            left.feature_type
            is not right.feature_type
        ):
            return False

        if (
            left.identifier_system
            is not right.identifier_system
        ):
            return False

        return (
            left.identity_key
            == right.identity_key
        )