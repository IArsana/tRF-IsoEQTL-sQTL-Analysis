"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/harmonization/disease.py

Description:
    Disease-context harmonization for M3.5.

    This module provides conservative cross-source harmonization of
    disease and cancer-context labels occurring in the pcatRFQTL
    datasets.

    Current primary disease:

        Prostate cancer
        TCGA code: PRAD

    Source representations may include values such as:

        PRAD
        PrCa
        prostate_cancer
        Prostate cancer
        prostate carcinoma

    M3.5 maps explicitly recognized synonyms to a canonical disease
    identity while preserving the original source value.

    Unknown or non-primary disease labels are preserved and are never
    automatically converted into prostate cancer.

    This module intentionally does not:

        - perform external ontology lookup
        - query EFO, MONDO, DOID, ICD, or NCIt
        - infer disease from genes or genomic loci
        - filter non-prostate-cancer records
        - modify source disease labels
        - infer cancer type from missing values

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


class DiseaseIdentityStatus(str, Enum):
    """Resolution status for one source disease label."""

    RESOLVED = "RESOLVED"
    PRESERVED = "PRESERVED"
    UNRESOLVED = "UNRESOLVED"


class DiseaseContextType(str, Enum):
    """Relationship of a disease label to the primary study disease."""

    PRIMARY_DISEASE = "PRIMARY_DISEASE"
    OTHER_DISEASE = "OTHER_DISEASE"
    UNKNOWN = "UNKNOWN"


class DiseaseIdentifierSystem(str, Enum):
    """Disease identifier namespaces recognized during M3.5."""

    TCGA = "TCGA"
    CANONICAL_INTERNAL = "CANONICAL_INTERNAL"
    SOURCE_LABEL = "SOURCE_LABEL"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class HarmonizedDisease:
    """
    Canonical representation of one disease-context label.

    Attributes:
        source_value:
            Original source disease/cancer label.

        normalized_source_value:
            Source label after representation-only normalization.

        canonical_disease_id:
            Internal canonical disease identifier.

        canonical_name:
            Human-readable canonical disease name.

        tcga_code:
            TCGA project/cancer code when explicitly known.

        identifier_system:
            Namespace supporting the canonical identity.

        context_type:
            Whether the record refers to the primary disease.

        identity_key:
            Stable disease identity key.

        usable:
            Whether this disease context may participate in canonical
            disease matching.

        is_primary_disease:
            Whether this represents prostate cancer.

        status:
            Harmonization status.

        note:
            Optional provenance information.
    """

    source_value: str | None
    normalized_source_value: str | None

    canonical_disease_id: str | None
    canonical_name: str | None
    tcga_code: str | None

    identifier_system: DiseaseIdentifierSystem
    context_type: DiseaseContextType

    identity_key: str | None

    usable: bool
    is_primary_disease: bool

    status: DiseaseIdentityStatus

    note: str | None = None

    def to_dict(
        self,
    ) -> dict[str, Any]:
        """Return JSON-serializable disease metadata."""

        result = asdict(
            self
        )

        result[
            "identifier_system"
        ] = self.identifier_system.value

        result[
            "context_type"
        ] = self.context_type.value

        result[
            "status"
        ] = self.status.value

        return result


class DiseaseHarmonizer:
    """
    Harmonize disease and cancer-context identifiers across sources.

    M3.5 currently resolves prostate-cancer synonyms required by the
    primary integration analysis. Other non-empty labels are preserved
    without asserting ontology equivalence.
    """

    PRIMARY_DISEASE_ID = "prostate_cancer"
    PRIMARY_DISEASE_NAME = "Prostate cancer"
    PRIMARY_TCGA_CODE = "PRAD"

    PRIMARY_ALIASES = frozenset(
        {
            "prad",
            "prca",
            "prostate cancer",
            "prostate carcinoma",
            "prostate_cancer",
            "prostate carcinoma cancer",
            "prostatic cancer",
            "prostatic carcinoma",
            "cancer of prostate",
            "cancer of the prostate",
            "malignant neoplasm of prostate",
            "malignant neoplasm of the prostate",
        }
    )

    SEPARATOR_PATTERN = re.compile(
        r"[\s_\-]+"
    )

    @staticmethod
    def normalize_source_value(
        value: Any,
    ) -> str | None:
        """
        Normalize disease label representation without changing meaning.

        The source value itself remains separately preserved.
        """

        if value is None:
            return None

        text = str(
            value
        ).strip()

        if not text:
            return None

        return text

    @classmethod
    def normalize_for_matching(
        cls,
        value: Any,
    ) -> str | None:
        """
        Normalize a disease label for conservative alias matching.

        Examples:
            "Prostate Cancer"  -> "prostate cancer"
            "prostate_cancer" -> "prostate cancer"
            "PRAD"            -> "prad"
        """

        text = cls.normalize_source_value(
            value
        )

        if text is None:
            return None

        normalized = text.casefold()

        normalized = (
            cls.SEPARATOR_PATTERN
            .sub(
                " ",
                normalized,
            )
        )

        normalized = " ".join(
            normalized.split()
        )

        return normalized

    @classmethod
    def is_primary_alias(
        cls,
        value: Any,
    ) -> bool:
        """Return whether a source label explicitly denotes prostate cancer."""

        normalized = cls.normalize_for_matching(
            value
        )

        if normalized is None:
            return False

        aliases = {
            cls.normalize_for_matching(
                alias
            )
            for alias in cls.PRIMARY_ALIASES
        }

        return (
            normalized
            in aliases
        )

    @classmethod
    def harmonize(
        cls,
        value: Any,
    ) -> HarmonizedDisease:
        """
        Harmonize one source disease-context label.

        Explicitly recognized prostate-cancer aliases are resolved to
        one canonical disease identity.

        Other non-empty values are preserved without claiming ontology
        equivalence.
        """

        source = cls.normalize_source_value(
            value
        )

        normalized = cls.normalize_for_matching(
            value
        )

        if source is None:
            return HarmonizedDisease(
                source_value=None,
                normalized_source_value=None,
                canonical_disease_id=None,
                canonical_name=None,
                tcga_code=None,
                identifier_system=(
                    DiseaseIdentifierSystem.UNKNOWN
                ),
                context_type=(
                    DiseaseContextType.UNKNOWN
                ),
                identity_key=None,
                usable=False,
                is_primary_disease=False,
                status=(
                    DiseaseIdentityStatus.UNRESOLVED
                ),
            )

        if cls.is_primary_alias(
            source
        ):
            return HarmonizedDisease(
                source_value=source,
                normalized_source_value=normalized,
                canonical_disease_id=(
                    cls.PRIMARY_DISEASE_ID
                ),
                canonical_name=(
                    cls.PRIMARY_DISEASE_NAME
                ),
                tcga_code=(
                    cls.PRIMARY_TCGA_CODE
                ),
                identifier_system=(
                    DiseaseIdentifierSystem
                    .CANONICAL_INTERNAL
                ),
                context_type=(
                    DiseaseContextType
                    .PRIMARY_DISEASE
                ),
                identity_key=(
                    "disease:"
                    f"{cls.PRIMARY_DISEASE_ID}"
                ),
                usable=True,
                is_primary_disease=True,
                status=(
                    DiseaseIdentityStatus.RESOLVED
                ),
            )

        return HarmonizedDisease(
            source_value=source,
            normalized_source_value=normalized,
            canonical_disease_id=None,
            canonical_name=None,
            tcga_code=None,
            identifier_system=(
                DiseaseIdentifierSystem
                .SOURCE_LABEL
            ),
            context_type=(
                DiseaseContextType
                .OTHER_DISEASE
            ),
            identity_key=(
                None
            ),
            usable=False,
            is_primary_disease=False,
            status=(
                DiseaseIdentityStatus
                .PRESERVED
            ),
            note=(
                "Non-primary disease label preserved without external "
                "ontology resolution."
            ),
        )

    @classmethod
    def diseases_equal(
        cls,
        left: HarmonizedDisease,
        right: HarmonizedDisease,
    ) -> bool:
        """Compare two canonical disease identities."""

        if not (
            left.usable
            and right.usable
        ):
            return False

        if (
            left.identity_key is None
            or right.identity_key is None
        ):
            return False

        return (
            left.identity_key
            == right.identity_key
        )

    @classmethod
    def both_primary_disease(
        cls,
        left: HarmonizedDisease,
        right: HarmonizedDisease,
    ) -> bool:
        """Return whether both records represent the primary disease."""

        return (
            left.is_primary_disease
            and right.is_primary_disease
        )