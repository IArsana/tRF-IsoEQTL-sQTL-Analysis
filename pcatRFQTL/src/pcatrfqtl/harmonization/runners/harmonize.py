"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/harmonization/runners/harmonize.py

Description:
    Dataset-level execution runner for M3 data harmonization.

    M3.6 applies the harmonization models developed in M3.1-M3.5 to
    standardized Parquet datasets generated during M2.

    Harmonization layers:

        - genome-build provenance
        - canonical genomic coordinates
        - QTL and tag/GWAS variant identity
        - biological feature identity
        - disease context

    Responsibilities:

        - read standardized Parquet datasets
        - preserve row cardinality
        - append harmonized metadata
        - preserve original standardized columns
        - distinguish QTL variants from LD tag/GWAS variants
        - write harmonized Parquet outputs
        - generate harmonization QC summaries

    This runner intentionally does not:

        - modify raw datasets
        - modify M2 standardized datasets
        - perform genome liftover
        - query external annotation databases
        - filter to prostate cancer
        - deduplicate records
        - resolve historical dbSNP aliases
        - perform final cross-source joins

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
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from pcatrfqtl.harmonization.disease import (
    DiseaseHarmonizer,
)
from pcatrfqtl.harmonization.features import (
    FeatureHarmonizer,
    FeatureType,
)
from pcatrfqtl.harmonization.variants import (
    VariantHarmonizer,
)
from pcatrfqtl.logging.logger import (
    get_logger,
)


logger = get_logger(__name__)


@dataclass(frozen=True)
class HarmonizationInputs:
    """
    Paths required for M3 harmonization.

    The standardized directory must correspond to the complete M2
    output tree.
    """

    standardized_directory: Path


class HarmonizationRunner:
    """Execute M3 harmonization over all standardized datasets."""

    # ==================================================================
    # Standardized input locations
    # ==================================================================

    TRFQTL_FILES = {
        "S2": "s2.parquet",
        "S10": "s10.parquet",
    }

    MORADI_QTL_FILES = {
        "s1_cis_intron": {
            "filename": "s1_cis_intron.parquet",
            "feature_type": FeatureType.INTRON_EVENT,
            "disease_context": "PrCa",
        },
        "s1_cis_exon": {
            "filename": "s1_cis_exon.parquet",
            "feature_type": FeatureType.EXON_EVENT,
            "disease_context": "PrCa",
        },
        "s1_cis_isoform": {
            "filename": "s1_cis_isoform.parquet",
            "feature_type": FeatureType.TRANSCRIPT_ISOFORM,
            "disease_context": "PrCa",
        },

        # Moradi S2 LD-linked cis QTL tables.
        "s2_cis_intron": {
            "filename": "s2_cis_intron.parquet",
            "feature_type": FeatureType.INTRON_EVENT,
            "disease_context": None,
        },
        "s2_cis_exon": {
            "filename": "s2_cis_exon.parquet",
            "feature_type": FeatureType.EXON_EVENT,
            "disease_context": None,
        },
        "s2_cis_isoform": {
            "filename": "s2_cis_isoform.parquet",
            "feature_type": FeatureType.TRANSCRIPT_ISOFORM,
            "disease_context": None,
        },

        "s3_trans_intron": {
            "filename": "s3_trans_intron.parquet",
            "feature_type": FeatureType.INTRON_EVENT,
            "disease_context": "PrCa",
        },
        "s3_trans_exon": {
            "filename": "s3_trans_exon.parquet",
            "feature_type": FeatureType.EXON_EVENT,
            "disease_context": "PrCa",
        },
        "s3_trans_isoform": {
            "filename": "s3_trans_isoform.parquet",
            "feature_type": FeatureType.TRANSCRIPT_ISOFORM,
            "disease_context": "PrCa",
        },

        # Moradi S4 source semantics are intentionally unusual:
        #
        # Trans-intron:
        #     mRNA_isoform contains transcript identifiers.
        #
        # Trans-exon:
        #     splicing_event contains exon identifiers.
        #
        # Trans-iso:
        #     splicing_event contains INT identifiers.
        #
        # These mappings therefore preserve the empirically verified
        # source semantics instead of trusting sheet names alone.
        "s4_trans_intron": {
            "filename": "s4_trans_intron.parquet",
            "feature_type": FeatureType.TRANSCRIPT_ISOFORM,
            "disease_context": None,
        },
        "s4_trans_exon": {
            "filename": "s4_trans_exon.parquet",
            "feature_type": FeatureType.EXON_EVENT,
            "disease_context": None,
        },
        "s4_trans_isoform": {
            "filename": "s4_trans_isoform.parquet",
            "feature_type": FeatureType.INTRON_EVENT,
            "disease_context": None,
        },
    }

    MORADI_DE_FILES = {
        "S5": {
            "filename": "s5.parquet",
            "feature_type": FeatureType.EXON_EVENT,
        },
        "S6": {
            "filename": "s6.parquet",
            "feature_type": FeatureType.EXON_EVENT,
        },
        "S7": {
            "filename": "s7.parquet",
            "feature_type": FeatureType.INTRON_EVENT,
        },
        "S8": {
            "filename": "s8.parquet",
            "feature_type": FeatureType.INTRON_EVENT,
        },
        "S9": {
            "filename": "s9.parquet",
            "feature_type": FeatureType.TRANSCRIPT_ISOFORM,
        },
        "S10": {
            "filename": "s10.parquet",
            "feature_type": FeatureType.TRANSCRIPT_ISOFORM,
        },
    }

    def __init__(
        self,
        inputs: HarmonizationInputs,
        output_directory: str | Path,
    ) -> None:
        """Initialize the M3 harmonization runner."""

        self.inputs = inputs

        self.output_directory = Path(
            output_directory
        )

    # ==================================================================
    # Filesystem helpers
    # ==================================================================

    @staticmethod
    def _ensure_file(
        path: Path,
    ) -> None:
        """Require a source file."""

        if not path.exists():
            raise FileNotFoundError(
                "Required harmonization input not found: "
                f"{path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                "Expected file but found: "
                f"{path}"
            )

    @staticmethod
    def _ensure_directory(
        path: Path,
    ) -> None:
        """Require a source directory."""

        if not path.exists():
            raise FileNotFoundError(
                "Required harmonization directory not found: "
                f"{path}"
            )

        if not path.is_dir():
            raise NotADirectoryError(
                "Expected directory but found: "
                f"{path}"
            )

    @staticmethod
    def _prepare_directory(
        path: Path,
        *,
        clear: bool = False,
    ) -> None:
        """Create an output directory."""

        if (
            clear
            and path.exists()
        ):
            shutil.rmtree(
                path
            )

        path.mkdir(
            parents=True,
            exist_ok=True,
        )

    @staticmethod
    def _write_parquet(
        dataframe: pd.DataFrame,
        path: Path,
    ) -> None:
        """Write harmonized output."""

        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        dataframe.to_parquet(
            path,
            index=False,
            engine="pyarrow",
        )

    # ==================================================================
    # Generic schema helpers
    # ==================================================================

    @classmethod
    def _value(
        cls,
        row: pd.Series,
        candidates: Iterable[str],
    ) -> Any:
        """
        Return the first non-missing value among candidate columns.

        Candidate ordering expresses semantic preference.

        A preferred column that exists but contains a missing scalar
        value must not block fallback to later standardized or raw
        provenance fields.

        Example:

            feature_id = NaN
            feature_raw = "INT1e+05"

        returns:

            "INT1e+05"
        """

        for column in candidates:

            if column not in row.index:
                continue

            value = row[
                column
            ]

            if value is None:
                continue

            if pd.api.types.is_scalar(
                value
            ):
                try:
                    if pd.isna(
                        value
                    ):
                        continue
                except (
                    TypeError,
                    ValueError,
                ):
                    pass

            return value

        return None

    @staticmethod
    def _verify_cardinality(
        *,
        dataset_name: str,
        input_rows: int,
        output_rows: int,
    ) -> None:
        """Require exact row-cardinality preservation."""

        if (
            input_rows
            != output_rows
        ):
            raise RuntimeError(
                f"{dataset_name} harmonization changed cardinality: "
                f"{input_rows} -> {output_rows}"
            )

    @staticmethod
    def _boolean_count(
        dataframe: pd.DataFrame,
        column: str,
    ) -> int:
        """Count explicit True values."""

        if (
            column
            not in dataframe.columns
        ):
            return 0

        return int(
            dataframe[
                column
            ]
            .fillna(False)
            .eq(True)
            .sum()
        )

    @staticmethod
    def _value_counts(
        dataframe: pd.DataFrame,
        column: str,
    ) -> dict[str, int]:
        """Return JSON-safe value counts."""

        if (
            column
            not in dataframe.columns
        ):
            return {}

        counts = (
            dataframe[
                column
            ]
            .fillna(
                "UNKNOWN"
            )
            .astype(
                str
            )
            .value_counts(
                dropna=False
            )
        )

        return {
            str(
                key
            ): int(
                value
            )
            for (
                key,
                value,
            ) in counts.items()
        }

    # ==================================================================
    # Flatten harmonization models
    # ==================================================================

    @staticmethod
    def _variant_columns(
        variant: Any,
    ) -> dict[str, Any]:
        """
        Flatten the primary/QTL HarmonizedVariant.

        harm_variant_* is retained as the general primary-variant
        namespace for backwards compatibility.
        """

        coordinate = (
            variant.coordinate
        )

        return {
            "harm_variant_source_id":
                variant.source_variant_id,

            "harm_variant_rsid":
                variant.canonical_rsid,

            "harm_variant_rsid_usable":
                variant.rsid_usable,

            "harm_chromosome":
                coordinate.chromosome,

            "harm_position":
                coordinate.position,

            "harm_coordinate_key":
                coordinate.coordinate_key,

            "harm_build_aware_coordinate_key":
                coordinate.build_aware_key,

            "harm_source_genome_build":
                coordinate.source_genome_build,

            "harm_source_assembly":
                coordinate.source_assembly,

            "harm_coordinate_usable":
                coordinate.coordinate_usable,

            "harm_coordinate_join_allowed":
                coordinate.coordinate_join_allowed,

            "harm_liftover_required":
                coordinate.liftover_required,

            "harm_liftover_performed":
                coordinate.liftover_performed,

            "harm_target_genome_build":
                coordinate.target_genome_build,

            "harm_target_coordinate_key":
                coordinate.target_coordinate_key,

            "harm_variant_identity_key":
                variant.identity_key,

            "harm_variant_identity_method":
                variant.identity_method.value,

            "harm_variant_status":
                variant.status.value,
        }

    @staticmethod
    def _prefixed_variant_columns(
        variant: Any,
        *,
        prefix: str,
    ) -> dict[str, Any]:
        """
        Flatten a variant using an explicit semantic namespace.

        Used for distinguishing Moradi QTL variants from LD-linked
        tag/GWAS variants.
        """

        coordinate = (
            variant.coordinate
        )

        return {
            f"{prefix}_source_id":
                variant.source_variant_id,

            f"{prefix}_rsid":
                variant.canonical_rsid,

            f"{prefix}_rsid_usable":
                variant.rsid_usable,

            f"{prefix}_chromosome":
                coordinate.chromosome,

            f"{prefix}_position":
                coordinate.position,

            f"{prefix}_coordinate_key":
                coordinate.coordinate_key,

            f"{prefix}_build_aware_coordinate_key":
                coordinate.build_aware_key,

            f"{prefix}_source_genome_build":
                coordinate.source_genome_build,

            f"{prefix}_source_assembly":
                coordinate.source_assembly,

            f"{prefix}_coordinate_usable":
                coordinate.coordinate_usable,

            f"{prefix}_coordinate_join_allowed":
                coordinate.coordinate_join_allowed,

            f"{prefix}_liftover_required":
                coordinate.liftover_required,

            f"{prefix}_liftover_performed":
                coordinate.liftover_performed,

            f"{prefix}_target_genome_build":
                coordinate.target_genome_build,

            f"{prefix}_target_coordinate_key":
                coordinate.target_coordinate_key,

            f"{prefix}_identity_key":
                variant.identity_key,

            f"{prefix}_identity_method":
                variant.identity_method.value,

            f"{prefix}_status":
                variant.status.value,
        }

    @staticmethod
    def _feature_columns(
        feature: Any,
    ) -> dict[str, Any]:
        """Flatten a HarmonizedFeature."""

        return {
            "harm_feature_source_id":
                feature.source_feature_id,

            "harm_feature_id":
                feature.canonical_feature_id,

            "harm_feature_type":
                feature.feature_type.value,

            "harm_feature_identifier_system":
                feature.identifier_system.value,

            "harm_feature_identity_key":
                feature.identity_key,

            "harm_feature_usable":
                feature.usable,

            "harm_feature_status":
                feature.status.value,

            "harm_feature_base_id":
                feature.base_identifier,

            "harm_feature_version":
                feature.version,

            "harm_feature_source_anomaly":
                feature.source_anomaly,
        }

    @staticmethod
    def _disease_columns(
        disease: Any,
        *,
        context_source: str,
    ) -> dict[str, Any]:
        """Flatten a HarmonizedDisease."""

        return {
            "harm_disease_source_value":
                disease.source_value,

            "harm_disease_normalized_value":
                disease.normalized_source_value,

            "harm_disease_id":
                disease.canonical_disease_id,

            "harm_disease_name":
                disease.canonical_name,

            "harm_disease_tcga_code":
                disease.tcga_code,

            "harm_disease_identity_key":
                disease.identity_key,

            "harm_disease_context_type":
                disease.context_type.value,

            "harm_disease_identifier_system":
                disease.identifier_system.value,

            "harm_disease_usable":
                disease.usable,

            "harm_is_primary_disease":
                disease.is_primary_disease,

            "harm_disease_status":
                disease.status.value,

            "harm_disease_context_source":
                context_source,
        }

    # ==================================================================
    # Generic row augmentation
    # ==================================================================

    @classmethod
    def _append_records(
        cls,
        dataframe: pd.DataFrame,
        records: list[
            dict[str, Any]
        ],
    ) -> pd.DataFrame:
        """Append harmonized columns without altering source columns."""

        if (
            len(
                dataframe
            )
            != len(
                records
            )
        ):
            raise RuntimeError(
                "Harmonization metadata cardinality does not match "
                "input DataFrame."
            )

        metadata = pd.DataFrame(
            records,
            index=dataframe.index,
        )

        duplicate_columns = (
            set(
                dataframe.columns
            )
            & set(
                metadata.columns
            )
        )

        if duplicate_columns:
            raise RuntimeError(
                "Harmonization attempted to overwrite existing "
                "columns: "
                f"{sorted(duplicate_columns)}"
            )

        return pd.concat(
            [
                dataframe,
                metadata,
            ],
            axis=1,
        )

    # ==================================================================
    # GWAS Catalog
    # ==================================================================

    def _harmonize_gwas_dataframe(
        self,
        dataframe: pd.DataFrame,
    ) -> pd.DataFrame:
        """Harmonize one standardized GWAS Catalog partition."""

        records: list[
            dict[str, Any]
        ] = []

        for _, row in dataframe.iterrows():

            variant = (
                VariantHarmonizer
                .harmonize(
                    dataset=(
                        "gwas_catalog"
                    ),

                    source_variant_id=self._value(
                        row,
                        (
                            "variant_descriptor_raw",
                            "variant_descriptor",
                            "SNPS",
                        ),
                    ),

                    rsid=self._value(
                        row,
                        (
                            "canonical_rsid",
                            "rsid",
                            "SNPS",
                        ),
                    ),

                    chromosome=self._value(
                        row,
                        (
                            "chromosome",
                            "chromosome_standardized",
                            "CHR_ID",
                        ),
                    ),

                    position=self._value(
                        row,
                        (
                            "position",
                            "position_standardized",
                            "CHR_POS",
                        ),
                    ),
                )
            )

            disease_value = (
                self._value(
                    row,
                    (
                        "disease_trait_raw",
                        "disease_trait",
                        "DISEASE/TRAIT",
                        "mapped_trait_raw",
                        "MAPPED_TRAIT",
                    ),
                )
            )

            disease = (
                DiseaseHarmonizer
                .harmonize(
                    disease_value
                )
            )

            record: dict[
                str,
                Any,
            ] = {}

            record.update(
                self._variant_columns(
                    variant
                )
            )

            record.update(
                self._disease_columns(
                    disease,
                    context_source="ROW",
                )
            )

            records.append(
                record
            )

        return self._append_records(
            dataframe,
            records,
        )

    def harmonize_gwas(
        self,
    ) -> dict[str, Any]:
        """Harmonize all standardized GWAS Catalog partitions."""

        source_directory = (
            self.inputs
            .standardized_directory
            / "gwas_catalog"
        )

        self._ensure_directory(
            source_directory
        )

        output_directory = (
            self.output_directory
            / "gwas_catalog"
        )

        self._prepare_directory(
            output_directory,
            clear=True,
        )

        parts = sorted(
            source_directory.glob(
                "part-*.parquet"
            )
        )

        if not parts:
            raise FileNotFoundError(
                "No GWAS Parquet parts found in "
                f"{source_directory}"
            )

        logger.info(
            "Starting GWAS harmonization: %d partitions",
            len(
                parts
            ),
        )

        total_input_rows = 0
        total_output_rows = 0

        rsid_usable_rows = 0
        primary_disease_rows = 0
        coordinate_join_allowed_rows = 0

        variant_status_counts: dict[
            str,
            int,
        ] = {}

        for index, source in enumerate(
            parts
        ):

            dataframe = (
                pd.read_parquet(
                    source
                )
            )

            input_rows = len(
                dataframe
            )

            harmonized = (
                self
                ._harmonize_gwas_dataframe(
                    dataframe
                )
            )

            self._verify_cardinality(
                dataset_name=(
                    f"GWAS {source.name}"
                ),
                input_rows=(
                    input_rows
                ),
                output_rows=len(
                    harmonized
                ),
            )

            output = (
                output_directory
                / source.name
            )

            self._write_parquet(
                harmonized,
                output,
            )

            total_input_rows += (
                input_rows
            )

            total_output_rows += len(
                harmonized
            )

            rsid_usable_rows += (
                self._boolean_count(
                    harmonized,
                    "harm_variant_rsid_usable",
                )
            )

            primary_disease_rows += (
                self._boolean_count(
                    harmonized,
                    "harm_is_primary_disease",
                )
            )

            coordinate_join_allowed_rows += (
                self._boolean_count(
                    harmonized,
                    "harm_coordinate_join_allowed",
                )
            )

            part_status_counts = (
                self._value_counts(
                    harmonized,
                    "harm_variant_status",
                )
            )

            for (
                status,
                count,
            ) in part_status_counts.items():

                variant_status_counts[
                    status
                ] = (
                    variant_status_counts
                    .get(
                        status,
                        0,
                    )
                    + count
                )

            logger.info(
                "GWAS harmonized part %d/%d: %d rows",
                index + 1,
                len(
                    parts
                ),
                len(
                    harmonized
                ),
            )

        return {
            "dataset":
                "gwas_catalog",

            "input_rows":
                total_input_rows,

            "output_rows":
                total_output_rows,

            "cardinality_preserved":
                (
                    total_input_rows
                    == total_output_rows
                ),

            "parts":
                len(
                    parts
                ),

            "rsid_usable_rows":
                rsid_usable_rows,

            "primary_disease_rows":
                primary_disease_rows,

            "coordinate_join_allowed_rows":
                coordinate_join_allowed_rows,

            "variant_status_counts":
                variant_status_counts,
        }

    # ==================================================================
    # Cancer-tRFQTL
    # ==================================================================

    def _harmonize_trfqtl_dataframe(
        self,
        dataframe: pd.DataFrame,
        *,
        table: str,
    ) -> pd.DataFrame:
        """Harmonize standardized Cancer-tRFQTL S2 or S10."""

        records: list[
            dict[str, Any]
        ] = []

        for _, row in dataframe.iterrows():

            variant = (
                VariantHarmonizer
                .harmonize(
                    dataset=(
                        "cancer_trfqtl"
                    ),

                    source_variant_id=self._value(
                        row,
                        (
                            "snp_id_raw",
                            "SNP ID",
                        ),
                    ),

                    rsid=self._value(
                        row,
                        (
                            "canonical_rsid",
                            "snp_id",
                            "SNP ID",
                        ),
                    ),

                    chromosome=self._value(
                        row,
                        (
                            "chromosome",
                            "chromosome_standardized",
                        ),
                    ),

                    position=self._value(
                        row,
                        (
                            "position",
                            "position_standardized",
                        ),
                    ),

                    coordinate=self._value(
                        row,
                        (
                            "coordinate_raw",
                            "coordinate",
                            "SNP position (hg19)",
                            "Position (hg19)",
                        ),
                    ),
                )
            )

            disease = (
                DiseaseHarmonizer
                .harmonize(
                    self._value(
                        row,
                        (
                            "cancer_type",
                            "cancer_type_raw",
                            "Cancer type",
                        ),
                    )
                )
            )

            record: dict[
                str,
                Any,
            ] = {}

            record.update(
                self._variant_columns(
                    variant
                )
            )

            record.update(
                self._disease_columns(
                    disease,
                    context_source="ROW",
                )
            )

            if (
                table
                == "S2"
            ):
                feature = (
                    FeatureHarmonizer
                    .harmonize(
                        self._value(
                            row,
                            (
                                "trf",
                                "trf_raw",
                                "tRF",
                            ),
                        ),
                        feature_type=(
                            FeatureType.TRF
                        ),
                    )
                )

                record.update(
                    self._feature_columns(
                        feature
                    )
                )

            records.append(
                record
            )

        return self._append_records(
            dataframe,
            records,
        )

    def harmonize_trfqtl(
        self,
    ) -> dict[str, Any]:
        """Harmonize standardized Cancer-tRFQTL datasets."""

        source_directory = (
            self.inputs
            .standardized_directory
            / "cancer_trfqtl"
        )

        self._ensure_directory(
            source_directory
        )

        output_directory = (
            self.output_directory
            / "cancer_trfqtl"
        )

        self._prepare_directory(
            output_directory,
            clear=True,
        )

        results: dict[
            str,
            Any,
        ] = {}

        for (
            table,
            filename,
        ) in self.TRFQTL_FILES.items():

            source = (
                source_directory
                / filename
            )

            self._ensure_file(
                source
            )

            logger.info(
                "Harmonizing Cancer-tRFQTL %s",
                table,
            )

            dataframe = (
                pd.read_parquet(
                    source
                )
            )

            harmonized = (
                self
                ._harmonize_trfqtl_dataframe(
                    dataframe,
                    table=table,
                )
            )

            self._verify_cardinality(
                dataset_name=(
                    f"Cancer-tRFQTL {table}"
                ),
                input_rows=len(
                    dataframe
                ),
                output_rows=len(
                    harmonized
                ),
            )

            output = (
                output_directory
                / filename
            )

            self._write_parquet(
                harmonized,
                output,
            )

            table_result: dict[
                str,
                Any,
            ] = {
                "input_rows":
                    int(
                        len(
                            dataframe
                        )
                    ),

                "output_rows":
                    int(
                        len(
                            harmonized
                        )
                    ),

                "cardinality_preserved":
                    (
                        len(
                            dataframe
                        )
                        == len(
                            harmonized
                        )
                    ),

                "primary_disease_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_is_primary_disease",
                    ),

                "variant_rsid_usable_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_variant_rsid_usable",
                    ),

                "coordinate_join_allowed_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_coordinate_join_allowed",
                    ),

                "variant_status_counts":
                    self._value_counts(
                        harmonized,
                        "harm_variant_status",
                    ),

                "disease_status_counts":
                    self._value_counts(
                        harmonized,
                        "harm_disease_status",
                    ),
            }

            if (
                "harm_feature_usable"
                in harmonized.columns
            ):
                table_result[
                    "feature_usable_rows"
                ] = (
                    self._boolean_count(
                        harmonized,
                        "harm_feature_usable",
                    )
                )

                table_result[
                    "feature_status_counts"
                ] = (
                    self._value_counts(
                        harmonized,
                        "harm_feature_status",
                    )
                )

            results[
                table
            ] = (
                table_result
            )

        return {
            "dataset":
                "cancer_trfqtl",

            "tables":
                results,
        }

    # ==================================================================
    # Moradi QTL
    # ==================================================================

    def _harmonize_moradi_qtl_dataframe(
        self,
        dataframe: pd.DataFrame,
        *,
        feature_type: FeatureType,
        implicit_disease: str | None,
    ) -> pd.DataFrame:
        """
        Harmonize one standardized Moradi QTL table.

        S1 and S3 contain a primary QTL variant.

        S2 and S4 additionally contain an LD-linked tag/GWAS variant.
        Both identities are retained explicitly.
        """

        records: list[
            dict[str, Any]
        ] = []

        for _, row in dataframe.iterrows():

            # ----------------------------------------------------------
            # Primary QTL variant
            # ----------------------------------------------------------

            qtl_variant = (
                VariantHarmonizer
                .harmonize(
                    dataset="moradi_qtl",

                    source_variant_id=self._value(
                        row,
                        (
                            "qtl_rsid_raw",
                            "snp_id_raw",
                            "variant_id_raw",
                            "SNP",
                            "T",
                        ),
                    ),

                    rsid=self._value(
                        row,
                        (
                            "qtl_rsid",
                            "canonical_rsid",
                            "rsid",
                            "snp_id",
                            "SNP",
                            "T",
                        ),
                    ),

                    chromosome=self._value(
                        row,
                        (
                            "qtl_chromosome",
                            "chromosome",
                        ),
                    ),

                    position=self._value(
                        row,
                        (
                            "qtl_position",
                            "position",
                        ),
                    ),

                    coordinate=self._value(
                        row,
                        (
                            "qtl_coordinate_raw",
                            "coordinate",
                            "coordinate_raw",
                            "SNP_pos",
                            "sQTL_SNP",
                        ),
                    ),
                )
            )

            # ----------------------------------------------------------
            # Optional LD tag/GWAS variant
            # ----------------------------------------------------------

            has_tag_variant = any(
                column
                in row.index
                for column in (
                    "tag_rsid_raw",
                    "tag_primary_rsid",
                    "tag_coordinate_raw",
                    "tag_chromosome",
                    "tag_position",
                )
            )

            tag_variant = None

            if has_tag_variant:
                tag_variant = (
                    VariantHarmonizer
                    .harmonize(
                        dataset="moradi_qtl",

                        source_variant_id=self._value(
                            row,
                            (
                                "tag_rsid_raw",
                                "tag_primary_rsid",
                            ),
                        ),

                        rsid=self._value(
                            row,
                            (
                                "tag_primary_rsid",
                                "tag_rsid_raw",
                            ),
                        ),

                        chromosome=self._value(
                            row,
                            (
                                "tag_chromosome",
                            ),
                        ),

                        position=self._value(
                            row,
                            (
                                "tag_position",
                            ),
                        ),

                        coordinate=self._value(
                            row,
                            (
                                "tag_coordinate_raw",
                            ),
                        ),
                    )
                )

            # ----------------------------------------------------------
            # Biological feature
            # ----------------------------------------------------------

            feature_value = (
                self._value(
                    row,
                    (
                        "feature_id",
                        "feature_raw",
                        "canonical_feature_id",
                        "feature_id_raw",
                        "splicing_event",
                        "mRNA_isoform",
                    ),
                )
            )

            feature = (
                FeatureHarmonizer
                .harmonize(
                    feature_value,
                    feature_type=(
                        feature_type
                    ),
                )
            )

            # ----------------------------------------------------------
            # Disease context
            # ----------------------------------------------------------

            disease_value = (
                self._value(
                    row,
                    (
                        "gwas_cancer",
                        "cancer_type",
                        "disease",
                    ),
                )
            )

            if (
                disease_value
                is None
            ):
                disease_value = (
                    implicit_disease
                )

                if (
                    implicit_disease
                    is not None
                ):
                    disease_source = (
                        "DATASET"
                    )
                else:
                    disease_source = (
                        "UNAVAILABLE"
                    )

            else:
                disease_source = (
                    "ROW"
                )

            disease = (
                DiseaseHarmonizer
                .harmonize(
                    disease_value
                )
            )

            # ----------------------------------------------------------
            # Assemble harmonized metadata
            # ----------------------------------------------------------

            record: dict[
                str,
                Any,
            ] = {}

            # General/legacy namespace = QTL variant.
            record.update(
                self._variant_columns(
                    qtl_variant
                )
            )

            # Explicit QTL namespace.
            record.update(
                self._prefixed_variant_columns(
                    qtl_variant,
                    prefix=(
                        "harm_qtl_variant"
                    ),
                )
            )

            # Explicit tag/GWAS variant namespace where applicable.
            if (
                tag_variant
                is not None
            ):
                record.update(
                    self._prefixed_variant_columns(
                        tag_variant,
                        prefix=(
                            "harm_tag_variant"
                        ),
                    )
                )

            record.update(
                self._feature_columns(
                    feature
                )
            )

            record.update(
                self._disease_columns(
                    disease,
                    context_source=(
                        disease_source
                    ),
                )
            )

            records.append(
                record
            )

        return self._append_records(
            dataframe,
            records,
        )

    def harmonize_moradi_qtl(
        self,
    ) -> dict[str, Any]:
        """Harmonize standardized Moradi S1-S4 QTL outputs."""

        source_directory = (
            self.inputs
            .standardized_directory
            / "moradi_qtl"
        )

        self._ensure_directory(
            source_directory
        )

        output_directory = (
            self.output_directory
            / "moradi_qtl"
        )

        self._prepare_directory(
            output_directory,
            clear=True,
        )

        results: dict[
            str,
            Any,
        ] = {}

        for (
            name,
            specification,
        ) in self.MORADI_QTL_FILES.items():

            source = (
                source_directory
                / specification[
                    "filename"
                ]
            )

            self._ensure_file(
                source
            )

            logger.info(
                "Harmonizing Moradi QTL: %s",
                name,
            )

            dataframe = (
                pd.read_parquet(
                    source
                )
            )

            harmonized = (
                self
                ._harmonize_moradi_qtl_dataframe(
                    dataframe,
                    feature_type=(
                        specification[
                            "feature_type"
                        ]
                    ),
                    implicit_disease=(
                        specification[
                            "disease_context"
                        ]
                    ),
                )
            )

            self._verify_cardinality(
                dataset_name=(
                    f"Moradi QTL {name}"
                ),
                input_rows=len(
                    dataframe
                ),
                output_rows=len(
                    harmonized
                ),
            )

            output = (
                output_directory
                / specification[
                    "filename"
                ]
            )

            self._write_parquet(
                harmonized,
                output,
            )

            table_result: dict[
                str,
                Any,
            ] = {
                "input_rows":
                    int(
                        len(
                            dataframe
                        )
                    ),

                "output_rows":
                    int(
                        len(
                            harmonized
                        )
                    ),

                "cardinality_preserved":
                    (
                        len(
                            dataframe
                        )
                        == len(
                            harmonized
                        )
                    ),

                # General alias retained for compatibility.
                "variant_rsid_usable_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_variant_rsid_usable",
                    ),

                # Explicit primary QTL identity.
                "qtl_variant_rsid_usable_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_qtl_variant_rsid_usable",
                    ),

                "qtl_coordinate_join_allowed_rows":
                    self._boolean_count(
                        harmonized,
                        (
                            "harm_qtl_variant_"
                            "coordinate_join_allowed"
                        ),
                    ),

                "qtl_variant_status_counts":
                    self._value_counts(
                        harmonized,
                        "harm_qtl_variant_status",
                    ),

                "feature_usable_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_feature_usable",
                    ),

                "feature_source_anomaly_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_feature_source_anomaly",
                    ),

                "feature_status_counts":
                    self._value_counts(
                        harmonized,
                        "harm_feature_status",
                    ),

                "primary_disease_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_is_primary_disease",
                    ),

                "disease_context_sources":
                    self._value_counts(
                        harmonized,
                        "harm_disease_context_source",
                    ),

                "disease_status_counts":
                    self._value_counts(
                        harmonized,
                        "harm_disease_status",
                    ),
            }

            if (
                "harm_tag_variant_rsid_usable"
                in harmonized.columns
            ):
                table_result[
                    "tag_variant_rsid_usable_rows"
                ] = (
                    self._boolean_count(
                        harmonized,
                        "harm_tag_variant_rsid_usable",
                    )
                )

                table_result[
                    "tag_coordinate_join_allowed_rows"
                ] = (
                    self._boolean_count(
                        harmonized,
                        (
                            "harm_tag_variant_"
                            "coordinate_join_allowed"
                        ),
                    )
                )

                table_result[
                    "tag_variant_status_counts"
                ] = (
                    self._value_counts(
                        harmonized,
                        "harm_tag_variant_status",
                    )
                )

            results[
                name
            ] = (
                table_result
            )

        return {
            "dataset":
                "moradi_qtl",

            "tables":
                results,
        }

    # ==================================================================
    # Moradi differential expression
    # ==================================================================

    def _harmonize_moradi_de_dataframe(
        self,
        dataframe: pd.DataFrame,
        *,
        feature_type: FeatureType,
    ) -> pd.DataFrame:
        """Harmonize one Moradi differential-expression table."""

        records: list[
            dict[str, Any]
        ] = []

        for _, row in dataframe.iterrows():

            feature_value = (
                self._value(
                    row,
                    (
                        "feature_id",
                        "feature_raw",
                        "canonical_feature_id",
                        "feature_id_raw",
                        "Unnamed: 0",
                    ),
                )
            )

            feature = (
                FeatureHarmonizer
                .harmonize(
                    feature_value,
                    feature_type=(
                        feature_type
                    ),
                )
            )

            disease = (
                DiseaseHarmonizer
                .harmonize(
                    "PrCa"
                )
            )

            record: dict[
                str,
                Any,
            ] = {}

            record.update(
                self._feature_columns(
                    feature
                )
            )

            record.update(
                self._disease_columns(
                    disease,
                    context_source="DATASET",
                )
            )

            records.append(
                record
            )

        return self._append_records(
            dataframe,
            records,
        )

    def harmonize_moradi_de(
        self,
    ) -> dict[str, Any]:
        """Harmonize Moradi differential-expression outputs."""

        source_directory = (
            self.inputs
            .standardized_directory
            / "moradi_de"
        )

        self._ensure_directory(
            source_directory
        )

        output_directory = (
            self.output_directory
            / "moradi_de"
        )

        self._prepare_directory(
            output_directory,
            clear=True,
        )

        results: dict[
            str,
            Any,
        ] = {}

        for (
            table,
            specification,
        ) in self.MORADI_DE_FILES.items():

            source = (
                source_directory
                / specification[
                    "filename"
                ]
            )

            self._ensure_file(
                source
            )

            logger.info(
                "Harmonizing Moradi DE: %s",
                table,
            )

            dataframe = (
                pd.read_parquet(
                    source
                )
            )

            harmonized = (
                self
                ._harmonize_moradi_de_dataframe(
                    dataframe,
                    feature_type=(
                        specification[
                            "feature_type"
                        ]
                    ),
                )
            )

            self._verify_cardinality(
                dataset_name=(
                    f"Moradi DE {table}"
                ),
                input_rows=len(
                    dataframe
                ),
                output_rows=len(
                    harmonized
                ),
            )

            output = (
                output_directory
                / specification[
                    "filename"
                ]
            )

            self._write_parquet(
                harmonized,
                output,
            )

            results[
                table
            ] = {
                "input_rows":
                    int(
                        len(
                            dataframe
                        )
                    ),

                "output_rows":
                    int(
                        len(
                            harmonized
                        )
                    ),

                "cardinality_preserved":
                    (
                        len(
                            dataframe
                        )
                        == len(
                            harmonized
                        )
                    ),

                "feature_usable_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_feature_usable",
                    ),

                "feature_source_anomaly_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_feature_source_anomaly",
                    ),

                "feature_status_counts":
                    self._value_counts(
                        harmonized,
                        "harm_feature_status",
                    ),

                "primary_disease_rows":
                    self._boolean_count(
                        harmonized,
                        "harm_is_primary_disease",
                    ),
            }

        return {
            "dataset":
                "moradi_de",

            "tables":
                results,
        }

    # ==================================================================
    # Complete M3 execution
    # ==================================================================

    def run_all(
        self,
    ) -> dict[str, Any]:
        """Execute complete M3 harmonized dataset generation."""

        self._ensure_directory(
            self.inputs
            .standardized_directory
        )

        self._prepare_directory(
            self.output_directory
        )

        logger.info(
            "Starting M3 — Data Harmonized."
        )

        report: dict[
            str,
            Any,
        ] = {
            "milestone":
                "M3",

            "stage":
                "data_harmonized",

            "inputs": {
                key: str(
                    value
                )
                for (
                    key,
                    value,
                ) in asdict(
                    self.inputs
                ).items()
            },

            "output_directory":
                str(
                    self.output_directory
                ),

            "policy": {
                "liftover_performed":
                    False,

                "disease_filtering_performed":
                    False,

                "deduplication_performed":
                    False,

                "cross_source_join_performed":
                    False,

                "row_cardinality_preserved":
                    True,

                "qtl_and_tag_variants_distinguished":
                    True,

                "source_anomalies_preserved":
                    True,
            },

            "datasets":
                {},
        }

        report[
            "datasets"
        ][
            "gwas_catalog"
        ] = (
            self.harmonize_gwas()
        )

        report[
            "datasets"
        ][
            "cancer_trfqtl"
        ] = (
            self.harmonize_trfqtl()
        )

        report[
            "datasets"
        ][
            "moradi_qtl"
        ] = (
            self.harmonize_moradi_qtl()
        )

        report[
            "datasets"
        ][
            "moradi_de"
        ] = (
            self.harmonize_moradi_de()
        )

        qc_directory = (
            self.output_directory
            / "qc"
        )

        self._prepare_directory(
            qc_directory
        )

        report_path = (
            qc_directory
            / "harmonization_summary.json"
        )

        report[
            "report_path"
        ] = str(
            report_path
        )

        with report_path.open(
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                report,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        logger.info(
            "M3 harmonization complete."
        )

        logger.info(
            "Harmonization report: %s",
            report_path,
        )

        return report