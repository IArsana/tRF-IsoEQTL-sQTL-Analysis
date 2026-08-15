"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/io/parquet.py

Description:
    Shared Parquet input/output utilities for the pcatRFQTL research
    pipeline.

    Pandas object columns may preserve heterogeneous source values such
    as strings, integers, and floating-point values within the same
    column. PyArrow requires a consistent physical representation when
    serializing such data to Parquet.

    This module provides a conservative preparation layer that:

        - preserves DataFrame row cardinality
        - preserves column names
        - preserves nested object columns
        - preserves homogeneous object columns
        - converts heterogeneous scalar object columns to pandas
          nullable string dtype before Parquet serialization
        - preserves missing values using pandas StringDtype
        - never mutates the caller's DataFrame in place

    This module intentionally does not:

        - alter scientific values
        - perform data standardization
        - perform identifier harmonization
        - flatten nested objects
        - infer biological meaning
        - change numerical columns that already have stable dtypes

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


class ParquetPreparationError(RuntimeError):
    """Raised when a DataFrame cannot be prepared safely for Parquet."""


def _is_nested_value(
    value: Any,
) -> bool:
    """
    Return whether a value represents a nested/container structure.

    Nested values are intentionally left untouched because converting
    them to strings would destroy their structured representation.
    """

    if isinstance(
        value,
        (
            list,
            tuple,
            dict,
            set,
        ),
    ):
        return True

    if isinstance(
        value,
        (
            str,
            bytes,
        ),
    ):
        return False

    shape = getattr(
        value,
        "shape",
        None,
    )

    if shape is not None:
        return True

    return False


def _python_scalar_types(
    series: pd.Series,
) -> set[type[Any]]:
    """
    Return Python types represented among non-missing scalar values.

    Nested/container values are excluded because they are handled
    independently.
    """

    value_types: set[
        type[Any]
    ] = set()

    for value in series:

        if value is None:
            continue

        if _is_nested_value(
            value
        ):
            continue

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

        value_types.add(
            type(
                value
            )
        )

    return value_types


def prepare_parquet_dataframe(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """
    Return a Parquet-safe copy of a DataFrame.

    Heterogeneous scalar ``object`` columns are converted to pandas
    nullable string dtype.

    Example:

        position_raw:
            "12345"
            67890.0
            None

    becomes:

        position_raw:
            "12345"
            "67890.0"
            <NA>

    Scientific/standardized numerical columns with stable numeric
    dtypes are not modified.

    Args:
        dataframe:
            DataFrame to prepare.

    Returns:
        A new DataFrame suitable for PyArrow Parquet serialization.

    Raises:
        TypeError:
            If ``dataframe`` is not a pandas DataFrame.

        ParquetPreparationError:
            If preparation unexpectedly alters row cardinality.
    """

    if not isinstance(
        dataframe,
        pd.DataFrame,
    ):
        raise TypeError(
            "prepare_parquet_dataframe expects a pandas DataFrame."
        )

    original_rows = len(
        dataframe
    )

    result = (
        dataframe.copy()
    )

    for column in result.columns:

        series = result[
            column
        ]

        if not pd.api.types.is_object_dtype(
            series.dtype
        ):
            continue

        non_missing = (
            series.dropna()
        )

        if non_missing.empty:
            continue

        contains_nested = any(
            _is_nested_value(
                value
            )
            for value in non_missing
        )

        if contains_nested:
            continue

        represented_types = (
            _python_scalar_types(
                non_missing
            )
        )

        if (
            len(
                represented_types
            )
            <= 1
        ):
            continue

        result[
            column
        ] = (
            series.astype(
                "string"
            )
        )

    if (
        len(
            result
        )
        != original_rows
    ):
        raise ParquetPreparationError(
            "Parquet preparation unexpectedly changed row cardinality: "
            f"{original_rows} -> {len(result)}"
        )

    return result


def write_parquet(
    dataframe: pd.DataFrame,
    path: str | Path,
    *,
    index: bool = False,
) -> Path:
    """
    Write a DataFrame to Parquet using conservative type preparation.

    Parent directories are created automatically.

    Args:
        dataframe:
            DataFrame to serialize.

        path:
            Destination Parquet file.

        index:
            Whether to persist the pandas index.

    Returns:
        Resolved destination Path.
    """

    output_path = Path(
        path
    )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    prepared = (
        prepare_parquet_dataframe(
            dataframe
        )
    )

    prepared.to_parquet(
        output_path,
        index=index,
        engine="pyarrow",
    )

    return output_path


def read_parquet(
    path: str | Path,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Read a Parquet file.

    This lightweight wrapper provides a common I/O entry point for the
    pipeline while leaving scientific interpretation to higher layers.
    """

    input_path = Path(
        path
    )

    if not input_path.exists():
        raise FileNotFoundError(
            f"Parquet file not found: {input_path}"
        )

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Expected Parquet file but found: {input_path}"
        )

    return pd.read_parquet(
        input_path,
        **kwargs,
    )