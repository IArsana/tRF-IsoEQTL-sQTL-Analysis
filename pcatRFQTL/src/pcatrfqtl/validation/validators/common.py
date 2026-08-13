"""
Common utilities for dataset-specific validators.

This module provides helper functions shared by validation modules.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

from typing import Iterable


def find_column(
    columns: Iterable[str],
    candidates: Iterable[str],
) -> str | None:
    """
    Find the first matching column from a list of candidates.

    Matching is case-insensitive and ignores surrounding whitespace.
    """

    normalized = {
        str(column).strip().lower(): str(column).strip()
        for column in columns
    }

    for candidate in candidates:
        match = normalized.get(
            str(candidate).strip().lower()
        )

        if match is not None:
            return match

    return None