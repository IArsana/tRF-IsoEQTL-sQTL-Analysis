"""
PCa-tRFQTL Research Pipeline
=============================

File:
    tests/test_cli.py

Description:
    Unit tests for the pcatRFQTL command-line interface.

    The CLI uses argparse subcommands. Tests therefore verify parser
    behavior without assuming that the root command can execute without
    a required subcommand.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import sys

import pytest

from pcatrfqtl.cli import main


def test_cli_requires_command(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Root CLI without a subcommand should exit with argparse code 2."""

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pcatrfqtl",
        ],
    )

    with pytest.raises(
        SystemExit
    ) as exc_info:
        main()

    assert (
        exc_info.value.code
        == 2
    )


def test_cli_help_runs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """CLI help should exit successfully."""

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pcatrfqtl",
            "--help",
        ],
    )

    with pytest.raises(
        SystemExit
    ) as exc_info:
        main()

    assert (
        exc_info.value.code
        == 0
    )