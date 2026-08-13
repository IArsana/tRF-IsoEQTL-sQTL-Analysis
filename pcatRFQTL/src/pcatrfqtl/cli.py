"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/cli.py

Description:
    Command-line interface for the PCa-tRFQTL research pipeline.

    This module exposes project-level commands for dataset inspection,
    data ingestion, validation, and downstream analysis.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pcatrfqtl.inspection.dataset_inspector import DatasetInspector
from pcatrfqtl.inspection.yaml_registry import DatasetRegistry
from pcatrfqtl.logging.logger import get_logger


logger = get_logger(__name__)


def inspect_data(config_path: Path) -> int:
    """
    Inspect all datasets registered in datasets.yaml.

    Parameters
    ----------
    config_path:
        Path to the datasets.yaml configuration file.

    Returns
    -------
    int
        Exit code.

        0:
            All datasets were successfully inspected.

        1:
            One or more datasets failed inspection.
    """

    logger.info("Starting dataset inspection.")
    logger.info("Dataset configuration: %s", config_path)

    try:
        registry = DatasetRegistry(config_path)
        logger.info(
            "Dataset registry loaded successfully."
        )

    except Exception:
        logger.exception(
            "Failed to load dataset registry: %s",
            config_path,
        )
        return 1

    try:
        inspector = DatasetInspector(registry)

        logger.info(
            "Inspecting registered datasets..."
        )

        results = inspector.inspect_all()

    except Exception:
        logger.exception(
            "Dataset inspection failed."
        )
        return 1

    output_path = (
        registry.project_root
        / "data"
        / "interim"
        / "qc"
        / "dataset_inspection.json"
    )

    try:
        inspector.save_json(
            results,
            output_path,
        )

        logger.info(
            "Inspection report saved: %s",
            output_path,
        )

    except Exception:
        logger.exception(
            "Failed to save inspection report: %s",
            output_path,
        )
        return 1

    total = len(results)

    successful = sum(
        result.readable and not result.errors
        for result in results
    )

    failed = total - successful

    logger.info(
        "Dataset inspection completed."
    )

    logger.info(
        "Datasets inspected: %d | Readable: %d | Failed: %d",
        total,
        successful,
        failed,
    )

    print("=" * 72)
    print("pcatRFQTL — Dataset Inspection")
    print("=" * 72)
    print(f"Datasets inspected : {total}")
    print(f"Readable           : {successful}")
    print(f"Failed             : {failed}")
    print(f"Report             : {output_path}")
    print("=" * 72)

    if failed:
        logger.warning(
            "%d dataset(s) failed inspection.",
            failed,
        )

    else:
        logger.info(
            "All registered datasets passed inspection."
        )

    return 0 if successful == total else 1


def build_parser() -> argparse.ArgumentParser:
    """
    Build the project CLI argument parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured command-line argument parser.
    """

    parser = argparse.ArgumentParser(
        prog="pcatrfqtl",
        description="pcatRFQTL analysis pipeline",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    inspect_parser = subparsers.add_parser(
        "inspect-data",
        help="Inspect all registered datasets.",
    )

    inspect_parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/datasets.yaml"),
        help="Path to datasets.yaml.",
    )

    return parser


def main() -> int:
    """
    Execute the command-line interface.

    Returns
    -------
    int
        CLI exit code.
    """

    parser = build_parser()
    args = parser.parse_args()

    logger.debug(
        "CLI command received: %s",
        args.command,
    )

    if args.command == "inspect-data":
        return inspect_data(args.config)

    parser.error(
        f"Unknown command: {args.command}"
    )

    return 1


if __name__ == "__main__":
    raise SystemExit(main())