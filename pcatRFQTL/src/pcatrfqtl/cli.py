"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/cli.py

Description:
    Command-line interface (CLI) entry point for the PCa-tRFQTL
    research pipeline.

    This module provides the primary interface for executing,
    configuring, and managing pipeline operations from the
    command line.

    The CLI will serve as the main entry point for future
    operations including data acquisition, validation,
    preprocessing, QTL integration, statistical analysis,
    colocalization, functional analysis, and validation.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

from pcatrfqtl.logging.logger import get_logger

logger = get_logger("pcatrfqtl")


def main() -> None:
    logger.info("PCa-tRFQTL pipeline initialized.")


if __name__ == "__main__":
    main()
