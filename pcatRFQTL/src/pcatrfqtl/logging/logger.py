"""
PCa-tRFQTL Research Pipeline
=============================

File:
    src/pcatrfqtl/logging/logger.py

Description:
    Provides centralized logging utilities for the PCa-tRFQTL
    research pipeline.

    This module is responsible for creating and configuring
    application loggers used throughout the pipeline. It provides
    consistent log formatting and logging levels for monitoring
    data acquisition, preprocessing, statistical analysis,
    integration, validation, and pipeline execution.

Author:
    I Putu Indra Arsana

Project:
    Integrative Analysis of Prostate Cancer Risk Variants,
    tRNA-Derived Fragment QTLs, and Transcript Isoform Regulation

License:
    MIT
"""

import logging
import sys


def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)

    if logger.handlers:
        return logger

    logger.setLevel(logging.INFO)

    handler = logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s")

    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger
