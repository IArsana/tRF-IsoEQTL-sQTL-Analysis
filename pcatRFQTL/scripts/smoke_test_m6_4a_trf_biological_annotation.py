"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m6_4a_trf_biological_annotation.py

Description:
    Smoke test for M6.4A Candidate/tRF Biological Annotation.

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
from pathlib import Path

from pcatrfqtl.analysis.m6.runners.annotate_trf_biology import (
    M64ATrfBiologicalAnnotationRunner,
)


ROOT = Path(__file__).resolve().parents[1]

CONFIG_PATH = (
    ROOT
    / "configs"
    / "m6_trf_biological_annotation.yaml"
)

OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "trf_biological_annotation"
)

QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "qc"
)


def main() -> None:

    runner = M64ATrfBiologicalAnnotationRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    summary = report["summary"]

    print()
    print("=" * 72)
    print("M6.4A RESULT SUMMARY")
    print("=" * 72)

    print(
        "Candidates assessed:                    "
        f"{summary['candidates_assessed']}"
    )

    print(
        "Candidates with tRF identity:           "
        f"{summary['candidates_with_trf_identity']}"
    )

    print(
        "Candidates with regulatory context:     "
        f"{summary['candidates_with_regulatory_feature_context']}"
    )

    print(
        "Candidates with resolved gene:          "
        f"{summary['candidates_with_resolved_gene']}"
    )

    print()

    for row in report[
        "candidate_trf_biological_annotation"
    ]:

        print(
            f"  Rank {row['priority_rank']}: "
            f"{row['lead_rsid']}"
        )

        print(
            "    tRF IDs:                 "
            f"{row['trf_ids']}"
        )

        print(
            "    Chromosome(s):            "
            f"{row['chromosomes']}"
        )

        print(
            "    tRF class(es):            "
            f"{row['trf_classes']}"
        )

        print(
            "    Parent tRNA(s):           "
            f"{row['parent_trnas']}"
        )

        print(
            "    Regulatory feature:       "
            f"{row['regulatory_feature_id']}"
        )

        print(
            "    RNA-processing context:   "
            f"{row['rna_processing_context']}"
        )

        print(
            "    Annotation status:        "
            f"{row['annotation_status']}"
        )

        print()

    print(
        "Next stage:"
    )

    print(
        f"  {summary['next_stage']}"
    )

    print()
    print("Full QC report:")
    print()

    print(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()