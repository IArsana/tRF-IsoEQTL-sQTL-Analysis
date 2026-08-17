"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m6_2b_regulatory_gene_annotation.py

Description:
    Smoke test for M6.2B source-native regulatory feature → gene annotation.

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

from pcatrfqtl.analysis.m6.runners.annotate_regulatory_genes import (
    M62BRegulatoryGeneAnnotationRunner,
)


ROOT = Path(
    __file__
).resolve().parents[
    1
]

CONFIG_PATH = (
    ROOT
    / "configs"
    / "m6_regulatory_gene_annotation.yaml"
)

OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "regulatory_gene_annotation"
)

QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m6"
    / "qc"
)


def main() -> None:

    print()
    print("=" * 72)
    print("M6.2B REGULATORY FEATURE → SOURCE GENE ANNOTATION")
    print("=" * 72)

    runner = M62BRegulatoryGeneAnnotationRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M6.2B RESULT SUMMARY")
    print("=" * 72)

    print(
        "Features assessed:             "
        f"{summary['features_assessed']}"
    )

    print(
        "Source genes resolved:         "
        f"{summary['features_source_gene_resolved']}"
    )

    print(
        "Source genes unresolved:       "
        f"{summary['features_source_gene_unresolved']}"
    )

    print(
        "Nearest-gene mapping:          "
        f"{summary['nearest_gene_mapping_performed']}"
    )

    print(
        "SNP-proximity mapping:         "
        f"{summary['snp_proximity_mapping_performed']}"
    )

    print(
        "External GTF mapping:          "
        f"{summary['external_gtf_mapping_performed']}"
    )

    print()
    print("Feature resolution:")
    print()

    for row in report[
        "feature_resolution"
    ]:

        print(
            f"  {row['regulatory_feature_id']}"
        )

        print(
            "    Parent gene:       "
            f"{row['parent_gene']}"
        )

        print(
            "    Source matches:    "
            f"{row['source_annotation_matches']}"
        )

        print(
            "    Status:            "
            f"{row['gene_annotation_status']}"
        )

        print(
            "    Continue M6.2C:    "
            f"{row['continue_to_m6_2c']}"
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

    print()
    print("=" * 72)
    print("M6.2B EXECUTION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()