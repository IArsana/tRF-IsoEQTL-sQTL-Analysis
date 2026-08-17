"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_4a_open_small_rna_discovery.py

Description:
    Smoke test for M7.4A Open/Public Small-RNA Dataset Discovery.

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

from pcatrfqtl.analysis.m7.runners.discover_open_small_rna_datasets import (
    M74AOpenSmallRnaDiscoveryRunner,
)


ROOT = Path(
    __file__
).resolve().parents[1]


CONFIG_PATH = (
    ROOT
    / "configs"
    / "m7_open_small_rna_discovery.yaml"
)


OUTPUT_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "open_small_rna_discovery"
)


QC_DIRECTORY = (
    ROOT
    / "data"
    / "processed"
    / "m7"
    / "qc"
)


def main() -> None:

    print()
    print("=" * 72)
    print("M7.4A OPEN/PUBLIC SMALL-RNA DATASET DISCOVERY")
    print("=" * 72)

    runner = M74AOpenSmallRnaDiscoveryRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
        output_directory=OUTPUT_DIRECTORY,
        qc_directory=QC_DIRECTORY,
    )

    report = runner.run()

    summary = report[
        "summary"
    ]

    assert (
        summary[
            "datasets_assessed"
        ]
        ==
        4
    )

    assert (
        summary[
            "datasets_discovery_eligible"
        ]
        >=
        1
    )

    assert (
        summary[
            "datasets_with_raw_sequence_route"
        ]
        >=
        1
    )

    assert (
        summary[
            "technical_eligibility_established"
        ]
        is False
    )

    assert (
        summary[
            "candidate_sequence_search_performed"
        ]
        is False
    )

    assert (
        summary[
            "trf_quantification_performed"
        ]
        is False
    )

    assert (
        summary[
            "clinical_association_performed"
        ]
        is False
    )

    assert (
        summary[
            "overall_discovery_status"
        ]
        ==
        "OPEN_SMALL_RNA_VALIDATION_DATASETS_IDENTIFIED"
    )

    print()
    print("Datasets assessed:             "
          f"{summary['datasets_assessed']}")

    print("Discovery eligible:            "
          f"{summary['datasets_discovery_eligible']}")

    print("High priority:                 "
          f"{summary['high_priority_datasets']}")

    print("Intermediate priority:         "
          f"{summary['intermediate_priority_datasets']}")

    print("Exploratory:                   "
          f"{summary['exploratory_datasets']}")

    print("Raw sequence route:            "
          f"{summary['datasets_with_raw_sequence_route']}")

    print("Reported tRNA context:         "
          f"{summary['datasets_with_trna_context']}")

    print("Recurrence metadata:           "
          f"{summary['datasets_with_recurrence_metadata']}")

    print()

    print("Leading dataset:")
    print(
        "  "
        f"{summary['leading_dataset_for_technical_audit']}"
    )

    print()
    print("Dataset discovery:")

    for row in report[
        "dataset_discovery"
    ]:

        print()
        print(
            f"  {row['dataset_id']} | "
            f"{row['priority_class']}"
        )

        print(
            "    samples:                  "
            f"{row['sample_count_series']}"
        )

        print(
            "    raw sequence route:       "
            f"{row['raw_sequence_route']}"
        )

        print(
            "    tRNA context:             "
            f"{row['trna_derived_rna_context_reported']}"
        )

        print(
            "    comparator:               "
            f"{row['tumor_comparator_available']}"
        )

        print(
            "    clinical metadata:        "
            f"{row['clinicopathologic_metadata_reported']}"
        )

        print(
            "    recurrence metadata:      "
            f"{row['recurrence_metadata_reported']}"
        )

        print(
            "    discovery eligible:       "
            f"{row['discovery_eligible']}"
        )

        print(
            "    technical eligibility:    "
            f"{row['technical_eligibility_established']}"
        )

        print(
            "    status:"
        )

        print(
            "      "
            f"{row['discovery_status']}"
        )

    print()
    print("Overall:")
    print(
        "  "
        f"{summary['overall_discovery_status']}"
    )

    print()
    print("Next stage:")
    print(
        "  "
        f"{summary['next_stage']}"
    )

    print()
    print("=" * 72)
    print("M7.4A SMOKE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    main()