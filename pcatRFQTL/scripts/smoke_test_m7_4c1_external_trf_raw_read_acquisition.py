"""
PCa-tRFQTL Research Pipeline
=============================

File:
    scripts/smoke_test_m7_4c1_external_trf_raw_read_acquisition.py

Description:
    Smoke test for M7.4C.1 External tRF Raw-Read Acquisition Audit.

    This test confirms that public SRA run metadata for GSE80400 can be
    acquired and converted into a raw-read acquisition manifest without
    downloading sequencing reads.

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

from pcatrfqtl.analysis.m7.runners.audit_external_trf_raw_read_acquisition import (
    M74C1ExternalTrfRawReadAcquisitionRunner,
)


ROOT = Path(
    __file__
).resolve().parents[
    1
]


CONFIG_PATH = (
    ROOT
    /
    "configs"
    /
    "m7_external_trf_raw_read_acquisition.yaml"
)


def validate_expected_state(
    report: dict,
) -> None:
    """Validate expected M7.4C.1 outcome."""

    assert (
        report[
            "milestone"
        ]
        ==
        "M7.4C.1"
    )

    assert (
        report[
            "stage"
        ]
        ==
        "external_trf_raw_read_acquisition_audit"
    )

    summary = report[
        "summary"
    ]

    assert (
        summary[
            "dataset_key"
        ]
        ==
        "m7_gse80400"
    )

    assert (
        summary[
            "geo_accession"
        ]
        ==
        "GSE80400"
    )

    assert (
        summary[
            "sra_study"
        ]
        ==
        "SRP073456"
    )

    assert (
        summary[
            "candidate_trf"
        ]
        ==
        "tRF-24-S3M8309N0Y"
    )

    assert (
        summary[
            "candidate_sequence"
        ]
        ==
        "GTAGTCGTGGCCGAGTGGTTAAGG"
    )

    assert (
        summary[
            "candidate_length_nt"
        ]
        ==
        24
    )

    assert (
        summary[
            "runs_identified"
        ]
        >=
        1
    )

    assert (
        summary[
            "total_spots"
        ]
        >
        0
    )

    assert (
        summary[
            "total_bases"
        ]
        >
        0
    )

    assert (
        summary[
            "public_runs_present"
        ]
        is True
    )

    assert (
        summary[
            "run_metadata_complete"
        ]
        is True
    )

    # Toolchain is intentionally NOT mandatory at M7.4C.1.
    assert isinstance(
        summary[
            "prefetch_available"
        ],
        bool,
    )

    assert isinstance(
        summary[
            "fasterq_dump_available"
        ],
        bool,
    )

    assert (
        summary[
            "raw_sequence_download_performed"
        ]
        is False
    )

    assert (
        summary[
            "fastq_generation_performed"
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
            "exact_match_counting_performed"
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
            "candidate_detection_claimed"
        ]
        is False
    )

    assert (
        summary[
            "candidate_absence_claimed"
        ]
        is False
    )

    assert (
        summary[
            "acquisition_manifest_ready"
        ]
        is True
    )

    assert (
        summary[
            "overall_status"
        ]
        ==
        "EXTERNAL_RAW_READ_ACQUISITION_READY"
    )

    assert (
        summary[
            "next_stage"
        ]
        ==
        "M7.4C.2_EXACT_CANDIDATE_SEQUENCE_COUNTING"
    )

    # ----------------------------------------------------------------------
    # Run-level safeguards
    # ----------------------------------------------------------------------

    runs = report[
        "run_inventory"
    ]

    manifest = report[
        "acquisition_manifest"
    ]

    assert (
        len(
            runs
        )
        ==
        summary[
            "runs_identified"
        ]
    )

    assert (
        len(
            manifest
        )
        ==
        summary[
            "runs_identified"
        ]
    )

    run_accessions = {
        row[
            "run_accession"
        ]
        for row in runs
    }

    assert (
        len(
            run_accessions
        )
        ==
        len(
            runs
        )
    )

    assert all(
        accession.startswith(
            "SRR"
        )
        for accession in run_accessions
    )

    for row in manifest:

        assert (
            row[
                "download_performed"
            ]
            is False
        )

        assert (
            row[
                "fastq_generated"
            ]
            is False
        )

        assert (
            row[
                "read_qc_performed"
            ]
            is False
        )

        assert (
            row[
                "adapter_trimmed"
            ]
            is False
        )

        assert (
            row[
                "candidate_sequence_search_performed"
            ]
            is False
        )

        assert (
            row[
                "candidate_exact_count_performed"
            ]
            is False
        )


def main() -> None:
    """Execute M7.4C.1 smoke test."""

    print()
    print("=" * 72)
    print("M7.4C.1 EXTERNAL tRF RAW-READ ACQUISITION AUDIT")
    print("=" * 72)

    runner = M74C1ExternalTrfRawReadAcquisitionRunner(
        project_root=ROOT,
        config_path=CONFIG_PATH,
    )

    report = runner.run()

    validate_expected_state(
        report
    )

    summary = report[
        "summary"
    ]

    print()
    print("=" * 72)
    print("M7.4C.1 RESULT SUMMARY")
    print("=" * 72)

    print(
        "Dataset:                         "
        f"{summary['geo_accession']}"
    )

    print(
        "SRA study:                       "
        f"{summary['sra_study']}"
    )

    print(
        "Runs identified:                 "
        f"{summary['runs_identified']}"
    )

    print(
        "Total spots:                     "
        f"{summary['total_spots']}"
    )

    print(
        "Total bases:                     "
        f"{summary['total_bases']}"
    )

    print()

    print(
        "Candidate:                       "
        f"{summary['candidate_trf']}"
    )

    print(
        "Sequence:                        "
        f"{summary['candidate_sequence']}"
    )

    print()

    print(
        "Public runs present:             "
        f"{summary['public_runs_present']}"
    )

    print(
        "Run metadata complete:           "
        f"{summary['run_metadata_complete']}"
    )

    print()

    print(
        "prefetch available:              "
        f"{summary['prefetch_available']}"
    )

    print(
        "fasterq-dump available:          "
        f"{summary['fasterq_dump_available']}"
    )

    print(
        "SRA Toolkit ready:               "
        f"{summary['sra_toolkit_ready']}"
    )

    print()

    print(
        "Raw read download performed:     "
        f"{summary['raw_sequence_download_performed']}"
    )

    print(
        "FASTQ generation performed:      "
        f"{summary['fastq_generation_performed']}"
    )

    print(
        "Candidate sequence search:       "
        f"{summary['candidate_sequence_search_performed']}"
    )

    print(
        "Exact match counting:            "
        f"{summary['exact_match_counting_performed']}"
    )

    print(
        "tRF quantification:              "
        f"{summary['trf_quantification_performed']}"
    )

    print()

    print(
        "Acquisition manifest ready:      "
        f"{summary['acquisition_manifest_ready']}"
    )

    print()

    print("Acquisition status:")
    print(
        "  "
        f"{summary['acquisition_status']}"
    )

    print()

    print("Overall:")
    print(
        "  "
        f"{summary['overall_status']}"
    )

    print()

    print("Next:")
    print(
        "  "
        f"{summary['next_stage']}"
    )

    print()
    print("Runs:")
    print()

    for row in report[
        "run_inventory"
    ]:

        print(
            f"  {row['run_accession']} | "
            f"{row['experiment_accession']} | "
            f"{row['spots']} spots | "
            f"{row['bases']} bases | "
            f"{row['library_layout']}"
        )

    print()
    print("=" * 72)
    print("M7.4C.1 SMOKE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    main()