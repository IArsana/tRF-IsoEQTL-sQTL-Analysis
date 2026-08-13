"""
Tests for validation result models.

Author:
    I Putu Indra Arsana
"""

from pcatrfqtl.validation.result import (
    ValidationCheck,
    ValidationReport,
    ValidationSeverity,
    ValidationStatus,
)


def test_validation_report_pass() -> None:
    """
    Test a passing validation report.
    """

    report = ValidationReport(
        dataset_id="test_dataset",
    )

    report.add_check(
        ValidationCheck(
            name="test",
            status=ValidationStatus.PASS,
            severity=ValidationSeverity.INFO,
            message="Validation passed.",
        )
    )

    assert report.status == ValidationStatus.PASS
    assert len(report.errors) == 0
    assert len(report.warnings) == 0


def test_validation_report_warning() -> None:
    """
    Test a warning validation report.
    """

    report = ValidationReport(
        dataset_id="test_dataset",
    )

    report.add_check(
        ValidationCheck(
            name="test",
            status=ValidationStatus.WARNING,
            severity=ValidationSeverity.WARNING,
            message="Validation warning.",
        )
    )

    assert report.status == ValidationStatus.WARNING
    assert len(report.warnings) == 1


def test_validation_report_fail() -> None:
    """
    Test a failed validation report.
    """

    report = ValidationReport(
        dataset_id="test_dataset",
    )

    report.add_check(
        ValidationCheck(
            name="test",
            status=ValidationStatus.FAIL,
            severity=ValidationSeverity.ERROR,
            message="Validation failed.",
        )
    )

    assert report.status == ValidationStatus.FAIL
    assert len(report.errors) == 1

def test_validation_report_statistics() -> None:
    """
    Test validation report statistics.
    """

    report = ValidationReport(
        dataset_id="test_dataset",
    )

    report.add_statistics(
        {
            "total_rows": 100,
            "fields": {
                "P-VALUE": {
                    "valid": 100,
                    "missing": 0,
                    "invalid": 0,
                }
            },
        }
    )

    result = report.to_dict()

    assert result["statistics"]["total_rows"] == 100
    assert (
        result["statistics"]["fields"]["P-VALUE"]["valid"]
        == 100
    )