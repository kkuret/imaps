"""Base validation."""
import datetime
import os

VALID_BED_EXTENSIONS = (
    ".bed",
    ".bed.gz",
)
VALID_BAM_EXTENSIONS = (
    ".bam",
    ".sam",
)


def validate_bed_file(fname, check_exist=False):
    """Validate BED file."""
    if not fname.endswith(VALID_BED_EXTENSIONS):
        raise ValueError(f"Bed file {fname} should have a valid bed extension.")

    if check_exist and not os.path.isfile(fname):
        raise ValueError(f"Bed file {fname} does not exist.")


def validate_bam_file(fname, check_exist=False):
    """Validate BAM file."""
    if not fname.endswith(VALID_BAM_EXTENSIONS):
        raise ValueError(f"Bam file {fname} should have a valid bam extension.")

    if check_exist and not os.path.isfile(fname):
        raise ValueError(f"Bam file {fname} does not exist.")


def validate_string(value, choices=None, allow_empty=False):
    """Validate string."""
    if not value and allow_empty:
        return

    if not isinstance(value, str):
        raise ValueError(f"Value {value} should be a string.")

    if choices and value not in choices:
        choices_ = ", ".join(choices)
        raise ValueError(f"Value {value} should be one of {choices_}.")


def validate_integer(value):
    """Validate integer."""
    if not isinstance(value, int):
        raise ValueError(f"Value {value} should be an integer.")


def validate_date(value, allow_empty=False):
    """Validate date format."""
    if not value and allow_empty:
        return

    try:
        datetime.datetime.strptime(value, "%Y-%m-%d")
    except ValueError:
        raise ValueError(f"Incorrect date format ({value}), should be YYYY-MM-DD.")
