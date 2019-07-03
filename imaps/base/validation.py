"""Base validation."""
import os

VALID_BED_EXTENSIONS = (
    '.bed',
    '.bed.gz',
)


def validate_bed_file(fname, check_exist=False):
    """Vaidate BED file."""
    if not fname.endswith(VALID_BED_EXTENSIONS):
        raise ValueError("Bed file {} should have a valid bed extesion.")

    if check_exist and not os.path.isfile(fname):
        raise ValueError("Bed file {} does not exist.")
