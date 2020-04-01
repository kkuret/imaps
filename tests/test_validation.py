"""Test validation functions."""
# pylint: disable=missing-docstring
from imaps.base.validation import (
    validate_bam_file,
    validate_bed_file,
    validate_date,
    validate_integer,
    validate_string,
)
from ngs_test_utils.testcase import NgsTestCase


class TestValidation(NgsTestCase):
    def test_validate_bed_file(self):
        message = "Bed file file.txt should have a valid bed extension."
        with self.assertRaisesRegex(ValueError, message):
            validate_bed_file("file.txt", check_exist=False)

        message = "Bed file file.bed does not exist."
        with self.assertRaisesRegex(ValueError, message):
            validate_bed_file("file.bed", check_exist=True)

        bed = self.make_bed(intervals=[["chr1", 10, 20, ".", 12, "+"]])
        validate_bed_file(bed, check_exist=True)

    def test_validate_bam_file(self):
        message = "Bam file file.txt should have a valid bam extension."
        with self.assertRaisesRegex(ValueError, message):
            validate_bam_file("file.txt", check_exist=False)

        message = "Bam file file.bam does not exist."
        with self.assertRaisesRegex(ValueError, message):
            validate_bam_file("file.bam", check_exist=True)

        bam = self.make_bam(chroms=[("chr1", 100)], segments=[dict(cigar=[(0, 75)])])
        validate_bam_file(bam, check_exist=True)

    def test_validate_string(self):
        message = "Value 123 should be a string."
        with self.assertRaisesRegex(ValueError, message):
            validate_string(123)

        message = "Value C should be one of A, B."
        with self.assertRaisesRegex(ValueError, message):
            validate_string("C", choices=["A", "B"])

        validate_string("A")
        validate_string("B", choices=["A", "B"])

        validate_string("", allow_empty=True)

    def test_validate_integer(self):
        message = "Value AAA should be an integer."
        with self.assertRaisesRegex(ValueError, message):
            validate_integer("AAA")

        validate_integer(123)

    def test_validate_date(self):
        message = "Incorrect date format \\(1.2.1990\\), should be YYYY-MM-DD."
        with self.assertRaisesRegex(ValueError, message):
            validate_date("1.2.1990")

        validate_date("1900-2-1")

        validate_date("", allow_empty=True)
