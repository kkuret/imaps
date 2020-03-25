"""Base imaps TestCase."""
# pylint: disable=protected-access

import os
import tempfile
import unittest

import pybedtools as pbt


class ImapsTestCase(unittest.TestCase):
    """Base class for imaps TestCases."""

    @staticmethod
    def get_filename(extension=None, directory=None):
        """Get availbale filename."""
        name = next(tempfile._get_candidate_names())

        if not directory:
            directory = tempfile._get_default_tempdir()
        if extension:
            name += '.' + extension

        return os.path.join(directory, name)

    def create_bed_from_list(self, intervals):
        """Create BED file from a list of lists."""
        fname = self.get_filename(extension='bed')
        pbt.BedTool(pbt.create_interval_from_list(item) for item in intervals).saveas(fname)
        return fname

    def assert_bed_equal(self, bedfile, expected):
        """Check if contents of a bedfile are equal to content of expected."""
        self.assertEqual(pbt.BedTool(bedfile).count(), len(expected))

        for real, expct in zip(pbt.BedTool(bedfile), expected):
            self.assertEqual(real.fields, expct)
