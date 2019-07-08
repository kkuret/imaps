"""Test example operation."""
# pylint: disable=missing-docstring

from imaps.operations.example import ExampleOperation

from .base import ImapsTestCase


class TestExampleOperation(ImapsTestCase):

    def test_run(self):
        sites = self.create_bed_from_list([
            ['chr1', '2', '3', '.', '9', '+'],
            ['chr1', '3', '4', '.', '3', '+'],

        ])
        outfile = ImapsTestCase.get_filename(extension='bed')

        ExampleOperation(sites, outfile).run()

        self.assert_bed_equal(outfile, [
            ['chr1', '2', '3', '.', '9', '+'],
        ])
