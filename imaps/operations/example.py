"""Example operation."""
import pybedtools as pbt
from imaps.base.operation import BaseOperation
from imaps.base.validation import validate_bed_file


class ExampleOperation(BaseOperation):
    """Simple example of an operaton."""

    def __init__(self, sites, outfile, threshold=5):
        """Initialize attributes.

        Parameters
        ----------
        sites : str
            Sites file (BED6 format).
        outfile : str
            Name of output file (BED6 format).
        threshold : int
            Number of sites to keep.

        """
        self.sites = sites
        self.outfile = outfile
        self.threshold = threshold

    def validate_inputs(self):
        """Validate inputs."""
        validate_bed_file(self.sites, check_exist=True)
        validate_bed_file(self.outfile)

    def main(self):
        """Filter out sites that have score lower than threshold."""
        filtered_sites = []
        for item in pbt.BedTool(self.sites):
            if int(item.score) >= self.threshold:
                filtered_sites.append(item)

        pbt.BedTool(filtered_sites).saveas(self.outfile)
