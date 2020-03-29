"""Extract sites from a bam file into a BED6 file."""
from collections import defaultdict

import pybedtools as pbt
import pysam
from imaps.base.operation import BaseOperation
from imaps.base.site import Site
from imaps.base.validation import validate_bam_file, validate_bed_file, validate_integer, validate_string


class Sites(BaseOperation):
    """Extract sites from a bam file into a BED6 file."""

    umi_separator = ":"

    def __init__(self, bam, output, quant="cDNA", multimax=1, group_by="start"):
        """
        Initialize attributes.

        :param bam: input BAM file
        :type bam: str
        :param output: Output BED6 file
        :type output: str
        :param quant: Report number of cDNA or reads
        :type quant: str
        :param multimax: Ignore reads, mapped more than multimax times
        :type multimax: int
        :param group_by: Assign score to start / middle / end
        :type group_by: str

        """
        # Inputs
        self.bam = bam
        self.output = output
        self.quant = quant
        self.multimax = int(multimax)
        self.group_by = group_by

        super().__init__()

        # Internal variables
        self._counts = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

    def validate_inputs(self):
        """Validate inputs."""
        validate_bam_file(self.bam, check_exist=True)
        validate_bed_file(self.output)
        validate_string(self.quant, choices=["cDNA", "reads"])
        validate_integer(self.multimax)
        validate_string(self.group_by, choices=["start", "middle", "end"])

    def main(self):
        """Root method of the analysis."""
        for read in pysam.AlignmentFile(self.bam, "rb", require_index=True):  # pylint: disable=no-member
            if not self.read_passes_filters(read):
                continue

            site = self.get_site(read)
            umi = self.get_umi(read)

            self._counts[site.chrom][site.pos][site.strand][umi].append(site.score)

        self.write_sites()

    def read_passes_filters(self, read):
        """Check if read passes all required filters."""
        if read.is_qcfail:
            return False
        if read.is_unmapped:
            return False
        if not read.has_tag("NH"):
            return False
        if read.get_tag("NH") > self.multimax:
            return False

        return True

    def get_site(self, read):
        """Get site."""
        chrom = read.reference_name
        positions = read.get_reference_positions()

        if read.is_reverse:
            strand = "-"
            if self.group_by == "start":
                pos = max(positions) + 1
            elif self.group_by == "middle":
                index = len(positions) // 2
                pos = positions[index]
            elif self.group_by == "end":
                pos = min(positions)
        else:
            strand = "+"
            if self.group_by == "start":
                pos = min(positions) - 1
            elif self.group_by == "middle":
                index = len(positions) // 2
                pos = positions[index]
            elif self.group_by == "end":
                pos = max(positions)

        score = 1 / read.get_tag("NH")

        return Site(chrom=chrom, strand=strand, pos=pos, score=score)

    def get_umi(self, read):
        """Get UMI from read."""
        return read.qname.split(self.umi_separator)[-1]

    def get_count(self, umis):
        """
        Get count on the basis of umis.

        Argument umis should be a dict where keys are UMI's and values are
        counts of reads with such UMI.
        """
        if self.quant == "reads":
            return sum([len(list_) for list_ in umis.values()])

        return sum([sum(list_) / len(list_) for list_ in umis.values()])

    def write_sites(self):
        """Write sites to file."""
        intervals = []
        for chrom, data1 in self._counts.items():
            for pos, data2 in data1.items():
                for strand, umis in data2.items():
                    # BED6 file needs an integer value for the "score" column
                    score = str(int(self.get_count(umis)))
                    intervals.append(pbt.create_interval_from_list([chrom, pos, pos + 1, ".", score, strand]))

        pbt.BedTool(interval for interval in intervals).sort().saveas(self.output)
