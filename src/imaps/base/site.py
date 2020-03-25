"""Site."""


class Site:  # pylint: disable=too-few-public-methods
    """Base class for CLIP site."""

    def __init__(self, chrom=None, strand=None, pos=None, score=None):
        """Initialize attributes."""
        self.chrom = chrom
        self.strand = strand
        self.pos = pos
        self.score = score
