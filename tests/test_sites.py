"""Test Sites operation."""
# pylint: disable=missing-docstring
import pysam
from imaps.operations.sites import Sites
from ngs_test_utils.testcase import NgsTestCase


class TestSites(NgsTestCase):
    def setUp(self):
        self.bed = self.get_filename(extension="bed")
        self.bam = self.make_bam(
            chroms=[("chr1", 1000)],
            segments=[
                dict(qname="r1:AAA", pos=100, cigar=[(0, 75)], tags={"NH": 1}),
                dict(qname="r2:CCC", pos=100, cigar=[(0, 75)], tags={"NH": 1}),
                dict(qname="r2:CCC", pos=100, cigar=[(0, 75)], tags={"NH": 1}),
                dict(qname="r3:CCC", pos=200, cigar=[(0, 75)], is_reverse=True, tags={"NH": 1}),
            ],
        )

    def test_filters(self):
        # We need sites object to run this method.
        sites = Sites(self.bam, self.bed)

        # QC fail
        read = self.make_bam_segment(cigar=[(0, 75)], is_qcfail=True)
        self.assertFalse(sites.read_passes_filters(read))

        # Unmapped
        read = self.make_bam_segment(cigar=[(0, 75)], is_unmapped=True)
        self.assertFalse(sites.read_passes_filters(read))

        # No NH tag
        read = self.make_bam_segment(cigar=[(0, 75)])
        self.assertFalse(sites.read_passes_filters(read))

        # NH > multimax
        read = self.make_bam_segment(cigar=[(0, 75)], tags={"NH": 2})
        self.assertFalse(sites.read_passes_filters(read))

        # All OK
        read = self.make_bam_segment(cigar=[(0, 75)], tags={"NH": 1})
        self.assertTrue(sites.read_passes_filters(read))

    def test_get_site(self):
        bam = self.make_bam(
            chroms=[("chr1", 1000)],
            segments=[
                dict(pos=40, cigar=[(0, 50)], tags={"NH": 2}),
                dict(pos=40, cigar=[(0, 50)], is_reverse=True, tags={"NH": 2}),
            ],
        )
        afile = pysam.AlignmentFile(bam, "rb")  # pylint: disable=no-member
        read_pos = next(afile)
        read_neg = next(afile)

        # Group by start
        sites = Sites(bam, self.bed, group_by="start")
        site_pos = sites.get_site(read_pos)
        self.assertEqual(site_pos.chrom, "chr1")
        self.assertEqual(site_pos.pos, 39)
        self.assertEqual(site_pos.strand, "+")
        self.assertEqual(site_pos.score, 0.5)
        site_neg = sites.get_site(read_neg)
        self.assertEqual(site_neg.chrom, "chr1")
        self.assertEqual(site_neg.pos, 90)
        self.assertEqual(site_neg.strand, "-")
        self.assertEqual(site_neg.score, 0.5)

        # Group by middle
        sites = Sites(bam, self.bed, group_by="middle")
        site_pos = sites.get_site(read_pos)
        site_neg = sites.get_site(read_neg)
        self.assertEqual(site_pos.pos, 65)
        self.assertEqual(site_neg.pos, 65)

        # Group by end
        sites = Sites(bam, self.bed, group_by="end")
        site_pos = sites.get_site(read_pos)
        site_neg = sites.get_site(read_neg)
        self.assertEqual(site_pos.pos, 89)
        self.assertEqual(site_neg.pos, 40)

    def test_get_umi(self):
        sites = Sites(self.bam, self.bed)

        read = self.make_bam_segment(qname="foo.bar:123:AAA", cigar=[(0, 75)])
        self.assertEqual(sites.get_umi(read), "AAA")

    def test_get_count(self):
        sites_cdna = Sites(self.bam, self.bed, quant="cDNA")
        sites_reads = Sites(self.bam, self.bed, quant="reads")

        # One read on one position uniquely mapped
        umis = {
            "AAA": [1],
        }
        self.assertEqual(sites_cdna.get_count(umis), 1)
        self.assertEqual(sites_reads.get_count(umis), 1)

        # One read on one position multi-mapped (2)
        umis = {
            "AAA": [0.5],
        }
        self.assertEqual(sites_cdna.get_count(umis), 0.5)
        self.assertEqual(sites_reads.get_count(umis), 1)

        # Two reads read on one position uniquely mapped
        umis = {
            "AAA": [1, 1],
        }
        self.assertEqual(sites_cdna.get_count(umis), 1)
        self.assertEqual(sites_reads.get_count(umis), 2)

        # Two reads read on one position, one of them multi-mapped (2)
        umis = {
            "AAA": [1, 0.5],
        }
        self.assertEqual(sites_cdna.get_count(umis), 0.75)
        self.assertEqual(sites_reads.get_count(umis), 2)

        # Complex
        umis = {
            "AAA": [1, 0.5],
            "BBB": [1],
        }
        self.assertEqual(sites_cdna.get_count(umis), 1.75)
        self.assertEqual(sites_reads.get_count(umis), 3)

    def test_run(self):
        Sites(self.bam, self.bed).run()

        self.assertEqual(
            self.tsv_to_list(self.bed), [["chr1", "99", "100", ".", "2", "+"], ["chr1", "275", "276", ".", "1", "-"]],
        )
