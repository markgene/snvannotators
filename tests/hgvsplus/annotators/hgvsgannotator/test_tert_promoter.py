"""Test HgvsGAnnotator class."""

import unittest

from hgvs.easy import parse

from snvannotators.hgvsplus.annotators.hgvsannotation import HgvsAnnotation
from snvannotators.hgvsplus.annotators.hgvsgannotator import HgvsGAnnotator
from snvannotators.hgvsplus.models import HgvsG, HgvsP, HgvsT


class HgvsGAnnotatorTertPromoterTestCase(unittest.TestCase):
    """Test HgvsGAnnotator class with TERT promoter."""

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        sequence_variant_g = parse("NC_000005.9:g.1295228G>A")
        hgvs_g = HgvsG.from_sequence_variant_g(sequence_variant_g=sequence_variant_g)
        alt_aln_method = "splign"
        tss_upstream_limit = 20000
        uncertain = False
        cls.hgvs_annotation = HgvsGAnnotator(
            hgvs_g=hgvs_g,
            alt_aln_method=alt_aln_method,
            tss_upstream_limit=tss_upstream_limit,
            uncertain=uncertain,
        ).annotate()
        cls.tx_ac = "NM_198253.2"

    def test_hgvs_annotation(self):
        self.assertEqual(self.hgvs_annotation, "")