"""Test HgvsParser."""

import unittest

from snvannotators.hgvsplus.models import HgvsC, HgvsG, HgvsP, HgvsT
from snvannotators.hgvsplus.parsers.hgvsparser import HgvsParser
from snvmodels.cpra import Cpra


class CpraGrch37AnnotatorArD891vTestCase(unittest.TestCase):
    """Test CpraGrch37Annotator class with AR D891V."""

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.hgvs_parser = HgvsParser()

    def test_hgvs_g(self):
        s = "NC_000012.11:g.25398266G>A"
        hgvs_g = self.hgvs_parser.parse(s=s)
        self.assertTrue(isinstance(hgvs_g, HgvsG))
        self.assertEqual(hgvs_g.format(), s)

    def test_hgvs_c(self):
        s = "NM_001369787.1(KRAS):c.53C>T"
        hgvs_c = self.hgvs_parser.parse(s=s)
        self.assertTrue(isinstance(hgvs_c, HgvsC))
        self.assertEqual(hgvs_c.format(), s)

    def test_hgvs_p(self):
        s = "NP_001356716.1:p.Ala18Val"
        hgvs_p = self.hgvs_parser.parse(s=s)
        self.assertTrue(isinstance(hgvs_p, HgvsP))
        self.assertEqual(hgvs_p.format(), s)
