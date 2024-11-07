"""Test HgvsG class."""

import unittest

from hgvs.easy import parse

from snvannotators.hgvsplus.models import HgvsG


class HgvsGTestCase(unittest.TestCase):
    """Test HgvsG class."""
    
    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
    
    def test_from_sequence_variant_g(self):
        sequence_variant_g = parse("NC_000005.9:g.1295228G>A")
        hgvs_g = HgvsG.from_sequence_variant_g(sequence_variant_g=sequence_variant_g)
        self.assertTrue(isinstance(hgvs_g, HgvsG))
