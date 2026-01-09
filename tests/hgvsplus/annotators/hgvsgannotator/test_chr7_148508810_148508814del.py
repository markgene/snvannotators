"""Test HgvsGAnnotator class with NC_000007.13:g.148508810_148508814del."""

import unittest

from hgvs.easy import parse

from snvannotators.hgvsplus.annotators.hgvsgannotator import HgvsGAnnotator
from snvannotators.hgvsplus.models import HgvsG

from .hgvsgannotatortesttemplate import HgvsGAnnotatorTestTemplate


class HgvsGAnnotatorChr7_148508809_GATGCT_G_TestCase(HgvsGAnnotatorTestTemplate, unittest.TestCase):
    """Test HgvsGAnnotator class with NC_000007.13:g.148508810_148508814del."""

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        sequence_variant_g = parse("NC_000007.13:g.148508810_148508814del")
        hgvs_g = HgvsG.from_sequence_variant_g(sequence_variant_g=sequence_variant_g)
        alt_aln_method = "splign"
        tss_upstream_limit = 20000
        uncertain = False
        cls.hgvs_annotation = HgvsGAnnotator(
            hgvs_g=hgvs_g,
            alt_aln_method=alt_aln_method,
            tss_upstream_limit=tss_upstream_limit,
            uncertain=uncertain,
            verbose=True,
        ).annotate()
        cls.hgvs_g_formatted = "NC_000007.13:g.148508810_148508814del"
        cls.hgvs_g_normalized_formatted = "NC_000007.13:g.148508811_148508815del"
        cls.tx_ac = "NM_004456.4"
        cls.hgvs_t_formmated = "NM_004456.4(EZH2):c.1852-3_1853del"
        cls.hgvs_p_formatted = "NP_004447.2:p.?"
