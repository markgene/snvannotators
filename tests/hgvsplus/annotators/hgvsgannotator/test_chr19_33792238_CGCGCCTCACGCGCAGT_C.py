"""Test HgvsGAnnotator class with NC_000019.9:g.33792238_33792242del."""

import unittest

from hgvs.easy import parse

from snvannotators.hgvsplus.annotators.hgvsgannotator import HgvsGAnnotator
from snvannotators.hgvsplus.models import HgvsG

from .hgvsgannotatortesttemplate import HgvsGAnnotatorTestTemplate


class HgvsGAnnotatorChr19_33792238_CGCGCCTCACGCGCAGT_C_TestCase(HgvsGAnnotatorTestTemplate, unittest.TestCase):
    """Test HgvsGAnnotator class with NC_000019.9:g.33792238_33792242del."""

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        sequence_variant_g = parse("NC_000019.9:g.33792239_33792254del")
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
        cls.hgvs_g_formatted = "NC_000019.9:g.33792239_33792254del"
        cls.hgvs_g_normalized_formatted = "NC_000019.9:g.33792239_33792254del"
        cls.tx_ac = "NM_004364.4"
        cls.hgvs_t_formmated = "NM_004364.4(CEBPA):c.1067_*5del"
        cls.hgvs_p_formatted = "NP_004355.2(CEBPA):p.Asn356SerfsTer61"
