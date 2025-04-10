"""Test HgvsGAnnotator class with TP53 c.754del."""

import unittest

from hgvs.easy import parse

from snvannotators.hgvsplus.annotators.hgvsgannotator import HgvsGAnnotator
from snvannotators.hgvsplus.models import HgvsG

from .hgvsgannotatortesttemplate import HgvsGAnnotatorTestTemplate


class HgvsGAnnotatorMetSpliceDonorSiteUpstreamV1001F1007DelTestCase(
    HgvsGAnnotatorTestTemplate, unittest.TestCase
):
    """Test HgvsGAnnotator class with in-frame deletion affecting the splice sites flanking MET exon 14 result."""

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        sequence_variant_g = parse("NC_000007.13:g.116412015_116412035del")
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
        cls.hgvs_g_formatted = "NC_000007.13:g.116412015_116412035del"
        cls.hgvs_g_normalized_formatted = "NC_000007.13:g.116412016_116412036del"
        cls.tx_ac = "NM_000245.3"
        cls.hgvs_t_formmated = "NM_000245.3(MET):c.3001_3021del"
        cls.hgvs_p_formatted = "NP_000236.2:p.Val1001_Phe1007del"
