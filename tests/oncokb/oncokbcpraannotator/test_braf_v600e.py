"""Test OncokbApi class."""

import unittest

from pyoncokb.oncokbapi import OncokbApi
from snvmodels.cpra import Cpra

from snvannotators.oncokb.oncokbcpraannotator import OncokbCpraAnnotator
from tests.testconfig import TestConfig

config = TestConfig()


class OncokbCpraAnnotatorBrafV600eTestCase(unittest.TestCase):
    """Test OncokbCpraAnnotator class with BRAF V600E."""

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        config = TestConfig()
        oncokb_auth = config.get_oncokb_authorization()
        oncokb_api = OncokbApi(auth=oncokb_auth)
        oncokb_cpra_annotator = OncokbCpraAnnotator(oncokb_api=oncokb_api)
        ref_genome = "GRCh37"
        cpra = Cpra(chrom="chr7", pos=140453136, ref="A", alt="T")
        cls.indicator_query_resp = oncokb_cpra_annotator.annotate(
            cpra=cpra, ref_genome=ref_genome
        )
        cls.query_alteration = "V600E"
        cls.query_entrez_gene_id = 673
        cls.query_hugo_symbol = "BRAF"
        cls.query_tumor_type = None
        cls.highest_diagnostic_implication_level = "LEVEL_Dx2"
        cls.highest_fda_level = "LEVEL_Fda2"
        cls.highest_prognostic_implication_level = None
        cls.highest_resistance_level = None
        cls.highest_sensitive_level = "LEVEL_1"
        cls.mutation_effect_known_effect = "Gain-of-function"
        cls.oncogenic = "Oncogenic"
        cls.tumor_type_summary = ""

    def test_allele_exist(self):
        self.assertTrue(self.indicator_query_resp.allele_exist)

    def test_query_alteration(self):
        self.assertEqual(
            self.indicator_query_resp.query.alteration, self.query_alteration
        )

    def test_query_entrez_gene_id(self):
        self.assertEqual(
            self.indicator_query_resp.query.entrez_gene_id, self.query_entrez_gene_id
        )

    def test_query_tumor_type(self):
        if self.query_tumor_type is None:
            self.assertIsNone(self.indicator_query_resp.query.tumor_type)
        else:
            self.assertEqual(
                self.indicator_query_resp.query.tumor_type,
                self.query_tumor_type,
            )

    def test_gene_exist(self):
        self.assertTrue(self.indicator_query_resp.gene_exist)

    def test_query_hugo_symbol(self):
        self.assertEqual(
            self.indicator_query_resp.query.hugo_symbol, self.query_hugo_symbol
        )

    def test_highest_diagnostic_implication_level(self):
        if self.highest_diagnostic_implication_level is None:
            self.assertIsNone(
                self.indicator_query_resp.highest_diagnostic_implication_level
            )
        else:
            self.assertEqual(
                self.indicator_query_resp.highest_diagnostic_implication_level,
                self.highest_diagnostic_implication_level,
            )

    def test_highest_fda_level(self):
        if self.highest_fda_level is None:
            self.assertIsNone(self.indicator_query_resp.highest_fda_level)
        else:
            self.assertEqual(
                self.indicator_query_resp.highest_fda_level, self.highest_fda_level
            )

    def test_highest_prognostic_implication_level(self):
        if self.highest_prognostic_implication_level is None:
            self.assertIsNone(
                self.indicator_query_resp.highest_prognostic_implication_level
            )
        else:
            self.assertEqual(
                self.indicator_query_resp.highest_prognostic_implication_level,
                self.highest_prognostic_implication_level,
            )

    def test_highest_resistance_level(self):
        if self.highest_resistance_level is None:
            self.assertIsNone(self.indicator_query_resp.highest_resistance_level)
        else:
            self.assertEqual(
                self.indicator_query_resp.highest_resistance_level,
                self.highest_resistance_level,
            )

    def test_highest_sensitive_level(self):
        if self.highest_sensitive_level is None:
            self.assertIsNone(self.indicator_query_resp.highest_sensitive_level)
        else:
            self.assertEqual(
                self.indicator_query_resp.highest_sensitive_level,
                self.highest_sensitive_level,
            )

    def test_hotspot(self):
        self.assertTrue(self.indicator_query_resp.hotspot)

    def test_mutation_effect_known_effect(self):
        self.assertEqual(
            self.indicator_query_resp.mutation_effect.known_effect,
            self.mutation_effect_known_effect,
        )

    def test_oncogenic(self):
        self.assertEqual(self.indicator_query_resp.oncogenic, self.oncogenic)

    def test_tumor_type_summary(self):
        self.assertEqual(
            self.indicator_query_resp.tumor_type_summary, self.tumor_type_summary
        )

    def test_variant_exist(self):
        self.assertTrue(self.indicator_query_resp.variant_exist)

    def test_vus(self):
        self.assertFalse(self.indicator_query_resp.vus)
