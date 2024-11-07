"""Test OncokbApi class."""

import unittest

from pyoncokb.oncokbapi import OncokbApi
from snvmodels.cpra import Cpra

from snvannotators.oncokb.oncokbcpragrch37annotator import OncokbCpraGrch37Annotator
from tests.testconfig import TestConfig

config = TestConfig()


class OncokbCpraGrch37AnnotatorTestCase(unittest.TestCase):
    """Test OncokbCpraGrch37Annotator class."""

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        config = TestConfig()
        oncokb_auth = config.get_oncokb_authorization()
        oncokb_api = OncokbApi(auth=oncokb_auth)
        cls.oncokb_cpra_grch37_annotator = OncokbCpraGrch37Annotator(
            oncokb_api=oncokb_api
        )

    def test_annotate_braf_v600e(self):
        """Test annotate method."""
        cpra = Cpra(chrom="chr7", pos=140453136, ref="A", alt="T")
        braf_v600e = self.oncokb_cpra_grch37_annotator.annotate(cpra=cpra)
        self.assertTrue(braf_v600e.allele_exist)
        self.assertEqual(braf_v600e.query.alteration, "V600E")
        self.assertEqual(braf_v600e.query.entrez_gene_id, 673)
        self.assertTrue(braf_v600e.gene_exist)
        self.assertEqual(braf_v600e.query.hugo_symbol, "BRAF")
        self.assertEqual(braf_v600e.highest_diagnostic_implication_level, "LEVEL_Dx2")
        self.assertEqual(braf_v600e.highest_fda_level, "LEVEL_Fda2")
        self.assertIsNone(braf_v600e.highest_prognostic_implication_level)
        self.assertIsNone(braf_v600e.highest_resistance_level)
        self.assertEqual(braf_v600e.highest_sensitive_level, "LEVEL_1")
        self.assertTrue(braf_v600e.hotspot)
        self.assertEqual(braf_v600e.mutation_effect.known_effect, "Gain-of-function")
        self.assertEqual(braf_v600e.oncogenic, "Oncogenic")
        self.assertIsNone(braf_v600e.query.tumor_type)
        self.assertEqual(braf_v600e.tumor_type_summary, "")
        self.assertTrue(braf_v600e.variant_exist)
        self.assertFalse(braf_v600e.vus)
