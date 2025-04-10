"""Test OncokbApi class."""

import unittest

from pyoncokb.oncokbapi import OncokbApi
from pyoncokb.models.indicatorqueryresp import IndicatorQueryResp
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

    def test_met_splice_donor_site_upstream_v1001_f1007_deletion(self):
        cpra = Cpra("chr7", 116412015, "TGTAGACTACCGAGCTACTTTT", "T")
        indicator_query_resp = self.oncokb_cpra_grch37_annotator.annotate(cpra=cpra)
        self.assertTrue(isinstance(indicator_query_resp, IndicatorQueryResp))
        self.assertIsNotNone(indicator_query_resp.query)
        self.assertEqual(indicator_query_resp.query.alteration, "V1001_F1007del")
        self.assertIsNotNone(indicator_query_resp.oncogenic)
        self.assertEqual(indicator_query_resp.oncogenic, "Likely Oncogenic")
