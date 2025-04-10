"""Test OncokbCpraAnnotator class with MET c.3028+1del."""

import unittest

from pyoncokb.oncokbapi import OncokbApi
from snvmodels.cpra import Cpra

from snvannotators.oncokb.oncokbcpragrch37annotator import OncokbCpraGrch37Annotator
from tests.testconfig import TestConfig

from .oncokbcpragrch37annotatortesttemplate import OncokbCpraGrch37AnnotatorTestTemplate

config = TestConfig()


class OncokbCpraGrch37AnnotatorMetSpliceDonorSiteC3028Plus1DelTestCase(
    OncokbCpraGrch37AnnotatorTestTemplate, unittest.TestCase
):
    """Test OncokbCpraGrch37Annotator class with MET c.3028+1del."""

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        config = TestConfig()
        oncokb_auth = config.get_oncokb_authorization()
        oncokb_api = OncokbApi(auth=oncokb_auth)
        oncokb_cpra_annotator = OncokbCpraGrch37Annotator(oncokb_api=oncokb_api)
        cpra = Cpra("chr7", 116412042, "AG", "A")
        cls.indicator_query_resp = oncokb_cpra_annotator.annotate(cpra=cpra)
        cls.query_alteration = "D1010fs"
        cls.query_entrez_gene_id = 4233
        cls.query_hugo_symbol = "MET"
        cls.query_tumor_type = None
        cls.highest_diagnostic_implication_level = None
        cls.highest_fda_level = None
        cls.highest_prognostic_implication_level = None
        cls.highest_resistance_level = None
        cls.highest_sensitive_level = None
        cls.mutation_effect_known_effect = "Unknown"
        cls.oncogenic = "Unknown"
        cls.allele_exist = False
        cls.variant_exist = True
        cls.hotspot = False
        cls.vus = True
        cls.tumor_type_summary = ""
