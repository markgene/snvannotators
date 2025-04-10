"""Test OncokbCpraAnnotator class with BRAF V600E."""

import unittest

from pyoncokb.oncokbapi import OncokbApi
from snvmodels.cpra import Cpra

from snvannotators.oncokb.oncokbcpragrch37annotator import OncokbCpraGrch37Annotator
from tests.testconfig import TestConfig

from .oncokbcpragrch37annotatortesttemplate import OncokbCpraGrch37AnnotatorTestTemplate

config = TestConfig()


class OncokbCpraGrch37AnnotatorBrafV600eTestCase(
    OncokbCpraGrch37AnnotatorTestTemplate, unittest.TestCase
):
    """Test OncokbCpraGrch37Annotator class with BRAF V600E."""

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        config = TestConfig()
        oncokb_auth = config.get_oncokb_authorization()
        oncokb_api = OncokbApi(auth=oncokb_auth)
        oncokb_cpra_annotator = OncokbCpraGrch37Annotator(oncokb_api=oncokb_api)
        cpra = Cpra(chrom="chr7", pos=140453136, ref="A", alt="T")
        cls.indicator_query_resp = oncokb_cpra_annotator.annotate(cpra=cpra)
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
        cls.allele_exist = True
        cls.variant_exist = True
        cls.hotspot = True
        cls.vus = False
        cls.tumor_type_summary = ""
