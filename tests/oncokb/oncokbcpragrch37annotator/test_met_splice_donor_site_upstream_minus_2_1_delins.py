"""Test OncokbCpraGrch37Annotator class with MET splice donor site mutation.

Deletion-insertion at position -2_-1 relative to intron 14, which is inside
exon 14 but right next to splice donor site.
"""

import unittest

from pyoncokb.oncokbapi import OncokbApi
from snvmodels.cpra import Cpra

from snvannotators.oncokb.oncokbcpragrch37annotator import OncokbCpraGrch37Annotator
from tests.testconfig import TestConfig

from .oncokbcpragrch37annotatortesttemplate import OncokbCpraGrch37AnnotatorTestTemplate

config = TestConfig()


class OncokbCpraGrch37AnnotatorMetSpliceDonorSiteUpstreamMinus2toMinus1DelinsTestCase(
    OncokbCpraGrch37AnnotatorTestTemplate, unittest.TestCase
):
    """Test OncokbCpraGrch37Annotator class with MET splice donor site mutation.

    Deletion-insertion at position -2_-1 relative to intron 14, which is inside
    exon 14 but right next to splice donor site.
    """

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        config = TestConfig()
        oncokb_auth = config.get_oncokb_authorization()
        oncokb_api = OncokbApi(auth=oncokb_auth)
        oncokb_cpra_annotator = OncokbCpraGrch37Annotator(oncokb_api=oncokb_api)
        cpra = Cpra("chr7", 116412042, "AG", "TT")
        cls.indicator_query_resp = oncokb_cpra_annotator.annotate(cpra=cpra)
        cls.query_alteration = "E1009_D1010delinsDY"
        cls.query_entrez_gene_id = 4233
        cls.query_hugo_symbol = "MET"
        cls.query_tumor_type = None
        cls.highest_diagnostic_implication_level = None
        cls.highest_fda_level = "LEVEL_Fda2"
        cls.highest_prognostic_implication_level = None
        cls.highest_resistance_level = None
        cls.highest_sensitive_level = "LEVEL_1"
        cls.mutation_effect_known_effect = "Likely Gain-of-function"
        cls.oncogenic = "Likely Oncogenic"
        cls.allele_exist = True
        cls.variant_exist = False
        cls.hotspot = False
        cls.vus = False
        cls.tumor_type_summary = ""
