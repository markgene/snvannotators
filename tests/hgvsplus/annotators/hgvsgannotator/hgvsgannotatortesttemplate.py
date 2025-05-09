"""Test template for HgvsGAnnotator class."""

# from hgvs.easy import parse

from snvannotators.hgvsplus.annotators.hgvsannotation import HgvsAnnotation
# from snvannotators.hgvsplus.annotators.hgvsgannotator import HgvsGAnnotator
from snvannotators.hgvsplus.models import HgvsG, HgvsP, HgvsT


class HgvsGAnnotatorTestTemplate:
    """Test template for HgvsGAnnotator class."""

    # @classmethod
    # def setUpClass(cls) -> None:
    #     super().setUpClass()
    #     sequence_variant_g = parse("NC_000005.9:g.1295228G>A")
    #     hgvs_g = HgvsG.from_sequence_variant_g(sequence_variant_g=sequence_variant_g)
    #     alt_aln_method = "splign"
    #     tss_upstream_limit = 20000
    #     uncertain = False
    #     cls.hgvs_annotation = HgvsGAnnotator(
    #         hgvs_g=hgvs_g,
    #         alt_aln_method=alt_aln_method,
    #         tss_upstream_limit=tss_upstream_limit,
    #         uncertain=uncertain,
    #     ).annotate()
    #     cls.hgvs_g_formatted = "NC_000005.9:g.1295228G>A"
    #     cls.hgvs_g_normalized_formatted = "NC_000005.9:g.1295228G>A"
    #     cls.tx_ac = "NM_198253.2"
    #     cls.hgvs_t_formmated = "NM_198253.2(TERT):c.-124C>T"
    #     cls.hgvs_p_formatted = "NP_937983.2(TERT):p.?"

    def test_hgvs_annotation(self):
        self.assertTrue(isinstance(self.hgvs_annotation, HgvsAnnotation))

    def test_hgvs_annotation_hgvs_g(self):
        self.assertTrue(isinstance(self.hgvs_annotation.hgvs_g, HgvsG))
        self.assertEqual(str(self.hgvs_annotation.hgvs_g), self.hgvs_g_formatted)

    def test_hgvs_annotation_hgvs_g_normalized(self):
        self.assertTrue(isinstance(self.hgvs_annotation.hgvs_g_normalized, HgvsG))
        self.assertEqual(
            str(self.hgvs_annotation.hgvs_g_normalized),
            self.hgvs_g_normalized_formatted,
        )

    def test_hgvs_annotation_hgvs_t(self):
        for hgvs_tp_annotation in self.hgvs_annotation.hgvs_tp_annotations:
            self.assertTrue(
                hgvs_tp_annotation.hgvs_t is None
                or isinstance(hgvs_tp_annotation.hgvs_t, HgvsT)
            )
            if hgvs_tp_annotation.tx_ac == self.tx_ac:
                self.assertEqual(str(hgvs_tp_annotation.hgvs_t), self.hgvs_t_formmated)
                return
        raise ValueError(f"transcript {self.tx_ac} is not found")

    def test_hgvs_annotation_hgvs_p(self):
        for hgvs_tp_annotation in self.hgvs_annotation.hgvs_tp_annotations:
            self.assertTrue(
                hgvs_tp_annotation.hgvs_p is None
                or isinstance(hgvs_tp_annotation.hgvs_p, HgvsP)
            )
            if hgvs_tp_annotation.tx_ac == self.tx_ac:
                self.assertEqual(str(hgvs_tp_annotation.hgvs_p), self.hgvs_p_formatted)
                return
        raise ValueError(f"transcript {self.tx_ac} is not found")
