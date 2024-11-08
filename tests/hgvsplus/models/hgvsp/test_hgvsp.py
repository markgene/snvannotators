"""Test HgvsP class."""

import unittest

from hgvs.easy import parse

from snvannotators.hgvsplus.models import HgvsP


class HgvsPTestCase(unittest.TestCase):
    """Test HgvsP class."""

    def test_missense(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_003997.1:p.(Trp24Cys)"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            HgvsP.MUTATION_TERM_LOOKUP["missense"],
        )

    def test_nonsense(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_004439.2:p.(Leu755Ter)"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            HgvsP.MUTATION_TERM_LOOKUP["nonsense"],
        )

    def test_synonymous(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_003997.1:p.Cys188="))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            HgvsP.MUTATION_TERM_LOOKUP["synonymous"],
        )

    def test_no_protein(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("LRG_199p1:p.0"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            HgvsP.MUTATION_TERM_LOOKUP["translation initiation codon: no protein"],
        )

    def test_tic_unknown(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("LRG_199p1:p.(Met1?)"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            HgvsP.MUTATION_TERM_LOOKUP["translation initiation codon: unknown"],
        )

    def test_inframe_deletion_one_amino_acid(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_003997.2:p.Val7del"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(), "in-frame deletion"
        )

    def test_inframe_duplication_several_amino_acids(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_003997.2:p.Lys23_Val25dup"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(), "in-frame duplication"
        )

    def test_inframe_stop_gain_insertion(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_004371.2:p.(Pro46_Asn47insSerSerTer)"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            "in-frame stop gain insertion",
        )

    def test_inframe_stop_gain_deletion_insertion(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_004371.2:p.(Asn47delinsSerSerTer)"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            "in-frame stop gain deletion-insertion",
        )

    def test_frameshift(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_000035.2:p.E124Gfs*34"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            HgvsP.MUTATION_TERM_LOOKUP["frameshift"],
        )

    def test_n_terminal_extension(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_003997.2:p.Met1ext-5"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            HgvsP.MUTATION_TERM_LOOKUP["N-terminal extension"],
        )

    def test_c_terminal_extension(self):
        hgvs_p = HgvsP.from_sequence_variant_p(parse("NP_003997.2:p.Ter110GlnextTer17"))
        self.assertEqual(
            hgvs_p.get_mutation_type_of_protein_impact(),
            HgvsP.MUTATION_TERM_LOOKUP["C-terminal extension"],
        )
