"""Knowledgebase."""

import unittest

from snvannotators.snvannotation.knowledgebase import Knowledgebase


class KnowledgebaseTestCase(unittest.TestCase):

    def test_valid_example_1(self):
        knowledgebase = Knowledgebase(
            name="Example 1", knowledgebase_type="variant", revision="v0.1"
        )
        self.assertTrue(isinstance(knowledgebase, Knowledgebase))
        self.assertTrue(knowledgebase.is_valid())
        

    def test_invalid_example_1(self):
        knowledgebase = Knowledgebase(
            name="Example 1", knowledgebase_type="invalid type", revision="v0.1"
        )
        self.assertTrue(isinstance(knowledgebase, Knowledgebase))
        self.assertFalse(knowledgebase.is_valid())
        with self.assertRaises(Exception) as context:
            knowledgebase.validate()
        self.assertTrue("knowledgebase_type invalid type must be" in str(context.exception))
