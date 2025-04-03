"""Knowledgebase."""

import unittest

from snvannotators.snvannotation.knowledgebase import Knowledgebase

# from snvannotators.snvannotation.knowledgebaseitem import KnowledgebaseItem


class KnowledgebaseTestCase(unittest.TestCase):

    def test_example_1(self):
        knowledgebase = Knowledgebase(
            name="Example 1", knowledgebase_type="variant", revision="v0.1"
        )
        self.assertTrue(isinstance(knowledgebase, Knowledgebase))

    def test_invalid_example_1(self):
        knowledgebase = Knowledgebase(
            name="Example 1", knowledgebase_type="invalid type", revision="v0.1"
        )
        self.assertFalse(knowledgebase.is_valid())
