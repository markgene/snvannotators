"""KnowledgebaseItem."""

import unittest

from snvannotators.snvannotation.knowledgebase import Knowledgebase
from snvannotators.snvannotation.knowledgebaseitem import KnowledgebaseItem


class KnowledgebaseItemTestCase(unittest.TestCase):

    def test_valid_example_1(self):
        knowledgebase = Knowledgebase(
            name="Example 1", knowledgebase_type="variant", revision="v0.1"
        )
        knowledgebase_item = KnowledgebaseItem(knowledgebase=knowledgebase, raw=1)
        self.assertTrue(isinstance(knowledgebase_item, KnowledgebaseItem))
        self.assertTrue(knowledgebase_item.is_valid())

    def test_invalid_example_1(self):
        knowledgebase = Knowledgebase(
            name="Example 1", knowledgebase_type="invalid type", revision="v0.1"
        )
        knowledgebase_item = KnowledgebaseItem(knowledgebase=knowledgebase, raw=1)
        self.assertTrue(isinstance(knowledgebase_item, KnowledgebaseItem))
        self.assertFalse(knowledgebase_item.is_valid())
        with self.assertRaises(Exception) as context:
            knowledgebase_item.validate()
        self.assertTrue("knowledgebase_type invalid type must be" in str(context.exception))
