"""Knowledgebase item."""

from dataclasses import dataclass
from typing import Any

from .knowledgebase import Knowledgebase


@dataclass
class KnowledgebaseItem:
    knowledgebase: Knowledgebase
    raw: Any

    def __post_init__(self):
        if not isinstance(self.knowledgebase, Knowledgebase):
            raise ValueError(
                f"knowledgebase {self.knowledgebase} must be a Knowledgebase object"
            )
