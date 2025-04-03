"""Knowledgebase."""

import logging

from dataclasses import dataclass
from typing import ClassVar, List

logger = logging.getLogger(__name__)


@dataclass
class Knowledgebase:
    name: str
    knowledgebase_type: str
    revision: str

    ALLOWED_KNOWLEDGEBASE_TYPES: ClassVar[List] = [
        "classification",
        "interpretation",
        "filter",
        "publication",
        "interp_assoc",
        "clinical_trial",
        "guideline",
        "clinical_evidence",
        "human_research_evidence",
        "variant",
    ]

    def __post_init__(self):
        if not isinstance(self.name, str):
            raise ValueError(f"name {self.name} must be a str")
        if not isinstance(self.knowledgebase_type, str):
            raise ValueError(
                f"knowledgebase_type {self.knowledgebase_type} must be a str"
            )
        if not isinstance(self.revision, str):
            raise ValueError(f"revision {self.revision} must be a str")

    def validate(self):
        if self.knowledgebase_type.lower() not in self.ALLOWED_KNOWLEDGEBASE_TYPES:
            s = ", ".join(self.ALLOWED_KNOWLEDGEBASE_TYPES)
            raise ValueError(
                f"knowledgebase_type {self.knowledgebase_type} must be one of {s}"
            )

    def is_valid(self) -> bool:
        if self.knowledgebase_type.lower() not in self.ALLOWED_KNOWLEDGEBASE_TYPES:
            s = ", ".join(self.ALLOWED_KNOWLEDGEBASE_TYPES)
            logger.warning(
                "knowledgebase_type %s must be one of %s", self.knowledgebase_type, s
            )
            return False
        return True
