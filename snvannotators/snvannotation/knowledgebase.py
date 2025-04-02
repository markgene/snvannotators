"""Knowledgebase."""

from dataclasses import dataclass
from typing import ClassVar, List


@dataclass
class Knowledgebase:
    name: str
    kb_type: str
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
        if not isinstance(self.kb_type, str):
            raise ValueError(f"kb_type {self.kb_type} must be a str")
        if not isinstance(self.revision, str):
            raise ValueError(f"revision {self.revision} must be a str")
        if self.kb_type.lower() not in self.ALLOWED_KNOWLEDGEBASE_TYPES:
            s = ", ".join(self.ALLOWED_KNOWLEDGEBASE_TYPES)
            raise ValueError(f"kb_type {self.kb_type} must be one of {s}")
