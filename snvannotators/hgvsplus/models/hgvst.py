"""Extend SequenceVariant class of c or n type."""

import logging

from hgvs.sequencevariant import SequenceVariant

logger = logging.getLogger(__name__)


class HgvsT(SequenceVariant):
    """Extend SequenceVariant class of c or n type."""

    EDIT_TYPE_LOOKUP = {
        "substitution": "sub",
        "deletion": "del",
        "duplication": "dup",
        "insertion": "ins",
        "inversion": "inv",
        "deletion-insertion": "delins",
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.type in ["c", "n"]

    @classmethod
    def from_sequence_variant_t(cls, sequence_variant_t: SequenceVariant):
        assert isinstance(sequence_variant_t, SequenceVariant)
        assert sequence_variant_t.type in ["c", "n"]
        return cls(
            ac=sequence_variant_t.ac,
            type=sequence_variant_t.type,
            posedit=sequence_variant_t.posedit,
        )

    def to_sequence_variant_t(self) -> SequenceVariant:
        sequence_variant_t = SequenceVariant(
            ac=self.ac, type=self.type, posedit=self.posedit
        )
        return sequence_variant_t

    def is_substitution(self) -> bool:
        """Is substitution?"""
        return self.posedit.edit.type == self.EDIT_TYPE_LOOKUP["substitution"]

    def is_deletion(self) -> bool:
        """Is deletion?"""
        return self.posedit.edit.type == self.EDIT_TYPE_LOOKUP["deletion"]

    def is_duplication(self) -> bool:
        """Is duplication?"""
        return self.posedit.edit.type == self.EDIT_TYPE_LOOKUP["duplication"]

    def is_insertion(self) -> bool:
        """Is insertion?"""
        return self.posedit.edit.type == self.EDIT_TYPE_LOOKUP["insertion"]

    def is_inversion(self) -> bool:
        """Is inversion?"""
        return self.posedit.edit.type == self.EDIT_TYPE_LOOKUP["inversion"]

    def is_deletion_insertion(self) -> bool:
        """Is deletion-insertion?"""
        return self.posedit.edit.type == self.EDIT_TYPE_LOOKUP["deletion-insertion"]
