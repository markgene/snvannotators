"""Extend SequenceVariant class of c or n type."""

import logging

from hgvs.easy import validate
from hgvs.exceptions import HGVSInvalidIntervalError
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

    def __init__(self, soft_validation: bool = True, *args, **kwargs):
        """init.

        :param soft_validation: only raise errors when they are not recognized, 
            defaults to True. If False, raise errors no matter they are recognized
            or not.
        :type soft_validation: bool, optional
        """
        super().__init__(*args, **kwargs)
        self.soft_validation = soft_validation
        assert self.type in ["c", "n"]
        self.is_valid()

    @classmethod
    def from_sequence_variant_t(
        cls, sequence_variant_t: SequenceVariant, soft_validation: bool = True
    ):
        assert isinstance(sequence_variant_t, SequenceVariant)
        assert sequence_variant_t.type in ["c", "n"]
        return cls(
            soft_validation=soft_validation,
            ac=sequence_variant_t.ac,
            type=sequence_variant_t.type,
            posedit=sequence_variant_t.posedit,
            gene=sequence_variant_t.gene,
        )

    def to_sequence_variant_t(self) -> SequenceVariant:
        sequence_variant_t = SequenceVariant(
            ac=self.ac, type=self.type, posedit=self.posedit, gene=self.gene
        )
        return sequence_variant_t

    def is_valid(self) -> bool:
        sequence_variant_t = self.to_sequence_variant_t()
        try:
            is_valid = validate(sequence_variant_t)
        except HGVSInvalidIntervalError as err:
            if "coordinate is out of bounds" in str(err):
                if self.soft_validation:
                    logger.warning(
                        "%s. The error is usually seen for promoter variant, e.g. c.-124"
                    )
                    is_valid = True
                else:
                    raise
            else:
                raise
        return is_valid

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
