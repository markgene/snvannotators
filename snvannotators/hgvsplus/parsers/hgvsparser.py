"""Parse."""

from hgvs.parser import Parser
from hgvs.sequencevariant import SequenceVariant

from snvannotators.hgvsplus.models import HgvsC, HgvsG, HgvsP, HgvsT


class HgvsParser:
    parser = Parser()

    instance = None

    def __new__(cls, *args, **kwargs):
        if cls.instance is None:
            cls.instance = super().__new__(cls)
        return cls.instance

    def parse(self, s: str, c_as_t: bool=True):
        sequence_variant = self.parse_sequence_variant(s=s)
        if sequence_variant.type == "g":
            return HgvsG.from_sequence_variant_g(sequence_variant_g=sequence_variant)
        if sequence_variant.type == "c":
            if c_as_t:
                return HgvsT.from_sequence_variant_t(sequence_variant_t=sequence_variant)
            else:
                return HgvsC.from_sequence_variant_c(sequence_variant_c=sequence_variant)
        if sequence_variant.type == "n":
            return HgvsT.from_sequence_variant_t(sequence_variant_t=sequence_variant)
        if sequence_variant.type == "p":
            return HgvsP.from_sequence_variant_p(sequence_variant_p=sequence_variant)
        raise ValueError(
            "fail to recognize the type {sequence_variant.type}. It must be one of g, c, n, p."
        )

    def parse_sequence_variant(self, s: str) -> SequenceVariant:
        return self.parser.parse(s)
