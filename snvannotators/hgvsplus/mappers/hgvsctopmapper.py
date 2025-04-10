"""Convert HgvsC to HgvsP.

Sometimes, the c_to_p method does not work as expected. For example, given
"NM_005343.2:c.2T>G", it will return a SequenceVariant object with its
posedit field of AARefAlt class, where a PosEdit class is expected.
"""

import logging
from typing import Optional

from hgvs.assemblymapper import AssemblyMapper
from hgvs.easy import hgvs_data_provider, parse
from hgvs.exceptions import HGVSInvalidIntervalError
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant

from snvannotators.hgvsplus.models import HgvsC, HgvsP
from snvannotators.hgvsplus.creators.hgvspunknownposedittranscriptaccessioncreator import (
    HgvsPUnknownPosEditTranscriptAccessionCreator,
)

logger = logging.getLogger(__name__)


class HgvsCToPMapper:
    """Convert SequenceVariant C to P type.

    Sometimes, the c_to_p method does not work as expected. For example, given
    "NM_005343.2:c.2T>G", it will return a SequenceVariant object with its
    posedit field of AARefAlt class, where a PosEdit class is expected.
    """

    def __init__(self, hgvs_c: HgvsC, assembly_name: str = "GRCh37"):
        self.hgvs_c = hgvs_c
        self.assembly_name = assembly_name
        self.assembly_mapper: Optional[AssemblyMapper] = None
        self.init_assembly_mapper()
        self.__post_init__()

    def __post_init__(self):
        if not isinstance(self.hgvs_c, HgvsC):
            raise ValueError(f"hgvs_c {self.hgvs_c} must be a HgvsC")
        if self.assembly_mapper is None:
            raise ValueError("assembly_mapper must be not None")
        if not isinstance(self.assembly_mapper, AssemblyMapper):
            raise ValueError("assembly_mapper must be an AssemblyMapper object")

    def init_assembly_mapper(self):
        assert self.assembly_name in ["hg19", "hg38", "GRCh37", "GRCh38"]
        if self.assembly_name == "hg19":
            assembly_name = "GRCh37"
        elif self.assembly_name == "hg38":
            assembly_name = "GRCh38"
        else:
            assembly_name = self.assembly_name
        self.assembly_mapper = AssemblyMapper(
            hgvs_data_provider, assembly_name=assembly_name
        )

    def map(self) -> HgvsP:
        """Map SequenceVariant C to P type.

        Sometimes, the c_to_p method does not work as expected. For example, given
        "NM_005343.2:c.2T>G", it will return a SequenceVariant object with its
        posedit field of AARefAlt class, where a PosEdit class is expected.
        """
        try:
            sequence_variant_p = self.get_sequence_variant_p()
        except HGVSInvalidIntervalError:
            hgvs_p = HgvsPUnknownPosEditTranscriptAccessionCreator(
                self.hgvs_c.ac
            ).create()
            return hgvs_p
        else:
            hgvs_p = HgvsP.from_sequence_variant_p(sequence_variant_p)
            if not isinstance(hgvs_p.posedit, PosEdit):
                logger.warning(
                    "posedit field is not a PosEdit object: %s. Trying to map.",
                    repr(hgvs_p.posedit),
                )
                hgvs_p_fixed = self.fix(sequence_variant_p=sequence_variant_p)
                if hgvs_p_fixed.posedit is not None and not isinstance(
                    hgvs_p_fixed.posedit, PosEdit
                ):
                    raise RuntimeError(f"fail to fix hgvs_p {hgvs_p}")
                else:
                    logger.warning("fix %s to %s", repr(hgvs_p), repr(hgvs_p_fixed))
                    return hgvs_p_fixed
            return hgvs_p

    def get_sequence_variant_p(self) -> SequenceVariant:
        try:
            sequence_variant_p = self.assembly_mapper.c_to_p(self.hgvs_c)
        except Exception as err:
            logger.info("fail to convert %s with GRCh37: %s", self.hgvs_c, err)
            sequence_variant_p = self.assembly_mapper.c_to_p(self.hgvs_c)
        return sequence_variant_p

    @staticmethod
    def fix(sequence_variant_p: SequenceVariant) -> HgvsP:
        sequence_variant_p_str = sequence_variant_p.format()
        sequence_variant_p_fixed = parse(sequence_variant_p_str)
        hgvs_p_fixed = HgvsP.from_sequence_variant_p(
            sequence_variant_p=sequence_variant_p_fixed
        )
        return hgvs_p_fixed
