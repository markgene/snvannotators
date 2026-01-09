"""Convert HgvsC to HgvsP.

Sometimes, the c_to_p method does not work as expected. For example, given
"NM_005343.2:c.2T>G", it will return a SequenceVariant object with its
posedit field of AARefAlt class, where a PosEdit class is expected.
"""

import logging

from hgvs.assemblymapper import AssemblyMapper
from hgvs.easy import am37, parse
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

    def __init__(
        self,
        hgvs_c: HgvsC,
        assembly_mapper: AssemblyMapper = am37,
        error_ok: bool = False,
        verbose: bool = False,
    ):
        self.hgvs_c = hgvs_c
        self.assembly_mapper = assembly_mapper
        self.error_ok = error_ok
        self.verbose = verbose
        self.__post_init__()

        # internal use only
        self.hgvs_p: HgvsP | None = None

    def __post_init__(self):
        if not isinstance(self.hgvs_c, HgvsC):
            raise ValueError(f"hgvs_c {self.hgvs_c} must be a HgvsC")
        if self.assembly_mapper is None:
            raise ValueError("assembly_mapper must be not None")
        if not isinstance(self.assembly_mapper, AssemblyMapper):
            raise ValueError("assembly_mapper must be an AssemblyMapper object")

    def map(self) -> HgvsP:
        """Map SequenceVariant C to P type.

        Sometimes, the c_to_p method does not work as expected. For example, given
        "NM_005343.2:c.2T>G", it will return a SequenceVariant object with its
        posedit field of AARefAlt class, where a PosEdit class is expected.
        """
        try:
            sequence_variant_p = self.get_sequence_variant_p()
        except HGVSInvalidIntervalError:
            self.hgvs_p = HgvsPUnknownPosEditTranscriptAccessionCreator(
                self.hgvs_c.ac
            ).create()
        else:
            self.hgvs_p = HgvsP.from_sequence_variant_p(sequence_variant_p)
            if not isinstance(self.hgvs_p.posedit, PosEdit):
                logger.warning(
                    "posedit field is not a PosEdit object: %s. Trying to map.",
                    repr(self.hgvs_p.posedit),
                )
                hgvs_p_fixed = self.fix(sequence_variant_p=sequence_variant_p)
                if hgvs_p_fixed.posedit is not None and not isinstance(
                    hgvs_p_fixed.posedit, PosEdit
                ):
                    raise RuntimeError(f"fail to fix hgvs_p {self.hgvs_p}")
                else:
                    logger.warning(
                        "fix %s to %s", repr(self.hgvs_p), repr(hgvs_p_fixed)
                    )
                    self.hgvs_p = hgvs_p_fixed
        self.update_hgvs_p_gene()
        return self.hgvs_p

    def update_hgvs_p_gene(self):
        """Update gene field of hgvs_p from hgvs_c."""
        if self.hgvs_p is None:
            if self.hgvs_c.gene is not None:
                if self.verbose:
                    logger.debug(
                        "hgvs_p is None, setting gene field from hgvs_c: %s",
                        self.hgvs_c.gene,
                    )
                    self.hgvs_p.gene = self.hgvs_c.gene
            else:
                logger.warning(
                    "hgvs_p gene field not updated because hgvs_c gene field is None"
                )
        else:
            if self.hgvs_c.gene is not None:
                if self.hgvs_p.gene != self.hgvs_c.gene:
                    logger.warning(
                        "updating hgvs_p gene field from %s to %s",
                        self.hgvs_p.gene,
                        self.hgvs_c.gene,
                    )

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
