"""Annotate Cpra of GRCh37."""

import logging
from typing import Dict, List

from hgvs.exceptions import HGVSDataNotAvailableError

from pyoncokb.models.indicatorqueryresp import IndicatorQueryResp
from pyoncokb.oncokbapi import OncokbApi

from snvmodels.converters.cpratocspragrch37converter import CpraToCspraGrch37Converter
from snvmodels.cpra.cpra import Cpra
from snvmodels.spra.cspra import Cspra
from transcriptfeatures.annotators.rangeannotators import (
    GenomicRange1BasedTranscriptFeatureAnnotator,
)
from transcriptfeatures.annotators.rangeannotators.transcriptfeaturerangeannotation import (
    TranscriptFeatureRangeAnnotation,
)
from transcriptfeatures.models.genomicpositions import GenomicRange1Based
from transcriptfeatures.transcriptfeaturescreators import TranscriptFeaturesCreator

from snvannotators.hgvsplus.annotators.hgvsannotation import HgvsAnnotation
from snvannotators.hgvsplus.annotators.hgvsgannotator import HgvsGAnnotator
from snvannotators.hgvsplus.creators.sequencevariantgcreator import (
    SequenceVariantGCreator,
)
from snvannotators.hgvsplus.models.hgvsg import HgvsG
from snvannotators.myvariant.annotators.myvariantcpraannotator import (
    MyvariantCpraAnnotator,
)
from snvannotators.myvariant.annotation import MyvariantAnnotation
from snvannotators.oncokb.oncokbcpragrch37annotator import OncokbCpraGrch37Annotator
from snvannotators.snvannotation import SnvAnnotation

cpra_to_cspra_converter = CpraToCspraGrch37Converter()
sequence_variant_g_creator = SequenceVariantGCreator()
logger = logging.getLogger(__name__)


class CpraGrch37Annotator:

    def __init__(
        self,
        cpra: Cpra,
        oncokb_api: OncokbApi,
        alt_aln_method: str = "splign",
        tss_upstream_limit: int = 20000,
        uncertain: bool = False,
        promoter_tss_upstream_offset: int = 1500,
        hgvs_g_to_tp_mapper_error_ok: bool = False,
        transcript_features_creator_error_ok: bool = True,
        error_ok: bool = True,
        verbose: bool = False,
    ):
        self.cpra = cpra
        self.oncokb_api = oncokb_api
        self.alt_aln_method = alt_aln_method
        self.tss_upstream_limit = tss_upstream_limit
        self.uncertain = uncertain
        self.promoter_tss_upstream_offset = promoter_tss_upstream_offset
        self.hgvs_g_to_tp_mapper_error_ok = hgvs_g_to_tp_mapper_error_ok
        self.transcript_features_creator_error_ok = transcript_features_creator_error_ok
        self.error_ok = error_ok
        self.verbose = verbose

        # lazy initialized attributes
        self.cspra = None
        self.genomic_range_1_based = None
        # initialize internal attributes
        self.cspra = self.get_cspra()
        self.genomic_range_1_based = self.get_genomic_range_1_based()
        self.oncokb_cpra_grch37_annotator = OncokbCpraGrch37Annotator(
            oncokb_api=oncokb_api
        )
        self.myvariant_cpra_annotator = MyvariantCpraAnnotator()
        self.indicator_query_resp: IndicatorQueryResp | None = None
        self.myvariant_annotation: MyvariantAnnotation | None = None
        self.hgvs_annotation: HgvsAnnotation | None = None
        self.transcript_feature_range_annotations: (
            List[TranscriptFeatureRangeAnnotation] | None
        ) = None

    def annotate(self) -> SnvAnnotation:
        cspra = self.get_cspra()
        indicator_query_resp = self.get_indicator_query_resp()
        myvariant_annotation = self.get_myvariant_annotation()
        hgvs_annotation = self.get_hgvs_annotation()
        transcript_feature_range_annotations = (
            self.get_transcript_feature_range_annotations()
        )
        meta = self.get_meta()
        snv_annotation = SnvAnnotation(
            snv=cspra,
            hgvs_annotation=hgvs_annotation,
            myvariant_annotation=myvariant_annotation,
            indicator_query_resp=indicator_query_resp,
            transcript_feature_range_annotations=transcript_feature_range_annotations,
            meta=meta,
        )
        return snv_annotation

    def get_indicator_query_resp(self) -> IndicatorQueryResp:
        """Get indicator query response from OncoKB."""
        if self.indicator_query_resp is None:
            self.indicator_query_resp = self.oncokb_cpra_grch37_annotator.annotate(
                cpra=self.cpra
            )
        return self.indicator_query_resp

    def get_myvariant_annotation(self) -> MyvariantAnnotation:
        if self.myvariant_annotation is None:
            self.myvariant_annotation = self.myvariant_cpra_annotator.annotate(
                cpra=self.cpra
            )
        return self.myvariant_annotation

    def get_hgvs_annotation(self) -> HgvsAnnotation:
        if self.hgvs_annotation is None:
            cspra = self.get_cspra()
            sequence_variant_g = sequence_variant_g_creator.create_from_spra(spra=cspra)
            hgvs_g = HgvsG.from_sequence_variant_g(
                sequence_variant_g=sequence_variant_g
            )
            hgvs_g_annotator = HgvsGAnnotator(
                hgvs_g=hgvs_g,
                alt_aln_method=self.alt_aln_method,
                tss_upstream_limit=self.tss_upstream_limit,
                uncertain=self.uncertain,
                hgvs_g_to_tp_mapper_error_ok=self.hgvs_g_to_tp_mapper_error_ok,
                error_ok=self.error_ok,
                verbose=self.verbose,
            )
            self.hgvs_annotation = hgvs_g_annotator.annotate()
        return self.hgvs_annotation

    def get_transcript_feature_range_annotations(
        self,
    ) -> List[TranscriptFeatureRangeAnnotation]:
        """Get transcript feature range annotations."""
        if self.transcript_feature_range_annotations is None:
            hgvs_annotation = self.get_hgvs_annotation()
            self.transcript_feature_range_annotations = (
                self.annotate_transcript_feature(hgvs_annotation=hgvs_annotation)
            )
        return self.transcript_feature_range_annotations

    def annotate_transcript_feature(
        self, hgvs_annotation: HgvsAnnotation
    ) -> List[TranscriptFeatureRangeAnnotation]:
        """Annotate transcript features for each transcript in hgvs_annotation."""
        genomic_range_1_based = self.get_genomic_range_1_based()
        transcript_feature_range_annotations = []
        for hgvs_tp_annotation in hgvs_annotation.hgvs_tp_annotations:
            try:
                transcript_features_creator = TranscriptFeaturesCreator(
                    tx_ac=hgvs_tp_annotation.tx_ac,
                    genome_ac=genomic_range_1_based.ac,
                    alt_aln_method=self.alt_aln_method,
                    promoter_tss_upstream_offset=self.promoter_tss_upstream_offset,
                )
                transcript_features = transcript_features_creator.create()
            except HGVSDataNotAvailableError as e:
                if self.transcript_features_creator_error_ok:
                    logger.error(
                        f"Data not available for transcript "
                        f"{hgvs_tp_annotation.tx_ac}: {e}"
                    )
                    transcript_features = None
                else:
                    raise e
            except Exception as e:
                if self.transcript_features_creator_error_ok:
                    logger.error(
                        f"Error creating transcript features for transcript "
                        f"{hgvs_tp_annotation.tx_ac}: {e}"
                    )
                    transcript_features = None
                else:
                    raise e

            if transcript_features is None:
                continue

            range_annotation = GenomicRange1BasedTranscriptFeatureAnnotator(
                genomic_range_1_based=genomic_range_1_based,
                transcript_features=transcript_features,
            ).annotate()
            transcript_feature_range_annotations.append(range_annotation)
        self.transcript_feature_range_annotations = transcript_feature_range_annotations
        return self.transcript_feature_range_annotations

    def get_genomic_range_1_based(self) -> GenomicRange1Based:
        if self.genomic_range_1_based is None:
            cspra = self.get_cspra()
            self.genomic_range_1_based = GenomicRange1Based(
                ac=cspra.ac, start=cspra.get_start_pos(), end=cspra.get_end_pos()
            )
        return self.genomic_range_1_based

    def get_cspra(self) -> Cspra:
        if self.cspra is None:
            self.cspra = cpra_to_cspra_converter.convert(cpra=self.cpra)
        return self.cspra

    def get_meta(self) -> Dict:
        meta = {
            "cpra": self.cpra,
            "cspra": self.get_cspra(),
            "alt_aln_method": self.alt_aln_method,
            "tss_upstream_limit": self.tss_upstream_limit,
            "uncertain": self.uncertain,
            "promoter_tss_upstream_offset": self.promoter_tss_upstream_offset,
        }
        return meta
