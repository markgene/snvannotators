"""Annotator."""

import myvariant
from ..myvariantannotation import MyvariantAnnotation


class MyvariantHgvsChrAnnotator:
    """MyVariant annotator for a variant in HGVS."""

    instance = None

    def __new__(cls, *args, **kwargs):
        if cls.instance is None:
            cls.instance = super().__new__(cls)
        return cls.instance

    def __init__(self):
        self.mv = myvariant.MyVariantInfo()

    def annotate(self, hgvs_chr: str) -> MyvariantAnnotation:
        """Annotate."""
        anno = self.mv.getvariant(hgvs_chr, fields="all", assembly="hg19")
        return MyvariantAnnotation(hgvs_chr=hgvs_chr, raw=anno)
