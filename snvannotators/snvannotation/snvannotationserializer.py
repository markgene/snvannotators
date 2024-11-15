"""Serialize SnvAnnotation.

Returns a dict instead of JSON string.
"""

from dataclasses import asdict
from datetime import datetime
from typing import Dict

from hgvs.sequencevariant import SequenceVariant

from snvannotators.hgvsplus.models import HgvsG, HgvsP, HgvsT
from .snvannotation import SnvAnnotation


def to_dict(
    obj,
    classkey=None,
    datetime_as_str: bool = False,
    sequence_variant_as_str: bool = True,
) -> dict:
    """Convert any object to dict.

    See https://stackoverflow.com/a/1118038/9721302
    """
    if isinstance(obj, dict):
        data = {}
        for k, v in obj.items():
            data[k] = to_dict(v, classkey, datetime_as_str=datetime_as_str)
        return data
    elif hasattr(obj, "_ast"):
        return to_dict(obj._ast(), datetime_as_str=datetime_as_str)
    elif hasattr(obj, "__iter__") and not isinstance(obj, str):
        return [to_dict(v, classkey, datetime_as_str=datetime_as_str) for v in obj]
    elif isinstance(obj, str):
        return str(obj)
    elif isinstance(obj, datetime):
        if datetime_as_str:
            return str(obj)
    elif (
        isinstance(obj, SequenceVariant)
        or isinstance(obj, HgvsG)
        or isinstance(obj, HgvsP)
        or isinstance(obj, HgvsT)
    ):
        if sequence_variant_as_str:
            return obj.format()
    elif hasattr(obj, "__dict__"):
        data = dict(
            [
                (key, to_dict(value, classkey, datetime_as_str=datetime_as_str))
                for key, value in obj.__dict__.items()
                if not callable(value) and not key.startswith("_")
            ]
        )
        if classkey is not None and hasattr(obj, "__class__"):
            data[classkey] = obj.__class__.__name__
        return data
    return obj


class SnvAnnotationSerializer:
    """Serialize SnvAnnotation.

    Returns a dict instead of JSON string.
    """
    
    def __init__(self, snv_annotation: SnvAnnotation):
        self.snv_annotation = snv_annotation

    def serialize(self) -> Dict:
        d = to_dict(
            asdict(self.snv_annotation),
            datetime_as_str=True,
            sequence_variant_as_str=True,
        )
        return d
