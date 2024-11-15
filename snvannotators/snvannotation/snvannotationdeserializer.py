"""Deserialize a serialized SnvAnnotation of dict."""

from typing import Any, Dict

from dacite import from_dict
from hgvs.sequencevariant import SequenceVariant

from snvannotators.hgvsplus.models import HgvsC, HgvsG, HgvsP, HgvsT
from snvannotators.hgvsplus.parsers.hgvsparser import HgvsParser

from .snvannotation import SnvAnnotation

hgvs_parser = HgvsParser()


def prepare(data) -> Any:
    sequence_variant_field_names = [
        "hgvs_g",
        "hgvs_g_normalized",
        "hgvs_c",
        "hgvs_t",
        "hgvs_p",
    ]
    if isinstance(data, dict):
        for key, value in data.items():
            if key in sequence_variant_field_names:
                if value is None:
                    data[key] = value
                elif isinstance(value, str):
                    if ":" in value:
                        data[key] = hgvs_parser.parse(s=value)
                    else:
                        data[key] = value
                elif (
                    isinstance(value, SequenceVariant)
                    or isinstance(value, HgvsC)
                    or isinstance(value, HgvsG)
                    or isinstance(value, HgvsP)
                    or isinstance(value, HgvsT)
                ):
                    data[key] = value
                else:
                    raise ValueError(
                        f"SequenceVariant field not have a string: {value}"
                    )
            else:
                data[key] = prepare(value)
        return data
    elif isinstance(data, list):
        return [prepare(elem) for elem in data]
    elif isinstance(data, tuple):
        return tuple([prepare(elem) for elem in data])
    else:
        return data


class SnvAnnotationDeserializer:
    instance = None

    def __new__(cls, *args, **kwargs):
        if cls.instance is None:
            cls.instance = super().__new__(cls)
        return cls.instance

    def deserialize(self, data: Dict) -> SnvAnnotation:
        prepared = prepare(data)
        snv_annotation = from_dict(data_class=SnvAnnotation, data=prepared)
        return snv_annotation
