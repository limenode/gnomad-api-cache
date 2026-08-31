"""Fetch and cache gnomAD annotations for VCF variants.

    from gnomad_api_cache import VariantCache, read_vcf

    with VariantCache("gnomad.sqlite") as cache:
        print(cache.fetch("cohort.vcf.gz"))
        cache.to_parquet("annotations.parquet")
"""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

__version__ = "0.2.0"

# Public name -> the module that defines it.
_LAZY_IMPORTS: dict[str, str] = {
    "BuildMismatchError": "gnomad_api_cache.adapters.vcf_adapter",
    "iter_variant_keys": "gnomad_api_cache.adapters.vcf_adapter",
    "read_variants": "gnomad_api_cache.adapters.vcf_adapter",
    "read_vcf": "gnomad_api_cache.adapters.vcf_adapter",
    "read_arrow": "gnomad_api_cache.adapters.record_adapter",
    "read_dicts": "gnomad_api_cache.adapters.record_adapter",
    "read_ids": "gnomad_api_cache.adapters.record_adapter",
    "read_tuples": "gnomad_api_cache.adapters.record_adapter",
    "to_variant_keys": "gnomad_api_cache.adapters.record_adapter",
    "VariantCache": "gnomad_api_cache.cache",
    "FetchSummary": "gnomad_api_cache.fetch",
    "fetch_into": "gnomad_api_cache.fetch",
    "InvalidVariantError": "gnomad_api_cache.keys",
    "VariantKey": "gnomad_api_cache.keys",
    "DEFAULT_DATASET": "gnomad_api_cache.query",
    "QUERY_VERSION": "gnomad_api_cache.query",
}

if TYPE_CHECKING:
    # Type checkers do not follow __getattr__ well enough to give these real
    # types, so state them here. Never executed at runtime.
    from gnomad_api_cache.adapters.record_adapter import (
        read_arrow,
        read_dicts,
        read_ids,
        read_tuples,
        to_variant_keys,
    )
    from gnomad_api_cache.adapters.vcf_adapter import (
        BuildMismatchError,
        iter_variant_keys,
        read_variants,
        read_vcf,
    )
    from gnomad_api_cache.cache import VariantCache
    from gnomad_api_cache.fetch import FetchSummary, fetch_into
    from gnomad_api_cache.keys import InvalidVariantError, VariantKey
    from gnomad_api_cache.query import DEFAULT_DATASET, QUERY_VERSION


def __getattr__(name: str) -> Any:
    module_name = _LAZY_IMPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(module_name), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted([*__all__, "__version__"])


__all__ = [
    "DEFAULT_DATASET",
    "QUERY_VERSION",
    "BuildMismatchError",
    "FetchSummary",
    "InvalidVariantError",
    "VariantCache",
    "VariantKey",
    "fetch_into",
    "iter_variant_keys",
    "read_arrow",
    "read_dicts",
    "read_ids",
    "read_tuples",
    "read_variants",
    "read_vcf",
    "to_variant_keys",
]
