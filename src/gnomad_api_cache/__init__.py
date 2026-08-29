"""Fetch and cache gnomAD annotations for VCF variants.

    from gnomad_api_cache import VariantCache, read_vcf

    with VariantCache("gnomad.sqlite") as cache:
        print(cache.fetch_vcf("cohort.vcf.gz"))
        cache.to_parquet("annotations.parquet")
"""

from __future__ import annotations

from gnomad_api_cache.adapters.vcf_adapter import (
    BuildMismatchError,
    iter_variant_keys,
    read_vcf,
)
from gnomad_api_cache.cache import VariantCache
from gnomad_api_cache.fetch import FetchSummary, fetch_gnomad
from gnomad_api_cache.keys import InvalidVariantError, VariantKey
from gnomad_api_cache.query import DEFAULT_DATASET, QUERY_VERSION

__version__ = "0.1.0"

__all__ = [
    "DEFAULT_DATASET",
    "QUERY_VERSION",
    "BuildMismatchError",
    "FetchSummary",
    "InvalidVariantError",
    "VariantCache",
    "VariantKey",
    "fetch_gnomad",
    "iter_variant_keys",
    "read_vcf",
]
