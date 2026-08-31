from __future__ import annotations

import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from typing import Any

import requests
from tqdm import tqdm

from gnomad_api_cache._utils import _chunked
from gnomad_api_cache.adapters.record_adapter import to_variant_keys
from gnomad_api_cache.cache import VariantCache
from gnomad_api_cache.client import post_gnomad
from gnomad_api_cache.keys import VariantKey
from gnomad_api_cache.query import (
    MAX_BATCH_SIZE,
    build_mitochondrial_query,
    build_variant_query,
    parse_batch_response,
)

# Documented limit is 10 requests per IP per 60s.
REQUEST_DELAY_SECONDS = 8

QueryBuilder = Callable[[list[str]], tuple[str, dict[str, str]]]


@dataclass(frozen=True)
class FetchSummary:
    """Fetch run summary. Returned by fetch_into and VariantCache.fetch."""

    requested: int = 0
    already_cached: int = 0
    fetched: int = 0
    not_found: int = 0
    errors: int = 0
    failed_batches: int = 0

    def __str__(self) -> str:
        return (
            f"{self.requested} requested, {self.already_cached} already cached, "
            f"{self.fetched} fetched, {self.not_found} not in gnomAD, "
            f"{self.errors} errored"
        )


def _fetch_batches(
    cache: VariantCache,
    keys: Sequence[VariantKey],
    build_query: QueryBuilder,
    progress: tqdm,
    delay: float,
) -> FetchSummary:
    """Query `keys` in batches, depositing each response as it arrives.

    Returns the tally for these batches only; the caller merges it."""
    fetched = not_found = errors = failed_batches = 0
    for batch in _chunked(keys, MAX_BATCH_SIZE):
        ids = [v.id for v in batch]
        try:
            query, variables = build_query(ids)
            payload = post_gnomad(query, variables)
            records = parse_batch_response(payload, ids)
            cache.put_many(records)
        except (requests.RequestException, KeyError) as e:
            cache.mark_errors(batch)
            errors += len(batch)
            failed_batches += 1
            tqdm.write(f"batch of {len(batch)} failed, recorded for retry: {e}")
        else:
            # A null record means gnomAD has no such variant, not a failure.
            found = sum(1 for record in records.values() if record is not None)
            fetched += found
            not_found += len(records) - found
        progress.update(len(batch))
        time.sleep(delay)
    return FetchSummary(
        fetched=fetched,
        not_found=not_found,
        errors=errors,
        failed_batches=failed_batches,
    )


def fetch_into(
    cache: VariantCache,
    variants: Any,
    retry_errors: bool = True,
    retry_not_found: bool = False,
    delay: float = REQUEST_DELAY_SECONDS,
) -> FetchSummary:
    """Populate an open cache with gnomAD records for `variants`.

    `variants` is anything record_adapter.to_variant_keys accepts: VariantKey
    objects, id strings, 4-tuples, row mappings, a dataframe, or a VCF path.
    Coercing here rather than in each caller means every entry point takes the
    same inputs.

    The cache is left open: it belongs to the caller, who may well want to
    export from it next.
    """
    variants = to_variant_keys(variants)
    variants_to_fetch = cache.needs_query(
        variants,
        retry_errors=retry_errors,
        retry_not_found=retry_not_found,
    )
    summary = FetchSummary(
        requested=len(variants),
        already_cached=len(variants) - len(variants_to_fetch),
    )

    if not variants_to_fetch:
        tqdm.write("All variants already cached.")
        return summary

    nuclear = [v for v in variants_to_fetch if not v.is_mitochondrial]
    mitochondrial = [v for v in variants_to_fetch if v.is_mitochondrial]

    with tqdm(
        total=len(variants_to_fetch), unit="variant", desc="Querying gnomAD"
    ) as progress:
        nuclear_tally = _fetch_batches(
            cache, nuclear, build_variant_query, progress, delay
        )
        mito_tally = _fetch_batches(
            cache, mitochondrial, build_mitochondrial_query, progress, delay
        )

    summary = replace(
        summary,
        fetched=nuclear_tally.fetched + mito_tally.fetched,
        not_found=nuclear_tally.not_found + mito_tally.not_found,
        errors=nuclear_tally.errors + mito_tally.errors,
        failed_batches=nuclear_tally.failed_batches + mito_tally.failed_batches,
    )

    if summary.failed_batches:
        tqdm.write(
            f"{summary.failed_batches} batch(es) failed and were recorded as "
            f"errors; re-run to retry them."
        )
    return summary
