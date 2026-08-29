import time
from collections.abc import Callable, Sequence

import requests
from tqdm import tqdm

from gnomad_api_caller.cache import VariantCache, _chunked
from gnomad_api_caller.keys import VariantKey
from gnomad_api_caller.query import (
    MAX_BATCH_SIZE,
    build_mitochondrial_query,
    build_variant_query,
    parse_batch_response,
    post_gnomad,
)

# Documented limit is 10 requests per IP per 60s.
REQUEST_DELAY_SECONDS = 6.1

QueryBuilder = Callable[[list[str]], tuple[str, dict[str, str]]]


def _fetch_batches(
    cache: VariantCache,
    keys: Sequence[VariantKey],
    build_query: QueryBuilder,
    progress: tqdm,
    delay: float,
) -> int:
    """Query `keys` in batches, depositing each response as it arrives.

    Returns the number of batches that failed."""
    failures = 0
    for batch in _chunked(keys, MAX_BATCH_SIZE):
        ids = [v.id for v in batch]
        try:
            query, variables = build_query(ids)
            payload = post_gnomad(query, variables)
            cache.put_many(parse_batch_response(payload, ids))
        except (requests.RequestException, KeyError) as e:
            cache.mark_errors(batch)
            failures += 1
            tqdm.write(f"batch of {len(batch)} failed, recorded for retry: {e}")
        progress.update(len(batch))
        time.sleep(delay)
    return failures


def fetch_gnomad(
    variants: list[VariantKey],
    cache_file: str,
    retry_errors: bool = True,
    retry_not_found: bool = False,
    delay: float = REQUEST_DELAY_SECONDS,
) -> None:
    """Populate the cache with gnomAD records for `variants`."""
    cache = VariantCache(cache_file)
    try:
        variants_to_fetch = cache.needs_query(
            variants,
            retry_errors=retry_errors,
            retry_not_found=retry_not_found,
        )
        nuclear = [v for v in variants_to_fetch if not v.is_mitochondrial]
        mitochondrial = [v for v in variants_to_fetch if v.is_mitochondrial]

        with tqdm(
            total=len(variants_to_fetch), unit="variant", desc="Querying gnomAD"
        ) as progress:
            failures = _fetch_batches(
                cache, nuclear, build_variant_query, progress, delay
            )
            failures += _fetch_batches(
                cache, mitochondrial, build_mitochondrial_query, progress, delay
            )

        if failures:
            tqdm.write(
                f"{failures} batch(es) failed and were recorded as errors; "
                f"re-run to retry them."
            )
    finally:
        cache.close()
