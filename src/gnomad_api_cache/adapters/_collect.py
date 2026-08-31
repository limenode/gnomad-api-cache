"""The row -> VariantKey loop."""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Callable, Iterable, Iterator
from typing import Any

from gnomad_api_cache.keys import InvalidVariantError, VariantKey

log = logging.getLogger(__name__)

ON_ERROR_SKIP = "skip"
ON_ERROR_RAISE = "raise"
ON_ERROR_CHOICES = (ON_ERROR_SKIP, ON_ERROR_RAISE)


def collect(
    rows: Iterable[Any],
    to_keys: Callable[[Any], VariantKey],
    *,
    on_error: str = ON_ERROR_SKIP,
    label: str | None = None,
) -> Iterator[VariantKey]:
    """Build keys from rows, tally and logging what could not be built successfully."""
    if on_error not in ON_ERROR_CHOICES:
        raise ValueError(
            f"on_error must be one of {ON_ERROR_CHOICES}, got {on_error!r}"
        )

    skipped: Counter[str] = Counter()
    yielded = 0

    for row in rows:
        try:
            key = to_keys(row)
        except InvalidVariantError as exc:
            if on_error == ON_ERROR_RAISE:
                raise
            skipped[exc.reason] += 1
            log.debug("skipping %r -- %s", row, exc)
            continue
        yielded += 1
        yield key

    prefix = f"{label}: " if label else ""
    if skipped:
        log.info(
            "%s%d keys, %d rows skipped (%s)",
            prefix,
            yielded,
            sum(skipped.values()),
            ", ".join(f"{k}={v}" for k, v in sorted(skipped.items())),
        )
    else:
        log.info("%s%d keys", prefix, yielded)
