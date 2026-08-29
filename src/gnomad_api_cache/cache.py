"""SQLite-backed cache of gnomAD records.

Stores each API response verbatim as compressed JSON, keyed by variant id and
dataset.

Reads behave like a dict:

    with VariantCache("gnomad.sqlite") as cache:
        record = cache["1-55039974-G-T"]     # -> dict, or None if not in gnomAD
        "1-55039974-G-T" in cache            # -> have we ever asked?
        len(cache)

Three columns carry the bookkeeping the fetch loop depends on:

  status         'found' | 'not_found' | 'error'. Caching a miss is what stops
                 absent variants from being re-queried on every run; 'error' is
                 kept distinct so transient failures can be retried without
                 being mistaken for real absences.
  query_version  which version of the GraphQL selection produced this row, so a
                 widened query re-fetches only the rows that predate it.
  fetched_at     provenance, and a basis for staleness policies later.
"""

from __future__ import annotations

import json
import sqlite3
import zlib
from collections.abc import Callable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, Self
from cyvcf2 import Variant

if TYPE_CHECKING:
    from gnomad_api_cache.fetch import FetchSummary

from gnomad_api_cache._utils import _chunked, _now
from gnomad_api_cache.keys import VariantKey
from gnomad_api_cache.query import DEFAULT_DATASET, QUERY_VERSION

# zlib level 6 measured ~8.1x on real gnomAD records (21.4 KB -> 2.6 KB mean)
# at ~0.3 ms/record. Level 9 buys 2% more for 2x the time.
COMPRESSION_LEVEL = 6

# Stay under SQLite's bound-parameter limit when building "IN (?, ?, ...)".
_SQL_CHUNK = 900

STATUS_FOUND = "found"
STATUS_NOT_FOUND = "not_found"
STATUS_ERROR = "error"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS variants (
    variant_id    TEXT    NOT NULL,
    dataset       TEXT    NOT NULL,
    status        TEXT    NOT NULL,
    data          BLOB,
    fetched_at    TEXT    NOT NULL,
    query_version INTEGER NOT NULL,
    PRIMARY KEY (variant_id, dataset)
);

CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""


class VariantCache(Mapping[str, "dict[str, Any] | None"]):
    def __init__(
        self,
        path: str | Path,
        dataset: str = DEFAULT_DATASET,
        query_version: int = QUERY_VERSION,
    ) -> None:
        self.path = str(path)
        self.dataset = dataset
        self.query_version = query_version
        self.db = sqlite3.connect(self.path)
        # WAL lets readers work while a fetch is writing. synchronous=NORMAL is
        # safe under WAL for this use: a crash can lose the last commit, which
        # costs at most one batch that the next run simply re-fetches.
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.execute("PRAGMA synchronous=NORMAL")
        self.db.executescript(_SCHEMA)
        self.db.commit()

    # --- lifecycle --------------------------------------------------------

    def close(self) -> None:
        self.db.commit()
        self.db.close()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def commit(self) -> None:
        self.db.commit()

    # --- what the fetch loop calls ---------------------------------------

    def needs_query(
        self,
        keys: Sequence[VariantKey],
        retry_errors: bool = True,
        retry_not_found: bool = False,
    ) -> list[VariantKey]:
        """Return the subset of `keys` that must be fetched from gnomAD.

        A key is skipped only when a row exists AND it was written by the
        current query_version AND its status is a real answer. Anything older
        or missing comes back for re-fetch.

        The two retry flags cover the two ways a key can be in the cache
        without usable data, which want opposite defaults:

          retry_errors     we never got an answer (network, bad response).
                           Usually transient, so retried by default -- leaving
                           these out would let a run report itself complete
                           while carrying silent holes.
          retry_not_found  gnomAD answered "no such variant". A real answer, so
                           skipped by default; set True to re-check after a
                           gnomAD release adds variants.

        Input order is preserved and duplicates are dropped, so the caller can
        pass a raw VCF-derived list straight in.
        """
        # Deduplicate up front
        unique: dict[str, VariantKey] = {}
        for key in keys:
            unique.setdefault(key.id, key)

        fresh: set[str] = set()
        ids = list(unique)
        for chunk in _chunked(ids, _SQL_CHUNK):
            placeholders = ",".join("?" * len(chunk))
            rows = self.db.execute(
                f"SELECT variant_id, status FROM variants "
                f"WHERE dataset = ? AND query_version >= ? "
                f"AND variant_id IN ({placeholders})",
                (self.dataset, self.query_version, *chunk),
            ).fetchall()
            for variant_id, status in rows:
                if retry_errors and status == STATUS_ERROR:
                    continue
                if retry_not_found and status == STATUS_NOT_FOUND:
                    continue
                fresh.add(variant_id)

        return [key for vid, key in unique.items() if vid not in fresh]

    def put(
        self,
        key: VariantKey | str,
        record: dict[str, Any] | None,
        status: str | None = None,
    ) -> None:
        """Stage one record. Call commit() (or put_many) to persist.

        `status` defaults to 'found' when a record is given and 'not_found'
        when it is None.
        """
        variant_id = key.id if isinstance(key, VariantKey) else key
        if status is None:
            status = STATUS_FOUND if record is not None else STATUS_NOT_FOUND
        blob = (
            zlib.compress(json.dumps(record).encode(), COMPRESSION_LEVEL)
            if record is not None
            else None
        )
        self.db.execute(
            "INSERT OR REPLACE INTO variants "
            "(variant_id, dataset, status, data, fetched_at, query_version) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (variant_id, self.dataset, status, blob, _now(), self.query_version),
        )

    def put_many(
        self,
        records: Mapping[str, dict[str, Any] | None],
        status: str | None = None,
    ) -> None:
        """Stage and commit a whole batch."""
        for variant_id, record in records.items():
            self.put(variant_id, record, status)
        self.db.commit()

    def mark_errors(self, keys: Sequence[VariantKey]) -> None:
        """Record a failed batch so a later run can retry just those keys."""
        for key in keys:
            self.put(key, None, STATUS_ERROR)
        self.db.commit()

    # --- dict-style reads -------------------------------------------------

    def __getitem__(self, variant_id: str) -> dict[str, Any] | None:
        row = self.db.execute(
            "SELECT status, data FROM variants WHERE variant_id = ? AND dataset = ?",
            (variant_id, self.dataset),
        ).fetchone()
        if row is None:
            raise KeyError(variant_id)
        status, blob = row
        if status != STATUS_FOUND or blob is None:
            return None
        return json.loads(zlib.decompress(blob))

    def __iter__(self) -> Iterator[str]:
        for (variant_id,) in self.db.execute(
            "SELECT variant_id FROM variants WHERE dataset = ?", (self.dataset,)
        ):
            yield variant_id

    def __len__(self) -> int:
        row = self.db.execute(
            "SELECT count(*) FROM variants WHERE dataset = ?", (self.dataset,)
        ).fetchone()
        return int(row[0])

    def many(self, variant_ids: Sequence[str]) -> dict[str, dict[str, Any] | None]:
        """Fetch many records in one round trip. Absent ids are simply omitted."""
        out: dict[str, dict[str, Any] | None] = {}
        for chunk in _chunked(list(variant_ids), _SQL_CHUNK):
            placeholders = ",".join("?" * len(chunk))
            for variant_id, status, blob in self.db.execute(
                f"SELECT variant_id, status, data FROM variants "
                f"WHERE dataset = ? AND variant_id IN ({placeholders})",
                (self.dataset, *chunk),
            ):
                out[variant_id] = (
                    json.loads(zlib.decompress(blob))
                    if status == STATUS_FOUND and blob is not None
                    else None
                )
        return out

    # --- export ----------------------------------------------------------
    #
    # Thin delegations to gnomad_api_cache.export, which holds the flattening
    # logic. Imported lazily so cache.py has no import-time dependency on the
    # export module (or, through it, on pyarrow).

    def to_json(
        self,
        path: str | Path,
        variant_ids: Sequence[str] | None = None,
        lines: bool = False,
        indent: int | None = None,
    ) -> int:
        """Write records verbatim as JSON. See export.to_json."""
        from gnomad_api_cache import export

        return export.to_json(self, path, variant_ids, lines, indent)

    def to_csv(
        self,
        path: str | Path,
        variant_ids: Sequence[str] | None = None,
        include_populations: bool = False,
    ) -> int:
        """Write flattened rows as CSV. See export.to_csv."""
        from gnomad_api_cache import export

        return export.to_csv(self, path, variant_ids, include_populations)

    def to_tsv(
        self,
        path: str | Path,
        variant_ids: Sequence[str] | None = None,
        include_populations: bool = False,
    ) -> int:
        """Write flattened rows as TSV. See export.to_tsv."""
        from gnomad_api_cache import export

        return export.to_tsv(self, path, variant_ids, include_populations)

    def to_parquet(
        self,
        path: str | Path,
        variant_ids: Sequence[str] | None = None,
        include_populations: bool = False,
    ) -> int:
        """Write flattened rows as Parquet. See export.to_parquet."""
        from gnomad_api_cache import export

        return export.to_parquet(self, path, variant_ids, include_populations)

    def status_counts(self) -> dict[str, int]:
        """Row counts by status -- a one-line health check on the cache."""
        return {
            status: count
            for status, count in self.db.execute(
                "SELECT status, count(*) FROM variants WHERE dataset = ? "
                "GROUP BY status",
                (self.dataset,),
            )
        }

    def fetch(self, variants: list[VariantKey], **kwargs) -> FetchSummary:
        """Fetch missing records into this cache. See fetch.fetch_into."""
        from gnomad_api_cache import fetch as _fetch

        return _fetch.fetch_into(self, variants, **kwargs)

    def fetch_vcf(
        self,
        vcf_path: str | Path,
        retry_errors: bool = True,
        retry_not_found: bool = False,
        filter_function: Callable[[Variant], bool] | None = None,
    ) -> FetchSummary:
        """Read a VCF and fetch missing records into this cache."""
        from gnomad_api_cache import fetch as _fetch
        from gnomad_api_cache.adapters import vcf_adapter

        return _fetch.fetch_into(
            self,
            vcf_adapter.read_vcf(
                vcf_path,
                filter_function=filter_function
            ),
            retry_errors=retry_errors,
            retry_not_found=retry_not_found,
        )
