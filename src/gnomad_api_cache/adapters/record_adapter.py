"""Build VariantKey objects from in-memory Python data.

Covers gnomAD id strings, 4-tuples, row mappings, and any dataframe exposing
the Arrow PyCapsule interface (pandas >= 2.2, polars >= 1.3, duckdb relations,
pyarrow tables). `to_variant_keys` handles all cases.

Each reader supplies callables to the collect() function in _collect.py.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping, Sequence
from itertools import chain
from pathlib import Path
from typing import Any

from gnomad_api_cache.adapters._collect import ON_ERROR_SKIP, collect
from gnomad_api_cache.keys import InvalidVariantError, VariantKey

FIELDS = ("chrom", "pos", "ref", "alt")

# Common headers for chrom, pos, ref, alt.
_FIELD_BY_ALIAS: Mapping[str, str] = {
    "chr": "chrom",
    "chrom": "chrom",
    "chromosome": "chrom",
    "contig": "chrom",
    "pos": "pos",
    "position": "pos",
    "start_1based": "pos",
    "ref": "ref",
    "reference": "ref",
    "ref_allele": "ref",
    "reference_allele": "ref",
    "alt": "alt",
    "alternate": "alt",
    "alt_allele": "alt",
    "alternate_allele": "alt",
}

# Interval-frame columns are 0-based half-open (BED, bioframe, pyranges) while
# gnomAD positions are 1-based. Treat them as an error to avoid off-by-one mistakes.
_INTERVAL_ALIASES = frozenset({"start", "chromstart", "end", "chromend", "stop"})

_VCF_SUFFIXES = (".vcf", ".vcf.gz", ".vcf.bgz", ".bcf", ".bcf.gz")


def _looks_like_vcf_path(text: str) -> bool:
    return text.lower().endswith(_VCF_SUFFIXES)


# --- row -> key functions ------------------------------------------------


def _key_from_id(row: Any) -> VariantKey:
    return VariantKey.from_id(str(row))


def _key_from_tuple(row: Any) -> VariantKey:
    if isinstance(row, (str, bytes)) or not isinstance(row, Sequence):
        raise InvalidVariantError(
            f"not a (chrom, pos, ref, alt) row: {row!r}", "malformed_row"
        )
    if len(row) != 4:
        raise InvalidVariantError(
            f"expected 4 fields (chrom, pos, ref, alt), got {len(row)}: {row!r}",
            "malformed_row",
        )
    chrom, pos, ref, alt = row
    return VariantKey.from_parts(str(chrom), pos, str(ref), str(alt))


def _key_from_mapping(
    row: Any,
    columns: Mapping[str, str] | None = None,
) -> VariantKey:
    if not isinstance(row, Mapping):
        raise InvalidVariantError(f"not a mapping: {row!r}", "malformed_row")

    if columns is not None:
        try:
            values = {field: row[columns[field]] for field in FIELDS}
        except KeyError as exc:
            raise InvalidVariantError(
                f"column {exc.args[0]!r} missing from row: {sorted(row)}",
                "missing_fields",
            ) from None
    else:
        # Lowercase the row's keys
        lowered = {str(key).lower(): value for key, value in row.items()}
        values = {}
        for alias, field in _FIELD_BY_ALIAS.items():
            if field not in values and alias in lowered:
                values[field] = lowered[alias]

        missing = [field for field in FIELDS if field not in values]
        if missing:
            if "pos" in missing and _INTERVAL_ALIASES & lowered.keys():
                raise InvalidVariantError(
                    "interval columns (start/end) are 0-based half-open while "
                    "gnomAD positions are 1-based; supply a 1-based 'pos' "
                    f"column instead: {sorted(row)}",
                    "interval_columns",
                )
            raise InvalidVariantError(
                f"no column for {missing} in row: {sorted(row)}", "missing_fields"
            )

    return VariantKey.from_parts(
        str(values["chrom"]), values["pos"], str(values["ref"]), str(values["alt"])
    )


# --- readers -------------------------------------------------------------


def iter_ids(
    variant_ids: Iterable[Any],
    *,
    on_error: str = ON_ERROR_SKIP,
    label: str | None = None,
) -> Iterator[VariantKey]:
    """Yield keys from "chrom-pos-ref-alt" (or colon-separated) id strings."""
    return collect(variant_ids, _key_from_id, on_error=on_error, label=label)


def iter_tuples(
    variants: Iterable[Any],
    *,
    on_error: str = ON_ERROR_SKIP,
    label: str | None = None,
) -> Iterator[VariantKey]:
    """Yield keys from (chrom, pos, ref, alt) rows, e.g. df.itertuples()."""
    return collect(variants, _key_from_tuple, on_error=on_error, label=label)


def iter_dicts(
    variants: Iterable[Any],
    *,
    columns: Mapping[str, str] | None = None,
    on_error: str = ON_ERROR_SKIP,
    label: str | None = None,
) -> Iterator[VariantKey]:
    """Yield keys from row mappings, e.g. df.to_dict("records").

    Column names are matched case-insensitively against a table of common
    spellings; `columns` overrides that with an explicit
    {"chrom": ..., "pos": ..., "ref": ..., "alt": ...} mapping.
    """
    return collect(
        variants,
        lambda row: _key_from_mapping(row, columns),
        on_error=on_error,
        label=label,
    )


def iter_arrow(
    source: Any,
    *,
    columns: Mapping[str, str] | None = None,
    on_error: str = ON_ERROR_SKIP,
    label: str | None = None,
) -> Iterator[VariantKey]:
    """Yield keys from any object implementing the Arrow PyCapsule interface."""
    import pyarrow

    def _rows() -> Iterator[dict[str, Any]]:
        # Stream batch by batch so a large frame is never fully materialized
        # as Python objects.
        reader = pyarrow.RecordBatchReader.from_stream(source)
        for batch in reader:
            yield from batch.to_pylist()

    return iter_dicts(_rows(), columns=columns, on_error=on_error, label=label)


def read_ids(variant_ids: Iterable[Any], **kwargs: Any) -> list[VariantKey]:
    """Eager convenience wrapper around iter_ids."""
    return list(iter_ids(variant_ids, **kwargs))


def read_tuples(variants: Iterable[Any], **kwargs: Any) -> list[VariantKey]:
    """Eager convenience wrapper around iter_tuples."""
    return list(iter_tuples(variants, **kwargs))


def read_dicts(variants: Iterable[Any], **kwargs: Any) -> list[VariantKey]:
    """Eager convenience wrapper around iter_dicts."""
    return list(iter_dicts(variants, **kwargs))


def read_arrow(source: Any, **kwargs: Any) -> list[VariantKey]:
    """Eager convenience wrapper around iter_arrow."""
    return list(iter_arrow(source, **kwargs))


# --- dispatch ------------------------------------------------------------


def to_variant_keys(
    source: Any,
    *,
    columns: Mapping[str, str] | None = None,
    on_error: str = ON_ERROR_SKIP,
) -> list[VariantKey]:
    """Coerce any supported input into a list of VariantKey.

    Accepts a VariantKey, a VCF path (str/Path ending in a VCF suffix), a
    single id string, a dataframe exposing the Arrow PyCapsule interface, a
    single row mapping, or an iterable of ids / 4-tuples / row mappings /
    cyvcf2 Variant records.
    """
    if isinstance(source, VariantKey):
        return [source]

    if isinstance(source, Path) or (
        isinstance(source, str) and _looks_like_vcf_path(source)
    ):
        from gnomad_api_cache.adapters.vcf_adapter import read_vcf

        return read_vcf(source, on_error=on_error)

    if isinstance(source, str):
        return [VariantKey.from_id(source)]

    if hasattr(source, "__arrow_c_stream__"):
        return read_arrow(source, columns=columns, on_error=on_error)

    # A dataframe predating the PyCapsule interface. Intercept it before the
    # generic iterable path, where iterating would yield column names.
    if hasattr(source, "to_dict") and hasattr(source, "columns"):
        return read_dicts(
            source.to_dict("records"), columns=columns, on_error=on_error
        )

    if isinstance(source, Mapping):
        return read_dicts([source], columns=columns, on_error=on_error)

    if not isinstance(source, Iterable):
        raise TypeError(f"cannot build variant keys from {type(source).__name__}")

    # Peek at the first row to pick a reader, then put it back.
    iterator = iter(source)
    try:
        first = next(iterator)
    except StopIteration:
        return []
    rows = chain([first], iterator)

    if isinstance(first, VariantKey):
        return [key for key in rows if isinstance(key, VariantKey)]
    if isinstance(first, (str, bytes)):
        return read_ids(rows, on_error=on_error)
    if isinstance(first, Mapping):
        return read_dicts(rows, columns=columns, on_error=on_error)
    if hasattr(first, "ALT") and hasattr(first, "CHROM"):
        from gnomad_api_cache.adapters.vcf_adapter import read_variants

        return read_variants(rows, on_error=on_error)
    return read_tuples(rows, on_error=on_error)
