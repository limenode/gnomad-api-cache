"""Export cached gnomAD records to Parquet, JSON, or delimited text.

Three formats, two scopes (whole cache, or a given list of variant ids):

    export.to_json(cache, "out.json")                  # nested, verbatim
    export.to_csv(cache, "out.csv", variant_ids=ids)   # flattened, one row each
    export.to_parquet(cache, "out.parquet")            # flattened, typed

JSON preserves the record exactly as gnomAD returned it. Parquet and CSV/TSV
are flattened to one row per variant, which drops the nested parts
(per-population breakdowns unless asked for, histograms, the full transcript
list) -- use JSON when you need those.

Mitochondrial records have a different shape from nuclear ones: no exome or
genome block, and heteroplasmy-aware counts instead of a single AC. They are
flattened into the same table with a `kind` column and `mito_`-prefixed
columns, so a mixed export loses nothing and stays self-describing.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from gnomad_api_cache.cache import VariantCache

# Predictors gnomAD v4 exposes; fixed so columns stay stable across exports
# even when a given variant is missing some.
IN_SILICO_PREDICTORS = (
    "cadd",
    "revel_max",
    "spliceai_ds_max",
    "pangolin_largest_ds",
    "phylop",
    "sift_max",
    "polyphen_max",
)

# gnomAD's top-level genetic ancestry groups.
ANCESTRY_GROUPS = (
    "afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining",
)

# exome and genome share a schema; joint has no `af` field, so it is derived.
SEQUENCING_TYPES = ("exome", "genome", "joint")


def _join(values: Any) -> str | None:
    """Collapse a list of scalars into one delimited cell, or None if empty."""
    if not values:
        return None
    return ",".join(str(v) for v in values)


def _pick_transcript(record: dict[str, Any]) -> dict[str, Any]:
    """The transcript a flat row should describe.

    Prefers MANE Select, then Ensembl canonical, then the first entry, matching
    how the gnomAD browser picks the consequence it shows.
    """
    consequences = record.get("transcript_consequences") or []
    for predicate in (
        lambda t: t.get("is_mane_select"),
        lambda t: t.get("is_canonical"),
    ):
        for transcript in consequences:
            if predicate(transcript):
                return transcript
    return consequences[0] if consequences else {}


def _sequencing_columns(
    block: dict[str, Any] | None,
    prefix: str,
    include_populations: bool,
) -> dict[str, Any]:
    """Frequency columns for one of exome / genome / joint."""
    block = block or {}
    ac, an = block.get("ac"), block.get("an")
    faf95 = block.get("faf95") or {}

    columns: dict[str, Any] = {
        f"{prefix}_ac": ac,
        f"{prefix}_an": an,
        # derive af from ac/an because joint has no af field
        f"{prefix}_af": (ac / an) if ac is not None and an else None,
        f"{prefix}_nhomalt": block.get("homozygote_count"),
        f"{prefix}_nhemi": block.get("hemizygote_count"),
        f"{prefix}_filters": _join(block.get("filters")),
        f"{prefix}_faf95": faf95.get("popmax"),
        f"{prefix}_faf95_population": faf95.get("popmax_population"),
    }

    if include_populations:
        by_id = {p.get("id"): p for p in block.get("populations") or []}
        for group in ANCESTRY_GROUPS:
            population = by_id.get(group) or {}
            p_ac, p_an = population.get("ac"), population.get("an")
            columns[f"{prefix}_ac_{group}"] = p_ac
            columns[f"{prefix}_an_{group}"] = p_an
            columns[f"{prefix}_af_{group}"] = (
                (p_ac / p_an) if p_ac is not None and p_an else None
            )
            columns[f"{prefix}_nhomalt_{group}"] = population.get("homozygote_count")

    return columns


def _flatten_mitochondrial(record: dict[str, Any]) -> dict[str, Any]:
    """Mito-specific columns. A mito call is homoplasmic or heteroplasmic, so
    there is no single AC to map onto the nuclear columns."""
    return {
        "mito_an": record.get("an"),
        "mito_ac_hom": record.get("ac_hom"),
        "mito_ac_het": record.get("ac_het"),
        "mito_max_heteroplasmy": record.get("max_heteroplasmy"),
        "mito_haplogroup_defining": record.get("haplogroup_defining"),
        "mito_mitotip_score": record.get("mitotip_score"),
        "mito_mitotip_prediction": record.get("mitotip_trna_prediction"),
        "mito_pon_mt_trna_prediction": record.get("pon_mt_trna_prediction"),
    }


def flatten_record(
    variant_id: str,
    record: dict[str, Any] | None,
    include_populations: bool = False,
) -> dict[str, Any]:
    """Reduce one cached record to a single flat row."""
    data = record or {}
    # Mito records carry no `chrom` and no `exome` block; the id is the only
    # place the contig appears.
    parts = variant_id.split("-")
    is_mito = parts[0] == "M"

    row: dict[str, Any] = {
        "variant_id": variant_id,
        "kind": "mitochondrial" if is_mito else "nuclear",
        "in_gnomad": record is not None,
        "chrom": data.get("chrom") or (parts[0] if len(parts) == 4 else None),
        "pos": data.get("pos")
        or (int(parts[1]) if len(parts) == 4 and parts[1].isdigit() else None),
        "ref": data.get("ref") or (parts[2] if len(parts) == 4 else None),
        "alt": data.get("alt") or (parts[3] if len(parts) == 4 else None),
        "rsid": (data.get("rsids") or [None])[0],
        "caid": data.get("caid"),
        "flags": _join(data.get("flags")),
    }

    row.update(_consequence_columns(_pick_transcript(data)))
    row["n_transcripts"] = len(data.get("transcript_consequences") or [])

    for prefix in SEQUENCING_TYPES:
        row.update(_sequencing_columns(data.get(prefix), prefix, include_populations))

    predictors = {
        p.get("id"): p.get("value") for p in data.get("in_silico_predictors") or []
    }
    row.update(
        {f"in_silico_{name}": predictors.get(name) for name in IN_SILICO_PREDICTORS}
    )

    coverage = data.get("coverage") or {}
    row["exome_mean_coverage"] = (coverage.get("exome") or {}).get("mean")
    row["genome_mean_coverage"] = (coverage.get("genome") or {}).get("mean")

    row.update(_flatten_mitochondrial(data if is_mito else {}))
    return row


def _consequence_columns(transcript: dict[str, Any]) -> dict[str, Any]:
    return {
        "gene_symbol": transcript.get("gene_symbol"),
        "gene_id": transcript.get("gene_id"),
        "transcript_id": transcript.get("transcript_id"),
        "consequence": transcript.get("major_consequence"),
        "hgvsc": transcript.get("hgvsc"),
        "hgvsp": transcript.get("hgvsp"),
        "lof": transcript.get("lof"),
        "lof_filter": transcript.get("lof_filter"),
        "lof_flags": transcript.get("lof_flags"),
        "polyphen": transcript.get("polyphen_prediction"),
        "sift": transcript.get("sift_prediction"),
        "n_transcripts": 0,
    }


def _select(
    cache: VariantCache,
    variant_ids: Sequence[str] | None,
) -> Iterator[tuple[str, dict[str, Any] | None]]:
    """Yield (variant_id, record) for the requested scope.

    Streams one row at a time so a whole-cache export does not need the entire
    cache resident in memory.
    """
    ids = list(cache) if variant_ids is None else list(variant_ids)
    for variant_id in ids:
        try:
            yield variant_id, cache[variant_id]
        except KeyError:
            # Asked for an id that was never queried: report it as a row with
            # no data rather than omitting it without explanation.
            yield variant_id, None


def iter_rows(
    cache: VariantCache,
    variant_ids: Sequence[str] | None = None,
    include_populations: bool = False,
) -> Iterator[dict[str, Any]]:
    """Flat rows for the requested scope, one dict per variant."""
    for variant_id, record in _select(cache, variant_ids):
        yield flatten_record(variant_id, record, include_populations)


def to_json(
    cache: VariantCache,
    path: str | Path,
    variant_ids: Sequence[str] | None = None,
    lines: bool = False,
    indent: int | None = None,
) -> int:
    """Write records verbatim, preserving all nesting.

    Default is one JSON object keyed by variant id. `lines=True` writes JSON
    Lines instead, which streams and is the better choice for a large cache --
    records average ~21 KB, so a 100k-variant dump is a couple of GB.
    """
    count = 0
    with open(path, "w") as handle:
        if lines:
            for variant_id, record in _select(cache, variant_ids):
                handle.write(
                    json.dumps({"variant_id": variant_id, "record": record}) + "\n"
                )
                count += 1
        else:
            payload = {vid: rec for vid, rec in _select(cache, variant_ids)}
            json.dump(payload, handle, indent=indent)
            count = len(payload)
    return count


def to_delimited(
    cache: VariantCache,
    path: str | Path,
    variant_ids: Sequence[str] | None = None,
    delimiter: str = ",",
    include_populations: bool = False,
) -> int:
    """Write flattened rows as delimited text.

    Columns are taken from the first row, so every row must share a schema --
    which is why flatten_record always emits the full column set, including for
    variants that are absent from gnomAD.
    """
    rows = iter_rows(cache, variant_ids, include_populations)
    first = next(rows, None)
    if first is None:
        raise ValueError("nothing to export")

    count = 0
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(first), delimiter=delimiter)
        writer.writeheader()
        writer.writerow(first)
        count = 1
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def to_csv(
    cache: VariantCache,
    path: str | Path,
    variant_ids: Sequence[str] | None = None,
    include_populations: bool = False,
) -> int:
    return to_delimited(cache, path, variant_ids, ",", include_populations)


def to_tsv(
    cache: VariantCache,
    path: str | Path,
    variant_ids: Sequence[str] | None = None,
    include_populations: bool = False,
) -> int:
    return to_delimited(cache, path, variant_ids, "\t", include_populations)


def to_parquet(
    cache: VariantCache,
    path: str | Path,
    variant_ids: Sequence[str] | None = None,
    include_populations: bool = False,
    compression: str = "zstd",
) -> int:
    """Write flattened rows as Parquet, ready for pandas/polars/duckdb.

    pyarrow is imported here rather than at module scope so that the JSON and
    CSV exports keep working without it installed.
    """
    import pyarrow
    import pyarrow.parquet

    rows = list(iter_rows(cache, variant_ids, include_populations))
    if not rows:
        raise ValueError("nothing to export")

    table = pyarrow.Table.from_pylist(rows)
    pyarrow.parquet.write_table(table, str(path), compression=compression)
    return len(rows)
