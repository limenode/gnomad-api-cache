# gnomad-api-cache

Fetch gnomAD annotations for the variants in a VCF, cache them in SQLite, and
export to various output formats (parquet, csv, tsv, json).

The public gnomAD API allows roughly 10 requests (of 25 variants each) per minute. 
Annotating a cohort twice — or exporting a second format — fetches from the locally 
built cache rather than requerying. 

Re-running into an existing cache fetches only what is missing.

## Install

```bash
pip install gnomad-api-cache

# or with uv
uv pip install gnomad-api-cache    # into the active environment
uv add gnomad-api-cache            # into a project
uv tool install gnomad-api-cache   # just the CLI
```

Requires Python 3.11 or newer.

## Usage

### Command line

Inputs are any VCF-like file (VCF, BCF, or bgzipped VCF) and a SQLite database to store the cache. The output retrieves all variants from the cache and writes them to a file in the specified format.
```bash
gnomad-api-cache -i cohort.vcf.gz -c gnomad.sqlite -o annotations.parquet
```

The output format is inferred from the file extension (`.csv`, `.tsv`,
`.parquet`, `.json`, `.jsonl`) and can be forced with `-f`:

```bash
gnomad-api-cache -i cohort.vcf.gz -c gnomad.sqlite -o annotations.txt -f tsv
```

Omit `-o` to populate the cache without writing a table:

```bash
gnomad-api-cache -i cohort.vcf.gz -c gnomad.sqlite
```

`python -m gnomad_api_cache` accepts the same arguments. Run
`gnomad-api-cache --help` for the full list, including `--dataset`,
`--include-populations`, and `--require-build`.

### Python

```python
from gnomad_api_cache import VariantCache

variant_ids = ["1-55051215-G-A", "M-8993-T-C", ...]

with VariantCache("gnomad.sqlite") as cache:
    summary = cache.fetch(variant_ids)
    print(summary)  # 1234 requested, 234 already cached, 1000 fetched, ...

    cache.to_parquet("annotations.parquet")
    cache.to_csv("annotations.csv")
```

`cache.fetch` takes in gnomAD id strings (hyphen- or colon-separated), 
`(chrom, pos, ref, alt)` tuples, row dicts, a dataframe, or a VCF path.

```python
with VariantCache("gnomad.sqlite") as cache:
    cache.fetch(df)  # any pandas/polars/duckdb frame with chrom/pos/ref/alt columns
    cache.fetch("cohort.vcf.gz")
```

Unusable rows are skipped and tallied; pass `on_error="raise"` to fail on the first one instead.

Reading the VCF yourself lets you filter before querying.
Define a function of the form `filter_function(cyvcf2.Variant) -> bool` 
and pass it to `read_vcf` to only fetch variants that pass the filter:

```python
from cyvcf2 import Variant
from gnomad_api_cache import VariantCache, read_vcf

def rare_only(variant: Variant) -> bool:
    af = variant.INFO.get("gnomad41_exome_AF", 0)
    return float(0 if af == "." else af) <= 0.05

with VariantCache("gnomad.sqlite") as cache:
    cache.fetch(read_vcf("cohort.vcf.gz", filter_function=rare_only))
```

An open cache is a read-only `Mapping` keyed by `chrom-pos-ref-alt`, so cached
records can be retrieved using dictionary-like syntax:

```python
record = cache["1-55051215-G-A"]  # dict, or None if gnomAD has no such variant
```

`to_variant_keys` builds those keys from any input `fetch` takes, so lookup keys
do not need to be manually constructed:

```python
from gnomad_api_cache import to_variant_keys

for key in to_variant_keys([("chr1", 55051215, "G", "A")]):
    record = cache[key.id]  # key.id == "1-55051215-G-A"
```

## Acknowledgement & Citation

This is an unofficial client. It queries the public gnomAD GraphQL API at
<https://gnomad.broadinstitute.org/api> and is not affiliated with or endorsed
by the Broad Institute or the gnomAD project. Please use the shared API
considerately — the default delay between requests is set to stay within the
documented rate limit.

If you use gnomAD data obtained through this tool, cite the current flagship
gnomAD paper. As of the latest release, it is as follows (v4 Preprint, Vancouver):

> Guez J, Goodrich JK, Moldovan MA, Chao KR, Kar P, Panchal R, Wilson MW, Laricchia KM, Rohlicek G, Biba D, Marten D. Integrating 730,947 exome sequences with clinical literature improves gene discovery. Medrxiv. 2026 Mar 25. <https://doi.org/10.64898/2026.03.23.26349081>

gnomAD data use terms: <https://gnomad.broadinstitute.org/terms>

## License

MIT — see [LICENSE](LICENSE).
