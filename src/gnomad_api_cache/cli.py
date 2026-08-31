"""Command-line interface for gnomad-api-cache."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from gnomad_api_cache import __version__
from gnomad_api_cache.query import DEFAULT_DATASET

if TYPE_CHECKING:
    from gnomad_api_cache.cache import VariantCache

# Output format per lowercased file extension. --format overrides this, and
# becomes mandatory when the extension is missing or not listed here.
EXTENSION_FORMATS = {
    ".csv": "csv",
    ".tsv": "tsv",
    ".tab": "tsv",
    ".parquet": "parquet",
    ".pq": "parquet",
    ".json": "json",
    ".jsonl": "jsonl",
    ".ndjson": "jsonl",
}

OUTPUT_FORMATS = ("csv", "tsv", "parquet", "json", "jsonl")

INPUT_FORMATS = ("vcf",)

# Sentinel for --require-build meaning "accept whatever the VCF declares".
ANY_BUILD = "any"

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_USAGE = 2
EXIT_INTERRUPTED = 130


def infer_format(path: Path) -> str | None:
    """Map an output path to an export format by extension, or None."""
    return EXTENSION_FORMATS.get(path.suffix.lower())


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser.

    Exposed separately from main() so documentation tooling and tests can
    introspect the interface without running anything.
    """
    parser = argparse.ArgumentParser(
        prog="gnomad-api-cache",
        description=(
            "Fetch gnomAD annotations for the variants in a VCF, cache them "
            "in SQLite, and optionally export a flattened table."
        ),
        epilog=(
            "examples:\n"
            "  gnomad-api-cache -i cohort.vcf.gz -c gnomad.sqlite "
            "-o out.parquet\n"
            "  gnomad-api-cache -i cohort.vcf.gz -c gnomad.sqlite "
            "-o out.txt -f tsv\n"
            "  gnomad-api-cache -i cohort.vcf.gz -c gnomad.sqlite\n"
            "\n"
            "The cache is the durable artifact: re-running against the same "
            "cache re-fetches\nonly what is missing, so exporting a second "
            "format costs no API calls."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        metavar="PATH",
        help="input variant file (VCF, optionally bgzipped)",
    )
    parser.add_argument(
        "-c",
        "--cache",
        required=True,
        type=Path,
        metavar="PATH",
        help="SQLite cache file; created if it does not exist",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="PATH",
        help=(
            "export destination; omit to populate the cache without "
            "writing a table"
        ),
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=OUTPUT_FORMATS,
        help=(
            "output format (default: inferred from the -o extension; "
            "required when the extension is unrecognised)"
        ),
    )
    parser.add_argument(
        "--input-format",
        choices=INPUT_FORMATS,
        default="vcf",
        help="input format (default: %(default)s)",
    )

    fetching = parser.add_argument_group("fetching")
    fetching.add_argument(
        "--dataset",
        default=DEFAULT_DATASET,
        metavar="NAME",
        help="gnomAD dataset id (default: %(default)s)",
    )
    fetching.add_argument(
        "--delay",
        type=float,
        default=None,
        metavar="SECONDS",
        help=(
            "seconds between API requests; the public API allows 10 per "
            "minute, so lower this only against a private instance"
        ),
    )
    fetching.add_argument(
        "--require-build",
        default="GRCh38",
        metavar="BUILD",
        help=(
            "reference build the VCF must declare, or 'any' to skip the "
            "check (default: %(default)s)"
        ),
    )
    fetching.add_argument(
        "--no-retry-errors",
        dest="retry_errors",
        action="store_false",
        help="leave previously errored variants alone instead of retrying",
    )
    fetching.add_argument(
        "--retry-not-found",
        action="store_true",
        help=(
            "re-query variants previously absent from gnomAD (useful only "
            "after a dataset release)"
        ),
    )

    export = parser.add_argument_group("export")
    export.add_argument(
        "--include-populations",
        action="store_true",
        help="add per-ancestry frequency columns to tabular output",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="suppress the run summary on stdout",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def _export(
    cache: VariantCache,
    path: Path,
    fmt: str,
    include_populations: bool,
) -> int:
    """Dispatch to the exporter for `fmt`. Returns the row/record count."""
    if fmt == "csv":
        return cache.to_csv(path, include_populations=include_populations)
    if fmt == "tsv":
        return cache.to_tsv(path, include_populations=include_populations)
    if fmt == "parquet":
        return cache.to_parquet(path, include_populations=include_populations)
    if fmt == "json":
        return cache.to_json(path)
    if fmt == "jsonl":
        return cache.to_json(path, lines=True)
    raise ValueError(f"unsupported output format: {fmt}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # Resolve the output format before fetching
    output_format: str | None = args.format
    if args.output is not None and output_format is None:
        output_format = infer_format(args.output)
        if output_format is None:
            parser.error(
                f"cannot infer an output format from '{args.output.name}'; "
                f"pass -f/--format with one of: {', '.join(OUTPUT_FORMATS)}"
            )

    if not args.input.exists():
        print(f"error: input file not found: {args.input}", file=sys.stderr)
        return EXIT_ERROR

    require_build = (
        None if args.require_build.lower() == ANY_BUILD else args.require_build
    )

    from gnomad_api_cache.adapters.vcf_adapter import BuildMismatchError, read_vcf
    from gnomad_api_cache.cache import VariantCache

    fetch_kwargs = {
        "retry_errors": args.retry_errors,
        "retry_not_found": args.retry_not_found,
    }
    if args.delay is not None:
        fetch_kwargs["delay"] = args.delay

    try:
        variants = read_vcf(args.input, require_build=require_build)
        if not variants:
            print(
                f"error: no usable variants in {args.input}", file=sys.stderr
            )
            return EXIT_ERROR

        with VariantCache(args.cache, dataset=args.dataset) as cache:
            summary = cache.fetch(variants, **fetch_kwargs)
            if not args.quiet:
                print(summary)

            if args.output is not None and output_format is not None:
                args.output.parent.mkdir(parents=True, exist_ok=True)
                count = _export(
                    cache,
                    args.output,
                    output_format,
                    args.include_populations,
                )
                if not args.quiet:
                    print(
                        f"wrote {count} records to {args.output} "
                        f"({output_format})"
                    )
    except BuildMismatchError as exc:
        print(f"error: {exc}", file=sys.stderr)
        print(
            "hint: pass --require-build any to skip the build check",
            file=sys.stderr,
        )
        return EXIT_ERROR
    except ImportError as exc:
        # Every export dependency ships with the package, so this means a
        # damaged environment rather than a missing optional install.
        print(f"error: {exc}", file=sys.stderr)
        print(
            "hint: the install looks incomplete; try "
            "pip install --force-reinstall gnomad-api-cache",
            file=sys.stderr,
        )
        return EXIT_ERROR
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_ERROR
    except KeyboardInterrupt:
        # The cache commits per batch, so an interrupted run keeps its work.
        print("\ninterrupted; cached progress is preserved", file=sys.stderr)
        return EXIT_INTERRUPTED

    return EXIT_OK
