"""Read a VCF, or already-loaded cyvcf2 records, into VariantKey objects."""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable, Iterator
from pathlib import Path

from cyvcf2 import VCF, Variant

from gnomad_api_cache.adapters._collect import ON_ERROR_SKIP, collect
from gnomad_api_cache.keys import InvalidVariantError, VariantKey

log = logging.getLogger(__name__)

# Use chr1 to detect build
CHR1_LENGTH_BY_BUILD = {248956422: "GRCh38", 249250621: "GRCh37"}

AltRow = tuple[Variant, str]


class BuildMismatchError(RuntimeError):
    """Raised when a VCF is not on the build the target gnomAD dataset uses."""


def detect_build(vcf: VCF) -> str | None:
    """Return "GRCh38", "GRCh37", or None if the header doesn't say."""
    lengths: dict[str, int] = dict(zip(vcf.seqnames, vcf.seqlens))
    chr1 = lengths.get("chr1") or lengths.get("1")
    if chr1 is None:
        return None
    return CHR1_LENGTH_BY_BUILD.get(chr1)


def is_normalized(vcf: VCF) -> bool:
    """True if the header shows bcftools norm was run."""
    return any(
        line.startswith("##bcftools_norm") for line in vcf.raw_header.split("\n")
    )


def _alt_syntax_reason(alt: str) -> str | None:
    """Name the VCF ALT spelling that cannot become a gnomAD ID, or None."""
    if alt.startswith("<"):
        return "symbolic"  # <DEL>, <NON_REF>, <*>
    if "[" in alt or "]" in alt:
        return "breakend"
    if alt == "*":
        return "spanning_deletion"
    return None


def _iter_alt_rows(
    records: Iterable[Variant],
    filter_function: Callable[[Variant], bool] | None,
) -> Iterator[AltRow]:
    """Fan records out to one (record, alt) row per ALT allele."""
    for record in records:
        if filter_function is not None and not filter_function(record):
            continue
        for alt in record.ALT:  # empty list when ALT is "."
            yield record, alt


def _key_from_alt_row(row: AltRow) -> VariantKey:
    record, alt = row
    reason = _alt_syntax_reason(alt)
    if reason is not None:
        raise InvalidVariantError(f"unusable ALT allele: {alt!r}", reason)
    return VariantKey.from_parts(record.CHROM, record.POS, record.REF, alt)


def iter_variant_keys_from_records(
    records: Iterable[Variant],
    *,
    filter_function: Callable[[Variant], bool] | None = None,
    on_error: str = ON_ERROR_SKIP,
    label: str | None = None,
) -> Iterator[VariantKey]:
    """Yield one VariantKey per ALT allele of already-loaded cyvcf2 records."""
    return collect(
        _iter_alt_rows(records, filter_function),
        _key_from_alt_row,
        on_error=on_error,
        label=label,
    )


def iter_variant_keys(
    path: str | Path,
    *,
    require_build: str | None = "GRCh38",
    warn_unnormalized: bool = True,
    filter_function: Callable[[Variant], bool] | None = None,
    on_error: str = ON_ERROR_SKIP,
) -> Iterator[VariantKey]:
    """Yield one VariantKey per ALT allele of every record in `path`."""
    vcf = VCF(str(path))
    try:
        build = detect_build(vcf)
        if require_build is not None:
            if build is None:
                log.warning(
                    "%s: could not determine build from header; assuming %s",
                    path,
                    require_build,
                )
            elif build != require_build:
                raise BuildMismatchError(
                    f"{path}: VCF is {build}, expected {require_build}. "
                    + "Lift over first, or target a matching gnomAD dataset."
                )

        if warn_unnormalized and not is_normalized(vcf):
            log.warning(
                "%s: no ##bcftools_norm header. Unnormalized indels will silently "
                + "miss in gnomAD. Run: bcftools norm -f <ref.fa> -m -any",
                path,
            )

        yield from iter_variant_keys_from_records(
            vcf,
            filter_function=filter_function,
            on_error=on_error,
            label=str(path),
        )
    finally:
        # try/finally because a caller may abandon this generator part-way.
        vcf.close()


def read_vcf(
    path: str | Path,
    *,
    require_build: str | None = "GRCh38",
    warn_unnormalized: bool = True,
    filter_function: Callable[[Variant], bool] | None = None,
    on_error: str = ON_ERROR_SKIP,
) -> list[VariantKey]:
    """Eager convenience wrapper around iter_variant_keys."""
    return list(
        iter_variant_keys(
            path,
            require_build=require_build,
            warn_unnormalized=warn_unnormalized,
            filter_function=filter_function,
            on_error=on_error,
        )
    )


def read_variants(
    variants: Iterable[Variant],
    *,
    filter_function: Callable[[Variant], bool] | None = None,
    on_error: str = ON_ERROR_SKIP,
) -> list[VariantKey]:
    """Eager convenience wrapper around iter_variant_keys_from_records."""
    return list(
        iter_variant_keys_from_records(
            variants,
            filter_function=filter_function,
            on_error=on_error,
        )
    )
