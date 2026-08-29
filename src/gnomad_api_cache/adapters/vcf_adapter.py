"""Read a VCF into VariantKey objects."""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Callable, Iterator
from pathlib import Path

from cyvcf2 import VCF, Variant

from gnomad_api_cache.keys import (
    CANONICAL_CHROMS,
    InvalidVariantError,
    VariantKey,
    normalize_chrom,
)

log = logging.getLogger(__name__)

# Use chr1 to detect build
CHR1_LENGTH_BY_BUILD = {248956422: "GRCh38", 249250621: "GRCh37"}


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


def _skip_reason(alt: str, chrom: str) -> str | None:
    """Return why this allele can't become a gnomAD ID, or None if it can."""
    if alt.startswith("<"):
        return "symbolic"  # <DEL>, <NON_REF>, <*>
    if "[" in alt or "]" in alt:
        return "breakend"
    if alt == "*":
        return "spanning_deletion"
    if normalize_chrom(chrom) not in CANONICAL_CHROMS:
        return "non_canonical_contig"  # *_random, chrUn_*, HLA-*
    return None


def iter_variant_keys(
    path: str | Path,
    *,
    require_build: str | None = "GRCh38",
    warn_unnormalized: bool = True,
    filter_function: Callable[[Variant], bool] | None = None,
) -> Iterator[VariantKey]:
    """Yield one VariantKey per ALT allele of every record in `path`."""
    vcf = VCF(str(path))

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

    skipped: Counter[str] = Counter()
    yielded = 0

    # main loop
    for record in vcf:
        if filter_function is not None and not filter_function(record):
            continue

        for alt in record.ALT:  # empty list when ALT is "."
            reason = _skip_reason(alt, record.CHROM)
            if reason is not None:
                skipped[reason] += 1
                continue
            try:
                key = VariantKey.from_parts(record.CHROM, record.POS, record.REF, alt)
            except InvalidVariantError as e:
                skipped["invalid"] += 1
                log.debug(
                    "skipping %s:%s %s>%s -- %s",
                    record.CHROM,
                    record.POS,
                    record.REF,
                    alt,
                    e,
                )
                continue
            yielded += 1
            yield key

    vcf.close()

    if skipped:
        log.info(
            "%s: %d keys, %d alleles skipped (%s)",
            path,
            yielded,
            sum(skipped.values()),
            ", ".join(f"{k}={v}" for k, v in sorted(skipped.items())),
        )
    else:
        log.info("%s: %d keys", path, yielded)


def read_vcf(
    path: str | Path,
    *,
    require_build: str | None = "GRCh38",
    warn_unnormalized: bool = True,
    filter_function: Callable[[Variant], bool] | None = None,
) -> list[VariantKey]:
    """Eager convenience wrapper around iter_variant_keys."""
    return list(
        iter_variant_keys(
            path,
            require_build=require_build,
            warn_unnormalized=warn_unnormalized,
            filter_function=filter_function
        )
    )
