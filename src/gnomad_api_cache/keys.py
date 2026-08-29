from __future__ import annotations

import re
from dataclasses import dataclass

CANONICAL_CHROMS = frozenset([str(i) for i in range(1, 23)] + ["X", "Y", "M"])

_ALLELE_RE = re.compile(r"^[ACGTN]+$")


class InvalidVariantError(ValueError):
    """Raised when a record cannot be expressed as a gnomAD variant ID."""


def normalize_chrom(chrom: str) -> str:
    """Canonicalize a contig name to the form gnomAD expects."""
    c = str(chrom).strip()

    # strip "chr" prefix
    if c.lower().startswith("chr"):
        c = c[3:]

    # convert "MT" to "M"
    c = c.upper()
    if c == "MT":
        c = "M"

    # remove leading zeros from numeric chromosomes
    if c.isdigit():
        c = str(int(c))  # "01" -> "1"

    return c


@dataclass(frozen=True)
class VariantKey:
    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def id(self) -> str:
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"

    @property
    def is_mitochondrial(self) -> bool:
        return self.chrom == "M"

    @classmethod
    def from_parts(cls, chrom: str, pos: int | str, ref: str, alt: str) -> VariantKey:
        """Build a key from raw record fields, canonicalizing and validating."""
        c = normalize_chrom(chrom)
        r = str(ref).strip().upper()
        a = str(alt).strip().upper()

        if c not in CANONICAL_CHROMS:
            raise InvalidVariantError(f"non-canonical contig: {chrom!r}")
        try:
            p = int(pos)
        except (TypeError, ValueError):
            raise InvalidVariantError(f"non-integer position: {pos!r}") from None
        if p < 1:
            raise InvalidVariantError(f"position must be 1-based: {pos!r}")
        if not _ALLELE_RE.match(r):
            raise InvalidVariantError(f"non-ACGTN ref allele: {ref!r}")
        if not _ALLELE_RE.match(a):
            raise InvalidVariantError(f"non-ACGTN alt allele: {alt!r}")

        return cls(chrom=c, pos=p, ref=r, alt=a)

    @classmethod
    def from_id(cls, variant_id: str) -> VariantKey:
        """Parse a "chrom-pos-ref-alt" string, e.g. from a text file or the cache."""
        parts = str(variant_id).strip().split("-")
        if len(parts) != 4:
            raise InvalidVariantError(f"expected chrom-pos-ref-alt, got {variant_id!r}")
        return cls.from_parts(*parts)
