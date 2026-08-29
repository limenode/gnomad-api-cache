"""Minimal local stubs for the cyvcf2 surface this project uses.

cyvcf2 is a compiled Cython extension that ships no py.typed and no .pyi, and
no stub package exists on PyPI or in typeshed. Rather than disable pyright's
unknown-type rules project-wide, we declare the members we touch. This covers
the parsing API plus the INFO/FORMAT accessors used for record filtering; add
to it when you reach for more of the API.

Upstream reference: https://brentp.github.io/cyvcf2/
"""

from collections.abc import Iterator
from typing import Any

class INFO:
    def get(self, key: str, default: Any = ...) -> Any: ...
    def __getitem__(self, key: str) -> Any: ...
    def __contains__(self, key: str) -> bool: ...

class Variant:
    CHROM: str
    POS: int
    REF: str
    # cyvcf2 splits multiallelic ALT into a list; empty when the VCF has "."
    ALT: list[str]
    # None when the VCF FILTER column is PASS, otherwise the filter string.
    FILTER: str | None
    ID: str | None
    QUAL: float | None
    INFO: INFO
    # 0-based half-open, unlike POS which is 1-based.
    start: int
    end: int
    is_snp: bool
    is_indel: bool
    var_type: str
    # numpy arrays; typed loosely since numpy stubs are not a dependency here.
    gt_types: Any
    gt_depths: Any
    gt_quals: Any
    genotypes: Any
    def format(self, field: str, vtype: str = ...) -> Any: ...

class VCF:
    def __init__(self, fname: str, /, *, gts012: bool = ..., lazy: bool = ...) -> None: ...
    def __iter__(self) -> Iterator[Variant]: ...
    def __next__(self) -> Variant: ...
    # Region query, e.g. vcf("chr1:1000-2000"); requires an index.
    def __call__(self, region: str) -> Iterator[Variant]: ...
    def close(self) -> None: ...
    @property
    def raw_header(self) -> str: ...
    @property
    def seqnames(self) -> list[str]: ...
    @property
    def seqlens(self) -> list[int]: ...
    @property
    def samples(self) -> list[str]: ...
    @property
    def num_records(self) -> int: ...
