"""GraphQL query construction for the gnomAD API."""

from __future__ import annotations

DEFAULT_DATASET = "gnomad_r4"

# Server-enforced: MAX_QUERY_COST is 25 and a root `variant` field costs 1.
MAX_BATCH_SIZE = 25

# Bump whenever the fragments below change, so the cache can tell which stored
# records predate the new fields and re-fetch only those.
QUERY_VERSION = 1


# --- Fragments ------------------------------------------------------------
#
# A fragment is a named, reusable selection set. `...VariantDetail` splices it
# in. Declaring it once is what keeps the request small no matter how many
# variants the batch carries.

# Frequency block shared by the `exome` and `genome` fields, which have the
# same type. `joint` is a *different* type with fewer fields, so it is spelled
# out separately in VariantDetail below.
SEQUENCING_DATA_FRAGMENT = """
fragment SeqData on VariantDetailsSequencingTypeData {
  ac
  an
  homozygote_count
  hemizygote_count
  filters
  flags
  faf95 { popmax popmax_population }
  faf99 { popmax popmax_population }
  populations { id ac an homozygote_count hemizygote_count }
  local_ancestry_populations { id ac an }
  age_distribution {
    het { bin_edges bin_freq n_smaller n_larger }
    hom { bin_edges bin_freq n_smaller n_larger }
  }
  quality_metrics {
    allele_balance { alt { bin_edges bin_freq n_smaller n_larger } }
    genotype_depth {
      all { bin_edges bin_freq n_smaller n_larger }
      alt { bin_edges bin_freq n_smaller n_larger }
    }
    genotype_quality {
      all { bin_edges bin_freq n_smaller n_larger }
      alt { bin_edges bin_freq n_smaller n_larger }
    }
    site_quality_metrics { metric value }
  }
}
"""

# Everything VariantDetails exposes that is worth caching. Deliberately omitted:
#   af, ac_hom, ac_hemi  -- deprecated; derive from ac/an and *_count
#   va                   -- GA4GH restatement of the frequency data already here;
#                           roughly doubles the payload for no new information
VARIANT_FRAGMENT = """
fragment VariantDetail on VariantDetails {
  variant_id
  reference_genome
  chrom
  pos
  ref
  alt
  caid
  rsids
  colocated_variants
  flags

  coverage {
    exome { mean median over_1 over_5 over_10 over_15 over_20 over_25 over_30 over_50 over_100 }
    genome { mean median over_1 over_5 over_10 over_15 over_20 over_25 over_30 over_50 over_100 }
  }

  exome { ...SeqData }
  genome { ...SeqData }

  joint {
    ac
    an
    homozygote_count
    hemizygote_count
    filters
    faf95 { popmax popmax_population }
    faf99 { popmax popmax_population }
    populations { id ac an homozygote_count hemizygote_count }
    age_distribution {
      het { bin_edges bin_freq n_smaller n_larger }
      hom { bin_edges bin_freq n_smaller n_larger }
    }
    quality_metrics { site_quality_metrics { metric value } }
    freq_comparison_stats {
      contingency_table_test { p_value odds_ratio }
      cochran_mantel_haenszel_test { chisq p_value }
      stat_union { p_value stat_test_name gen_ancs }
    }
  }

  transcript_consequences {
    consequence_terms
    domains
    gene_id
    gene_version
    gene_symbol
    hgvs
    hgvsc
    hgvsp
    is_canonical
    is_mane_select
    is_mane_select_version
    lof
    lof_flags
    lof_filter
    major_consequence
    polyphen_prediction
    sift_prediction
    refseq_id
    refseq_version
    transcript_id
    transcript_version
  }

  in_silico_predictors { id value flags }
  lof_curations { gene_id gene_version gene_symbol verdict flags project }
  multi_nucleotide_variants {
    combined_variant_id
    changes_amino_acids
    n_individuals
    other_constituent_snvs
  }
  non_coding_constraint { chrom start stop element_id possible observed expected oe z }
  vrs {
    _id
    type
    location {
      _id
      type
      sequence_id
      interval { type start { type value } end { type value } }
    }
    state { type sequence }
  }
}
"""

# Mitochondrial variants live behind a different root field with a different
# type. `variant(...)` returns null for them, which would otherwise be cached
# as "not in gnomAD". Note the heteroplasmy-aware counts: there is no single
# `ac`, because a mitochondrial call is homoplasmic or heteroplasmic.
MITOCHONDRIAL_FRAGMENT = """
fragment MitoDetail on MitochondrialVariantDetails {
  variant_id
  reference_genome
  pos
  ref
  alt
  rsids
  flags
  filters
  an
  ac_hom
  ac_het
  ac_hom_mnv
  excluded_ac
  max_heteroplasmy
  heteroplasmy_distribution { bin_edges bin_freq n_smaller n_larger }
  haplogroup_defining
  haplogroups { id an ac_het ac_hom faf faf_hom }
  populations { id an ac_het ac_hom heteroplasmy_distribution { bin_edges bin_freq } }
  mitotip_score
  mitotip_trna_prediction
  pon_ml_probability_of_pathogenicity
  pon_mt_trna_prediction
  age_distribution {
    het { bin_edges bin_freq n_smaller n_larger }
    hom { bin_edges bin_freq n_smaller n_larger }
  }
  site_quality_metrics { name value }
  genotype_quality_metrics { name all { bin_edges bin_freq } alt { bin_edges bin_freq } }
  genotype_quality_filters { name filtered { bin_edges bin_freq } }
  transcript_consequences {
    consequence_terms
    gene_id
    gene_symbol
    hgvsc
    hgvsp
    is_canonical
    is_mane_select
    lof
    lof_flags
    lof_filter
    major_consequence
    transcript_id
  }
}
"""


def alias_for(index: int) -> str:
    """Alias used for the Nth variant in a batch.

    Every field in a GraphQL response is keyed by its name, so 25 `variant`
    fields would collide. Aliasing renames each one in the response.
    """
    return f"v{index}"


def build_variant_query(
    variant_ids: list[str],
    dataset: str = DEFAULT_DATASET,
) -> tuple[str, dict[str, str]]:
    """Build one batched query for up to MAX_BATCH_SIZE nuclear variants.

    Returns (query_text, variables), to be POSTed as
    {"query": query_text, "variables": variables}.
    """
    if not variant_ids:
        raise ValueError("no variant ids given")
    if len(variant_ids) > MAX_BATCH_SIZE:
        raise ValueError(
            f"{len(variant_ids)} variants exceeds the server's per-request " +
            f"cost ceiling of {MAX_BATCH_SIZE}"
        )

    # One "$vN: String!" per variant, plus the dataset. Declaring the dataset as
    # a variable keeps the query text identical for every batch of a given size,
    # which is friendlier to the server's response cache.
    declarations = ", ".join(
        ["$dataset: DatasetId!"]
        + [f"${alias_for(i)}: String!" for i in range(len(variant_ids))]
    )
    selections = "\n".join(
        f"  {alias_for(i)}: variant(variantId: ${alias_for(i)}, dataset: $dataset)" +
        f" {{ ...VariantDetail }}"  # noqa: F541
        for i in range(len(variant_ids))
    )

    query = (
        f"query VariantBatch({declarations}) {{\n{selections}\n}}\n"
        f"{VARIANT_FRAGMENT}{SEQUENCING_DATA_FRAGMENT}"
    )
    variables: dict[str, str] = {"dataset": dataset}
    for i, vid in enumerate(variant_ids):
        variables[alias_for(i)] = vid

    return query, variables


def build_mitochondrial_query(
    variant_ids: list[str],
    dataset: str = DEFAULT_DATASET,
) -> tuple[str, dict[str, str]]:
    """Same as build_variant_query, for chrM variants."""
    if not variant_ids:
        raise ValueError("no variant ids given")
    if len(variant_ids) > MAX_BATCH_SIZE:
        raise ValueError(
            f"{len(variant_ids)} variants exceeds the server's per-request " +
            f"cost ceiling of {MAX_BATCH_SIZE}"
        )

    declarations = ", ".join(
        ["$dataset: DatasetId!"]
        + [f"${alias_for(i)}: String!" for i in range(len(variant_ids))]
    )
    selections = "\n".join(
        f"  {alias_for(i)}: mitochondrial_variant" +
        f"(variant_id: ${alias_for(i)}, dataset: $dataset) {{ ...MitoDetail }}"
        for i in range(len(variant_ids))
    )

    query = (
        f"query MitoBatch({declarations}) {{\n{selections}\n}}\n"
        f"{MITOCHONDRIAL_FRAGMENT}"
    )
    variables: dict[str, str] = {"dataset": dataset}
    for i, vid in enumerate(variant_ids):
        variables[alias_for(i)] = vid

    return query, variables


def parse_batch_response(
    payload: dict,
    variant_ids: list[str],
) -> dict[str, dict | None]:
    """Map a response body back onto the variant ids that produced it.

    Returns {variant_id: record_or_None}, where None means the server had no
    such variant.

    Raises KeyError if the body has no `data` at all, which indicates a
    transport- or query-level failure the caller should retry rather than
    interpret as 25 absent variants.
    """
    data = payload.get("data")
    if data is None:
        raise KeyError("response contained no 'data'")

    return {
        vid: data.get(alias_for(i))
        for i, vid in enumerate(variant_ids)
    }
