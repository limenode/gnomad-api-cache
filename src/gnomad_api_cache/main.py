from __future__ import annotations

from cyvcf2 import Variant

from gnomad_api_cache import VariantCache

if __name__ == "__main__":
    vcf_file = "/lab01/Projects/Lionel_Projects/dukeData/vcf/duke_quartet.annovar_slivar.fix.vcf.gz"
    cache_file = "/lab01/Projects/Lionel_Projects/gnomad_api_caller/data/gnomad.sqlite"

    def filter_function_1(variant: Variant) -> bool:
        val = variant.INFO.get("gnomad41_exome_AF", 0)
        if val == ".":
            val = 0.0
        val = float(val)

        return val <= 0.05

    with VariantCache(cache_file) as cache:
        print(cache.fetch_vcf(vcf_file, filter_function=filter_function_1))
        cache.to_parquet(
            "/lab01/Projects/Lionel_Projects/gnomad_api_caller/data/gnomad.parquet"
        )
        cache.to_csv(
            "/lab01/Projects/Lionel_Projects/gnomad_api_caller/data/gnomad.csv"
        )
