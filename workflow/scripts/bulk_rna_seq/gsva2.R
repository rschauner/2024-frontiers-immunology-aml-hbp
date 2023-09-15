
box::use(dplyr)

lazyLoad("data/cached_objects_from_pitzer/GSVA DGE object_4f5313a141c6deea3bfc02d341f41dc4")

info <- readRDS(snakemake@input[["dge"]])$samples

library(Biobase)

BinExpressionLevels <- function(ratio, percent, high, low) {
    ratio_ord <- ratio[order(ratio)]
    n_samples <- floor(length(ratio_ord) * percent)
    dplyr$case_when(
        ratio > rev(ratio_ord)[n_samples] ~ high,
        ratio < ratio_ord[n_samples] ~ low,
        .default = glue$glue("middle{(1 - percent * 2) * 100}_pct")
    )
}

phenoData(gsva_eset) <- phenoData(gsva_eset) |>
    as("data.frame") |>
    dplyr$mutate(
        OGT_OGA_ratio_bin = info$OGT_OGA_ratio_bin,
        OGT = info$OGT,
        OGA = info$OGA,
        OGT_bin = info$OGT_bin,
        OGA_bin = info$OGA_bin
    ) |>
    as("AnnotatedDataFrame")

saveRDS(gsva_eset, snakemake@output[["dge"]])