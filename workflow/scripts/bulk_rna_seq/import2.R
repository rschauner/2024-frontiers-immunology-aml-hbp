
box::use(dplyr, glue)
lazyLoad("data/cached_objects_from_pitzer/Make DGE List_a9c21d1667770a5650264e11dc418228")

percent <- snakemake@params[["percent"]]
percent <- ifelse(percent > 1, percent / 100, percent)

BinExpressionLevels <- function(ratio, percent, high, low) {
    ratio_ord <- ratio[order(ratio)]
    n_samples <- floor(length(ratio_ord) * percent)
    dplyr$case_when(
        ratio > rev(ratio_ord)[n_samples] ~ high,
        ratio < ratio_ord[n_samples] ~ low,
        .default = glue$glue("middle{(1 - percent * 2) * 100}_pct")
    )
}

dge$samples <- dplyr$mutate(
    dge$samples,
    OGT_OGA_ratio = OGT / OGA,
    OGT_OGA_ratio_bin = BinExpressionLevels(OGT_OGA_ratio, percent, "OGT", "OGA"),
    OGT_bin = BinExpressionLevels(OGT, percent, "OGT_high", "OGT_low"),
    OGA_bin = BinExpressionLevels(OGA, percent, "OGA_high", "OGA_low")
)

saveRDS(dge, snakemake@output[["dge"]])