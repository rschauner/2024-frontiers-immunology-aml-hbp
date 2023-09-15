
box::use(
    Seurat[FindAllMarkers],
    glue[glue],
    readr[write_tsv],
    tibble[remove_rownames],
    dplyr[relocate]
)

seurat <- readRDS(snakemake@input[["seurat"]])
assay <- snakemake@params[["assay"]]

if (snakemake@params[["group_by"]] == "clusters") {
    snakemake@params[["group_by"]] <- glue("{assay}_clusters")
}


res <- FindAllMarkers(
    seurat,
    features = c("OGT", "OGA"),
    logfc.threshold = 0,
    test.use = "LR",
    assay = assay,
    group.by = snakemake@params[["group_by"]],
    latent.vars = c("S.Score", "G2M.Score")
)
if (nrow(res) > 0) {
    res <- res |>
        relocate(gene, cluster) |>
        remove_rownames()
}

write_tsv(res, file = snakemake@output[["tsv"]])
