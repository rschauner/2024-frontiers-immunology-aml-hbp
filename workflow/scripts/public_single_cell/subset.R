
seurat <- readRDS(snakemake@input[["seurat"]])

seurat <- subset(seurat, subset = day != "Post-treatment" | sample_origin == "Healthy")

saveRDS(seurat, snakemake@output[["diagnostic"]])

seurat2 <- subset(seurat, subset = collapsed_cell_types %in% c("Myeloid Progenitor (HD)", "LSC-like", "AML Blast"))

saveRDS(seurat2, snakemake@output[["cell_types_subset"]])