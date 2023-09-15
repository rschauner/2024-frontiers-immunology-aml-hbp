
box::use(
    Seurat[NormalizeData, FindVariableFeatures, ScaleData, SCTransform]
)

seurat <- readRDS(snakemake@input[["seurat"]])

verbose <- TRUE

seurat <- NormalizeData(seurat, verbose = verbose)
seurat <- FindVariableFeatures(seurat, verbose = verbose)
seurat <- ScaleData(seurat, verbose = verbose)

if (snakemake@params[["SCTransform"]] == "yes") {
    seurat <- SCTransform(seurat, return.only.var.genes = FALSE, verbose = verbose)
}

saveRDS(seurat, snakemake@output[["seurat"]])