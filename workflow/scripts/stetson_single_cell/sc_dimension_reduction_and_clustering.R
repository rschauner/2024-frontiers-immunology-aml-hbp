
box::use(
    Seurat,
    so = SeuratObject[`DefaultAssay<-`],
    glue[glue],
    purrr
)

seurat <- readRDS(snakemake@input[["seurat"]])

seurat_split <- Seurat$SplitObject(seurat, split.by = "patient_id")

seurat_split <- purrr$map(seurat_split, Seurat$SCTransform, return.only.var.genes = FALSE, verbose = FALSE)
features <- Seurat$SelectIntegrationFeatures(object.list = seurat_split, nfeatures = length(rownames(seurat)), verbose = FALSE)
seurat_split <- Seurat$PrepSCTIntegration(object.list = seurat_split, anchor.features = features, verbose = FALSE)

anchors <- Seurat$FindIntegrationAnchors(object.list = seurat_split, normalization.method = "SCT", anchor.features = features, verbose = FALSE)
integrated <- Seurat$IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 50, verbose = FALSE)

seurat[["integrated"]] <- integrated[["integrated"]]

graphs <- c("nn", "snn")

for (assay in so$Assays(seurat)) {
    DefaultAssay(seurat) <- assay

    seurat <- Seurat$RunPCA(seurat, reduction.name = glue("{assay}_pca"), reduction.key = glue("{assay}PC_"), verbose = FALSE)
    seurat <- Seurat$FindNeighbors(seurat, reduction = glue("{assay}_pca"), graph.name = glue("{assay}_{graphs}"), verbose = FALSE)
    seurat <- Seurat$FindClusters(seurat, graph.name = glue("{assay}_snn"), verbose = FALSE)
    seurat[[glue("{assay}_clusters")]] <- seurat[["seurat_clusters"]]

    seurat <- Seurat$RunUMAP(
        seurat,
        dims = 1:20,
        reduction = glue("{assay}_pca"),
        return.model = TRUE,
        reduction.name = glue("{assay}_umap"),
        reduction.key = glue("{assay}UMAP_"),
        verbose = FALSE
    )
}

saveRDS(seurat, snakemake@output[["seurat"]])