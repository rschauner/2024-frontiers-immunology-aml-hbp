box::use(
    SeuratObject[`DefaultAssay<-`],
    Seurat,
    Cairo,
    dplyr,
    ggplot2,
    glue,
    tidyselect,
    future
)

seurat <- readRDS(snakemake@input[["seurat"]])

DrawGroupedUMAP <- function(seurat, group_by, reduction) {
    Seurat$DimPlot(seurat, group.by = group_by, reduction = reduction)
}
DrawFeatureUMAP <- function(seurat, features, reduction) {
    Seurat$FeaturePlot(seurat, features = features, combine = FALSE, order = TRUE, reduction = reduction)
}
DrawDotPlot <- function(seurat, group_by, scale = TRUE) {
    Seurat$DotPlot(
        seurat,
        features = c("OGT", "OGA"),
        group.by = group_by,
        #scale.min = 20, scale.max = 60,
        col.min = -2, col.max = 2,
        scale = scale
    )
}
DrawBarGraph <- function(seurat, assay) {
    cols_use <- c("patient_id", glue$glue("{assay}_clusters"))
    seurat[[cols_use]] |>
        dplyr$group_by(dplyr$across(tidyselect$all_of(cols_use))) |>
        dplyr$summarise(n = dplyr$n()) |>
        dplyr$mutate(f = n / sum(n) * 100) |>
        ggplot2$ggplot(mapping = ggplot2$aes_string(x = "patient_id", fill = glue$glue("{assay}_clusters"), y = "f")) +
        ggplot2$geom_col() +
        ggplot2$theme_classic()
}

DrawClusterHeatmap <- function(seurat, assay) {
    res <- Seurat$FindAllMarkers(
        seurat,
        #test.use = "roc",
        assay = assay,
        only.pos = TRUE,
        group.by = "integrated_clusters"
    )
    features <- res |>
        dplyr$group_by(cluster) |>
        dplyr$slice_head(n = 15) |>
        dplyr$pull(gene)
    if (assay == "SCT") seurat <- Seurat$GetResidual(seurat, features = features)
    Seurat$DoHeatmap(seurat, features = features, raster = FALSE)
}

if (snakemake@params[["type"]] == "Grouped UMAP") {
    plot <- DrawGroupedUMAP(seurat, group_by = snakemake@params[["group_by"]], reduction = snakemake@params[["reduction"]])
} else if (snakemake@params[["type"]] == "Dot Plot") {
    DefaultAssay(seurat) <- snakemake@params[["fea_assay"]]
    plot <- list(
        DrawDotPlot(seurat, snakemake@params[["group_by"]]),
        DrawDotPlot(seurat, snakemake@params[["group_by"]], scale = FALSE)
    )
} else if (snakemake@params[["type"]] == "Bar Plot") {
    assay <- snakemake@params[["assay"]]
    plot <- DrawBarGraph(seurat, snakemake@params[["assay"]])
} else if (snakemake@params[["type"]] == "Feature UMAP") {
    DefaultAssay(seurat) <- snakemake@params[["assay"]]
    plot <- DrawFeatureUMAP(seurat, snakemake@params[["features"]], snakemake@params[["reduction"]])
} else if (snakemake@params[["type"]] == "Cluster Heatmap") {
    future$plan("multicore", workers = future$availableCores())
    DefaultAssay(seurat) <- snakemake@params[["assay"]]
    plot <- DrawClusterHeatmap(seurat, snakemake@params[["assay"]])
} else {
    stop("Plot type not recognized")
}

Cairo$CairoPDF(snakemake@output[["pdf"]], dpi = 600, width = 180, height = 180)
plot
graphics.off()