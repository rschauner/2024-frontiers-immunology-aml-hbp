suppressPackageStartupMessages(
    box::use(
        Seurat,
        so=SeuratObject[`DefaultAssay<-`],
        patchwork,
        glue[glue]
    )
)

seurat <- readRDS(snakemake@input[["seurat"]])

seurat[["cell_types_origin"]] <- paste0(seurat$CellType, "_", seurat$sample_origin)

group <- snakemake@params[["group"]]
DefaultAssay(seurat) <- assay <- snakemake@params[["assay"]]
wildcard <- snakemake@params[["subset"]]
title <- glue("Assay: {assay} | Group: {group} | Subset: {wildcard}")

features <- snakemake@params[["features"]]


p1 <- Seurat$VlnPlot(seurat, features = features, group.by = group, pt.size = 0) +
    patchwork$plot_annotation(title = title) +
    Seurat$NoLegend()

p2 <- Seurat$DotPlot(
    seurat,
    features = features,
    group.by = group,
    scale.min = 20, scale.max = 60,
    col.min = -2, col.max = 2
) +
    patchwork$plot_annotation(title = title)

p3 <- Seurat$DotPlot(
    seurat,
    features = features,
    group.by = group,
    col.min = -2, col.max = 2,
    scale = FALSE
) +
    patchwork$plot_annotation(title = title)

pdf(snakemake@output[["pdf"]], width = 8, height = 10)
p1
p2
p3
graphics.off()
