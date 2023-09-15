
ChangeFeatureNames <- function(object, old, new) {
    assay <- object[["RNA"]]
    for (slot_use in c("counts", "data", "scale.data", "meta.features")) {
        x <- slot(assay, slot_use)
        rownames(x)[rownames(x) == old] <- new
        slot(assay, slot_use) <- x
    }
    if (old %in% assay@var.features) {
        box::use(stringr)
        slot(assay, "var.features") <- stringr$str_replace(slot(assay, "var.features"), old, new)
    }
    object[["RNA"]] <- assay
    return(object)
}

box::use(
    Seurat,
    so=SeuratObject,
    stringr,
    purrr,
    dplyr,
    tidyr,
    data.table[fread]
)

# data: data/GSE116256_RAW/GSM3587923_AML1012-D0.dem.txt.gz
#  ann: data/GSE116256_RAW/GSM3587924_AML1012-D0.anno.txt.gz
files <- list.files(snakemake@input[[1]], full.names = TRUE)

data <- list(
    ann = files[stringr$str_detect(files, "anno")],
    gene = files[stringr$str_detect(files, "dem")]
)

data[["name"]] <- data$gene |>
    basename() |>
    stringr$str_remove("GSM[[:digit:]]{7}_") |>
    stringr$str_remove("\\.[[:lower:]|\\.]+")

data <- purrr$list_transpose(data)

seurat <- list()
for (i in seq_along(data)) {
    ann <- as.data.frame(fread(data[[i]][["ann"]]))
    rownames(ann) <- ann[ , "Cell", drop = TRUE]
    ann[["sample_id"]] <- data[[i]][["name"]]

    x <- fread(data[[i]][["gene"]])
    genes <- x[["Gene"]]
    x[ , "Gene"] <- NULL
    x <- as.matrix(x)
    rownames(x) <- genes

    seurat[[i]] <- suppressWarnings(so$CreateSeuratObject(counts = x, meta.data = ann))
}

seurat <- merge(seurat[[1]], seurat[2:length(seurat)])
seurat <- ChangeFeatureNames(seurat, old = "MGEA5", new = "OGA")
progenetor_cells <- c("GMP", "HSC", "Prog")

md <- seurat[[]] |>
    tidyr$separate(orig.ident, c("pat", "day"), sep = "-") |>
    dplyr$mutate(
        sample_origin = dplyr$case_when(
            stringr$str_detect(sample_id,"^BM") ~ "Healthy",
            stringr$str_detect(sample_id,"^AML") ~ "AML",
            .default = "cell_line"
        ),
        CellType2 = dplyr$if_else(sample_origin == "Healthy", sample_origin, CellType),
        collapsed_cell_types = dplyr$case_when(
            CellType %in% progenetor_cells & sample_origin == "Healthy" ~ "Myeloid Progenitor (HD)",
            CellType %in% progenetor_cells ~ "Myeloid Progenitor (TME)",
            CellType %in% paste0(progenetor_cells, "-like") ~ "LSC-like",
            PredictionRefined == "malignant" ~ "AML Blast",
            sample_origin == "Healthy" ~ "Bone Marrow (HD)",
            .default = "Bone Marrow (TME)"
        ),
        day = dplyr$if_else(day == "D0", "Diagnostic", "Post-treatment"),
        day = dplyr$if_else(sample_origin == "Healthy", NA_character_, day),
        .keep = "none"
    )

seurat <- Seurat$AddMetaData(seurat, md)
seurat <- subset(seurat, subset = sample_origin != "cell_line")
saveRDS(seurat, snakemake@output[["seurat"]])