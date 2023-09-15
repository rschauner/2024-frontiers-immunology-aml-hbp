
box::use(
    purrr[map],
    stringr,
    SeuratObject
)

ChangeFeatureNames <- function(object, old, new) {
    assay <- object[["RNA"]]
    for (slot_use in c("counts", "data", "scale.data", "meta.features")) {
        x <- slot(assay, slot_use)
        rownames(x)[rownames(x) == old] <- new
        slot(assay, slot_use) <- x
    }
    if (old %in% assay@var.features) {
        slot(assay, "var.features") <- stringr$str_replace(slot(assay, "var.features"), old, new)
    }
    object[["RNA"]] <- assay
    return(object)
}

seurat <- readRDS(snakemake@input[["seurat"]])
seurat <- map(seurat, ChangeFeatureNames, old = "MGEA5", new = "OGA")
seurat_merge <- merge(seurat[[1]], seurat[2:5])

seurat_merge[["patient_id"]] <- seurat_merge[["orig.ident"]]

counts <- SeuratObject$GetAssayData(seurat_merge, slot = "counts")
seurat_merge[["RNA"]] <- SeuratObject$CreateAssayObject(counts = counts[!stringr$str_detect(rownames(counts), "^ERCC-"), ])

saveRDS(seurat_merge, snakemake@output[["seurat"]])