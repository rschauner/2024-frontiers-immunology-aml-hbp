

ExtractAssayData <- function(seurat) {
    box::use(
        so=SeuratObject
    )

    assay <- so$GetAssayData(seurat, slot = "data")
    assay <- assay[Matrix::rowSums(assay) != 0, ]

    return(
         list(
             matrix = MatrixExtra::t_shallow(assay),
             genes = as.list(rownames(assay))
         )
    )
}

CreateMAGICAssay <- function(data) {
    box::use(
        so = SeuratObject
    )
    data <- t(data)
    so$CreateAssayObject(data = data)
}

seurat <- readRDS(snakemake@input[["seurat"]])

box::use(
    readr[read_csv, write_csv],
    arrow[read_feather, write_feather],
    Matrix[writeMM],
    data.table[fwrite]
)

if (!is.null(snakemake@input[["magic"]])) {
    data <- read_feather(snakemake@input[["magic"]])
    rownames(data) <- colnames(seurat)
    assay <- CreateMAGICAssay(data)
    seurat[["MAGIC"]] <- assay
    saveRDS(seurat, snakemake@output[["seurat"]])
} else {
    out <- ExtractAssayData(seurat)
    writeMM(out[["matrix"]], snakemake@output[["data"]])
    fwrite(out[["genes"]], snakemake@output[["genes"]])
}
