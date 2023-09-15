
box::use(edgeR)

dge <- readRDS(snakemake@input[["dge"]])
dge <- edgeR$calcNormFactors(dge)

saveRDS(dge, snakemake@output[["dge"]])