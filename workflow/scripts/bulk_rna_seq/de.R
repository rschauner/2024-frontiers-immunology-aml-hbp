options(box.path = here::here())

box::use(
    limma,
    data.table[fwrite],
    glue[glue],
    lib = workflow/scripts/bulk_rna_seq/contrasts
)

dge <- readRDS(snakemake@input[["dge"]])

model <- lib$GetContrasts(snakemake@params[["contrast"]], dge)

pdf(snakemake@output[["voom"]])
v <- limma$voom(dge, model[[1]], plot = TRUE)
graphics.off()

vfit <- limma$lmFit(v, model[[1]])

dir.create(snakemake@output[['de']])
for (coef in colnames(model[[2]])) {
    tsv_file <- glue("{snakemake@output[['de']]}/{coef}.tsv")
    pdf_file <- glue("{snakemake@output[['de']]}/{coef}.pdf")

    cfit <- limma$contrasts.fit(vfit, contrasts = model[[2]])
    cfit <- limma$eBayes(cfit)

    top_table <- limma$topTable(cfit, coef = coef, n = Inf)
    top_table$gene <- rownames(top_table)
    fwrite(top_table, tsv_file)

    pdf(file = pdf_file)
    par(mfrow=c(1,2))
    limma$plotMA(cfit, coef = coef)
    limma$volcanoplot(cfit, coef = coef)
    graphics.off()
}
