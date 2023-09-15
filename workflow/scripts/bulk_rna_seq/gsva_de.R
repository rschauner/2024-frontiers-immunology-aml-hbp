options(box.path = here::here())
box::use(
    Biobase,
    stringr[str_replace_all, str_remove, str_to_title],
    limma,
    data.table[...],
    glue[glue],
    readxl,
    here[here],
    lib = workflow/scripts/bulk_rna_seq/contrasts
)

DrawAnnotatedGeneSetHeatmap <- function(data, gene_set_df, ann_colors) {
    box::use(
        hm = ComplexHeatmap,
        Biobase[assayData, phenoData],
        circlize,
        stringr[str_replace_all, str_remove, str_to_title]
    )
    FixRowNames <- function(x) {
        x |>
            str_replace_all("_", " ") |>
            str_remove("HALLMARK |GOBP |GOMF |GOCC ") |>
            str_to_title()
    }

    pathways <- gene_set_df$gene_set
    rownames(gene_set_df) <- pathways
    gene_set_df$gene_set <- NULL

    rownames(gene_set_df) <- FixRowNames(pathways)

    x <- assayData(data)$exprs[pathways, ]
    rownames(x) <- FixRowNames(rownames(x))
    cann <- hm$HeatmapAnnotation(
        df = phenoData(data)[, c("age_group", "OGN", "OGT", "OGA", "OGT_OGA_ratio_bin", "OGT_bin", "OGA_bin")]@data,
        col = ann_colors
    )
    if ("adj.P.Val" %in% colnames(gene_set_df)) {
        p_vals <- gene_set_df[["adj.P.Val"]]
        p_vals <- fcase(
            p_vals > 0.05, "",
            p_vals > 0.01, "*",
            p_vals > 0.001, "**",
            default = "***"
        )
        gene_set_df[["adj.P.Val"]] <- NULL
        rann <- hm$HeatmapAnnotation(
            p_val = hm$anno_text(p_vals),
            df = as.data.frame(gene_set_df),
            which = "row"
        )
    } else {
        rann <- hm$HeatmapAnnotation(
            df = as.data.frame(gene_set_df),
            which = "row"
        )
    }


    bound <- c(min = min(x), max = max(x))
    bound_max <- max(abs(bound))
    h <- hm$Heatmap(
        x,
        show_column_names = FALSE,
        show_row_names = TRUE,
        cluster_rows = FALSE,
        name = "Enrichment Score",
        left_annotation = rann,
        top_annotation = cann,
        col = circlize$colorRamp2(c(-bound_max, 0, bound_max), c("blue", "black", "red"))
    )
    hm$draw(h)
}


dge <- readRDS(snakemake@input[["dge"]])
gene_sets_plot <- readxl$read_excel(here("annotated_gene_sets.xlsx"))
setDT(gene_sets_plot)
setkey(gene_sets_plot, gene_set)

model <- lib$GetContrasts(snakemake@params[["contrast"]], dge)

vfit <- limma$lmFit(dge, model[[1]])
md_col_use <- names(attr(model[[1]], "contrasts"))

dir.create(snakemake@output[['de']])

for (coef in colnames(model[[2]])) {
    tsv_file <- glue("{snakemake@output[['de']]}/{coef}.tsv")
    pdf_file <- glue("{snakemake@output[['de']]}/{coef}.pdf")

    cfit <- limma$contrasts.fit(vfit, contrasts = model[[2]])
    cfit <- limma$eBayes(cfit)

    vals_use <- rownames(model[[2]])[model[[2]][ , coef] != 0] |> str_remove(md_col_use)

    hm_data <- dge[, Biobase$phenoData(dge)@data[[md_col_use]] %in% vals_use]

    row_labels <- rownames(Biobase$assayData(hm_data)$exprs) |>
        str_replace_all("_", " ") |>
        str_remove("HALLMARK |GOBP ") |>
        str_to_title()
    names(row_labels) <- rownames(Biobase$assayData(hm_data)$exprs)

    top_table <- limma$topTable(cfit, coef = coef, n = Inf)
    top_table$gene <- rownames(top_table)
    fwrite(top_table, tsv_file)

    setDT(top_table)
    setkey(top_table, "gene")

    gene_sets_plot_p_val <- top_table[gene_sets_plot, list(gene_set, db, category, adj.P.Val)]
    setorder(gene_sets_plot_p_val, category)

    pdf(file = pdf_file, width = 14)
    old_par <- par(mfrow = c(1,2))
    limma$plotMA(cfit, coef = coef)
    limma$volcanoplot(cfit, coef = coef)

    par(old_par)
    DrawAnnotatedGeneSetHeatmap(hm_data, gene_sets_plot_p_val, lib$AnnotationColors())
    graphics.off()
}
