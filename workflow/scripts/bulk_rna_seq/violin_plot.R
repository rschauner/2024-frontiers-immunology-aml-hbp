

MakeViolinPlot <- function(object, genes, contrast_mat = NULL, contrast = NULL, palette = NULL, diseases = NULL, viz_age = FALSE) {
    box::use(
        data.table[...],
        ggplot2
    )
    if (!is.null(contrast_mat)) {
        contrast <- ifelse(is.null(contrast), 1, contrast)
        cols_use <- rownames(contrast_mat)[contrast_mat[, contrast] != 0]
        design <- object$design[, colnames(object$design) %in% cols_use]
        samples_use <- rownames(design)[rowSums(design) != 0]
        object <- object[, colnames(object) %in% samples_use]
    }

    if ("EList" %in% class(object)) {
        x <- object$E
        ann <- object$targets
    } else if ("DGEList" %in% class(object)) {
        x <- object$counts
        ann <- object$samples
    } else {
        return(invisible(NULL))
    }

    x <- t(x[genes, ])
    x <- as.data.table(x, keep.rownames = ".id")
    ann <- as.data.table(ann)
    dt <- melt(x[ann, on = ".id"], measure.vars = genes)

    if (!is.null(diseases)) {
        dt <- dt[simple_disease %in% diseases, ]
        if (viz_age) {
            dt[simple_disease == "AML", simple_disease := age_group]
        }
    }

    ggplot2$ggplot(data = dt, mapping = ggplot2$aes(x = simple_disease, y = value, fill = simple_disease)) +
    ggplot2$geom_violin() +
    ggplot2$facet_wrap(~ variable, scales = "free_y", ncol = 4) +
    ggplot2$theme_classic() +
    ggplot2$geom_jitter(size = 0.2) +
    ggplot2$stat_summary(fun.data = stat_n_values, geom = "text", hjust = 0.5, vjust = 0.9) +
    ggplot2$labs(x = NULL, y = "Normalized Expression", legend = NULL) +
    ggplot2$scale_fill_manual(values = palette)
}

stat_n_values <- function(y) {
    box::use(
        glue[glue]
    )
    y_max <- max(y)
    return(
        data.frame(
            y = 1.1 * y_max,
            label = glue("n = {length(y)}")
        )
    )
}

options(box.path = here::here())

box::use(lib = workflow/scripts/bulk_rna_seq/contrasts[AnnotationColors])

v <- readRDS(snakemake@input[["dge"]])

hbp_genes <- c("NAGK", "GNPNAT1", "PGM3", "UAP1", "OGT", "OGA", "GFPT2", "GFPT1")

pdf(snakemake@output[["pdf"]])
MakeViolinPlot(v, genes = hbp_genes, palette = AnnotationColors()$simple_disease, diseases = c("AML", "Whole_Blood"))
MakeViolinPlot(v, genes = hbp_genes, palette = AnnotationColors()$simple_disease, diseases = c("AML", "Whole_Blood"), viz_age = TRUE)
graphics.off()