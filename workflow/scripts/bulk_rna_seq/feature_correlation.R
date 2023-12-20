
options(box.path = here::here())
box::use(
    tibble,
    stringr,
    broom,
    stats,
    readr,
    dplyr,
    purrr,
    ggplot2,
    tidyr,
    forcats,
    Cairo[CairoPDF],
    limma,
    lib = workflow/scripts/bulk_rna_seq/contrasts
)

dge <- readRDS(snakemake@input[["dge"]])
v <- limma$voom(dge, model.matrix(~ 0 + dataset, data = dge$samples), plot = FALSE)
data <- v$E[c("OGA", "OGT"), ] |> t() |> tibble$as_tibble()
data$group <- dge$samples$simple_disease

data <- data |>
    dplyr$filter(group %in% c("AML", "Whole_Blood")) |>
    dplyr$mutate(
        group = group |> forcats$as_factor() |> forcats$fct_relevel("Whole_Blood", "AML")
    ) |>
    dplyr$group_by(group)

data |>
    dplyr$summarise(
        cov = cov(OGA, OGT),
        cor = list(broom$tidy(cor.test(OGA, OGT)))
    ) |>
    tidyr$unnest(cor) |>
    readr$write_tsv(file = snakemake@output[["cor"]])

lm(OGT ~ OGA * group, data = data) |>
    broom$tidy() |>
    readr$write_tsv(file = snakemake@output[["lm"]])

CairoPDF(snakemake@output[["pdf"]])
ggplot2$ggplot(data = data, mapping = ggplot2$aes(x = OGA, y = OGT, color = group)) +
    ggplot2$geom_point(alpha = 0.2) +
    ggplot2$theme_bw() +
    ggplot2$theme(aspect.ratio = 1) +
    ggplot2$scale_color_manual(values = lib$AnnotationColors()$simple_disease)
graphics.off()

CairoPDF(snakemake@output[["pdf_lm"]])
ggplot2$ggplot(data = data, mapping = ggplot2$aes(x = OGA, y = OGT, color = group)) +
    ggplot2$geom_point(alpha = 0.2) +
    ggplot2$theme_bw() +
    ggplot2$scale_color_manual(values = lib$AnnotationColors()$simple_disease) +
    ggplot2$theme(aspect.ratio = 1) +
    ggplot2$geom_smooth(method = "lm")
graphics.off()