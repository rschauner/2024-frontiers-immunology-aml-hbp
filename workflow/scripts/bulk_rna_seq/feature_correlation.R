

box::use(
    tibble,
    stringr,
    broom,
    stats,
    readr,
    dplyr,
    purrr
)

dge <- readRDS(snakemake@input[["dge"]])
data <- dge$counts[c("OGA", "OGT"), ] |> t() |> as.data.frame()

fns <- list(
    cor = function(x) broom$tidy(cor.test(x[[1]], x[[2]])),
    cov = function(x) stats$cov(x[[1]], x[[2]])
)

purrr$map_dfc(fns, ~ .x(data)) |>
    t() |>
    as.data.frame() |>
    tibble$rownames_to_column() |>
    dplyr$mutate(
        rowname = stringr$str_pad(rowname, max(stringr$str_width(rowname))),
        V1 = stringr$str_pad(V1, max(stringr$str_width(V1)), side = "right")
    ) |>
    readr$write_tsv(file = snakemake@output[[1]], col_names = FALSE)
