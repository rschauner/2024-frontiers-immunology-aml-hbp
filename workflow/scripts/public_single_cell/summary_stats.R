
box::use(
    dplyr,
    readr,
    tidyr,
    tidyselect
)

seurat <- readRDS(snakemake@input[["seurat"]])

seurat[[]] |>
    dplyr$mutate(
        CellType = dplyr$if_else(CellType == "", "Uncalled", CellType)
    ) |>
    dplyr$group_by(CellType, sample_origin, day) |>
    dplyr$summarise(n = dplyr$n()) |>
    tidyr$pivot_wider(names_from = CellType, values_from = n) |>
    dplyr$mutate(
        dplyr$across(
            .cols = tidyselect$where(is.integer),
            .fns = ~ dplyr$if_else(is.na(.x), 0, as.integer(.x))
        ),
        total = rowSums(dplyr$pick(tidyselect$where(is.double)))
    ) |>
    readr$write_tsv(snakemake@output[["tsv"]])