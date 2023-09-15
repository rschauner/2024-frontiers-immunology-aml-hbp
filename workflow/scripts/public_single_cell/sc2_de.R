
box::use(
    Seurat,
    so = SeuratObject[`DefaultAssay<-`],
    dplyr,
    readr,
    stringr,
    tidyr,
    purrr,
    tibble,
    stats
)

seurat <- readRDS(snakemake@input[["seurat"]])

DefaultAssay(seurat) <- assay <- snakemake@params[["assay"]]

cell_types <- c("HSC", "GMP", "ProMono", "Mono")
seurat[["cell_type_origin"]] <- paste(seurat$CellType, seurat$sample_origin, sep = "_")

# Create Argument tibble for CellType and Healthy v AML contrasts
de_args1 <- tibble$tibble(
    CellType = 1,
    cell_type_origin = 1,
    ident1 = cell_types
) |>
    tidyr$pivot_longer(cols = !ident1, names_to = "group") |>
    dplyr$mutate(
        ident2 = dplyr$if_else(group == "CellType", paste0(ident1, "-like"), paste0(ident1, "_AML")),
        cell_type = ident1,
        ident1 = dplyr$if_else(group == "cell_type_origin", paste0(ident1, "_Healthy"), ident1)
    ) |>
    dplyr$select(-value)

# Create Argument tibble for the dollapsed cell types DE
collapsed_cell_types <- unique(seurat$collapsed_cell_types)
de_args2 <- tidyr$expand_grid(ident1 = collapsed_cell_types, ident2 = collapsed_cell_types) |>
    dplyr$filter(ident1 != ident2) |>
    dplyr$mutate(group = "collapsed_cell_types")

dplyr$bind_rows(de_args1, de_args2) |>
    dplyr$mutate(
        res = purrr$pmap(
            list(ident1, ident2, group),
            ~ Seurat$FindMarkers(
                seurat,
                features = snakemake@params[["features"]],
                logfc.threshold = 0,
                min.pct = 0,
                ident.1 = ..1,
                ident.2 = ..2,
                group.by = ..3,
                test.use = "LR"
            ) |>
                tibble$rownames_to_column(var = "gene")
        )
    ) |>
    tidyr$unnest(res) |>
    dplyr$relocate(group, cell_type, .before = 1) |>
    dplyr$group_by(group, ident1, ident2) |>
    dplyr$mutate(
        ident1 = dplyr$if_else(group == "cell_type_origin", "Healthy", ident1),
        ident2 = dplyr$if_else(group == "cell_type_origin", "AML", ident2),
        group = dplyr$case_when(
            group == "cell_type_origin" ~ cell_type,
            group == "CellType" ~ "AML",
            group == "collapsed_cell_types" ~ "Annotated Cell Groups"
        ),
        p_val_adj = stats$p.adjust(p_val, method = "holm")
    ) |>
    dplyr$select(-cell_type) |>
    readr$write_tsv(snakemake@output[["tsv"]])





