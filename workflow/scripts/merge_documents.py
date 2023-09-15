import pathlib
from typing import Callable, List
import pandas as pd
import click
import PyPDF2


def get_import_method_from_extension(file_type: str) -> Callable:
    match file_type:
        case ".csv":
            return pd.read_csv
        case ".tsv":
            return lambda f: pd.read_csv(f, sep="\t")
        case ".pickle":
            return pd.read_pickle
        case ".feather":
            return pd.read_feather
        case ".parquet":
            return pd.read_parquet
        case ".txt":
            raise ValueError(
                "File Types that end in txt are ambiguous. "
                "Please provide --file-type or use a different file extension"
            )
        case _:
            raise ValueError(f"{file_type} is an unknown file extension")


@click.command(help="Merge multiple flat table files into one")
@click.argument(
    "files",
    nargs=-1,
    type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path),
)
@click.option("--output", help="file to write outputs", type=click.File("wb"))
@click.option(
    "--file-type",
    help="File type for input files. This will be automatically detected per file but can be set",
    type=click.Choice(
        ["csv", "tsv", "pickle", "feather", "parquet"], case_sensitive=False
    ),
)
def dataframe(files: List[pathlib.Path], output: click.File, file_type: str) -> None:
    files_work = [f for f in files if f.stat().st_size != 0]

    if len(files_work) == 0:
        with open(output, mode="x") as f:
            f.close()
    else:
        if file_type:
            ftypes = [f".{file_type}"] * len(files_work)
        else:
            ftypes = [x.suffix for x in files_work]

        import_functions = [get_import_method_from_extension(x) for x in ftypes]
        dfs = []

        for f, fun in zip(files_work, import_functions, strict=True):
            dfs.append(fun(f))

        res: pd.DataFrame = pd.concat([x for x in dfs if x is not None])
        res.to_csv(output, sep="\t", index=False)


@click.command(help="Merge multiple PDF files into one")
@click.argument(
    "file",
    nargs=-1,
    type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path),
)
@click.option(
    "--output",
    help="File to output to",
    type=click.Path(exists=False, dir_okay=False, path_type=pathlib.Path),
)
def pdf(file, output) -> None:
    pdfs: List[click.Path] = [f for f in file if f.stat().st_size != 0]

    # if len(pdfs) == 0:
    #     with open(output, mode="x") as f:
    #         f.close()
    # else:
    #     merger = PdfFileMerger()
    #     for f in pdfs:
    #         merger.append(f)
    #     merger.write(output)
    #     merger.close()

    if len(pdfs) == 0:
        with output.open(mode="x") as f:
            writer = PyPDF2.PdfWriter()
            writer.add_page(
                PyPDF2.PageObject.create_blank_page(width=7 * 72, height=7 * 72)
            )

    else:
        with PyPDF2.PdfFileMerger() as merger:
            for f in pdfs:
                merger.append(f)
            merger.write(output)


@click.group(help="Merge flat tables or pdfs into one document")
def cli():
    pass


cli.add_command(dataframe)
cli.add_command(pdf)

if __name__ == "__main__":
    cli()
