import multiprocessing
import pathlib
import pickle
from typing import Tuple, Callable, Any

import click
import magic
import scipy
import pandas as pd
import numpy as np


def load_magic(file: pathlib.Path) -> magic.MAGIC:
    try:
        with file.open("rb") as f:
            fit = pickle.load(f)

        if not isinstance(fit, magic.MAGIC):
            raise TypeError("Loaded object is not of type magic.MAGIC")

        return fit
    except (FileNotFoundError, IOError) as e:
        click.echo(f"Error: File {file} not found or could not be opened.")
        raise e
    except pickle.PickleError as e:
        click.echo(f"Error: Pickle error while loading file {file}.")
        raise e


def fit_model(df: pd.DataFrame, threads: int = 1, solver: str = "exact") -> magic.MAGIC:
    return magic.MAGIC(n_jobs=threads, solver=solver, t="auto").fit(df)


def transform_model(fit: magic.MAGIC, genes: Tuple[str]) -> pd.DataFrame:
    return fit.transform(genes=genes)


def get_import_method_from_extension(file: pathlib.Path) -> Callable:
    match file.suffix:
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
        case ".mtx":
            return scipy.io.mmread
        case ".txt":
            raise ValueError(
                "File Types that end in txt are ambiguous. "
                "Please use a different file extension"
            )
        case _:
            raise ValueError(f"{file.suffix} is an unknown file extension")


def get_export_method_from_extension(df: pd.DataFrame, file: pathlib.Path) -> None:
    file_type = file.suffix
    match file_type:
        case ".csv":
            df.to_csv(file)
        case ".tsv":
            df.to_csv(file, sep="\t")
        case ".pickle":
            df.to_pickle(file)
        case ".feather":
            df.to_feather(file)
        case ".parquet":
            df.to_parquet(file)
        case ".txt":
            raise ValueError(
                "File Types that end in txt are ambiguous. "
                "Please use a different file extension"
            )
        case _:
            raise ValueError(f"{file_type} is an unknown file extension")


def read_data(file: pathlib.Path):
    import_method = get_import_method_from_extension(file)
    return import_method(file)


def save_model(fit: Any, output: pathlib.Path) -> None:
    output.write_bytes(pickle.dumps(fit))


def get_gene_indices(genes_names, genes):
    return [i for i, name in enumerate(genes_names) if name in genes]


@click.command()
@click.argument(
    "file", type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path)
)
@click.argument(
    "output",
    type=click.Path(
        exists=False, dir_okay=False, path_type=pathlib.Path, writable=True
    ),
)
@click.option(
    "--threads",
    type=click.IntRange(1, multiprocessing.cpu_count(), clamp=True),
    default=1,
    help="Number of threads to use for MAGIC",
)
@click.option(
    "--solver",
    type=click.Choice(["exact", "approximate"]),
    default="exact",
    help="Which method to use during computation",
)
def fit(file: pathlib.Path, output: click.Path, threads: int, solver: str) -> None:
    """
    --threads and --solver will be saved in OUTPUT for calls to transform
    """
    df = read_data(file)
    fit = fit_model(df, threads=threads, solver=solver)
    save_model(fit, output)


@click.command()
@click.argument(
    "file", type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path)
)
@click.argument(
    "output",
    type=click.Path(
        exists=False, dir_okay=False, path_type=pathlib.Path, writable=True
    ),
)
@click.argument("genes", nargs=-1)
@click.option(
    "--genes-path",
    type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path),
    help="Optional if input to MAGIC was a dataframe",
)
def transform(
    file: pathlib.Path, output: click.Path, genes: Tuple[str], genes_path: pathlib.Path
) -> None:
    """
    OUTPUT a csv file

    GENES a list of genes, if empty all genes are used
    """
    fit = load_magic(file)
    genes_use = genes or "all_cells"

    if not isinstance(fit.X, pd.DataFrame) and genes_use != "all_cells":
        if genes_path is None:
            raise click.MissingParameter(
                "--genes-path is required if a matrix input was used and genes have been provided"
            )
        genes_names = read_data(genes_path)
        genes_use = get_gene_indices(genes_names, genes)

    data = transform_model(fit, genes_use)

    if isinstance(data, np.ndarray):
        data = pd.DataFrame(data, columns=list(genes))

    get_export_method_from_extension(data, output)


@click.group(help="Run magic-impute")
@click.version_option("v1.0")
def cli():
    pass


cli.add_command(transform)
cli.add_command(fit)


if __name__ == "__main__":
    cli()
