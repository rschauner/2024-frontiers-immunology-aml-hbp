# Single Cell and Bulk Expression Analyses Identifies Enhanced Hexosamine Biosynthetic Pathway and O-GlcNAcylation in Acute Myeloid Leukemia

## Running Analysis

The workflow uses [Snakemake](https://snakemake.github.io/) and [mamba](https://mamba.readthedocs.io/en/latest/index.html). 
To get started, install mamba by running the following code in the terminal (MacOS or Linux only):

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh

mamba env create -f enviornment.yml
```

*MacOS: The analysis requires the Intel versions of R and R packages which will be handled by the `enviornment.yml` file.*

After [downloading the data](##data-download), run the analysis with the following code:
```bash
snakemake \
    --use-conda
    --cores 4 # can use as many as you would like
```

## Data Download

### Data from St. Jude

St. Jude Cloud was used to acquire the RNA sequencing data (accession numbers: SJC-DS-1013 and SJC-DS-1009). 
All feature count files for diagnostic samples were collected. 
You must request access to these files from [St. Jude Cloud](https://platform.stjude.cloud/data/cohorts). 

### All other data

In the conda env defined at `enviornment.yml`, run the following command: 

```bash
chmod +x dataset_download.sh
./dataset_download.sh
```
