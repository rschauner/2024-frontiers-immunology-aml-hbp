#!/bin/bash

cd data || exit

echo "Downloading RNA Data from the NIH"
echo "This will take some time"
mkdir "nih_counts"
mkdir "clinical_info"
gdc-client download -m "gdc_rna_manifest.txt" --dir nih_counts > /dev/null
gdc-client download -m "tcga_target_clinical_data_manifest.txt" --dir clinical_info > /dev/null

echo "Downloading data from GTEx"
mkdir "GTEx"
curl -o \
    GTEx/gene_reads_2017-06-05_v8_whole_blood.gct.gz \
    https://storage.googleapis.com/gtex_analysis_v8/rna_seq_gene_reads/gene_reads_2017-06-05_v8_whole_blood.gct.gz

# TODO: Downloading St. Jude data

echo "Downloading GSE198919"
curl -O https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE198919&format=file&file=GSE198919%5FRaw%5FCounts%5Ftable%5FRNAseq%2Etxt%2Egz

echo "Downloading GSE116256"
curl -O https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE116256&format=file
tar -xf GSE116256_RAW.tar

echo "Downloading GSE126068"
mkdir "GSE126068"

for ((  i = 1;  i <= 5;  i++ )); do
    curl -o \
    "GSE126068/GSE126068_Patient${i}_RawCounts.csv.gz" \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126068&format=file&file=GSE126068%5FPatient${i}%5FRawCounts.csv%2Egz"
done

chdir ..