# dbSNP to Parquet Converter

A Snakemake pipeline to download dbSNP data and convert it to Parquet format for efficient querying.

## Overview

This pipeline:
1. Downloads dbSNP data (GRCh38 build 156)
2. Filters for SNVs only
3. Converts chromosome contigs to standard naming
4. Splits data by chromosome 
5. Creates Parquet lookup tables with RSID mappings

## Requirements

- bcftools
- Python with cyvcf2, pandas, pyarrow, polars
- Snakemake

## Usage

1. Configure resources in `config.yaml`
2. Run the pipeline:
   ```bash
   snakemake --cluster "sbatch -p {resources.partition} --mem={resources.mem} -t {resources.time} -c {threads}" -j 50
   ```

## Output

- `output/dbSNP_156.bcf` - Full filtered BCF file
- `output/dbSNP_156.chr*.bcf` - Per-chromosome BCF files  
- `output/dbSNP_156.chr*.lookup.parquet` - Per-chromosome RSID lookup tables