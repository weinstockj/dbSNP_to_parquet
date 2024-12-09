import cyvcf2 as cy 
import pandas as pd 
import pyarrow as pa
import pyarrow.parquet as pq
import polars as pl
import logging
import gc
import sys

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO, stream=sys.stdout)

input_bcf = sys.argv[1]
output_file = sys.argv[2]

count = 0
BATCH_SIZE = 500000

schema = pa.schema([
        pa.field("RSID", pa.string()),
        pa.field("ID", pa.string())
    ])

writer = pq.ParquetWriter(output_file, schema)

for variant in cy.VCF(input_bcf):

    if not len(variant.ALT) >= 1:
        logging.warning(f"Skipping variant {variant.ID} as it has no ALT alleles variant.ALT = {variant.ALT}, variant.end = {variant.end}")
        continue

    row = {
        "RSID"    : variant.ID,
        "ID"   : f"{variant.CHROM}_{variant.end}_{variant.REF}_{variant.ALT[0]}",
    }

    row_df = pl.DataFrame(row)

    if count == 0:

        df = row_df.clone()

    if (count % BATCH_SIZE == 0) & (count > 0):

        df.vstack(row_df, in_place=True)

        batch = pa.RecordBatch.from_pandas(df.to_pandas())

        logging.info(f"At count = {count}, now writing batch to {output_file}")

        writer.write_batch(batch)

        df = row_df.clone() # empty dataframe
    else:

        df.vstack(row_df, in_place=True)

    count += 1

    if count % 100000 == 0:
        logging.info(f"Processed {count} variants")
        gc.collect()

batch = pa.RecordBatch.from_pandas(df.to_pandas())
logging.info(f"At count = {count}, now writing last batch to {output_file}")

writer.write_batch(batch)

writer.close()
