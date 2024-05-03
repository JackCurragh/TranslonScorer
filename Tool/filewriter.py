'''Script for writing files that need to be stored during analysis'''
import os
import polars as pl

def saveorfsandexons(orf_df, exon_df):
    if not os.path.exists("data/files"):
            os.mkdir("data/files")
    orf_df.write_csv("data/files/annotated_orfs.csv")
    exon_df = exon_df.with_columns(pl.col("exon_start", "exon_stop", "tran_coord_start", "tran_coord_stop")
                                        .apply(lambda x: ','.join(map(str, x))))
    exon_df.write_csv("data/files/exons.csv")

