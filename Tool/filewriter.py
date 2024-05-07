'''Script for writing files that need to be stored during analysis'''
import os
import polars as pl

def saveorfsandexons(orf_df, exon_df):
   '''
   docstring
   '''
   if not os.path.exists("data/files"):
            os.mkdir("data/files")
   orf_df.write_csv("data/files/annotated_orfs.csv")

   fixchromcol = exon_df.with_columns(pl.col("chr").apply(lambda x: x[0]))

   exon_df = fixchromcol.with_columns(pl.col("start", "stop", "tran_start", "tran_stop")
                                        .apply(lambda x: ','.join(map(str, x))))
   exon_df.write_csv("data/files/exons.csv")