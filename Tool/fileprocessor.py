'''This script contains functions to transform data frames into different file types'''
import polars as pl
import os

def asitecalc(df):
    '''
    docstring
    '''
    #GROUP TO CALCULATE A-SITE
    y = []
    for value, data in df.group_by("length"):
        x = df.filter(pl.col("length") == value).with_columns((pl.col("pos") + offsets[value]).alias("A-site"))
        y.append(x)
    df_asite = pl.concat(y)

    #GROUP ON A-SITE
    df_final = df_asite.group_by("chr", "A-site").agg(pl.col("count").sum())
    df_tobed = df_final.with_columns((pl.col("A-site") + 1).alias("end"))
    df_bed = df_tobed.sort(["chr","A-site"]).select(["chr", "A-site", "end", "count"])
    
    return (df_bed)

def dftobed(df):
    '''
    docstring
    '''
    df_filter = df.select(pl.col("qname","rname","pos","end"))
    df_filtered = df.with_columns((pl.col("end") - pl.col("pos")).alias("length")).select(pl.all().exclude("end"))
    
    
    #split specifically on _x -> implement more!!!
    split_func = lambda s: int(s.split("_x")[1])
    df_namesplit = df_filtered.with_columns(
            pl.col("qname").apply(split_func)).rename({"qname":"count", "rname":"chr"})
    
    #OFFSETS -> DIFFERENT FUNCTION!!
    length_values = df_namesplit.get_column("length")
    offsets = {
            i: 15 for i in length_values
    }
    
    bed = asitecalc(df_namesplit)

    return bed

def bedtobigwig(file, chromsize):
    '''
    docstring
    '''
    os.system(f"bedGraphToBigWig {file} {chromsize} file.bw")