
import polars as pl
   # Explode exon start and stop lists to individual rows
    exploded_exon_df = exon_df.explode(["start", "stop", "tran_start", "tran_stop"])
    exploded_cds_df = cds_df.explode(['start', 'stop'])
    # Join the DataFrames on the transcript ID
    combined_df = exploded_cds_df.join(exploded_exon_df, on="tran_id", how="inner")
    print(combined_df)
    # Calculate transcript-level start and stop coordinates
    combined_df = combined_df.with_columns([
        (pl.when((pl.col("start") >= pl.col("start_right")) & (pl.col("start") <= pl.col("stop_right")))
         .then(pl.col("tran_start") + (pl.col("start") - pl.col("start_right")))
         .otherwise(None)).alias("tran_start_cd"),
        
        (pl.when((pl.col("stop") >= pl.col("start_right")) & (pl.col("stop") <= pl.col("stop_right")))
         .then(pl.col("tran_stop") - (pl.col("stop_right") - pl.col("stop")))
         .otherwise(None)).alias("tran_stop_cd")
    ])
    print(combined_df)
    # Drop rows with None values in calculated columns
    combined_df = combined_df.drop_nulls(["tran_start_cd", "tran_stop_cd"])
    
    # Select and rename relevant columns
    result_df = combined_df.select([
        pl.col("tran_id"),
        pl.col("start"),
        pl.col("stop"),
        pl.col("tran_start_cd").alias("tran_start"),
        pl.col("tran_stop_cd").alias("tran_stop")
    ])
    print(result_df.filter(pl.col('tran_id') == 'ENST00000439326.8'))
    return result_df