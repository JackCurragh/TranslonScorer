import polars as pl

def extract_transcript_id(attr_str):
    '''
    docstring
    '''
    for attr in attr_str.split(";"):
        if attr.startswith("Parent=transcript:") \
                or attr.startswith("ID=transcript:"):
            return attr.split(":")[1]
        elif attr.startswith("transcript_id="):
            return attr.split("=")[1]
        elif attr.startswith(" transcript_id "):
            return attr.split(" ")[2].replace('"', "")
    return pl.nan

def getexons_and_cds(annotation_file):
    '''
    docstring
    '''
    # Read the annotation file into a DataFrame
    df = pl.read_csv(annotation_file, separator='\t', ignore_errors=True,
                     has_header=False, truncate_ragged_lines=True,
                     comment_prefix = "#").select(["column_1", "column_2", "column_3","column_4", "column_5", "column_9"]).rename({"column_1":"chr",
                     "column_2":"source", "column_3":"type", "column_4":"start", "column_5":"stop", "column_9":"attributes"})
    # Filter to retain only coding regions (CDS)
    coding_regions = df.filter((pl.col("type") == "CDS"))
    
    codinginfo = coding_regions.with_columns(pl.col("attributes")
                                    .apply(lambda attributes: extract_transcript_id(attributes)
                                    ).alias("tran_id")).select(pl.all().exclude("attributes", "source"))

    groupedcds = codinginfo.group_by("tran_id").agg(
                 pl.col("start").alias("cds_start"),
                 pl.col("stop").alias("cds_stop"))                                 
    #Getting exons
    exon_regions = df.filter((pl.col("type") == "exon"))

    exoninfo = exon_regions.with_columns(pl.col("attributes")
                                    .apply(lambda attributes: extract_transcript_id(attributes)
                                    ).alias("tran_id")).select(pl.all().exclude("attributes", "source"))

    groupedexons = exoninfo.group_by("tran_id").agg(
                 pl.col("start").alias("exon_start"),
                 pl.col("stop").alias("exon_stop"))

    exon_coords = exontranscriptcoords(groupedexons, "exon_start", "exon_stop")
    
    cds_coords = cdstranscriptcoords(groupedcds, exon_coords)
    print(cds_coords)
    
    return cds_coords

def exontranscriptcoords(df: pl.DataFrame, start_column: str, end_column: str) -> pl.DataFrame:
    '''
    docstring
    '''
    # Initialize new columns
    new_column_1 = []
    new_column_2 = []
    # Iterate over rows
    for i in range(len(df)):
        start_values = df[start_column][i]
        end_values = df[end_column][i]
        new_start_values = []  # Starting value is 0
        new_stop_values = []
        for j in range(len(start_values)):
            if j == 0:
                new_start = 0
            else:
                new_start = new_stop_values[j-1] + 1
            # Calculate stop coordinate
            stop_coordinate = end_values[j] - start_values[j]
            new_start_values.append(new_start)
            new_stop_values.append(new_start + stop_coordinate)
        new_column_1.append(new_start_values)
        new_column_2.append(new_stop_values)
    # Add new columns to the dataframe
    df = df.with_columns((pl.Series(new_column_1)).alias("tran_coord_start"))
    df = df.with_columns((pl.Series(new_column_2)).alias("tran_coord_stop"))
    return df
             

def cdstranscriptcoords(cds_df, exon_df):
    '''
    docstring
    '''
    # Initialize new columns to store calculated transcript coordinates
    coord_start_list = []
    coord_stop_list = []

    # Iterate over each row in the cds DataFrame
    for i in range(len(cds_df)):
        # Get the transcript present in the current row of the cds DataFrame
        transcript_id = cds_df['tran_id'][i]

        # Find the corresponding row in the exon DataFrame with the same transcript_id
        exon_row = exon_df.filter(exon_df['tran_id'] == transcript_id).select(pl.all())
        
        if exon_row is not None:
            # Get start and stop values from cds DataFrame for the current row
            cds_start_values = cds_df['cds_start'].gather(i)
            cds_stop_values = cds_df['cds_stop'].gather(i)
            print(cds_start_values)
            # Get corresponding exon start and stop values from exon DataFrame
            exon_start_values = exon_row['exon_start']
            exon_stop_values = exon_row['exon_stop']
            print(exon_start_values)
            # Calculate the differences between cds_start and exon_start, and cds_stop and exon_stop for each value
            start_diff = [cds_start - exon_start for cds_start, exon_start in zip(cds_start_values, exon_start_values)]
            stop_diff = [cds_stop - exon_stop for cds_stop, exon_stop in zip(cds_stop_values, exon_stop_values)]

            # Calculate transcript coordinates by adding differences to exon transcript coordinates for each value
            coord_start_values = [tran_start + diff for tran_start, diff in zip(exon_row['tran_coord_start'], start_diff)]
            coord_stop_values = [tran_stop + diff for tran_stop, diff in zip(exon_row['tran_coord_stop'], stop_diff)]

            coord_start_list.append(coord_start_values)
            coord_stop_list.append(coord_stop_values)
        else:
            # If the transcript is not found in the exon DataFrame, append empty lists
            coord_start_list.append([])
            coord_stop_list.append([])

    # Add the transcript coordinates lists as new columns to the cds DataFrame
    cds_df = cds_df.drop('coord_start')
    cds_df = cds_df.drop('coord_stop')
    cds_df = cds_df.with_columns((pl.Series(coord_start_list)).alias("coord_start"))
    cds_df = cds_df.with_columns((pl.Series(coord_stop_list)).alias("coord_stop"))

    return cds_df






   





