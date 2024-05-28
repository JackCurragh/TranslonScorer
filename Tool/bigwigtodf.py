"""This script contains functions to  and calculate the transcriptomic coordinates"""

import polars as pl
import pyBigWig as bw

from scoring import sru_score, calculate_scores
def checkchromnot(bigwig, annot):
    """
    Checks if the chromosome notation in a BigWig file is identical to the notation in an annotation file.

    Parameters:
    - bigwig (str): Path to the BigWig file.
    - annot (str): Path to the annotation file.

    Returns:
    - True if the chromosome notation is identical.
    
    Raises:
    - Exception: if the chromosome notation in the BigWig file is not identical to the notation in the annotation file.

    This function compares the chromosome notation in the BigWig file with the notation in the annotation file.
    If the notations are identical, it returns True. Otherwise, it raises an exception with an informative error message.
    """
    if bigwig == annot:
        return True
    else:
        raise Exception(
            "Chromosome notation in Bigwig file has to be identical to the notation in the annotation file!"
        )


def transcriptreads(bwfile, exon_df):
    """
    Converts a BigWig file to a DataFrame based on provided exon annotation.

    Parameters:
    - bigwig (str): Path to the BigWig file.
    - exon (str): Path to the exon annotation file.

    Returns:
    - df_tran (DataFrame): DataFrame containing transcript information derived from the BigWig file.

    Raises:
    - Exception: if a BigWig file is not provided.

    This function reads a BigWig file and an exon annotation file. It performs various operations to extract transcript information
    from the BigWig file based on the exon coordinates. The resulting transcript information is stored in a DataFrame named `df_tran`.
    The DataFrame includes columns for transcript ID, transcript start and stop coordinates, and counts.

    The function first checks if the given file is a BigWig file. It then reads the exon annotation file and extracts necessary
    information, such as the chromosome notation. It ensures that the chromosome notation in the BigWig file matches the exon annotation.
    Then, it iterates over each exon in the annotation file and retrieves intervals from the BigWig file that correspond to the exon
    coordinates. It calculates the transcript start and stop coordinates for each interval and stores the information in corresponding lists.
    Finally, it constructs the `df_tran` DataFrame using the extracted transcript information and returns it.
    """
    read_list=[]
    if exon_df["chr"][0] == "M":
        return pl.DataFrame()
    exon_pairs = zip(exon_df["start"][0], exon_df["stop"][0])
    transcript_pairs = list(
        zip(exon_df["tran_start"][0], exon_df["tran_stop"][0])
    )
    for idx, exon in enumerate(exon_pairs):
        if int(exon[0]) == int(exon[1]):
            continue
        intervals = bwfile.intervals(
            str(exon_df["chr"][0]), int(exon[0]), int(exon[1])
        )
        # Filter out "None" type intervals
        if intervals is not None:
            for interval in intervals:
                # Tran start coordinate
                diff_start = interval[0] - int(exon[0])
                transtart = int(transcript_pairs[idx][0]) + diff_start
                # Tran stop coordinate
                diff_stop = int(exon[1]) - interval[1]
                transtop = int(transcript_pairs[idx][1]) - diff_stop
                # Counts
                counts = interval[2]
                #Add info to list
                dataframe_dict={
                    'tran_start':transtart,
                    'tran_stop':transtop,
                    'counts':counts
                }
                if dataframe_dict:
                    read_list.append(dataframe_dict)
    if read_list:
        read_df = pl.from_dicts(read_list)
        return read_df
    else:
        return pl.DataFrame()

def scoredf(df, tran_reads, sru_range, typeorf):
    '''
    docstring
    '''
    
    df = df.with_columns([
        pl.struct(['start'])\
        .apply(lambda x: sru_score(x['start'], tran_reads, sru_range,0))\
        .alias('rise_up')
        ])
    df = df.with_columns([
        pl.struct(['stop'])\
        .apply(lambda x: sru_score(x['stop'], tran_reads, sru_range,1))\
        .alias('step_down')
        ])
    df = df.with_columns([
        pl.struct(['start', 'stop'])\
        .apply(lambda x: calculate_scores(x['start'], x['stop'], tran_reads))\
        .alias('list_scores')
        ])
    df = df.with_columns(
        (pl.col('list_scores').apply(lambda x: x[0]).alias('hrf')),
        (pl.col('list_scores').apply(lambda x: x[1]).alias('avg')),
        (pl.col('list_scores').apply(lambda x: x[2]).alias('nzc')))
    
    df_dict = df.with_columns(
        score=pl.sum_horizontal('rise_up','step_down','hrf','avg','nzc'))\
        .select(pl.all().exclude('list_scores')).to_dict(as_series=False)
    return df_dict


def scoring(bigwig, exon, orfs, sru_range=15):
    '''
    docstring
    '''
   
    bwfile = bw.open(bigwig)
    if bwfile.isBigWig():
        exon_df = pl.read_csv(exon, has_header=True, separator=",")
        exon_df = exon_df.with_columns(
            pl.col("start", "stop", "tran_start", "tran_stop").apply(lambda x: x.split(",")),
            pl.col('chr').apply(lambda x: x.split('r')[1]))

        orf_df = pl.read_csv(orfs, has_header=True, separator=',')
        counter=0
        orfscores = []
        for tran in orf_df['tran_id'].unique():
            if counter % 1000 == 0:
                print(f'{counter} transcripts done', sep='\r')
            exons = exon_df.filter(pl.col('tran_id') == tran)
            orfs = orf_df.filter(pl.col('tran_id') == tran)
            tran_reads = transcriptreads(bwfile, exons)
            if not tran_reads.is_empty():
                for typeorf in orfs['type'].unique():
                    orfs_filtered = orfs.filter(pl.col('type') == typeorf)
                    #score orfs based on type
                    orf_dict = scoredf(orfs_filtered, tran_reads, sru_range, typeorf)
                    orfscores.append(orf_dict)
            counter+=1
        orfscores_df = pl.from_dicts(orfscores).explode('tran_id','start','stop','type',
            'length','startorf','stoporf','rise_up',
            'step_down','hrf','avg','nzc','score')
        return orfscores_df
    else:
        raise Exception("Must provide a bigwig file to convert")
