"""This script contains functions to  and calculate the transcriptomic coordinates"""

import polars as pl
import pyBigWig as bw

from scoring import sru_score, calculate_scores

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
    if exon_df["chr"][0] == "chrM":
        return pl.DataFrame()
    exon_pairs = zip(exon_df["start"][0], exon_df["stop"][0])
    transcript_pairs = list(
        zip(exon_df["tran_start"][0], exon_df["tran_stop"][0])
    )
    for idx, exon in enumerate(exon_pairs):
        if int(exon[0]) == int(exon[1]):
            continue
        
        if str(exon_df["chr"][0]) not in bwfile.chroms().keys():
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

def oldscoring(df, tran_reads, sru_range, typeorf):
    """
    Applies specific scoring calculations to a DataFrame based on the type of ORF.

    This function modifies the given DataFrame `df` by adding new columns based on the 
    specified type of ORF (`uoORF`, `doORF`, or other). It calculates scores using 
    provided `tran_reads` and `sru_range`, and then computes a final score.

    Parameters:
    df (pl.DataFrame): The input Data frame containing 'start' and 'stop' columns.
    tran_reads (df): Data frame containing reads on transcriptomic level required for scoring functions.
    sru_range (int): Range parameter required for SRU scoring functions.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.

    Returns:
    dict: A dictionary representation of the modified DataFrame.

    Notes:
    - For 'uoORF', the function calculates the 'rise_up' score based on the 'start' column.
    - For 'doORF', the function calculates the 'step_down' score based on the 'stop' column.
    - For other types, it calculates a list of scores and extracts 'hrf', 'avg', and 'nzc' from it.
    - The final score is a sum of 'rise_up', 'step_down', 'hrf', 'avg', and 'nzc' columns.
    """
    if typeorf == 'uoORF':
        df = df.with_columns([
            pl.struct(['start'])\
            .apply(lambda x: sru_score(x['start'], tran_reads, sru_range,0))\
            .alias('rise_up')
            ])
    elif typeorf == 'doORF':
        df = df.with_columns([
            pl.struct(['stop'])\
            .apply(lambda x: sru_score(x['stop'], tran_reads, sru_range,1))\
            .alias('step_down')
            ])
    else:
        df = df.with_columns((
        pl.struct(['start', 'stop'])\
        .apply(lambda x: calculate_scores(x['start'], x['stop'], tran_reads))\
        .alias('list_scores')
        ))
    df = df.with_columns(
        (pl.col('list_scores').apply(lambda x: x[0]).alias('hrf')),
        (pl.col('list_scores').apply(lambda x: x[1]).alias('avg')),
        (pl.col('list_scores').apply(lambda x: x[2]).alias('nzc')))
    
    df = df.with_columns(
        score=pl.sum_horizontal('rise_up','step_down','hrf','avg','nzc'))\
        .select(pl.all().exclude('list_scores')).to_dict(as_series=False)
    return df

def newscoring(df, tran_reads, sru_range, typeorf, scoredict):
    """
    Updates a scoring dictionary with 'rise_up' and 'step_down' scores based on the DataFrame and ORF type.

    This function calculates and updates the `scoredict` with 'rise_up' and 'step_down' scores
    based on unique 'start' and 'stop' values from the DataFrame `df`. The scores are determined
    using the `sru_score` function, and the type of ORF (`typeorf`) dictates which scores are computed.

    Parameters:
    df (pl.DataFrame or pl.Series): The input DataFrame or Series containing 'start' and 'stop' columns.
    tran_reads (df): Data frame containing reads on transcriptomic level required for scoring functions.
    sru_range (int): Range parameter required for SRU scoring functions.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.
    scoredict (dict): Dictionary to store the computed 'rise_up' and 'step_down' scores.

    Returns:
    dict: The updated scoring dictionary with 'rise_up' and 'step_down' scores.

    Notes:
    - For types other than 'doORF', it calculates 'rise_up' scores based on unique 'start' values.
    - For types other than 'uoORF', it calculates 'step_down' scores based on unique 'stop' values.
    - The `scoredict` is updated with these scores, where keys are unique 'start' or 'stop' values
      and values are the corresponding scores from the `sru_score` function.
    """
    if not typeorf == 'doORF':
        if type(df) != type(pl.Series()):
            startvalues = df.get_column('start')\
                            .unique()
        else: startvalues = df
        startscore = startvalues\
                    .map_elements(lambda x: sru_score(x, tran_reads, sru_range,0))\
                    .to_list()
        for i, value in enumerate(startscore):
            scoredict['rise_up'].update({startvalues[i]:startscore[i]})

    if not typeorf == 'uoORF':
        if type(df) != type(pl.Series()):
            stopvalues = df.get_column('stop')\
                           .unique()
        else: stopvalues = df
        stopscore = stopvalues\
                    .map_elements(lambda x: sru_score(x, tran_reads, sru_range,1))\
                    .to_list()
        for i,value in enumerate(stopscore):
            scoredict['step_down'].update({stopvalues[i]:stopscore[i]})
    return scoredict

def globalscores(df, tran_reads, typeorf):
    """
    Computes global scores for a DataFrame based on 'start' and 'stop' values and the type of ORF.

    This function modifies the given DataFrame `df` by calculating a list of scores from 'start' 
    and 'stop' values and adds columns for specific scores (Highes reading frame/hrf, Average/avg, Non-zero coverage/nzc). 
    It then computes a final score by summing up relevant columns based on the type of ORF (`typeorf`).

    Parameters:
    df (pl.DataFrame): The input DataFrame containing 'start' and 'stop' columns.
    tran_reads (df): Data frame containing reads on transcriptomic level required for scoring functions.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.

    Returns:
    dict: A dictionary representation of the modified DataFrame with computed scores.

    Notes:
    - Computes 'hrf', 'avg', and 'nzc' scores from 'start' and 'stop' columns using `calculate_scores`.
    - For 'doORF', the final score is the sum of 'step_down', 'hrf', 'avg', and 'nzc' columns.
    - For 'uoORF', the final score is the sum of 'rise_up', 'hrf', 'avg', and 'nzc' columns.
    - For other types, the final score is the sum of 'rise_up', 'step_down', 'hrf', 'avg', and 'nzc' columns.
    """
    df = df.with_columns((
        pl.struct(['start', 'stop'])\
        .apply(lambda x: calculate_scores(x['start'], x['stop'], tran_reads))\
        .alias('list_scores')
        ))
    df = df.with_columns(
        (pl.col('list_scores').apply(lambda x: x[0]).alias('hrf')),
        (pl.col('list_scores').apply(lambda x: x[1]).alias('avg')),
        (pl.col('list_scores').apply(lambda x: x[2]).alias('nzc')))\
        .select(pl.all().exclude('list_scores'))
    
    if typeorf == 'doORF':
        df = df.with_columns(
            score=pl.sum_horizontal('step_down','hrf','avg','nzc')
            ).to_dict(as_series=False)
    elif typeorf == 'uoORF':
        df = df.with_columns(
            score=pl.sum_horizontal('rise_up','hrf','avg','nzc')
            ).to_dict(as_series=False)
    else:
        df = df.with_columns(
            score=pl.sum_horizontal('rise_up','step_down','hrf','avg','nzc')
            ).to_dict(as_series=False)
    return df

def existingscore(df, typeorf, scoredict):
    """
    Filters the DataFrame to exclude rows with existing scores in the scoring dictionary.

    This function filters out 'start' and 'stop' values from the DataFrame `df` that already have
    corresponding scores in the `scoredict`. The filtering behavior depends on the type of ORF (`typeorf`).

    Parameters:
    df (pl.DataFrame): The input DataFrame containing 'start' and 'stop' columns.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.
    scoredict (dict): Dictionary containing the existing scores for 'rise_up' and 'step_down'.

    Returns:
    pl.DataFrame or pl.Series: A filtered DataFrame or Series excluding rows with existing scores.

    Notes:
    - For 'uoORF', filters out 'start' values that exist in `scoredict['rise_up']`.
    - For 'doORF', filters out 'stop' values that exist in `scoredict['step_down']`.
    - For other types, filters out rows where either 'start' is in `scoredict['rise_up']` or 
      'stop' is in `scoredict['step_down']`, and retains rows where at least one of these conditions is met.
    """
    if typeorf == 'uoORF':
        df = df['start'].map_elements(lambda x: x if not x in scoredict['rise_up'] else None).drop_nulls()
        return df
    elif typeorf == 'doORF':
        df = df['stop'].map_elements(lambda x: x if not x in scoredict['step_down'] else None).drop_nulls()
        return df
    else:
        df = df.select(['start', 'stop']).with_columns(
            (pl.col('start')
             .apply(lambda x: x if not x in scoredict['rise_up'] else None)
             .alias('in_ru')),
            (pl.col('stop')
             .apply(lambda x: x if not x in scoredict['step_down'] else None)
             .alias('in_sd'))
        )
        df = df.filter((pl.col("in_ru").is_not_null()) | (pl.col("in_sd").is_not_null()))\
               .select(['start','stop'])
        return df

def assigningscore(df, scoredict, typeorf):
    """
    Assigns scores from a dictionary to a DataFrame based on the type of ORF.

    This function updates the DataFrame `df` by assigning scores from the `scoredict` to the 
    'rise_up' and 'step_down' columns based on the 'start' and 'stop' values. The type of ORF (`typeorf`)
    determines which scores are assigned.

    Parameters:
    df (pl.DataFrame): The input DataFrame containing 'start' and 'stop' columns.
    scoredict (dict): Dictionary containing the scores for 'rise_up' and 'step_down'.
    typeorf (str): Type of ORF, can be 'uoORF', 'doORF', or any other value for different processing.

    Returns:
    pl.DataFrame: The modified DataFrame with assigned scores.

    Notes:
    - For 'uoORF', assigns 'rise_up' scores from `scoredict` based on 'start' values and sets 'step_down' to 0.0.
    - For 'doORF', assigns 'step_down' scores from `scoredict` based on 'stop' values and sets 'rise_up' to 0.0.
    - For other types, assigns both 'rise_up' and 'step_down' scores from `scoredict` based on 'start' and 'stop' values.
    """
    if typeorf == 'uoORF':
        df = df.with_columns(
                (pl.col('start')
                .apply(lambda x: scoredict['rise_up'][x])
                .alias('rise_up')),
                (pl.lit(0.0)
                 .alias('step_down')))
    
    elif typeorf == 'doORF':
        df = df.with_columns(
                (pl.col('stop')
                .apply(lambda x: scoredict['step_down'][x])
                .alias('step_down')),
                (pl.lit(0.0)
                 .alias('rise_up'))
                )
    else:
        df = df.with_columns(
                (pl.col('start')
                .apply(lambda x: scoredict['rise_up'][x])
                .alias('rise_up')),
                (pl.col('stop')
                .apply(lambda x: scoredict['step_down'][x])
                .alias('step_down'))
                )
    return df

def scoring(bigwig, exon, orfs, old_scoring, sru_range):
    '''
    Perform scoring on ORFs based on transcript data from a BigWig file and exon annotations.

    Parameters:
    - bigwig (str): Path to the BigWig file containing transcript data.
    - exon (str): Path to the CSV file containing exon annotations.
    - orfs (str): Path to the CSV file containing ORF annotations.
    - old_scoring (bool): Flag indicating whether to use the old scoring method.
    - sru_range (int): Range parameter for scoring.

    Returns:
    - DataFrame: A Pandas DataFrame containing scored ORFs with additional metrics.

    This function reads transcript and ORF annotations from provided files, 
    computes scores based on transcript reads obtained from the BigWig file,
    and optionally applies scoring methods based on the `old_scoring` flag.

    If `old_scoring` is True, the function uses the `oldscoring` function to score ORFs,
    appending the results to `orfscores`. If False, it first checks existing scores 
    using `existingscore`, updates the scores with `newscoring`, and then assigns scores 
    to each ORF using `assigningscore`. Additionally, global scores are computed using 
    `globalscores`.

    The final scored ORFs are returned as a flattened DataFrame (`orfscores_df`) where each 
    row represents an individual ORF with associated scoring metrics.

    Note: This function assumes the existence of helper functions like `transcriptreads`, 
    `oldscoring`, `existingscore`, `newscoring`, `assigningscore`, and `globalscores`.
    '''
    bwfile = bw.open(bigwig)
    if bwfile.isBigWig():
        exon_df = pl.read_csv(exon, has_header=True, separator=",")
        exon_df = exon_df.with_columns(
            pl.col("start", "stop", "tran_start", "tran_stop").apply(lambda x: x.split(",")))
        orf_df = pl.read_csv(orfs, has_header=True, separator=',')

        counter=0
        orfscores = []
        for tran in orf_df['tran_id'].unique():
            scoredict = {'rise_up':{},'step_down':{}}
            
            if counter % 1000 == 0:
                print('\r' + f'{counter} transcripts scored', end='')
            
            exons = exon_df.filter(pl.col('tran_id') == tran)
            orfs = orf_df.filter(pl.col('tran_id') == tran)
            
            if exons.is_empty():
                continue
            
            tran_reads = transcriptreads(bwfile, exons)
            if not tran_reads.is_empty():
                for typeorf in orfs['type'].unique():
                    orfs_filtered = orfs.filter(pl.col('type') == typeorf)
                    
                    if old_scoring:
                        orfs_filtered = oldscoring(orfs_filtered, tran_reads, sru_range, typeorf)
                        orfscores.append(orfs_filtered)
                    else:
                        emptyscore_df = existingscore(orfs_filtered, typeorf, scoredict)
                        if not emptyscore_df.is_empty():
                            scoredict = newscoring(emptyscore_df, tran_reads, sru_range, typeorf, scoredict)
                        #run one apply to get rise up and step down scores from updated dictionary
                        orfs_filtered = assigningscore(orfs_filtered, scoredict, typeorf)
                        orfs_filtered = globalscores(orfs_filtered, tran_reads, typeorf)
                        orfscores.append(orfs_filtered)
            counter+=1
        orfscores_df = pl.from_dicts(orfscores).explode(pl.all())
        print('\n')
        return orfscores_df
    
    else:
        raise Exception("Must provide a bigwig file to convert")
