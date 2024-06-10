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
    if exon_df["chr"][0] == "M":
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
    '''
    docstring
    '''
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
    '''
    docstring
    '''
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
    '''
    docstring
    '''
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
    '''
    docstring
    '''
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
    '''
    docstring
    '''
    if typeorf == 'uoORF':
        df = df.with_columns(
                (pl.col('start')
                .apply(lambda x: scoredict['rise_up'][x])
                .alias('rise_up')))
    
    elif typeorf == 'doORF':
        df = df.with_columns(
                (pl.col('stop')
                .apply(lambda x: scoredict['step_down'][x])
                .alias('step_down'))
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
            scoredict = {'rise_up':{},'step_down':{}}
            
            if counter % 1000 == 0:
                print('\r' + f'{counter} transcripts done', end='')
            
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
        orfscores_df = pl.DataFrame(orfscores)
        print(orfscores_df)
        return orfscores_df
    
    else:
        raise Exception("Must provide a bigwig file to convert")
