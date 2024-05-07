'''This script contains functions to  and calculate the transcriptomic coordinates'''
import polars as pl
import pyBigWig as bw


from findexonscds import gettranscriptcoords

def checkchromnot(bigwig, annot):
    '''
    This function checks whether the chromosome notation of the Bigwig file is identical to the notation used in the annotation file
    '''
    if bigwig == annot:
        return True
    else:
        raise Exception("Chromosome notation in Bigwig file has to be identical to the notation in the annotation file!")

def bigwigtodf(bigwig, exon):
    '''
    docstring
    '''
    bwfile = bw.open(bigwig)
    chromcheck = list(bwfile.chroms().keys())[0]
    tran_start = []
    tran_stop = []
    if bwfile.isBigWig():
        exon_df = pl.read_csv(exon, has_header= True, separator = ",")
        exon_not = exon_df.get_column("chr").to_list()[0]
        checkchromnot(int(chromcheck), exon_not)
        exon_df = exon_df.with_columns(pl.col("start", "stop", "tran_start", "tran_stop").apply(lambda x: x.split(",")))
        for i in range(len(exon_df)):
            exon_pairs = zip(exon_df["start"][i], exon_df["stop"][i])
            transcript_pairs = list(zip(exon_df["tran_start"][i], exon_df["tran_stop"][i]))
            for idx, exon in enumerate(exon_pairs):
                intervals = bwfile.intervals(str(exon_df["chr"][i]), int(exon[0]), int(exon[1]))
                if intervals is not None:
                    for interval in intervals:
                        print(exon, transcript_pairs[idx], interval)
                        diff_start = interval[0] - int(exon[0])
                        transtart = int(transcript_pairs[idx][0]) + diff_start
                        tran_start.append(transtart)
                        diff_stop = int(exon[1]) - interval[1]
                        transtop = int(transcript_pairs[idx][1]) - diff_stop
                        tran_stop.append(transtop)
        df_tran = pl.DataFrame(tran_start).rename({"column_0":"tran_start"})
        df_tran = df_tran.with_columns((pl.Series(tran_stop)).alias("tran_stop"))
        print(exon, df_tran)
        
        

    else: raise Exception("Must provide a bigwig file to convert")