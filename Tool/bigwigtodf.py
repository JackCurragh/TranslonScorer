"""This script contains functions to  and calculate the transcriptomic coordinates"""

import polars as pl
import pyBigWig as bw

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


def bigwigtodf(bigwig, exon):
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
    bwfile = bw.open(bigwig)
    chromcheck = list(bwfile.chroms().keys())[0]
    tran_id = []
    tran_start = []
    tran_stop = []
    counts = []
    if bwfile.isBigWig():
        exon_df = pl.read_csv(exon, has_header=True, separator=",")
        #Function commented since it doesn't work
        #exon_not = exon_df.get_column("chr").to_list()[0]
        #checkchromnot(int(chromcheck), exon_not)
        exon_df = exon_df.with_columns(
            pl.col("start", "stop", "tran_start", "tran_stop").apply(
                lambda x: x.split(",")
            )
        )
        for i in range(len(exon_df)):
            if exon_df["chr"][i] == "M":
                print("M chrom")
            else: 
                exon_pairs = zip(exon_df["start"][i], exon_df["stop"][i])
                transcript_pairs = list(
                    zip(exon_df["tran_start"][i], exon_df["tran_stop"][i])
                )
                for idx, exon in enumerate(exon_pairs):
                    if int(exon[0]) == int(exon[1]):
                        print("skipped")
                    else:
                        intervals = bwfile.intervals(
                            str(exon_df["chr"][i]), int(exon[0]), int(exon[1])
                        )
                    # Filter out "None" type intervals
                        if intervals is not None:
                            for interval in intervals:
                                tran_id.append(exon_df["tran_id"][i])
                                # Tran start coordinate
                                diff_start = interval[0] - int(exon[0])
                                transtart = int(transcript_pairs[idx][0]) + diff_start
                                tran_start.append(transtart)
                                # Tran stop coordinate
                                diff_stop = int(exon[1]) - interval[1]
                                transtop = int(transcript_pairs[idx][1]) - diff_stop
                                tran_stop.append(transtop)
                                # Counts
                                counts.append(interval[2])
        df_tran = pl.DataFrame(tran_start).rename({"column_0": "tran_start"})
        df_tran = df_tran.with_columns(
            (pl.Series(tran_stop)).alias("tran_stop"),
            (pl.Series(counts)).alias("counts"),
            (pl.Series(tran_id)).alias("tran_id"),
        )
        df_tran = df_tran.select(["tran_id", "tran_start", "tran_stop", "counts"])
        return df_tran
    else:
        raise Exception("Must provide a bigwig file to convert")
