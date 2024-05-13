import pysam
import polars as pl
import oxbow as ox

def readbam(bampath):
    """
    Reads a given BAM file, extracts relevant information, and returns it as a DataFrame.

    Parameters:
    - bampath (str): Path to the BAM file to be processed.

    Returns:
    - df (DataFrame): Polars DataFrame containing the extracted information from the BAM file.

    This function indexes the BAM file using pysam, reads the indexed file using ox.read_bam, 
    and then reads the data into a DataFrame using pl.read_ipc. The DataFrame containing the 
    relevant information extracted from the BAM file is returned for further processing.
    """
    pysam.index(bampath)
    bamfile = ox.read_bam(bampath)
    df = pl.read_ipc(bamfile)
    return df


"""def readofst(ofstpath):
        docstring

        """
