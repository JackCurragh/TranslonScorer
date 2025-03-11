import pysam
import polars as pl
import oxbow as ox
from ..utils.logging import log_info

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
    log_info("BAM file indexed successfully")
    bamfile = ox.read_bam(bampath)
    log_info("BAM file read successfully")
    df = pl.read_ipc(bamfile)
    log_info("DataFrame created successfully")
    return df


"""def readofst(ofstpath):
        docstring

        """
