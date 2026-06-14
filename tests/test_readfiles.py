import os

import pytest

pysam = pytest.importorskip("pysam")
ox = pytest.importorskip("oxbow")
import polars as pl

# Legacy test built around a local BAM fixture that is not shipped (data/ is
# gitignored). Skip the whole module when the fixture is absent rather than
# invoking readbam() at import time and erroring during collection.
if not os.path.exists("data/SRR11005879.bam"):
    pytest.skip(
        "data/SRR11005879.bam fixture not present", allow_module_level=True
    )


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
    print(df)
    return df

readbam('data/SRR11005879.bam')

def test_answer():
    assert readbam("data/SRR11005879.bam") == 1