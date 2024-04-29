import pysam
import polars as pl
import oxbow as ox

def readbam(bampath):
    '''
    docstring
    
    '''
    pysam.index(bampath)     
    bamfile = ox.read_bam(bampath)
    df = pl.read_ipc(bamfile)

    return df

'''def readofst(ofstpath):
        docstring

        '''