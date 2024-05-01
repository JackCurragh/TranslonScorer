import pybedtools as pybed
import shutil
from Bio import SeqIO
import polars as pl


from orffinder import find_orfs
from findexonscds import getexons_and_cds

def gettranscripts(seq, annotation, outfile="transcripts.fa"):
    '''
    docstring
    '''
    annot = pybed.BedTool(annotation)
    
    fasta = pybed.BedTool(seq)
    transcripts = annot.filter(lambda feature: feature[2] == 'transcript')
    annot = transcripts.sequence(fi = fasta, name = True, s = True)
    
    shutil.copyfile(annot.seqfn, outfile)
    return outfile

def matchcoordinates(annotation, df):
    '''
    docstring
    '''
    cds_df = getexons_and_cds(annotation)

    idcol = df.get_column("tran_id")
    for orfid in idcol:
        exonregion = exon_df.filter(pl.col("tran_id") == orfid)
        print(exonregion)

    return exonregion


def preporfs(transcript, starts, stops, minlength, maxlength):
    '''
    docstring
    '''
    df = pl.DataFrame()
    counter = 0
    with open(transcript) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if counter < 100: 
                orfs = find_orfs(record.seq, starts, stops, minlength, maxlength)
                tran_id = str(record.id).split("|")[0]          
                if orfs:
                    for start, stop, length, orf in orfs:
                        df_toadd = pl.DataFrame({"tran_id": [tran_id],"pos": [start], "end": [stop], "length": [length], "startorf": [orf[:3]], "stoporf": [orf[-3:]]})
                        df = pl.concat([df, df_toadd])
                        counter = counter +1
        return df
