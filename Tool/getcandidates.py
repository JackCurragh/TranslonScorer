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

def orfrelativeposition(annotation, df):
    '''
    docstring
    '''
    cds_df = getexons_and_cds(annotation)
    filtered_df = cds_df.select(pl.all().exclude("cds_start", "cds_stop"))

    idcol = df.select(pl.col("tran_id", "pos","end"))
    
    for orfid in idcol["tran_id"].unique():
        cdsregion = filtered_df.filter(pl.col("tran_id") == orfid).select(pl.all())
        orfregion = idcol.filter(pl.col("tran_id") == orfid).select(pl.all())
        if not cdsregion.is_empty():
            cdspair = [cdsregion['cds_tran_start'][0], cdsregion['cds_tran_stop'][0]]
            orfpair = list(zip(orfregion['pos'], orfregion['end']))
            for orf in orfpair:
                if orf[0] < cdspair[0] and orf[1] < cdspair[0]:
                    print(orf, cdspair, "upstream ORF detected")
                elif orf[0] > cdspair[1] and orf[1] > cdspair[1]:
                    print(orf, cdspair, "downstream ORF detected")
                elif orf[0] < cdspair[0] and orf[1] < cdspair[1] and orf[1] > cdspair[0]:
                    print(orf, cdspair, "Overlapping upstream ORF detected")
                elif orf[0] < cdspair[1] and orf[0] > cdspair[0] and orf[1] > cdspair[1]:
                    print(orf, cdspair, "Overlapping downstream ORF detected")
                elif orf[0] >= cdspair[0] and orf[1] <= cdspair[1]:
                    print(orf, cdspair, "Internal ORF detected")
        else: print(f"non coding transcript: {orfid}")      
            

    return idcol


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
