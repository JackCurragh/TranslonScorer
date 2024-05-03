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

    cds_df, exon_coords = getexons_and_cds(annotation, list(df["tran_id"].unique()))
    filtered_df = cds_df.select(pl.all().exclude("cds_start", "cds_stop"))
    
   
    orftype = []
    for orfid in sorted(df["tran_id"].unique()):
        cdsregion = filtered_df.filter(pl.col("tran_id") == orfid).select(pl.all())
        orfregion = df.filter(pl.col("tran_id") == orfid).select(pl.all())
        orfpair = list(zip(orfregion['pos'], orfregion['end']))
        if not cdsregion.is_empty():
            if cdsregion['cds_tran_start'][0] > cdsregion['cds_tran_stop'][0]:
                cdspair = [cdsregion['cds_tran_stop'][0], cdsregion['cds_tran_start'][0]]
            else : cdspair = [cdsregion['cds_tran_start'][0], cdsregion['cds_tran_stop'][0]]
            for orf in orfpair:
                if orf[0] < cdspair[0] and orf[1] < cdspair[0]:
                    orftype.append("uORF")
                elif orf[0] == cdspair[0] and orf[1] == cdspair[1]:
                    orftype.append("CDS")
                elif orf[0] > cdspair[1] and orf[1] > cdspair[1]:
                    orftype.append("dORF")
                elif orf[0] < cdspair[0] and orf[1] < cdspair[1] and orf[1] > cdspair[0]:
                    orftype.append("uoORF")
                elif orf[0] <= cdspair[1] and orf[0] >= cdspair[0] and orf[1] > cdspair[1]:
                    orftype.append("doORF")
                elif orf[0] >= cdspair[0] and orf[1] <= cdspair[1]:
                    orftype.append("iORF")
                elif orf[0] < cdspair[0] and orf[1] > cdspair[1]:
                    orftype.append("eoORF")
                elif orf[0] < cdspair[0] and orf[1] == cdspair[1]:
                    orftype.append("extORF")        
                else: 
                    print("unexpected", orf, cdspair)
                    orftype.append("Unexpected")
        else: 
            for orf in orfpair:
                orftype.append("Non Coding")      
    df = df.with_columns((pl.Series(orftype)).alias("type"))
    return df, exon_coords


def preporfs(transcript, starts, stops, minlength, maxlength):
    '''
    docstring
    '''
    df = pl.DataFrame()
    
    with open(transcript) as handle:
        for record in SeqIO.parse(handle, "fasta"): 
                orfs = find_orfs(record.seq, starts, stops, minlength, maxlength)
                tran_id = str(record.id).split("|")[0]          
                if orfs:
                    for start, stop, length, orf in orfs:
                        df_toadd = pl.DataFrame({"tran_id": [tran_id],"pos": [start], "end": [stop], "length": [length], "startorf": [str(orf[:3])], "stoporf": [str(orf[-3:])]})
                        df = pl.concat([df, df_toadd])
        df = df.sort("tran_id")
        return df
