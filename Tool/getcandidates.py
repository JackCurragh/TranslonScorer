import pybedtools as pybed
import shutil
from Bio import SeqIO
import polars as pl


from orffinder import find_orfs

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

##def findcoordinates(record_id, annotation):
    '''
    docstring
    '''
#annot = pybed.BedTool(annotation)









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
                if orfs:
                    for start, stop, length, orf in orfs:
                        df_toadd = pl.DataFrame({"pos": [start], "end": [stop], "length": [length], "startorf": [orf[:3]], "stoporf": [orf[-3:]]})
                        df = pl.concat([df, df_toadd])
            counter = counter + 1
        print(df)
#address problem with data frame! position|end and length do not get along!