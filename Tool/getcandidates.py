from Bio import SeqIO
import pyranges as pr
import polars as pl
import ahocorasick

from orffinder import find_orfs
from findexonscds import getexons_and_cds


def gettranscripts(seq, annotation, outfile="data/transcripts.fa"):
    """
    Extracts transcript sequences from a genome annotation file.

    This function takes a genome sequence file and a genome annotation file,
    extracts transcript sequences from the annotation file, and saves them
    to a FASTA file.

    Parameters:
        seq (str): Path to the genome sequence file in FASTA format.
        annotation (str): Path to the genome annotation file in BED/GFF/GTF format.
        outfile (str): Path to save the output transcript sequences in FASTA format.
                       Default is 'transcripts.fa'.

    Returns:
        str: Path to the output FASTA file containing transcript sequences.

    Example:
        output_file = gettranscripts("genome.fa", "annotation.gff", outfile="transcripts.fa")
    """
    ann = pr.read_gtf(annotation)
    transcripts = ann[ann.Feature == 'exon']
    tran_seq = transcripts.get_transcript_sequence(transcript_id='transcript_id', path=seq)
    with open(outfile, 'w') as fw:
        for index, id, seq in tran_seq.itertuples():
            fw.write(f'>{id}\n{seq}\n')
            
    return outfile

def classify_orf(row):
    
    if row['stop'] < row['tran_start']:
        return "uORF"
    elif row['start'] == row['tran_start'] and row['stop'] == row['tran_stop']:
        return "CDS"
    elif row['start'] > row['tran_stop']:
        return "dORF"
    elif row['start'] < row['tran_start'] and row['stop'] >= row['tran_start']:
        return "uoORF"
    elif row['start'] <= row['tran_stop'] and row['stop'] > row['tran_stop']:
        return "doORF"
    elif row['start'] >= row['tran_start'] and row['stop'] <= row['tran_stop']:
        return "iORF"
    elif row['start'] < row['tran_start'] and row['stop'] > row['tran_stop']:
        return "eoORF"
    elif row['start'] < row['tran_start'] and row['stop'] == row['tran_stop']:
        return "extORF"
    else:
        print("unexpected", row['start'], row['stop'], row['tran_start'], row['tran_stop'])
        return "Unexpected"

def orfrelativeposition(annotation, df, exon_df):
    '''
    Determines the relative position of ORFs to coding sequences (CDS).

    This function takes a genome annotation file and a DataFrame containing ORF coordinates,
    and determines the relative position of each ORF with respect to coding sequences (CDS).
    It classifies each ORF into different categories based on its relationship with CDS.

    Parameters:
        annotation (str): Path to the genome annotation file in BED/GFF/GTF format.
        df (polars.DataFrame): DataFrame containing ORF coordinates. It must have columns
                               'tran_id', 'pos', and 'end' representing transcript ID, start
                               position, and end position of each ORF respectively.

    Returns:
        tuple: A tuple containing two polars DataFrames:
               - The first DataFrame contains ORF coordinates with an additional column
                 'type' indicating the relative position of each ORF to CDS.
               - The second DataFrame contains exon coordinates.

    Example:
        orf_df, exon_coords = orfrelativeposition("annotation.gff", orf_df)
    '''
    orflist=[]
    cds_df, exon_coords = getexons_and_cds(annotation, list(df["tran_id"].unique()))
    print('typing ORFS')

    tranids = list(cds_df['tran_id'].unique())
    #TYPING ORFS
    codingorfs = df.with_columns(shared=pl.col('tran_id').is_in(tranids))\
                .filter(pl.col('shared') == True)\
                .select(pl.all().exclude('shared'))
    codingorfs = codingorfs.join(cds_df, on='tran_id')
    codingorfs = codingorfs.with_columns(
                pl.struct(['start', 'stop', 'tran_start', 'tran_stop'])
                .apply(lambda row: classify_orf(row))
                .alias('type'))\
                .select(pl.all().exclude('tran_start', 'tran_stop'))\
                .to_dict(as_series=False)
    orflist.append(codingorfs)
    #NON CODING ORFS
    noncodingorfs = df.with_columns(shared=pl.col('tran_id').is_in(tranids))\
                .filter(pl.col('shared') == False)\
                .select(pl.all().exclude('shared'))
    noncodingorfs = noncodingorfs.with_columns(type=pl.lit('Non Coding')).to_dict(as_series=False)
    orflist.append(noncodingorfs)
    #MAKE ONE DF
    df = pl.from_dicts(orflist).explode('tran_id', 'start', 'stop', 'length', 'startorf', 'stoporf', 'type')
    print('done typing')
    return df, exon_coords


def create_automaton(codons):
    automaton = ahocorasick.Automaton()

    for idx, key in enumerate(codons):
        automaton.add_word(key, (idx, key))
    automaton.make_automaton()
    return automaton

def preporfs(transcript, starts, stops, minlength, maxlength):
    """
    Extracts potential ORFs (Open Reading Frames) from transcript sequences.

    This function reads transcript sequences from a FASTA file, identifies potential
    ORFs within the sequences based on start and stop codons, and returns a DataFrame
    containing information about the identified ORFs.

    Parameters:
        transcript (str): Path to the FASTA file containing transcript sequences.
        starts (list): List of start codons to search for (e.g., ['ATG']).
        stops (list): List of stop codons to search for (e.g., ['TAA', 'TAG', 'TGA']).
        minlength (int): Minimum length of ORFs to consider.
        maxlength (int): Maximum length of ORFs to consider.

    Returns:
        polars.DataFrame: DataFrame containing information about identified ORFs.
                          It has columns 'tran_id', 'pos', 'end', 'length', 'startorf',
                          and 'stoporf', indicating transcript ID, start position, end
                          position, length, start codon, and stop codon of each ORF respectively.

    Example:
        orf_df = preporfs("transcripts.fasta", ['ATG'], ['TAA', 'TAG', 'TGA'], 50, 500)
    """
    startautomaton = create_automaton(starts)
    stopautomaton = create_automaton(stops)
    dict_list=[]
    # COUNTER!!!!!!!!!!
    counter = 0
    with open(transcript) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if counter % 20000 == 0:
                print(f"read {counter} Transcripts")
            tran_id = str(record.id).split("|")[0]
            append_list = find_orfs(str(record.seq), tran_id, startautomaton, stopautomaton, minlength, maxlength)
            dict_list.extend(append_list)
            counter = counter + 1
        df = pl.from_dicts(dict_list)
        df = df.sort("tran_id")
        print(df.filter(pl.col('tran_id') == 'ENST00000379265.5'))
        return df
