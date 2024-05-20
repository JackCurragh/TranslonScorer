from Bio import SeqIO
import pyranges as pr
import polars as pl
import ahocorasick

from orffinder import find_orfs
from findexonscds import getexons_and_cds


def gettranscripts(seq, annotation, outfile="transcripts.fa"):
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
    transcripts = ann[ann.Feature == 'transcript']
    tran_seq = transcripts.get_transcript_sequence(transcript_id='transcript_id', path=seq)
    with open(outfile, 'w') as fw:
        for index, id, seq in tran_seq.itertuples():
            fw.write(f'>{id}\n{seq}\n')
            
    return outfile

def orfrelativeposition(annotation, df):
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
    cds_df, exon_coords = getexons_and_cds(annotation, list(df["tran_id"].unique()))
    filtered_df = cds_df.select(pl.all().exclude("start", "stop"))
    orftype = []
    for orfid in sorted(df["tran_id"].unique()):
        cdsregion = filtered_df.filter(pl.col("tran_id") == orfid)
        orfregion = df.filter(pl.col("tran_id") == orfid)
        orfpair = list(zip(orfregion["pos"], orfregion["end"]))
        if not cdsregion.is_empty():
            if cdsregion["tran_start"][0] > cdsregion["tran_stop"][0]:
                cdspair = [cdsregion["tran_stop"][0], cdsregion["tran_start"][0]]
            else:
                cdspair = [cdsregion["tran_start"][0], cdsregion["tran_stop"][0]]
            for orf in orfpair:
                if orf[0] < cdspair[0] and orf[1] < cdspair[0]:
                    orftype.append("uORF")
                elif orf[0] == cdspair[0] and orf[1] == cdspair[1]:
                    orftype.append("CDS")
                elif orf[0] > cdspair[1] and orf[1] > cdspair[1]:
                    orftype.append("dORF")
                elif (
                    orf[0] < cdspair[0] and orf[1] < cdspair[1] and orf[1] > cdspair[0]
                ):
                    orftype.append("uoORF")
                elif (
                    orf[0] <= cdspair[1]
                    and orf[0] >= cdspair[0]
                    and orf[1] > cdspair[1]
                ):
                    orftype.append("doORF")
                elif (
                    orf[0] < cdspair[0]
                    and orf[1] <= cdspair[0]
                ):
                    orftype.append("uoORF")
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
        return df
