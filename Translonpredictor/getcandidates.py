from Bio import SeqIO
import pyranges as pr
import polars as pl
import ahocorasick

from orffinder import find_orfs
from findexonscds import getexons_and_cds


def gettranscripts(seq, annotation, outfilename):
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
    ann = pr.read_gtf(annotation, ignore_bad=True)
    transcripts = ann[ann.Feature == "exon"]
    tran_seq = transcripts.get_transcript_sequence(
        transcript_id="transcript_id", path=seq
    )
    with open(f"{outfilename}_transcripts.fa", "w") as fw:
        for index, id, seq in tran_seq.itertuples():
            fw.write(f">{id}\n{seq}\n")

    return f"{outfilename}_transcripts.fa"


def classify_orf(row):
    """
    Classify an ORF (Open Reading Frame) based on its relative position to transcript start and stop sites.

    Parameters:
    - row (Series or dict-like): A Pandas Series or dictionary-like object containing ORF information,
                                including 'start', 'stop', 'tran_start', and 'tran_stop' values.

    Returns:
    - str: A string indicating the classification of the ORF based on its relative position:
           - "uORF": Upstream ORF (stop < tran_start)
           - "CDS": Coding Sequence (start == tran_start and stop == tran_stop)
           - "dORF": Downstream ORF (start > tran_stop)
           - "uoORF": Upstream Overlapping ORF (start < tran_start and stop >= tran_start)
           - "doORF": Downstream Overlapping ORF (start <= tran_stop and stop > tran_stop)
           - "iORF": Internal ORF (start >= tran_start and stop <= tran_stop)
           - "eoORF": Encapsulated Overlapping ORF (start < tran_start and stop > tran_stop)
           - "extORF": Extended ORF (start < tran_start and stop == tran_stop)
           - "Unexpected": Indicates unexpected conditions where none of the above criteria are met.

    This function categorizes an ORF based on its positional relationship with the transcript start (`tran_start`)
    and stop (`tran_stop`) sites. It evaluates the relative positions of 'start' and 'stop' compared to
    'tran_start' and 'tran_stop' to determine the appropriate classification.

    If the relative position does not match any expected categories, an "Unexpected" classification is printed
    with the values of 'start', 'stop', 'tran_start', and 'tran_stop', and "Unexpected" is returned.

    Note: This function assumes the input `row` contains numerical values for 'start', 'stop', 'tran_start',
    and 'tran_stop', typically retrieved from a Pandas DataFrame or similar data structure.
    """
    if row["stop"] < row["tran_start"]:
        return "uORF"
    elif row["start"] == row["tran_start"] and row["stop"] == row["tran_stop"]:
        return "CDS"
    elif row["start"] > row["tran_stop"]:
        return "dORF"
    elif row["start"] < row["tran_start"] and row["stop"] >= row["tran_start"]:
        return "uoORF"
    elif row["start"] <= row["tran_stop"] and row["stop"] > row["tran_stop"]:
        return "doORF"
    elif row["start"] >= row["tran_start"] and row["stop"] <= row["tran_stop"]:
        return "iORF"
    elif row["start"] < row["tran_start"] and row["stop"] > row["tran_stop"]:
        return "eoORF"
    elif row["start"] < row["tran_start"] and row["stop"] == row["tran_stop"]:
        return "extORF"
    else:
        print(
            "unexpected", row["start"], row["stop"], row["tran_start"], row["tran_stop"]
        )
        return "Unexpected"


def orfrelativeposition(annotation, df, cds_df):
    """
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
    """
    orflist = []
    if not "cdsdf" in globals():
        cds_df, exon_coords = getexons_and_cds(annotation, list(df["tran_id"].unique()))

    print("Typing ORFS")
    tranids = list(cds_df["tran_id"].unique())
    # TYPING ORFS
    codingorfs = (
        df.with_columns(shared=pl.col("tran_id").is_in(tranids))
        .filter(pl.col("shared") == True)
        .select(pl.all().exclude("shared"))
    )

    codingorfs = codingorfs.join(cds_df, on="tran_id")
    codingorfs = (
        codingorfs.with_columns(
            pl.struct(["start", "stop", "tran_start", "tran_stop"])
            .apply(lambda row: classify_orf(row))
            .alias("type")
        )
        .select(pl.all().exclude("tran_start", "tran_stop"))
        .to_dict(as_series=False)
    )
    orflist.append(codingorfs)

    # NON CODING ORFS
    noncodingorfs = (
        df.with_columns(shared=pl.col("tran_id").is_in(tranids))
        .filter(pl.col("shared") == False)
        .select(pl.all().exclude("shared"))
    )

    noncodingorfs = noncodingorfs.with_columns(type=pl.lit("Non Coding")).to_dict(
        as_series=False
    )
    orflist.append(noncodingorfs)
    # MAKE ONE DF
    df = pl.from_dicts(orflist).explode(
        "tran_id", "start", "stop", "length", "startorf", "stoporf", "type"
    )

    cdslist = df.filter(pl.col("type") == "CDS")
    cdslist = list(cdslist["tran_id"].unique())

    return df, exon_coords


def create_automaton(codons):
    """
    Create an Aho-Corasick automaton for efficient substring search.

    Parameters:
    - codons (list): A list of strings representing codons to be added to the automaton.

    Returns:
    - ahocorasick.Automaton: An Aho-Corasick automaton initialized with the provided codons.

    This function initializes an Aho-Corasick automaton, a data structure used for efficient
    substring matching. It adds each codon from the `codons` list to the automaton and assigns
    each codon a unique index and key. After adding all codons, the automaton is finalized with
    `make_automaton()`.

    The resulting `automaton` object is returned, which can be used to quickly search for any
    of the added codons within a given text or sequence.

    Note: This function assumes the use of the `ahocorasick` library, which provides an efficient
    implementation of the Aho-Corasick algorithm for pattern matching.
    """
    automaton = ahocorasick.Automaton()

    for idx, key in enumerate(codons):
        automaton.add_word(key, (idx, key))
    automaton.make_automaton()
    return automaton


def preporfs(transcript, starts, stops, minlength, maxlength):
    """
    Predict ORFs (Open Reading Frames) from transcript sequences using start and stop codon patterns.

    Parameters:
    - transcript (str): Path to a FASTA file containing transcript sequences.
    - starts (list): List of strings representing start codon patterns (e.g., ['ATG', 'GTG', 'TTG']).
    - stops (list): List of strings representing stop codon patterns (e.g., ['TAA', 'TAG', 'TGA']).
    - minlength (int): Minimum length threshold for predicted ORFs.
    - maxlength (int): Maximum length threshold for predicted ORFs.

    Returns:
    - DataFrame: A Pandas DataFrame containing predicted ORFs for each transcript sequence.

    This function reads transcript sequences from the provided FASTA file (`transcript`). For each
    transcript sequence, it creates Aho-Corasick automata (`startautomaton` and `stopautomaton`)
    using the provided lists of start and stop codon patterns (`starts` and `stops`).

    It iterates through each transcript sequence, identifies ORFs using the `find_orfs` function,
    and appends the results to `dict_list`. After processing all transcripts, `dict_list` is converted
    to a Pandas DataFrame (`df`) where each dictionary represents a row of ORF predictions.

    The DataFrame `df` is sorted by `tran_id` and returned as the final output.

    Note: This function assumes the use of the Biopython library (`SeqIO` for parsing FASTA files)
    and a custom function `find_orfs` for ORF prediction based on Aho-Corasick automata.

    Example usage:
    ```python
    transcript = 'transcript.fasta'
    starts = ['ATG', 'GTG', 'TTG']
    stops = ['TAA', 'TAG', 'TGA']
    minlength = 50
    maxlength = 5000
    predicted_orfs = preporfs(transcript, starts, stops, minlength, maxlength)
    ```
    """
    dict_list = []
    # COUNTER!!!!!!!!!!
    counter = 0
    with open(transcript) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            startautomaton = create_automaton(starts)
            stopautomaton = create_automaton(stops)
            if counter % 20000 == 0:
                print("\r" + f"Read {counter} transcripts", end="")
            tran_id = str(record.id).split("|")[0]
            append_list = find_orfs(
                str(record.seq),
                tran_id,
                startautomaton,
                stopautomaton,
                minlength,
                maxlength,
            )
            dict_list.extend(append_list)
            counter = counter + 1
        df = pl.from_dicts(dict_list)
        df = df.sort("tran_id")
        print("\n")
        return df
