import orfipy_core as orf
import re
import time
import ahocorasick
from Bio import SeqIO
import polars as pl


def create_automaton(codons):
    automaton = ahocorasick.Automaton()

    for idx, key in enumerate(codons):
        automaton.add_word(key, (idx, key))
    automaton.make_automaton()
    return automaton

def find_all_positions(sequence, automaton):
    frames = {0:[], 1:[], 2:[]}
    codons = {}
    for i, (order, codon) in automaton.iter(sequence):
        frames[(i-2)%3 ].append(i)
        codons[i] = codon
    return frames, codons




def find_orfs(
    sequence,
    tran_id,
    startautomaton,
    stopautomaton,
    minlength=0,
    maxlength=1000000
):

    orf_list=[]
    startpositions, start_codons = find_all_positions(sequence, startautomaton)
    stoppositions, stop_codons = find_all_positions(sequence, stopautomaton)

   
    for frame, startpositions in startpositions.items():
        for position in startpositions:
            valid_stops = [i for i in stoppositions[frame] if i > position]
            if valid_stops:
                stopposition = min(valid_stops)
                stopcodon = stop_codons[stopposition]
            else:
                stopposition = len(sequence)
                stopcodon = sequence[-3:]
            orf_data= {
                "tran_id":tran_id,
                "pos":position-1,
                "end": stopposition+1,
                "length":stopposition+1 - position-1,
                "startorf":start_codons[position],
                "stoporf": stopcodon,
                }
            if orf_data['length'] < maxlength and orf_data['length'] > minlength:
                orf_list.append(orf_data)
            print(orf_list)
    return orf_list

#test data
sequence = "ATGCATGACTAGCATCAGCATCAGCATATGGACGAATTAAGTTAA"
tran_id = 'ENST1'
start_codons=["ATG"]
stop_codons=["TAA", "TAG", "TGA"]
startautomaton = create_automaton(start_codons)
stopautomaton = create_automaton(stop_codons)

find_orfs(sequence, tran_id, startautomaton, stopautomaton)

def preporfs(transcript, starts, stops, minlength, maxlength):
    startautomaton = create_automaton(starts)
    stopautomaton = create_automaton(stops)
    dict_list=[]
    # COUNTER!!!!!!!!!!
    counter = 0
    with open(transcript) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if counter % 100000 == 0:
                print(f"read {counter} Transcripts")
            tran_id = str(record.id).split("|")[0]
            append_list = find_orfs(str(record.seq), tran_id, startautomaton, stopautomaton)
            dict_list.extend(append_list)
            counter = counter + 1
        df = pl.from_dicts(dict_list)
        df = df.sort("tran_id")
        print(df)
        return df


#preporfs('data/gencode.v45.transcripts.fa', start_codons,stop_codons, minlength=0, maxlength=100000)

data = [(10,27,18,'ATGAAATGGACGAATTAA'), (15,32,18,'ATGGACGAATTAAGTTAA')]
#test orffinder
def test_orffinder():
    assert find_orfs(sequence) == data