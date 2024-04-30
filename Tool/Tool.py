#/usr/bin/python3
"""Script to specify arguments on CLI"""
import os
import click
import polars as pl

from readfiles import readbam
from fileprocessor import dftobed, bedtobigwig
from getcandidates import gettranscripts, preporfs

@click.group()
def function():
    pass

@function.command()
@click.argument("input")
@click.option("--chromsize", "-c", help="Provide a file containing the chromosome sizes")
def cli(input, chromsize):
    location = os.getcwd() + '/' + input
    #if file is provided
    if  os.path.isfile(location):
        #read in bam file
        df = readbam(location)
        #calculate asite + converting to BedGraph
        beddf = dftobed(df)
        
        bedfile = beddf.write_csv("file.bedGraph", separator="\t")

        #Converting Bedgrapgh to Bigwig format
        bedtobigwig(bedfile, chromsize)
    
    #if folder with BAM files is provided
    elif os.path.isdir(location):
        list_dir = os.listdir(location)
        y = []
        filename = os.path.basename(location)
        for file in list_dir:
            bamcheck = file.split(".", 1)
            if bamcheck[1] == "bam":
                location = location + '/' + file
                #read in bam file
                df = readbam(location)
                #calculate asite and converting dataframe to BEDGRAPH format
                asitedf = dftobed(df)
                y.append(asitedf)
            else: click.echo("No .bam files found in folder")
        df_merged = pl.concat(y)
        df_merged.write_csv(f"{filename}.bedGraph", separator="\t")



@function.command()
@click.option("--seq", "-s", help="Provide a file containing the genomic sequence (.fa)")
@click.option("--tran", "-t", help="Provide a file containing the transcript sequences (.fa)")
@click.option("--ann", "-a", help="Provide a file containing the annotation (.gtf)")
@click.option("--starts", "-sta", default="ATG", help="Provide a list of start codons")
@click.option("--stops", "-stp", default="TAA, TAG, TGA", help="Provide a list of stop codons")
@click.option("--minlen", "-min", default=0, help="Provide the minimum length")
@click.option("--maxlen", "-max", default=1000000, help="Provide the maximum length")

def orfprep(seq, ann, tran, starts, stops, minlen, maxlen):
    '''
    docstring
    '''
    if seq and ann:
        transcript = gettranscripts(seq, ann)
    elif tran:
        transcript = tran
    else: raise Exception("Must provide genomic sequence and annotation or transcript sequences")

    preporfs(transcript, starts.split(","), stops.split(","), minlen, maxlen)

    
if __name__ == "__main__":
    function()
