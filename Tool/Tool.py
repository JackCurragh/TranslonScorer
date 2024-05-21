# /usr/bin/python3
"""Script to specify arguments on CLI"""
import os
import click
import polars as pl

from readfiles import readbam
from fileprocessor import dftobed, bedtobigwig
from getcandidates import gettranscripts, preporfs, orfrelativeposition
from filewriter import saveorfsandexons
from bigwigtodf import bigwigtodf
from plotting import scoreandplot

@click.group()
def function():
    pass

@function.command()
@click.option("--bam", "-b", help="Provide a bam file")
@click.option("--chromsize", "-c", help="Provide a file containing the chromosome sizes")
@click.option("--seq", "-s", help="Provide a file containing the genomic sequence (.fa)")
@click.option("--tran", "-t", help="Provide a file containing the transcript sequences (.fa)")
@click.option("--ann", "-a", help="Provide a file containing the annotation (.gtf)")
@click.option("--starts", "-sta", default="ATG", help="Provide a list of start codons")
@click.option("--stops", "-stp", default="TAA, TAG, TGA", help="Provide a list of stop codons")
@click.option("--minlen", "-min", default=0, help="Provide the minimum length")
@click.option("--maxlen", "-max", default=1000000, help="Provide the maximum length")
@click.option("--bigwig", "-bw", help="Provide a Bigwig file to convert")
@click.option("--exon", "-ex", help="Provide a file containing exon positions")
@click.option("--bedfile", "-bw", help="Provide a Bigwig file to convert")
@click.option("--orfs", "-of", help="Provide a file containing annotated orfs")
@click.option("--range_param", "-rp",
             help="Provide an integer that indicates the range in which a plot will be constructed \
                 around the relative start position")
@click.option("--sru_range", "-sru", help="Provide an integer that indicates the range for \
             the Start Rise Up score. This sets the amount of nucleotides before and after \
             the stop codon will regarded when calculating")


def tool(bam, bedfile, chromsize, bigwig, seq, tran, ann, starts, stops, minlen, maxlen,
         exon, orfs, range_param, sru_range):
    """
    docstring
    """
    if bam or chromsize or bedfile:
        if bam and chromsize and ann:
            print('bam+chromsize')
            #READ IN BAM START
            location = os.getcwd() + "/" + bam
            # if file is provided
            if os.path.isfile(location):
                # read in bam file
                df = readbam(location)
                # calculate asite + converting to BedGraph
                beddf, exondf = dftobed(df, ann)

                if not os.path.exists("data/file.bedGraph"):
                    bedfile = beddf.write_csv("file.bedGraph", separator="\t")

                    # Converting Bedgrapgh to Bigwig format
                    bigwig = bedtobigwig(bedfile, chromsize)

        elif bedfile and chromsize:
            print('bedfile, chromsize')
            bigwig = bedtobigwig(bedfile, chromsize)

        elif bigwig:
            print('bigwig')
        else: raise Exception(
            "Must provide valuable input. The options are the following:\n \
                1. Bam file (.bam) + file containing chromosome information\n \
                2. BedGraph (.bedGraph) file + file containing chromosome information\n \
                3. Bigwig (.bw) file"
            )
        
    #ORFPREP START
    if ann:
        if seq and ann:
            print('getting transcript')
            transcript = gettranscripts(seq, ann)
        elif tran and ann:
            print('transcript set')
            transcript = tran
        print('getting orfs')
        orfdf = preporfs(transcript, starts.split(", "), stops.split(", "), minlen, maxlen)
        exondf = 0
        if exondf == 0:
            exondf = pl.DataFrame()
        orf_ann_df, exon_df = orfrelativeposition(ann, orfdf, exondf)   
        orfs, exon = saveorfsandexons(orf_ann_df, exon_df)
        
        print('bigwigtodf')
        bwtrancoords = bigwigtodf(bigwig, exon)
        bwtrancoords.write_csv("data/files/bw_tran.csv")


    elif orfs and exon:
        print('orfs and exon')
        print('bigwigtodf')
        bwtrancoords = bigwigtodf(bigwig, exon)
        bwtrancoords.write_csv("data/files/bw_tran.csv")
        
        print('scoring and plotting')
        scoreandplot(orfs, "data/files/bw_tran.csv" ,range_param, sru_range)
    
    else:
        raise Exception(
            "Must provide valuable for ORFS. The options are the following:\n \
                1. A file containing FASTA sequences (.fa) + annotation file (.gtf)\n \
                2. A file containing transcript sequences (.fa) + annotation file (.gtf)\n \
                3. A file containing annotated ORFS (.csv) + file containg exon information (.csv)"
        )



############################################################################################################
    # if folder with BAM files is provided
    '''elif os.path.isdir(location):
        list_dir = os.listdir(location)
        y = []
        filename = os.path.basename(location)
        for file in list_dir:
            bamcheck = file.split(".", 1)
            if bamcheck[1] == "bam":
                location = location + "/" + file
                # read in bam file
                df = readbam(location)
                # calculate asite and converting dataframe to BEDGRAPH format
                asitedf = dftobed(df)
                y.append(asitedf)
            else:
                click.echo("No .bam files found in folder")
        df_merged = pl.concat(y)
        df_merged.write_csv(f"{filename}.bedGraph", separator="\t")
        bedtobigwig(bedfile, chromsize)'''

if __name__ == "__main__":
    function()
