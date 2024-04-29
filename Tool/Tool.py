#/usr/bin/python3
"""Script to specify arguments on CLI"""
import os
import click
import polars as pl

from readfiles import readbam
from fileprocessor import dftobed, bedtobigwig


@click.command()
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
        #foldername = os.path.
        for file in list_dir:
            bamcheck = file.split(".", 1)
            if bamcheck[1] == "bam":
                location = location + '/' + file
                df = readbam(location)
                #calculate asite and converting dataframe to BEDGRAPH format
                asitedf = dftobed(df)
                y.append(asitedf)
            else: click.echo("No .bam files found in folder")
        df_merged = pl.concat(y)
        df_merged.write_csv(f"{filename}.bedGraph", separator="\t")

               
               

        
    
    click.echo()
if __name__ == "__main__":
    cli()
