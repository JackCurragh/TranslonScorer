"""Script to specify arguments on CLI"""
import os
import click
import pysam
import oxbow as ox
import polars as pl
import time 
@click.command()
@click.argument("input")
def cli(input):
    location = os.getcwd() + '/' + input
    #if BAM file is provided
    if  os.path.isfile(location) == True:
        #index file before reading
        pysam.index(location)     
        bamfile = ox.read_bam(location)
        df = pl.read_ipc(bamfile)
        df.head(5)

    #if folder with BAM files is provided
    elif os.path.isdir(location) == True:
        list_dir = os.listdir(location)
        for file in list_dir:
            bamcheck = file.split(".", 1) 
            if bamcheck[1] == "bam":
                location = location + '/' + file
                pysam.index(location)     
                bamfile = pysam.AlignmentFile(location, "rb")
                click.echo("Reading for file: {}".format(os.path.basename(location)))
                for read in bamfile.fetch():
                    click.echo(read)
                bamfile.close()
            else: click.echo("No .bam files found in folder")
    
    click.echo()
if __name__ == "__main__":
    time.sleep(100)
    cli()
