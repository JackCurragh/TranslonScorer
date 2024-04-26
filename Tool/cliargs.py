#/usr/bin/python3
"""Script to specify arguments on CLI"""
import os
import click
import pysam
import oxbow as ox
import polars as pl

@click.command()
@click.argument("input")
def cli(input):
    location = os.getcwd() + '/' + input
    #if BAM file is provided
    if  os.path.isfile(location):
        #index file before reading
        pysam.index(location)     
        bamfile = ox.read_bam(location)
        df = pl.read_ipc(bamfile)

####################################################################################################################################
        #Code to implement in new script -> transforming columns

        df_select = df.select(pl.col("qname","rname","pos","end"))
        df_filtered = df_select.with_columns((pl.col("end") - pl.col("pos")).alias("length")).select(pl.all().exclude("end"))
        # Define a function to split the string and return the parts separately
        split_func = lambda s: int(s.split("_x")[1])
        # Apply the split function to the "qname" column and create two new columns
        df_namesplit = df_filtered.with_columns(
            pl.col("qname").apply(split_func)).rename({"qname":"count"})

        #OFFSETS
        length_values = df_namesplit.get_columns()[3]
        offsets = {
            i: 15 for i in length_values
        }
        #GROUP TO CALC A-SITE
        y = []
        for value, data in df_namesplit.group_by("length"):
            x = df_namesplit.filter(pl.col("length") == value).with_columns((pl.col("length") + offsets[value]).alias("A-site"))
            y.append(x)
    
        df_asite = pl.concat(y)

        #GROUP ON A-SITE!
        df_final = df_asite.group_by("rname", "A-site").agg(pl.col("count").sum())
        
        df_final.write_csv("out.csv")
        print(df_final)
####################################################################################################################################


    #if folder with BAM files is provided
    elif os.path.isdir(location):
        list_dir = os.listdir(location)
        for file in list_dir:
            bamcheck = file.split(".", 1) 
            if bamcheck[1] == "bam":
                location = location + '/' + file
                pysam.index(location)     
                bamfile = ox.read_bam(location)
                df = pl.read_ipc(bamfile)




            else: click.echo("No .bam files found in folder")
    
    click.echo()
if __name__ == "__main__":
    cli()
