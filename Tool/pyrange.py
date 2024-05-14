'''This is a temporary script to try and implement the pyrange package into the script'''

import pyranges as pr
import polars as pl

def pyrange(annotation, tran=[]):
    
    ann = pr.read_gtf(annotation)
    ann = ann.drop(['Source', 'Score', 'havana_transcript', 'exon_number', 'hgnc_id', 'havana_gene',
                    'ont', 'protein_id', 'ccdsid', 'artif_dupl', 'Frame', 'tag', 'level',
                    'gene_id', 'gene_name','gene_type', 'transcript_type', 'transcript_name',
                    'transcript_support_level', 'exon_id'], axis=1)

    if tran:
        ann = ann[ann['transcript_id'].isin(tran)]

    # Getting CDS
    coding_regions = ann[ann.Feature == 'CDS']

    groupedcds = coding_regions.groupby('transcript_id').agg(lambda x: list(x))
    print(groupedcds)
    
    #groupedcds = (
    #    coding_regions.group_by("tran_id")
    #    .agg(
    #        pl.col("start").alias("start"), pl.col("stop").alias("stop"), pl.col("chr")
    #    )
    #    .select(["chr", "tran_id", "start", "stop"])
    #)

    # Getting exons
    exon_regions = ann[ann.Feature == 'exon']

    pos_exons = exon_regions[exon_regions.Strand == '+']
    neg_exons = exon_regions[exon_regions.Strand == '-']
    


    pos_exons, neg_exons = procesexons(exon_regions)

    exon_coords_plus = exontranscriptcoords(pos_exons, posstrand=True)
    # column names switched to calculate inverse of positions for negative strands
    exon_coords_neg = exontranscriptcoords(neg_exons, posstrand=False)

    exon_coords = pl.concat([exon_coords_plus, exon_coords_neg]).select(
        pl.all().exclude("strand")
    )
    cds_coords = gettranscriptcoords(groupedcds, exon_coords)

    return cds_coords, exon_coords

pyrange("data/gencode.v45.annotation.gtf")