



def process_bam_workflow(bam_path, annotation_file):
    """Process BAM file and get exons/CDS."""
    log_info('Processing BAM file...')
    # bam_df = bam.readbam(bam_path)
    
    # # Get exons and CDS
    # cds_df, exon_df = bam.getexons_and_cds(annotation_file)
    
    # # Process BAM data
    # bam_type, _ = bam.detect_bam_type(bam_df, exon_df)
    # if bam_type == 'genomic':
    #     # First convert genomic coordinates to transcript coordinates
    #     bam_df = bam.bamtranscript(bam_df, exon_df)
    #     # Then calculate positions relative to CDS
    #     bam_df = bam.process_transcriptomic_bam(bam_df, cds_df)
    # else:
    #     # For transcriptomic BAM, just calculate CDS positions
    #     bam_df = bam.process_transcriptomic_bam(bam_df, cds_df)
    
    # return bam_df, exon_df


def find_orfs_workflow(sequence_file, annotation_file, start_codons, stop_codons, min_len, max_len):
    """Find ORFs in transcript sequences."""
    log_info('Finding ORFs...')
    # orf_df = orffinder.preporfs(sequence_file, start_codons.split(","), stop_codons.split(","), min_len, max_len)
    
    # log_info('Determining relative position of ORFs to CDS...')
    # # Determine the relative position of ORFs to CDS
    # orf_df, exon_coords = orfrelativeposition(annotation_file, orf_df)
    
    # return orf_df, exon_coords


def score_orfs_workflow(bigwig_paths, exon_df, orf_df, scoring_method, sru_range, stranded):
    """Score ORFs using BigWig data."""
    log_info('Scoring ORFs...')
    # scored_orfs = bigwig.scoring(bigwig_paths, exon_df, orf_df, scoring_method == 'classic', sru_range, stranded)
    
    # return scored_orfs


def plot_workflow(scored_orfs, bigwig_path, exon_df, plot_range, outfile):
    """Generate plots for top 10 ORFs."""
    log_info('Generating plots...')
    # plots.plottop10(scored_orfs, bigwig_path, exon_df, plot_range, outfile)






def all_workflow(bam_path, bigwig_path, forward_bigwig, reverse_bigwig, chromsizes, annotation, sequence, start_codons, stop_codons, min_len, max_len, scoring_method, sru_range, plot_range, outfile):

    log_info("Starting pipeline...")


    # # Validate inputs
    # if not bam_path and not bigwig_path and not (forward_bigwig and reverse_bigwig):
    #     raise click.BadParameter(
    #         "Either BAM file (-b) or BigWig file (-bw) or forward/reverse BigWig files must be provided"
    #     )
    
    # if bam_path and not chromsizes:
    #     raise click.BadParameter("Chromosome sizes file (-c) is required when processing BAM files")
    
    # # Handle strand-specific bigwig inputs
    # if forward_bigwig and reverse_bigwig:
    #     bigwig_paths = {'forward': forward_bigwig, 'reverse': reverse_bigwig}
    #     stranded = True
    # elif bigwig_path:
    #     bigwig_paths = bigwig_path
    # else:
    #     bigwig_paths = None
    
    # # Initialize exon_df
    # exon_df = None

    # # Determine pipeline stages
    # if bam_path:
    #     if not bigwig_paths:
    #         log_info("BAM file provided without BigWig. Will process BAM to generate coverage.")
    #         location = os.path.abspath(bam_path)
    #         if not os.path.isfile(location):
    #             raise click.BadParameter(f"BAM file not found: {bam_path}")
            
    #         # Process BAM and get annotations
    #         bam_df, exon_df = process_bam_file(location, annotation)
            
    #         # Calculate A-site positions
    #         offsets = coordinates.change_point_analysis(bam_df)
            
    #         # Process by strand if requested
    #         if stranded:
    #             # Split BAM by strand
    #             fwd_bam = bam_df.filter(pl.col("strand") == "+")
    #             rev_bam = bam_df.filter(pl.col("strand") == "-")
                
    #             # Process forward strand
    #             fwd_bed = bed.asitecalc(fwd_bam, offsets)
    #             fwd_bedgraph = f"{outfile}_fwd.bedGraph"
    #             fwd_bed.write_csv(fwd_bedgraph, separator="\t", include_header=False)
    #             fwd_bigwig = f"{outfile}_fwd.bw"
    #             bed.bedtobigwig(fwd_bedgraph, chromsizes, f"{outfile}_fwd")
                
    #             # Process reverse strand
    #             rev_bed = bed.asitecalc(rev_bam, offsets)
    #             rev_bedgraph = f"{outfile}_rev.bedGraph"
    #             rev_bed.write_csv(rev_bedgraph, separator="\t", include_header=False)
    #             rev_bigwig = f"{outfile}_rev.bw"
    #             bed.bedtobigwig(rev_bedgraph, chromsizes, f"{outfile}_rev")
                
    #             # Set up paths for scoring
    #             bigwig_paths = {'forward': f"{outfile}_fwd.bw", 'reverse': f"{outfile}_rev.bw"}
    #         else:
    #             # Process combined strands (original behavior)
    #             bed_df = bed.asitecalc(bam_df, offsets)
    #             bedgraph_path = f"{outfile}.bedGraph"
    #             bed_df.write_csv(bedgraph_path, separator="\t", include_header=False)
    #             bigwig_path = f"{outfile}.bw"
    #             bed.bedtobigwig(bedgraph_path, chromsizes, outfile)
    #             bigwig_paths = bigwig_path
    #     else:
    #         log_info("Both BAM and BigWig provided. Using BigWig directly.")
    
    # # Ensure exon_df is set when using BigWig directly
    # if exon_df is None:
    #     log_info("Loading exon data from annotation file...")
    #     cds_df, exon_df = bam.getexons_and_cds(annotation)
    
    # # Ensure transcriptomic input for ORF finding
    # log_info("Ensuring transcriptomic input for ORF finding...")
    # transcript_fasta = sequence
    # if not sequence.endswith('_transcripts.fa'):
    #     log_info("Generating transcript sequences from genomic FASTA and GTF annotation...")
    #     transcript_fasta = coordinates.gettranscripts(sequence, annotation, outfile)
    
    # # Find ORFs
    # log_info("Finding ORFs...")
    # orf_df = orffinder.preporfs(transcript_fasta, start_codons.split(","), stop_codons.split(","), min_len, max_len)
    
    # log_info("Determining relative position of ORFs to CDS...")
    # # Determine the relative position of ORFs to CDS
    # orf_df, exon_coords = orfrelativeposition(annotation, orf_df, cds_df)

    # # Score ORFs using bigWig data
    # log_info("Scoring ORFs...")
    # scored_orfs = bigwig.scoring(bigwig_paths, exon_df, orf_df, scoring_method == 'classic', sru_range, stranded)
    
    # scored_orfs.write_csv(f"{outfile}_orfs_scored.csv")
    
    # # Generate plots
    # log_info("Generating plots...")
    # if stranded and isinstance(bigwig_paths, dict):
    #     # Use forward strand for plotting if we have strand-specific data
    #     plot_bigwig = bigwig_paths['forward']
    # else:
    #     plot_bigwig = bigwig_paths
        
    # plots.plottop10(scored_orfs, plot_bigwig, exon_df, plot_range, outfile)
    
    # log_info("Pipeline completed successfully!")