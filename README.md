# TranslonScorer

## Overview

`TranslonScorer` is a command-line tool for scoring potential translational events (translons) from Ribo-seq data. The tool processes Ribo-seq BAM files, identifies potential ORFs from transcript sequences, and scores them based on ribosome coverage patterns. It provides detailed scoring metrics that allow users to apply their own thresholds for translon classification based on their required stringency. The tool supports multiple input file formats and generates output files in several formats including `.bedGraph`, `.bw`, `.html`, and `.csv`.

## Installation

To use this tool, you need to have Python and the necessary dependencies installed. Install the required packages using `pip`:

```sh
pip install TranslonScorer@git+https://github.com/JackCurragh/TranslonScorer#egg=TRANSLONSCORER
```

## Usage

The tool provides several workflows for different analysis needs:

### 1. Complete Pipeline (Recommended)

Run the entire pipeline end-to-end with a single command:

```sh
translonscorer all \
    -b ribo.bam \
    -c chrom.sizes \
    -s genome.fa \
    -a anno.gtf \
    -o output
```

This will:
1. Process the Ribo-seq BAM file
2. Extract transcripts
3. Find and score potential ORFs
4. Generate visualization reports

### 2. Individual Steps

For more control, you can run each step separately:

#### Process BAM Files
```sh
translonscorer process-bam \
    -b ribo.bam \
    -c chrom.sizes \
    -a anno.gtf \
    -o output
```

#### Find and Score ORFs
```sh
translonscorer find-orfs \
    -s genome.fa \
    -a anno.gtf \
    -bw coverage.bw \
    -o output
```

#### Score Existing ORFs
```sh
translonscorer score-orfs \
    -f orfs.csv \
    -bw coverage.bw \
    -e exons.csv \
    -o output
```

#### Generate Visualization Report
```sh
translonscorer plot \
    -s scored_orfs.csv \
    -bw coverage.bw \
    -e exons.csv \
    -o report
```

## Command Options

### Common Options
- `-o, --outfile`: Base name for output files (required)

### Process BAM Command
- `-b, --bam`: Input BAM file from Ribo-seq data (required)
- `-c, --chromsizes`: Chromosome sizes file (required)
- `-a, --annotation`: GTF annotation file (required)
- `-off, --offsets`: File containing read length-specific offsets for A-site calculation

### Find ORFs Command
- `-s, --sequence`: Input FASTA file (genomic or transcriptomic) (required)
- `-a, --annotation`: GTF annotation file (required)
- `-bw, --bigwig`: BigWig file containing Ribo-seq coverage (required)
- `--start-codons`: Comma-separated list of start codons (default: ATG)
- `--stop-codons`: Comma-separated list of stop codons (default: TAA,TAG,TGA)
- `--min-len`: Minimum ORF length in nucleotides (default: 0)
- `--max-len`: Maximum ORF length in nucleotides (default: 1000000)
- `--sru-range`: Nucleotide range for Start Rise Up score calculation (default: 15)
- `--scoring-method`: Scoring algorithm to use (classic/modern, default: modern)

### Score ORFs Command
- `-f, --orfs`: CSV file containing pre-annotated ORFs (required)
- `-bw, --bigwig`: BigWig file containing Ribo-seq coverage (required)
- `-e, --exons`: CSV file containing exon positions (required)
- `--scoring-method`: Scoring algorithm to use (classic/modern, default: modern)
- `--sru-range`: Nucleotide range for Start Rise Up score calculation (default: 15)

### Plot Command
- `-s, --scored-orfs`: CSV file containing scored ORFs (required)
- `-bw, --bigwig`: BigWig file containing Ribo-seq coverage (required)
- `-e, --exons`: CSV file containing exon positions (required)
- `--plot-range`: Plot range around start position (default: 30)

## Output Files

The tool generates several output files depending on the command used:

### Process BAM
- `{outfile}.bedGraph`: Coverage in bedGraph format
- `{outfile}.bw`: Coverage in bigWig format

### Find ORFs / Score ORFs
- `{outfile}_orfs_scored.csv`: Scored ORFs with multiple metrics for classification
- `{outfile}_report.html`: Visualization report showing score distributions and features

## Interpreting Results

The tool provides several scoring metrics for each potential ORF:
- Start Rise Up (SRU) score: Measures ribosome accumulation at start codons
- High Read Frame (HRF) score: Quantifies reading frame preference
- Average coverage
- Non-Zero Codon ratio

Users should determine appropriate score thresholds based on their specific requirements and experimental context. The visualization report includes score distributions to help inform threshold selection.

## Error Handling

The tool requires specific combinations of input files to function correctly. If the necessary files are not provided, it will raise an exception with guidance on the required files.

Please ensure that:
1. All required input files exist and are readable
2. File formats match the expected types (BAM, FASTA, GTF, etc.)
3. Chromosome notation is consistent across annotation and coverage files

## Contributing

Contributions are welcome. Please fork the repository and submit a pull request.

## License

This project is licensed under the MIT License.
```