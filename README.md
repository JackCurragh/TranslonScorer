# Translon Scorer

## Overview

`TranslonScorer` is a command-line tool for translon calling. The process consists of processing Ribo-seq BAM files, extracting and scoring ORFs from transcript sequences based on the annotation and codons provided. It supports multiple input file formats and generates output files in several formats including `.bedGraph`, `.bw`, `.html`, and `.csv`.

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
translonpredictor all \
    -b ribo.bam \
    -c chrom.sizes \
    -s genome.fa \
    -a anno.gtf \
    -o output
```

This will:
1. Process the Ribo-seq BAM file
2. Extract transcripts
3. Find and score ORFs
4. Generate visualization reports

### 2. Individual Steps

For more control, you can run each step separately:

#### Process BAM Files
```sh
translonpredictor process-bam \
    -b ribo.bam \
    -c chrom.sizes \
    -a anno.gtf \
    -o output
```

#### Find and Score ORFs
```sh
translonpredictor find-orfs \
    -s genome.fa \
    -a anno.gtf \
    -bw coverage.bw \
    -o output
```

#### Score Existing ORFs
```sh
translonpredictor score-orfs \
    -f orfs.csv \
    -bw coverage.bw \
    -e exons.csv \
    -o output
```

#### Generate Visualization Report
```sh
translonpredictor plot \
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
- `{outfile}_orfs_scored.csv`: Scored ORFs
- `{outfile}_report.html`: Visualization report

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
