# Translon Predictor

## Overview

`translonpredictor` is a command-line tool for translon calling. The process consists of processing BAM files, extracting and scoring ORFs from transcript sequences based on the annotation and codons provided. It supports multiple input file formats and generates output files in several formats including `.bedGraph`, `.bw`, `.html`, and `.csv`.

## Installation

To use this tool, you need to have Python and the necessary dependencies installed. Install the required packages using `pip`:

```sh
pip install click pandas pyBigWig
```
## Usage
The tool is invoked using the `translonpredictor` command. Below are the options and their descriptions:

```sh
Usage: translonpredictor [OPTIONS]

Options:
  -b, --bam TEXT              Provide a BAM file
  -c, --chromsize TEXT        Provide a file containing the chromosome sizes
  -s, --seq TEXT              Provide a file containing the genomic sequence (.fa)
  -t, --tran TEXT             Provide a file containing the transcript sequences (.fa)
  -a, --ann TEXT              Provide a file containing the annotation (.gtf)
  -sta, --starts TEXT         Provide a list of start codons (default: "ATG")
  -stp, --stops TEXT          Provide a list of stop codons (default: "TAA,TAG,TGA")
  -min, --minlen INTEGER      Provide the minimum length (default: 0)
  -max, --maxlen INTEGER      Provide the maximum length (default: 1000000)
  -bw, --bigwig TEXT          Provide a Bigwig file to convert
  -ex, --exon TEXT            Provide a file containing exon positions
  -bw, --bedfile TEXT         Provide a Bigwig file to convert
  -of, --orfs TEXT            Provide a file containing annotated ORFs
  -rp, --range_param INTEGER  Provide an integer for the plot range around the relative start position (default: 30)
  -sru, --sru_range INTEGER   Provide an integer for the Start Rise Up score range (default: 15)
  -ofs, --offsets TEXT        Provide a file containing offset parameters
  -s, --scoretype BOOLEAN     Select the scoring algorithm (default: False for old scoring algorithm)
  -pf, --plotfile TEXT        Provide a '.csv' file containing scored ORFs to use for plotting
  -ofn, --outfilename TEXT    Provide a name for the output files

  --help                      Show this message and exit.
```
## Examples
### Processing BAM File
To process a BAM file and generate the necessary outputs:
```sh
translonpredictor --bam example.bam --chromsize chrom.sizes --ann annotations.gtf --outfilename output_name
```
### Extracting and scoring ORFs from transcript sequences.
To extract and score ORFs from transcript sequences:

```sh
translonpredictor --seq genome.fa --ann annotations.gtf --outfilename output_name
```
The file containing transcript sequences can also be provided:
```sh
translonpredictor --tran transcripts.fa --ann annotations.gtf --outfilename output_name
```
### Generating a Report from a Plot File
To generate a report using a previously scored ORFs file:

```sh
translonpredictor --plotfile scored_orfs.csv --bigwig example.bw --exon exon_positions.csv --outfilename output_name
```
## Output Files
The tool generates several output files depending on the provided inputs:

.bedGraph files containing bedGraph formatted data.
.bw BigWig files.
.csv files with scored ORFs.
.html report containing translon information.

## Error Handling
The tool requires specific combinations of input files to function correctly. If the necessary files are not provided, it will raise an exception with guidance on the required files.

If you only wish to process a BAM file, an error will still pop up for the ORF handling. In this case the error can be ignored.

Please make sure to chromosome notation is the same for the annotation and the bigwig file.
## Contributing
Contributions are welcome. Please fork the repository and submit a pull request.

## License
This project is licensed under the MIT License.
```r
This README file provides an overview, installation instructions, usage examples, output desc
```