#!/bin/bash

# Exit on any error
set -e

echo "Starting TranslonScorer test script..."

# Function to check if last command was successful
check_status() {
    if [ $? -eq 0 ]; then
        echo "✓ $1 successful"
    else
        echo "✗ Error during $1"
        exit 1
    fi
}

# Function to check if a file exists and is readable
check_file() {
    if [ -r "$1" ]; then
        echo "✓ Found $2: $1"
        return 0
    else
        echo "✗ $2 not found: $1"
        return 1
    fi
}

# Activate micromamba environment
echo "Activating TranslonScorer environment..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate TranslonScorer
check_status "Environment activation"

# Check for uncommitted changes
echo "Checking git status..."
if [[ -n $(git status -s) ]]; then
    echo "Uncommitted changes found. Committing changes..."
    
    # Stage and commit changes with detailed message
    git add .
    git commit -m "refactor: make pipeline more flexible

- Update test script to handle both BAM and BigWig inputs
- Add automatic input detection and pipeline selection
- Improve error handling and file validation
- Add progress tracking and status checks"
    check_status "Git commit"
    
    # Push changes
    echo "Pushing changes to remote..."
    git push
    check_status "Git push"
else
    echo "No changes to commit"
fi

# Ensure we're starting with a clean slate
echo "Uninstalling any existing versions..."
pip uninstall -y TranslonScorer || true
check_status "Uninstall"

# First, test the development version
echo "Testing development version installation..."
pip install -e .
check_status "Development installation"

echo "Testing development version execution..."
python -c "from TranslonScorer.core import scoring, coordinates; from TranslonScorer.file_handlers import bam, bed, bigwig; print('Import test successful')"
check_status "Development version test"

# Clean up development installation
pip uninstall -y TranslonScorer
check_status "Development version cleanup"

# Now test as a user would experience it
echo "Testing user installation from GitHub..."
pip install --no-cache-dir --force-reinstall  git+https://github.com/JackCurragh/TranslonScorer
check_status "User installation"

echo "Testing user installation execution..."
python -c "import TranslonScorer; print('Package path:', TranslonScorer.__path__[0]); import os; print('Files:', os.listdir(TranslonScorer.__path__[0]))"

python -c "from TranslonScorer.core import scoring, coordinates; from TranslonScorer.file_handlers import bam, bed, bigwig; print('Import test successful')"
check_status "User version test"

# If all tests pass, proceed with the actual run
echo "All installation tests passed. Proceeding with analysis..."

# Define file paths
BAM_FILE=~/Processed_test/star_align/bam/SRR25602018.Aligned.sortedByCoord.out.bam
BIGWIG_FILE=test_output/test.bw
CHROM_SIZES=data/chrom.sizes
GENOME_FA=data/genome.fa
GTF_FILE=data/MANE.gtf
OUTPUT_PREFIX=test_output/test_orfs
FW_BIGWIG=data/all_forward.bigWig
RV_BIGWIG=data/all_reverse.bigWig

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
mkdir -p "$OUTPUT_DIR"

# Check which files exist and determine pipeline path
echo "Checking available input files..."

# BAM_EXISTS=0
# BIGWIG_EXISTS=0
# check_file "$BAM_FILE" "BAM file" && BAM_EXISTS=1
# check_file "$BIGWIG_FILE" "BigWig file" && BIGWIG_EXISTS=1
# check_file "$GENOME_FA" "Genome FASTA" || exit 1
# check_file "$GTF_FILE" "GTF annotation" || exit 1
# check_file "$FW_BIGWIG" "Forward BigWig" && FW_BIGWIG_EXISTS=1
# check_file "$RV_BIGWIG" "Reverse BigWig" && RV_BIGWIG_EXISTS=1

PROFILE=1 memray run /Users/jackt/mamba/envs/TranslonScorer/bin/translonscorer \
--forward_bigwig $FW_BIGWIG \
--reverse_bigwig $RV_BIGWIG \
--sequence "$GENOME_FA" \
--annotation "$GTF_FILE" \
--output "$OUTPUT_PREFIX" \
--sru-range 15 \
--stranded
# if [ $FW_BIGWIG_EXISTS -eq 1 ] & [ $RV_BIGWIG_EXISTS -eq 1 ]; then
#     echo "Forward and reverse BigWig file found, using direct ORF finding path..."
#     echo "Running TranslonScorer with memory monitoring..."
#     PROFILE=1 memray run /Users/jackt/mamba/envs/TranslonScorer/bin/translonscorer all \
#     --forward_bigwig $FW_BIGWIG \
#     --reverse_bigwig $RV_BIGWIG \
#     --sequence "$GENOME_FA" \
#     --annotation "$GTF_FILE" \
#     --outfile "$OUTPUT_PREFIX" \
#     --scoring-method modern \
#     --sru-range 15 \
#     --stranded

# elif [ $BIGWIG_EXISTS -eq 1 ]; then
#     echo "BigWig file found, using direct ORF finding path..."
#     echo "Running TranslonScorer with memory monitoring..."
#     PROFILE=1 memray run /Users/jackt/mamba/envs/TranslonScorer/bin/translonscorer all \
#         --sequence "$GENOME_FA" \
#         --annotation "$GTF_FILE" \
#         --bigwig_path "$BIGWIG_FILE" \
#         --outfile "$OUTPUT_PREFIX" \
#         --scoring-method modern \
#         --sru-range 15
# elif [ $BAM_EXISTS -eq 1 ]; then
#     if ! check_file "$CHROM_SIZES" "Chromosome sizes"; then
#         echo "Error: Chromosome sizes file required for BAM processing"
#         exit 1
#     fi
#     echo "BAM file found, using full pipeline path..."
#     echo "Running TranslonScorer with memory monitoring..."
#     PROFILE=1 translonscorer all \
#         --bam_path "$BAM_FILE" \
#         --chromsizes "$CHROM_SIZES" \
#         --sequence "$GENOME_FA" \
#         --annotation "$GTF_FILE" \
#         --outfile "$OUTPUT_PREFIX" \
#         --scoring-method modern \
#         --sru-range 15
# else
#     echo "Error: Neither BAM nor BigWig file found"
#     exit 1
# fi

echo "Script completed successfully!"
echo "Output files can be found with prefix: $OUTPUT_PREFIX"

# Deactivate the environment
micromamba deactivate 