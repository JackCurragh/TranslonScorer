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

# Check for uncommitted changes
echo "Checking git status..."
if [[ -n $(git status -s) ]]; then
    echo "Uncommitted changes found. Please enter a commit message:"
    read -r commit_msg
    
    # Stage and commit changes
    git add .
    git commit -m "$commit_msg"
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
pip uninstall -y TRANSLONSCORER || true
pip uninstall -y translonscorer || true
check_status "Uninstall"

# First, test the development version
echo "Testing development version installation..."
pip install -e .
check_status "Development installation"

echo "Testing development version execution..."
translonscorer --help
check_status "Development version test"

# Clean up development installation
pip uninstall -y TranslonScorer
check_status "Development version cleanup"

# Now test as a user would experience it
echo "Testing user installation from GitHub..."
pip install git+https://github.com/JackCurragh/TranslonScorer
check_status "User installation"

echo "Testing user installation execution..."
translonscorer --help
check_status "User version test"

# If all tests pass, proceed with the actual run
echo "All installation tests passed. Proceeding with analysis..."

# Check if required files exist
BAM_FILE=~/Processed_test/star_align/bam/SRR25602018.Aligned.sortedByCoord.out.bam
CHROM_SIZES=TranslonScorer/data/chrom.sizes
GENOME_FA=TranslonScorer/data/genome.fa
GTF_FILE=TranslonScorer/data/MANE.gtf
OUTPUT_PREFIX=TranslonScorer/data/test

# Check input files
for file in "$BAM_FILE" "$CHROM_SIZES" "$GENOME_FA" "$GTF_FILE"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required file not found: $file"
        exit 1
    fi
done

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
mkdir -p "$OUTPUT_DIR"

echo "Running TranslonScorer..."
translonscorer all \
    -b "$BAM_FILE" \
    -c "$CHROM_SIZES" \
    -s "$GENOME_FA" \
    -a "$GTF_FILE" \
    -o "$OUTPUT_PREFIX"
check_status "TranslonScorer execution"

echo "Script completed successfully!"
echo "Output files can be found with prefix: $OUTPUT_PREFIX" 