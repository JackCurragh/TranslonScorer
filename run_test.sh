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

# Activate micromamba environment
echo "Activating TranslonScorer environment..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate TranslonScorer
check_status "Environment activation"

# Check for uncommitted changes
echo "Checking git status..."
if [[ -n $(git status -s) ]]; then
    echo "Uncommitted changes found. Committing changes..."
    
    # Stage and commit changes with generic message
    git add .
    git commit -m "Update: Automated test commit"
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
pip install git+https://github.com/JackCurragh/TranslonScorer-1
check_status "User installation"

echo "Testing user installation execution..."
python -c "from TranslonScorer.core import scoring, coordinates; from TranslonScorer.file_handlers import bam, bed, bigwig; print('Import test successful')"
check_status "User version test"

# If all tests pass, proceed with the actual run
echo "All installation tests passed. Proceeding with analysis..."

# Check if required files exist
BAM_FILE=~/Processed_test/star_align/bam/SRR25602018.Aligned.sortedByCoord.out.bam
CHROM_SIZES=data/chrom.sizes
GENOME_FA=data/genome.fa
GTF_FILE=data/MANE.gtf
OUTPUT_PREFIX=test_output/test

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

# Run the test script that exercises the new implementation
python -c "
from TranslonScorer.file_handlers import bam, bed, bigwig
from TranslonScorer.core import scoring, coordinates
import polars as pl

# Process BAM file
print('Processing BAM file...')
bam_df = pl.read_csv('$BAM_FILE')
exon_df = pl.read_csv('$GTF_FILE')

# Get exons and CDS
cds_df, exon_df = bam.getexons_and_cds('$GTF_FILE')

# Process BAM data
bam_type, _ = bam.detect_bam_type(bam_df, exon_df)
if bam_type == 'genomic':
    bam_df = bam.bamtranscript(bam_df, exon_df)
else:
    bam_df = bam.process_transcriptomic_bam(bam_df, cds_df)

# Calculate A-site positions
offsets = coordinates.change_point_analysis(bam_df)
bed_df = bed.asitecalc(bam_df, offsets)

# Convert to BigWig
bed.bedtobigwig('${OUTPUT_PREFIX}.bedGraph', '$CHROM_SIZES', '${OUTPUT_PREFIX}')

# Score ORFs
orfs_df = bigwig.scoring('${OUTPUT_PREFIX}.bw', exon_df, cds_df, False, 50)

# Save results
bed.saveorfsandexons(orfs_df, exon_df, '${OUTPUT_PREFIX}')

print('Analysis complete!')
"
check_status "TranslonScorer execution"

echo "Script completed successfully!"
echo "Output files can be found with prefix: $OUTPUT_PREFIX"

# Deactivate the environment
micromamba deactivate 