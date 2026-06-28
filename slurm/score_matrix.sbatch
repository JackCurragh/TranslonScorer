#!/bin/bash
#SBATCH --job-name=ts-score-matrix
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=ts-score-matrix-%j.log

source /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/translonscorer_fresh/.venv/bin/activate

translonscorer score-matrix \
  --events-dir /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/translonscorer_6k/events \
  --matrix-dir /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/6k_runs/global_partitioned \
  --gtf /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/gencode.v45.annotation.gtf \
  --psite-index /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/translonscorer_6k/psite_index_filtered \
  --store-dir /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/translonscorer_6k/scores_pif \
  --data-version pif_v1 \
  --annotation-version gencode_v45 \
  --n-workers 32
