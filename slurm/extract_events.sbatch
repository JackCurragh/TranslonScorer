#!/bin/bash
#SBATCH --job-name=ts-extract-events
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=ts-extract-events-%j.log

source /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/translonscorer_fresh/.venv/bin/activate

translonscorer extract-events \
  --bigbed /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/translonscorer_6k/input/TransCODE_GENCODE_v45_RiboSeqORFs_TransCODE_GENCODEv45_RiboseqORF_catalog_comprehensive.bigBed \
  --context-gtf /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/gencode.v45.annotation.gtf \
  --out-dir /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/translonscorer_6k/events
