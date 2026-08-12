#!/bin/sh
#SBATCH --job-name="af3"
#SBATCH --partition=chorus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1024GB
#SBATCH --gpus=2


# Do not change these paths:

# path to model parameters inside the container
export MODEL_PATH=/root/public_databases/parameter_models/
export DOWNLOAD_DIR=/shared/biodata/alphafold3
SIF=/app/software/AlphaFold/containers/alpafold3.sif

#### User defined locations, please modify these:

## example input.json file can be found at
## https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#installation-and-running-your-first-prediction

## starting point from scicomp webpage:
# export BASE=$HOME/alphafold3_data
# export OUTPUT_DIR=$BASE/output
# export JSON_PATH=$BASE/example_input.json

## try this:
export BASE=/fh/fast/malik_h/user/jayoung/git_more_repos/MalikLab_bioinformaticsResources/janets_NOTES_forMyself/proteomics/alphafold/alphafold_3_scicomp_test/
export OUTPUT_DIR=$BASE/fold_output
# fold_input.json is the example input.json file can be found at
# https://github.com/google-deepmind/alphafold3?tab=readme-ov-file#installation-and-running-your-first-prediction
export JSON_PATH=$BASE/fold_input.json


mkdir -p $OUTPUT_DIR


apptainer exec \
  --nv \
  --bind $DOWNLOAD_DIR:/root/public_databases \
  --bind $BASE:$BASE \
  $SIF \
  python /app/alphafold/run_alphafold.py \
  --json_path=$JSON_PATH \
  --model_dir=$MODEL_PATH \
  --output_dir=$OUTPUT_DIR \
  --db_dir=/root/public_databases
