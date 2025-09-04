#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-04:00:00
#SBATCH --array=[0-1]%2 # how many jobs and batch size
#SBATCH --output=slurm.%J.out
#SBATCH --job-name="fastqDump"

## This script uses two tools from NCBI's SRA-Toolkit (prefetch and fastq-dump) to download public fastq.gz files from NCBI.  

## Run this script by typing:
# sbatch fastqDump_sbatch.sh

## fasterq-dump - possibly a better way...
# There is a newer, faster tool called fasterq-dump. I haven't used it much because for a while there were problems that meant we could not use it on gizmo/rhino (info here - https://sciwiki.fredhutch.org/scicompannounce/2023-09-01-fasterq-dump-issue/). But it sounds like that's no longer a problem
# more information on fasterq dump from NCBI: https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump


source /app/lmod/lmod/init/profile
module load SRA-Toolkit/3.1.1-gompi-2023b 

## define all accessions that we want
ACCESSIONS=(SRR9290530 SRR9290503)

## choose single SRA accession, one for each slurm job, and put its name in a variable
SINGLE_ACCESSION="${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}"

## figure out what the .sra intermediate file will be called
SRAFILE="${SINGLE_ACCESSION}/${SINGLE_ACCESSION}.sra"


## run prefetch, if the sra file does not already exist
echo "My accession: ${SINGLE_ACCESSION}"
if test -f "${SRAFILE}"; then
    echo "$SRAFILE already exists"
else
    echo ""
    echo "running prefetch"
    prefetch ${SINGLE_ACCESSION}
fi

## run fastq-dump 
echo ""
echo "running fastq-dump"
fastq-dump --split-spot --split-3 --gzip  ${SRAFILE}
