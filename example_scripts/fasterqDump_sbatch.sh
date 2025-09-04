#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-04:00:00
#SBATCH --array=[0-1]%2 # how many jobs and batch size
#SBATCH --output=slurm.%J.out
#SBATCH --job-name="fasterqDump"

## This script uses two tools from NCBI's SRA-Toolkit (prefetch and fasterq-dump) to download public fastq.gz files from NCBI.  

## Run this script by typing:
# sbatch fasterqDump_sbatch.sh

## More information on fasterq dump from NCBI: https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

## For a while there were problems that meant we could not use fasterq-dump on gizmo/rhino (info here - https://sciwiki.fredhutch.org/scicompannounce/2023-09-01-fasterq-dump-issue/). But it sounds like that's no longer a problem.  If problems recur, see script called fastqDump_sbatch.old.sh that uses the older fastq-dump algorithm instead of fasterq-dump


## load SRA-Toolkit module
source /app/lmod/lmod/init/profile
module load SRA-Toolkit/3.1.1-gompi-2023b 

## define all accessions that we want
ACCESSIONS=(SRR9290530 SRR9290503)

## choose single SRA accession, one for each slurm job, and put its name in a variable
SINGLE_ACCESSION="${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}"

## figure out what the .sra intermediate file will be called
SRAFILE="${SINGLE_ACCESSION}/${SINGLE_ACCESSION}.sra"


##### run prefetch, if the sra file does not already exist
echo ""
echo "### My accession: ${SINGLE_ACCESSION}"
if test -f "${SRAFILE}"; then
    echo "$SRAFILE already exists"
else
    echo ""
    echo "# Running prefetch"
    prefetch ${SINGLE_ACCESSION}
fi

##### run fasterq-dump 
echo ""
echo "# Running fasterq-dump"
fasterq-dump -e 4 ${SRAFILE}
# -e is num threads. if you change it, you should also adjust --cpus-per-task at the top of the script
# --split-3 is the default

##### compress output
R1_FILE="${SINGLE_ACCESSION}_1.fastq"
R2_FILE="${SINGLE_ACCESSION}_2.fastq"
UNPAIRED_FILE="${SINGLE_ACCESSION}.fastq"

echo ""
echo "# Compressing output"
if test -f "${R1_FILE}"; then
    gzip "${R1_FILE}"
fi
if test -f "${R2_FILE}"; then
    gzip "${R2_FILE}"
fi
if test -f "${UNPAIRED_FILE}"; then
    gzip "${UNPAIRED_FILE}"
fi


##### clean up - remove the dir containing the .sra file
echo ""
echo "# Cleaning up"
rm -r ${SINGLE_ACCESSION}