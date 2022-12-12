#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-08:00:00
#SBATCH --array=[0-2]%3 
#SBATCH --output=slurm.%J.out
#SBATCH --job-name="testSlurm"

#### load module(s)
source /app/lmod/lmod/init/profile
module load SAMtools/1.11-GCC-10.2.0

#### define sample names
SAMPLE_IDS=(sample1 sample2 sample3)

#### choose sample for each job
SINGLE_ID="${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}"

#### specify input/output file names
BAM_FILE="${SINGLE_ID}.bam"
INDEX_FILE="${SINGLE_ID}.bam.bai"
STATS_FILE="${SINGLE_ID}.bwa.flagstats"

#### run commands
samtools index ${BAM_FILE}
samtools flagstat ${BAM_FILE} > ${STATS_FILE}

