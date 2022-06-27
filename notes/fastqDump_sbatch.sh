#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-04:00:00
#SBATCH --array=[0-8]%9 # how many jobs and batch size
#SBATCH --output=slurm.%J.out
#SBATCH --job-name="fastqDump"

ACCESSIONS=(SRR2927735 SRR2927736 SRR2927737 SRR2927738 SRR2927739 SRR2927740 SRR2927741 SRR2927742 SRR2927743)

CACHE="/fh/fast/malik_h/grp/public_databases/NCBI/SRA/cache/sra"
SINGLE_ACCESSION="${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}"
SRAFILE="${CACHE}/${ACCESSIONS[$SLURM_ARRAY_TASK_ID]}.sra"

echo "My accession: ${SINGLE_ACCESSION}"
if test -f "${SRAFILE}"; then
    echo "$SRAFILE already exists"
else
    echo ""
    echo "running prefetch"
    prefetch ${SINGLE_ACCESSION}
fi

echo ""
echo "running fastq-dump"
fastq-dump --split-spot --split-e --gzip  ${SRAFILE}
