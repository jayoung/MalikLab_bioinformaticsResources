#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-08:00:00
#SBATCH --array=[0-2]%3        ### in square parentheses, we specify a range of numerical indices - the final index will be one fewer than the number of jobs, because we start counting at 0.  After the %, we specify the batch size (i.e. max number of sbatch jobs that should run simultaneously). When there aren't many jobs, you can probably run them all at once. When you have a lot of jobs, or if each job is intensive, you might want to run only a limited number at once, e.g. 20 or 40.
#SBATCH --output=slurm.%J.out
#SBATCH --job-name="bwa_slurm_rescue_selection"

#### load modules
source /app/lmod/lmod/init/profile
module load SAMtools/1.11-GCC-10.2.0
module load BWA/0.7.17-GCCcore-11.2.0

#### specify an array of sample names we'll be working on one by one
SAMPLE_IDS=(sample1 sample2 sample3)

#### set up a variable that represents INDIVIDUAL sample names
# $SLURM_ARRAY_TASK_ID is a magical variable - it derives from the "--array=[0-2]%3" line above, each job getting an index, in this case 0, 1 and 2.   We use that index to pull out an individual sample name from SAMPLE_IDS (in bash we start counting from 0)
SINGLE_ID="${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}"

REF_GENOME="/fh/fast/malik_h/grp/public_databases/UCSC/fly_Aug2014_dm6/dm6.fa_bwaFormat/dm6.fa"

#### specify input file names - filenames contain the sample names
R1_FILE="../${SINGLE_ID}_R1.fq.gz"
R2_FILE="../${SINGLE_ID}_R2.fq.gz"

#### specify output file names - filenames contain the sample names
BAM_FILE="${SINGLE_ID}.bwa.bam"
COUNTS_FILE="${SINGLE_ID}.bwa.bam.counts"
LOGS_FILE="${SINGLE_ID}.bwa.logs.txt"
FINAL_OUTPUT_FILE="${SINGLE_ID}.bwa.bam.flagstats"

#### now actually run the code
echo "My sample: ${SINGLE_ID}"

# if the final output file already exists, we don't need to do anything
if test -f "${FINAL_OUTPUT_FILE}"; then
    echo "$FINAL_OUTPUT_FILE already exists"
else
    echo ""
    echo "running bwa"
    # run bwa mem to map paired reads to ref genome (output = sam format), use samtools view to convert sam format to bam format, and use samtools sort to sort by position in the reference genome. Save any screen output to the log file.
    (bwa mem -t 4 ${REF_GENOME} ${R1_FILE} ${R2_FILE} | samtools view -@ 4 -Sb -  | samtools sort -@ 4 -O bam > ${BAM_FILE}) 2>> ${LOGS_FILE}
    # index the resulting sorted bam file
    samtools index ${BAM_FILE}
    # use samtools flagstat to see how well the reads mapped
    samtools flagstat ${BAM_FILE} > ${FINAL_OUTPUT_FILE}
fi

# clean up (not really needed, but why not)
module purge
