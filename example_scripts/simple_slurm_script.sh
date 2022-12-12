#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-08:00:00
#SBATCH --array=[0-2]%3        ### in square parentheses, we specify a range of numerical indices - the final index will be one fewer than the number of jobs, because we start counting at 0.  After the %, we specify the batch size (i.e. max number of sbatch jobs that should run simultaneously). When there aren't many jobs, you can probably run them all at once. When you have a lot of jobs, or if each job is intensive, you might want to run only a limited number at once, e.g. 20 or 40.
#SBATCH --output=slurm.%J.out
#SBATCH --job-name="testSlurm"

#### load modules
# if we're using any of the software installed by Hutch scientific computing, we'll need to load the relevant modules. 
# First we run the 'source' line to set up the computer to work with modules
source /app/lmod/lmod/init/profile
# Then we actually load any modules we need, just like we would on the command line:
module load SAMtools/1.11-GCC-10.2.0

#### We set up an array of sample names (SAMPLE_IDS) we'll be working on one by one. Use regular parentheses, separate values by spaces, no quotes needed.
SAMPLE_IDS=(sample1 sample2 sample3)

#### We set up a variable called SINGLE_ID that represent one sample name. This helps us run the same task(s) for each sample
# you'll see we select those SINGLE_IDs from the SAMPLE_IDS array
# to perform the selection (using square brackets, []), we use a special variable called $SLURM_ARRAY_TASK_ID.  It derives from the "--array=[0-2]%3" statement on line 5. Each job is given a numerical index, in this case from 0 to 2.   We use that index to pull out an individual sample name from SAMPLE_IDS (in bash we start counting from 0)
SINGLE_ID="${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}"


#### specify input file names - filenames contain the sample names
BAM_FILE="${SINGLE_ID}.bam"

#### specify output file names - filenames contain the sample names
INDEX_FILE="${SINGLE_ID}.bam.bai"
STATS_FILE="${SINGLE_ID}.bwa.flagstats"

#### now actually run the code. Echo statements help us track what's going on - they'll appear in the slurm.JOBID.out files
echo "My sample: ${SINGLE_ID}"

# for long-running programs we might want to test whether the final output file already exists, in which case we don't need to do anything:
if test -f "${STATS_FILE}"; then
    echo "$STATS_FILE already exists"
else
    echo "running index and flagstat"
    samtools index ${BAM_FILE}
    samtools flagstat ${BAM_FILE} > ${STATS_FILE}
fi

# clean up (not really needed, but why not)
module purge
