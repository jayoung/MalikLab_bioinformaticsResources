#!/bin/bash
source /app/lmod/lmod/init/profile

module purge
module load BBMap/38.91-GCC-10.2.0

## declare an array variable
declare -a FILES=(sample1.fastq.gz
                  sample2.fastq.gz)

## now loop through the above array
for i in "${FILES[@]}"
do
   echo "$i"
   OUT="$i.testformat.txt"
   testformat.sh -in=$i > $OUT
done

module purge
