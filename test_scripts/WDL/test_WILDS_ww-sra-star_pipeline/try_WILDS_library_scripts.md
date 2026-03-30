# Try using WDL scripts to download and map reads

March 2026

Using `ww-sra-star.wdl` to download reads from SRA and then map using STAR

##  First attempt: doesn't work with single end reads

March 17, 2026 (and before)

I tried using WDL to download AND map reads, using the provided script ww-sra-star.wdl, but it doesn't seem to work well for single ends. Submitted an issue: https://github.com/getwilds/wilds-wdl-library/issues/281

Files are here:

```
~/public_databases/NCBI/SRA/data/mammalian_expression_profiles/mouse/mouse_CardosoMoreira_expressionSurvey/try_WILDS_library_scripts/old_WDLs
```

## Second attempt, after Taylor fixed the scripts - works!

He said "As such, we've made some updates in #290 to enable this functionality, and I think we've got it working! Would you mind downloading a copy of the ww-sra-star script from that branch and trying it out again with ERR2588371?"

Get the new scripts and try again - work in `git_branch_fix-sra-star-jyoung` subdir

```
cd ~/public_databases/NCBI/SRA/data/mammalian_expression_profiles/mouse/mouse_CardosoMoreira_expressionSurvey/try_WILDS_library_scripts/git_branch_fix-sra-star-jyoung/try1

wget https://raw.githubusercontent.com/getwilds/wilds-wdl-library/refs/heads/fix-sra-star-jyoung/pipelines/ww-sra-star/ww-sra-star.wdl

# the originals of these came from the github last week but I can tell they haven't been changed since
cp ../old_WDLs/small_mouse_test/cromwell-options.json .
cp ../old_WDLs/small_mouse_test/inputs.json  .

# I upped the max retries to 2

# submitted using PROOF workbench
# ID = 10b8f256-8daa-483b-af88-072f85be206d
# working dir = /hpc/temp/malik_h/user/jayoung/cromwell-scratch/sra_star/10b8f256-8daa-483b-af88-072f85be206d
```

Seems like it worked.

Outputs were copied to the `final_workflow_outputs_dir` I specified in `cromwell-options.json`, which is `.`

Check a couple of things
- flagstats on the bams
- readcounts on the fastq.gz

Check flagstat
```
run_SamtoolsFlagstats.pl *bam
```
Compare with older results - roughly the same, within the boundaries of what I'd expect given that I changed the parameters in my STAR mapping
```
more ERR2588370*stats 
more ~/public_databases/NCBI/SRA/data/mammalian_expression_profiles/mouse/mouse_CardosoMoreira_expressionSurvey/STAR_mm39_maxMultiHits1/ERR2588370_STAR/*stats 

more ERR2588371*stats 
more ~/public_databases/NCBI/SRA/data/mammalian_expression_profiles/mouse/mouse_CardosoMoreira_expressionSurvey/STAR_mm39_maxMultiHits1/ERR2588371_STAR/*stats
```


Check readcounts
```
# because one sample was cached and the other wasn't, these are differently nested:
cp  /hpc/temp/malik_h/user/jayoung/cromwell-scratch/sra_star/10b8f256-8daa-483b-af88-072f85be206d/call-fastqdump/*/*/*gz .
cp  /hpc/temp/malik_h/user/jayoung/cromwell-scratch/sra_star/10b8f256-8daa-483b-af88-072f85be206d/call-fastqdump/*/*/*/*gz .

countFastqReads.pl *gz
```

Compare with readcounts.txt I got the other way - looks good.
```
grep 'ERR258837[01]' ~/public_databases/NCBI/SRA/data/mammalian_expression_profiles/mouse/mouse_CardosoMoreira_expressionSurvey/fastqFiles/readcounts.txt 
ERR2588370.fastq.gz             15986818
ERR2588371.fastq.gz             41188903
```

I removed the output files:

```
rm ER* readcounts.txt
```