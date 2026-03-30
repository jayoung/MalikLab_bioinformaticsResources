# janet-learning-WDL
notes on learning WDL/Cromwell

# Useful links

## Tutorials, documentation

[OpenWDL](https://openwdl.org) has links to tutorials

[Cromwell docs](https://cromwell.readthedocs.io/en/stable/) - Amy says some stuff there is out of date

Learning WDL [readthedocs](https://wdl-docs.readthedocs.io/en/stable/) (recommended by Amy Dec 2022)

WDL [version 1.0 specs](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md)

[JSON format](https://en.wikipedia.org/wiki/JSON#Syntax) specification

More tutorials listed (here)(https://sciwiki.fredhutch.org/compdemos/Cromwell/#additional-wdl-resources)

Dave Tang's [Learning WDL](https://davetang.org/muse/2020/01/09/learning-wdl/) blog post

[WDL Best Practices / Style Guide](https://gist.github.com/scottfrazer/aa4ab1945a6a4c331211)

## Hutch WDL tutorials (and other repos)

Starting a Cromwell server, using Shiny app, etc: [FH_WDL101](https://hutchdatascience.org/FH_WDL101_Cromwell/introduction.html) 

Hutch [Cromwell server configuration repo](https://github.com/FredHutch/diy-cromwell-server) (NOT intended for me to mess with). Used in cromwell.sh for server setup.

Amy's [example/test workflows](https://github.com/FredHutch/wdl-test-workflows/).  I have a clone that I edited for clarity, in `~/FH_fast_storage/cromwell-home/janet-learning-WDL/wdl-test-workflows`

More [example workflows](https://github.com/FredHutch/reproducible-workflows) from Amy.  I cloned it in `~/FH_fast_storage/cromwell-home/janet-learning-WDL`

Learning WDL [FH_WDL102](https://hutchdatascience.orghttps://hutchdatascience.org/FH_WDL102_Workflows/) 

## Hutch WDL-related links/locations

Shiny job submission/monitoring [dashboard](https://cromwellapp.fredhutch.org) 

My cromwell home dir: `~/FH_fast_storage/cromwell-home`

Cromwell scratch dir: `/fh/scratch/delete90/malik_h/jayoung`

My current Cromwell server details will be here: `~/FH_fast_storage/cromwell-home/README.md`


### Places to copy WDLs/tasks from:

Broad's [WARP repository](https://github.com/broadinstitute/warp): WDL Analysis Research Pipelines (WARP) repository is a collection of cloud-optimized pipelines for processing biological data from the Broad Institute Data Sciences Platform and collaborators.

[Dockstore](https://www.dockstore.org/search?descriptorType=WDL&entryType=workflows&searchMode=files)

Search the Hutch github repo with 'wdl'

## Things to check out

[pipeline-builder](https://github.com/epam/pipeline-builder)  visualizes pipelines, maybe even creates WDL?

WDL [functions](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#standard-library)

# To do

finish going through tutorials

understand better when Cromwell uses cached results or not


# My learning progress

I think I understand how to run WDLs on our cluster: I worked through the [Hutch WDL101 course](https://hutchdatascience.org/FH_WDL101_Cromwell/index.html) in detail, providing lots of feedback

I want to look at [OpenWDL tutorials](https://github.com/openwdl/learn-wdl) too

My first real WDL is working (`~/FH_fast_storage/bat_reproduction/wdl_scripts/dnaseq_fq_to_vcf.wdl`)

It's going well.  One thing I'm confused about is when it decides to re-run tasks versus reuse existing results. It's rerunning more often than I would think.  Is it to do with the anonymized scatters?  i.e. some tasks get called ScatterAt93_15 rather than call-TaskName

xxx think about how to handle copying output files out of scratch

xxx I want to try jobname as a runtime variable

# WDL-related habits I want to think about

How to handle copying output files out of scratch, given that directory structure is not what I'm used to. See how it looks once I get my bat-Orr-batch2 WDL running

Where to keep the wdl+jsons - I think IN the project repo.  One JSON per analysis batch?

How to split work - one giant WDL script, or break it up into separate ones?

Naming workflows - the name used in the WDL file's workflow block is used elsewhere too, e.g. it's what the output dir in /fh/scratch is called, and is found in the JSON input files. To me it makes sense to use the same name for the wdl script file, too. I probably want to try NOT to reuse the workflow block name in different projects, as the output dirs will get muddled up and they'll try to reuse each others' results.


# WDL syntax questions
- indents don't matter: is that correct?  is that true for WDL files AND json files?

# WDL syntax notes

The order of `call` statements in the `workflow{}` block does NOT matter.  the WDL interpreter magically figures out which tasks depend on each other and runs any that it can run at appropriate times.

## Command block
The `command` block can be enclosed with `<<< >>>` or `{ }`.   If the command itself might use `{}` (e.g. perl or python stuff) then the `<<< >>>` notation is better.   When we use `<<< >>>`, variables used in the command block are in the form `~{var}`.  If we use the {} format, variables could also be specified as `${}`. This all seems confusing: I think I will stick to this format: `command <<< ~{var} >>>`.

## String manipulation

[Basename](https://wdl-docs.readthedocs.io/en/stable/WDL/basename/) gets a file name without the path.  Simplest usage:
```
File input_file = "/Users/chris/input.bam"
String base = basename(input_file)
    # Result: input.bam
```
Can also strip off a specified extension:
```
File input_file = "/Users/chris/input.bam"
String stripped = basename(input_file, ".bam") 
    # Result: input
```

## Input data files (e.g fq.gz)

Specifying input data files that are actually links: I thought I might do this so that my fastq file aliasing would be used.  But it seems like Cromwell follows the link through to its source and uses the original filename instead. So my trick of using a link to make a shorter file alias is not useful.

## Using directories as intermediate results files

Turns out you can't pass around directories between tasks (at least not right now, in version 1.0 of WDL+Cromwell) - you have to tar.gz them. This seems like an unnecessary step when running on a local cluster (but sounds like there's a good reason to do it when running on the cloud - hard to pass around directories, easier to pass around files)

See `~/FH_fast_storage/cromwell-home/janet-workflow-tests/janet_caching_question/janet_caching_question.md`

## Specifying workflow inputs

Specify inputs using:
- JSON file
- tab-delimited files: see [this example](https://github.com/FredHutch/wdl-test-workflows/tree/main/localBatchFileScatter)
- hard-coded within the WDL script.

When I provide inputs in a JSON file, I need to tell the WDL what format each one is in, by including it in the `input{}` block like this:
```
input {
    Array[String] allSamples
    Map[String, Array[String]] samplesToPairs
    Map[String, Map[String,File]] pairsToFq
}
```

Hard-coded inputs/variables do NOT go in the input{} block, and they look like this:
```
Array[String] allSamples = ["s1", "s2"]
Map[String, Array[String]] samplesToPairs = {
   "s1": ["s1_pair1"],
   "s2": ["s2_pair1","s2_pair2"]
}
```

JSON inputs: it seems like it's OK to include a variable in the JSON file that's does not have a format declared in the WDL file, either directly or via a `struct` (structure). The variable names in `struct` DO need to match the names in the JSON, but they do NOT need to be in the same order.

JSON format: there's no way to include true comments, but we can use something like this so that comments masquerade as data:
```
  "##_COMMENT1": "INPUT BAM",
```

To merge two or more json files from the command line, we use [`jq`](https://stedolan.github.io/jq)
```
jq -s '.[0] * .[1]' file1 file2 > out.json
jq -s '.[0] * .[1] * .[2]' file1 file2 file3 > out.json
```

## Scattering

Scattering over a map can be done - it is described [here](https://bioinformatics.stackexchange.com/questions/16100/extracting-wdl-map-keys-as-a-task) and [here](https://github.com/openwdl/wdl/issues/106#issuecomment-356047538). It may or may not be possible in v1.0 of WDL - there were some changes in v1.1.

## Map structures

Nested map structures might be possible, or might not - see [here](https://bioinformatics.stackexchange.com/questions/14673/reading-nested-map-data-structures-in-wdl)

## Outputs

To get output files properly listed, if I don't know what they're called ahead of time, I may at some point need to use a trick involving glob - see [example](https://github.com/FredHutch/tg-wdl-cellRanger/blob/main/cellRanger.wdl): `Array[File] outputDir = glob("./countrun/outs/*")`


## modularizing code

See example [here](https://github.com/theiagen/terra_utilities/blob/main/workflows/wf_cat_column.wdl)

Can put task code blocks in a separate file for reuse, and import. e.g.
```
import "../tasks/task_file_handling.wdl" as file_handling
```
after which a task called `cat_files` (found in `task_file_handling.wdl`) is available in the importing workflow via the name `file_handling.cat_files`

# Misc notes

WDL arrays are 0-indexed

If the cromwell server times out while a workflow is running, jobs do continue.  "Note: when servers go down, all jobs they'd sent will continue.  When you start up a server the next time using the same database, the new server will pick up whereever the previous workflows left off."

Unix `tree` command: for cromwell dirs, it's often useful to only show some number of levels
```
tree -L 5 
```

## VScode extensions

(installed these on my laptop but maybe not the work desktop)
- WDL DevTools 
- WDL syntax highlighter 
- Prettify Json 
- JSON Tools 

## temporary dirs used by Cromwell

When running through cromwell, the tmpdir is set in each script, and is always a subdir of wherever cromwell-executions lives (in my case this is on `/fh/scratch`). 

When running a plain sbatch job, the tmpdir is elsewhere, e.g. `/loc/scratch/7617088`.   

I want to deliberately set the tmp dir for picard.jar jobs (e.g. MarkDuplicates), because they seem to run a lot faster using `/loc/scratch/` dirs than they do using subdirs of `/fh/scratch`.  I set that in two places: for java using `-Djava.io.tmpdir`, and using the `-TMP_DIR` argument to `picard.jar`
```
java -Xmx10g -Djava.io. -Djava.io.tmpdir=$SCRATCH_LOCAL -jar $EBROOTPICARD/picard.jar MarkDuplicates -TMP_DIR $SCRATCH_LOCAL
```

## job submission at the Hutch

For a more sophisticated way to interact with a Cromwell server: "The Cromwell server has an API so the app and the R package are just wrappers for talking to the API. if you start up a server, and then in your browser go to http://gizmok1:55555 (or whatever), you'll see the SWagger UI and you can futz around all you like with the endpoint definitions. Then you can use any old command line tool to send workflows or do any of this stuff."  
I THINK you would use that swagger UI to produce an examplecurl command (it also submits the workflow). Then maybe on another occasion, you could copy and modify that curl command, maybe playing with additional options listed in the swaggerUI, and then just copy-paste the curl command onto the rhino/gizmo command line to submit an additional job?


# womtool.jar

## Validation
```
module load cromwell/84

# with a single json
java -jar $EBROOTCROMWELL/womtool.jar validate --inputs dnaseq_fq_to_vcf.consolidatedInputs.v2.small_tinyRegions.json dnaseq_fq_to_vcf.wdl 

# I don't (yet) know how to do this when I have subworkflow dependencies (this workflow requires dnaseq_fq_to_vcf.skeleton.subworkflows_bundle.zip)
java -jar $EBROOTCROMWELL/womtool.jar validate --inputs dnaseq_fq_to_vcf_skeleton.consolidatedInputs.v2.small_tinyRegions.json dnaseq_fq_to_vcf_skeleton.wdl 

# maybe I validate individual subworkflows, but in this case it cannot find the inputs (not surprising, as they pass through from the main workflow)
java -jar $EBROOTCROMWELL/womtool.jar validate --inputs dnaseq_fq_to_vcf_skeleton.consolidatedInputs.v2.small_tinyRegions.json dnaseq_fq_to_vcf.skeleton.subworkflow.eachPair.wdl 

module purge
```

## workflow graph

To plot structure of a WDL workflow as a graph:
```
module load cromwell/84
java -jar $EBROOTCROMWELL/womtool.jar graph dnaseq_fq_to_vcf.wdl | dot -Tpng > dnaseq_fq_to_vcf.graph.png
module purge
```
or a much more detailed version (in my case, far too detailed to be useful)
```
java -jar $EBROOTCROMWELL/womtool.jar womgraph dnaseq_fq_to_vcf.wdl | dot -Tpng > dnaseq_fq_to_vcf.graphAll.png
```

# When does a job rerun?

It was able to reuse old results when the only change to a task was to add a memory request to a runtime block (`memory: "12GB"`).  

This [page](https://cromwell.readthedocs.io/en/stable/cromwell_features/CallCaching/) explains that while changes in some runtime variables WOULD trigger a task to be rerun (ContinueOnReturnCode, Docker, FailOnStderr), other runtime variables (including `memory`, `cpu`, and `disks`) can be changed and the cache can still be used.

Interesting - I copy-pasted a task from one workflow to another.  Cromwell was able to figure out that I had run that task in an identical way in a totally different workflow, and it re-used the results.
```
task hostname {
  command {
    echo $(hostname)
    echo "output of which codeml"
    which codeml
    echo "done"
  }
  output {
    File out = stdout()
  }
}
```
That brings up a question - Cromwell server is storing a lot of old results. Does the search through the database of cached results ever get so slow that it's impractical?  Should I be clearing the database occasionally?

xxxx I want to do a test.  let's say I have a pipeline, I run it on two samples.  I then edit the json to add a third sample but I don't change anything else. does cromwell rerun the first two samples, or not?  is the answer the same regardless of sample order in the json?


# Failure Modes

Documented [here](https://cromwell.readthedocs.io/en/stable/execution/ExecutionTwists/): 

"Cromwell supports two failure modes, which specify how Cromwell behaves when a job fails during the execution of a workflow.
`NoNewCalls` (default)
Cromwell does not start any new call as soon as a job fails. Cromwell will still monitor the rest of the jobs until they complete (successfully or not).
`ContinueWhilePossible`
Cromwell attempts to run as many jobs as possible until no more can be started. When all running jobs are complete, the workflow fails."

# Memory requests

On our cluster, [slurm does NOT let us control memory allocations](https://sciwiki.fredhutch.org/scicomputing/compute_jobs/#memory):  

"Currently memory (or RAM) is not scheduled by Slurm. This means that requesting memory has little effect on the job or its available resources. Memory is currently only advisory: Slurm will only ensure that the node allocated to the job has more memory installed than the amount requested by the job- it does not look at memory availability or what is consumed by yours or other jobs on the node.

When your job needs “a lot” of memory use CPUs as a proxy for the memory you expect to be needed. If you think your job will need more than 4GB of memory, request one CPU for every 4GB required. For example, if you think your job will need 6GB of RAM, you would request 2 CPUs (adjust upward when the desired memory isn’t a multiple of four).

If you still want to add a memory request, use the --mem option. This option takes an argument: a number indicating the amount of memory required on the node. The default unit is megabytes- to specify the unit, append K, M, G, or T for kilobytes, megabytes, gigabytes, or terabytes."

So in WDL I will specify memory AND cpu.  
```
    memory: "12GB"
    cpu: 3    ##  memory ~ 4 x cpu
```
Memory will be ignored (unless it's super-big and would affect which nodes I might be allocated to) - it's CPUs that will be more likely to play into whether memory is actually available.


xxx for gatk:  what's the difference between -Xmx10g and -Xms10g ?

The flag Xmx specifies the maximum memory allocation pool for a Java Virtual Machine (JVM), while Xms specifies the initial memory allocation pool. The Xms flag has no default value, and Xmx typically has a default value of 256 MB. A common use for these flags is when you encounter a java.lang.OutOfMemoryError.

This means that your JVM will be started with Xms amount of memory and will be able to use a maximum of Xmx amount of memory. For example, starting a JVM like below will start it with 256 MB of memory and will allow the process to use up to 2048 MB of memory:

java -Xms256m -Xmx2048m

xxx where's that page that tells which runtime variables make cromwell think it needs to rerun a job?

# Notes on conversations, tutorials, etc

## Video "Deciphering a mystery workflow written in WDL" (Broad)

[youtube video](https://terra.bio/deciphering-a-mystery-workflow-written-in-wdl/)

look for two things in the code: 
- workflow block (links the tasks)
- task block - contains 'command block'

womtools (a java script) can analyze a WDL and make a DOT file, and graphviz (online tool, or can install it) can visualize DOT files. In those graphs:
- ovals = call a task
- hexagon = operater (modify/control execution of task)
- boxes = scope of control/operator
- dashed lines = conditional (we don't always do it)
- arrows = input-to-output relationstip

can nest wdl files (workflows within workflows) and the graph can handle that.  Makes double-walled ovals for steps that are entire workflows


## Conversation with Amy Oct 14 2022

This was motivated by me trying to understand how to run the example workflows, and failing.

The Shiny app's 'validate workflow': useful!  Can validate with or without the input json files (use `cat` if there's >1 input json)

When submitting jobs via the Shiny app: order of input json files does not matter.

On the 'track jobs' tab of the shiny app, when a job is selected, the 'workflow specific job information' are shows the jobs spun off by the workflow.  

If the workflow entirely failed to run, then copy the workflow ID, go back to the 'submit jobs' tab, and put that ID into the 'troubleshoot workflow' box.  Work from the bottom up of that output to see where the errors were.

If the workflow did start running but did not complete, there may be some helpful info on which steps failed in that 'track jobs' tab, or it might be easier to go to the workflow output dir and look at the stderr etc files there

Broad uses the google cloud - some of their WDL files have a few google-specific things in them

If some tasks of a workflow fail, resubmit it and Cromwell SHOULD figure out which things to rerun and which to copy from before. Default we think is to make soft links? Cache or not?  

Hutch github - Amy's workflows all have 'wdl' in their name.  Amy also made a bunch of docker containers that are set up for common tasks - useful for reproducibility or for when I transition to using AWS.

Workflow options - see examples folder. E.g. caching or not.

Vocab:  when we scatter, each job is a 'shard'

In WDL code:   
`~{variable}` is a WDL variable (as opposed to a shell variable)

see 'Output copying' in [Cromwell docs](https://cromwell.readthedocs.io/en/stable/) for how to copy the files we want as final output from /fh/scratch to a specified place in /fh/fast

The workflow block has an `output` section at the end that specifies which files are final output.  Each task also has an `output` section.

## Conversation with Amy Dec 14 2022

### Rerunning tasks
When does cromwell decide to re-run a task?  
- it checks the last modified time and file path of any input files - if those have changed it will re-run.
- it checks the shell script that would be produced to run a task: if this differs at all from the previous run, we rerun it (even if the only difference is number of cores)
- if any upstream tasks changed, ALL downstream tasks will be rerun too.

It does NOT usually look at md5 checksums to really really make sure a file is identical.  There is a way to make it do that, but it makes everything very slow.



### How to move/copy output files to a more sensible place. 

Several options: 
- Add a 'cp' task at the very end of the workflow
- Use the table of output file names I can get from the Shiny app (or via the fh.wdlR package) as input for a shell/perl/R script to copy files locally.
- Use the [workflow options json](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/).   One of the many things I can do with this file is to specify a place to copy the final output files to. The usual behavior is to copy the whole crazy nested structure of the workflow directory - this guarantees we never overwrite files if they have the same names as one another. There is some sort of option to un-nest, but you have to be very careful about file naming.  


### Tips for developing WDL code

Can run a task the first time without any output block. This is especially useful if we're running a program and we can't remember all the names of the output files.

WDL uses error codes returned by each task to determine completion.  Knowing that could help me figure out how to use error codes in my perl (etc) scripts to make the whole pipeline stop if there was an error.

Inputs:  Amy's wdl scripts tend to specify every single input file needed for a task, even if the command-line for the task only takes one of them. Example:  to map reads using `bwa`, we only need to specify the reference genome fasta file name, but bwa expects a bunch of index files to be present too.  If we are only ever running on /fh/fast it is unnecessary to specify all of those files, because cromwell does not need to copy the input files anywhere in order to get the task to run.  However, if we are running using docker and/or AWS, the files DO need to be copied over, so we DO need to specify every single one.

### Scattering tasks

Should be able to do a nested scatter (e.g. scatter over samples, then scatter each sample over regions)

In my case (bat SNP calling) I will want to hard-code the ~100 bed file names into a json file to specify them as inputs. Amy says it's unlikely I can use their sequential numbering to help me.


# using sub-workflows

If I want to include additional WDL files, e.g. if I have split some tasks into sub-workflows, I need to supply the extra WDL files to the Cromwell server when I submit the job. I have two options.

## Option 1: use https

I make the sub.wdl files web-accessible (e.g. on github), and import the raw version, like this:
```
import "https://raw.githubusercontent.com/jayoung/janet-learning-WDL/main/scripts_for_import/dnaseq_fq_to_vcf.skeleton.subworkflow.eachPair.wdl" as sub_eachPair
```
Sometimes I might want to specify a particular commit. For example, I was having trouble seeing the updated version of a wdl soon after I synced to git, because github only updates the raw cache every 5 minutes. See [this post](https://stackoverflow.com/questions/62785962/get-raw-file-from-github-without-waiting-for-5-minute-cache-update). I do it like this:
```
import "https://raw.githubusercontent.com/jayoung/janet-learning-WDL/dfed3d18be8c7150d01f420b60312470e8b3749c/scripts_for_import/dnaseq_fq_to_vcf.skeleton.subworkflow.eachPair.wdl" as sub_eachPair
```

## Option 2 - make a zip bundle

I put the extra WDLs into a zip bundle like this:
```
zip subwdls sub1.wdl sub2.wdl
```
And I can submit ththe zip bundle to Cromwell via `fh.wdlR` using the `Dependencies` option:
```
thisJob <- cromwellSubmitBatch(WDL = "my_workflow.wdl", 
                               Params = "my_inputs.json", 
                               Dependencies = "subwdls.zip")
```


## April 4th, 2024: DaSL PROOF demo session

PROOF (Production Research On-ramp for Optimization and Feasibility) - new tool to run WDLs at the Hutch

i.e. the Cromwell server.  Can interact via:
- web shiny app [proof.fredhutch.org](proof.fredhutch.org)
- an R package called `proofr` (with rcromwell)

New [guide](https://hutchdatascience.org/WDL_Workflows_Guide/index.html) to writing WDLs 

The app:
- first step is always to log in (button on top right)
- for file uploads, can drag-drop (or browse)
- more info [here](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/)

My questions, maybe
- on submitting a workflow we can give it nicknames. where do they appear after the workflow is running?
- would seem natural to have a troubleshoot/cancel link from the 'track jobs' tab, rather than having to copy-paste the ID to the 'troubleshoot' tab

Resources for pre-made WDLs:
- DASL's [WILDS github collection](https://github.com/orgs/getwilds/repositories?) - any repo whose name begins WW is WDL-based
- [BioWDL](https://biowdl.github.io)

WILDS also has a [library of docker images](https://github.com/getwilds/wilds-docker-library) that can be used by WDLs

