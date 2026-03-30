# Running WDLs at the Hutch 

## Current notes

(Notes from March 2026) 

We now use the [PROOF workbench](https://pwb.fredhutch.org/) to spin up a cromwell server, and to validate/run/monitor WDL scripts 

With default configuration:

```
Scratch directory: /hpc/temp/malik_h/user/jayoung/cromwell-scratch 
Workflow log directory: /home/jayoung/.cromwell/workflow-logs 
```

Hutch documentation: [Getting Started with PROOF Workbench](https://sciwiki.fredhutch.org/datademos/pwb-tutorial/)


A workflow run will typically have (at least) three files:
- the WDL script (e.g. `ww-sra-star.wdl`)
- a json file specifying input files, parameters, etc (e.g. `inputs.json`)
- a json file specifying options for this run of the workflow, e.g. `cromwell-options.json`. This file includes some useful parameters:
    - `final_workflow_outputs_dir`
    - sets workflow default for `maxRetries` (I think it can be overridden for individual tasks)



## Old way to run WDLS at the Hutch

See `~/FH_fast_storage/cromwell-home/janet-learning-WDL`

Pre-2025ish:  the way to run WDLs at the Hutch was by starting your own Cromwell server on the command line, and interacting with it using the [PROOF shiny app](https://proof.fredhutch.org/)

There were also some R package called `fh.wdlR` that allowed interaction with a running cromwell server and I think I was using it to help me copy pipeline outputs. 

Before, I found it frustrating to copy final pipeline outputs to a useful location. Now it seems much more doable, by setting the `final_workflow_outputs_dir` parameter in `cromwell-options.json`

# notes files

xxx to do after reorge
- get rid of my janet-learning-WDL repo
- any other repos to get rid of?
- any links to that repo?

look at xxx and see if there's anything useful to keep 

