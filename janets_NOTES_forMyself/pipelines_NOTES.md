# Makefiles

[Makefiles tutorial](http://byronjsmith.com/make-bml/) 

Makefiles are an old-school way to make a pipeline.  Originally intended for building software, but they work for pipelines too: they are able to only run parts of the pipeline for files that have been modified

This is an example of a "rule" in a Makefile (the whole thing is a rule):
```
isles.png: isles.dat
	./plotcount.py isles.dat isles.png
```
First line:   `target: input`  

Remaining line(s): commands used in that step of the pipeline, after a TAB character (NOT spaces - make is fussy)

The "target" is `isles.png` and is the name of the output file created by the command(s) specified in the rule.  Then there's a colon, and the names of any input file(s) needed by the rule

The command `make` looks for a file called `./Makefile` and by default runs the first rule in the file, which by convention often describes a target called 'all'.  You can specify specific rules to be run by supplying the target name e.g. `make isles.dat` or `make clean`.

Makefiles that have several rules often can be thought of as a dependency graph.

`make` knows how to interpret the Makefile dependency graph to see that some of the jobs can be done in parallel! Run using `make --jobs` to allow it to paralellize.

Scripts should be specified as prerequesites in a Makefile as well as data and outputs, so that things get re-run when a script is rerun.


# Nextflow

See `~/FH_fast_storage/general_notes/computing/nextflow_test/nextflow_NOTES_2022.md` and other files there

See also `~/FH_fast_storage/getSameSpeciesGenomicEquivalent` where I'm trying to use nextflow in a real project

# WDL

See `~/FH_fast_storage/cromwell-home/janet-learning-WDL`