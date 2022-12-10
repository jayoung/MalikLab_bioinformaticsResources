# Useful links

Shiny job submission/monitoring [dashboard](https://cromwellapp.fredhutch.org) (current db ID `gizmok122:46707`)

Amy's [original instructions](https://github.com/FredHutch/diy-cromwell-server) from before Oct 2022 on starting a Cromwell server, and example jobs

Amy's [newer instructions](https://hutchdatascience.org/FH_WDL101_Cromwell/introduction.html) from Nov 2022 

My cromwell home dir: `~/FH_fast_storage/cromwell-home`

Cromwell scratch dir: `/fh/scratch/delete90/malik_h/jayoung`

[OpenWDL](https://openwdl.org) has links to tutorials

[Cromwell docs](https://cromwell.readthedocs.io/en/stable/) - Amy says some stuff there is out of date

Learning WDL [readthedocs](https://wdl-docs.readthedocs.io/en/stable/) (recommended by Amy Dec 2022)

# progress

I think I undertand how to run WDLs on our cluster: I worked through the [Hutch WDL101 course](https://hutchdatascience.org/FH_WDL101_Cromwell/index.html) in detail, providing lots of feedback

Next I need to learn WDL.  I got up to [here](Learning WDL [readthedocs](https://wdl-docs.readthedocs.io/en/stable/WDL/execute/) in the readthedocs tutorial


# questions for amy

- how to run workflows
- how to troubleshoot - where are the errors?
- SNP call design
- how does cromwell keep track of when to reuse results versus when to rerun


# habits
copying output files

where to keep the wdl+jsons

to me it makes sense to use the same name for the wdl script file as for the workflow block.  The directory structure in /fh/scratch uses the workflow block name 

# syntax questions
- indents don't matter: is that correct?  is that true for WDL files AND json files?
- what's the syntax for the command block?  in places I see
```
command <<<

>>>
```
and in places I see 
```
command {

}
```