# Fred Hutch compute resources

Much more documentation [here](https://sciwiki.fredhutch.org/scicomputing/comp_index/). These are just a few notes to help me get you up to speed.

If you’re working remotely, you’ll need to install and use [VPN software](https://centernet.fredhutch.org/cn/u/center-it/help-desk/vpn.html) to connect to Hutch computing resources.


## Sources of help

There's a Slack workspace called [FH-BCR](https://fhbig.slack.com/), with several helpful channels (e.g. R, python, Q&A). Search to see if someone asked a related question before, or ask your own question - it's a very helpful community.

Scientific computing staff are helpful - you can contact them by email (scicomp@fredhutch.org) and they monitor Slack for questions.


# The compute cluster

The cluster consists of many "nodes" (one node = one computer), each named `gizmoXXX` (XXX = a letter and a number). 

Other computers named `rhino1`, `rhino2` etc control access to FHCRC's compute cluster. ([technical details](https://sciwiki.fredhutch.org/scicomputing/compute_platforms/))

We should not use the rhino computers to do any significant computing - they simply provide access to the gizmo nodes, where the real work happens. Be careful - rhino will allow you to run stuff, but the systems admin people don't like it. It's good to keep scicomp people happy - they are very helpful to us, and they know what they're doing.

There are two ways to work on the cluster:
- "interactively", where we type commands, stuff happens, we type more commands
- "batch" mode (more advanced), where we might want to run longer-running command or set of commands without paying attention to it. We can submit batch jobs to the cluster from a rhino computer using the `sbatch` or `srun` commands (more below).

Each node has multiple "CPUs" (also known as "threads") (central processing unit): some nodes have 4 CPUs, others 6, 14, 24, 36 CPUs. Often we just use one CPU per program we're running, but some programs (e.g. `blastn`) are written so they can use >1 CPU at once ("parallel processing"). We will probably be sharing a node with other people using it at the same time: try to stick to the number of CPUs you requested, and don't exceed that. 

Cluster specifications, last time I looked: (up-to-date version should be [here](https://sciwiki.fredhutch.org/scicomputing/compute_platforms/)):  
There are four classes of gizmo nodes (e.g. gizmok27 is an K class node).  
- G class (18 computers) each = 6 CPUs, 256Gb RAM
- H class (3 computers) each = 14 CPUs, 768Gb RAM
- J class (42 computers) each = 24 CPUs, 384Gb RAM
- K class (170 computers) each = 36 CPUs, 768Gb RAM
 -(D class nodes are obsolete: they had 12 CPUs, 48Gb RAM) 

 

# Getting interactive access to a node

See the [documentation](https://sciwiki.fredhutch.org/scicomputing/access_methods/)

To get interactive access to a node and do some work, there are several steps:

1. log on to rhino:
- start up a Terminal window from your Mac
- figure out whether your unix-style user name on your Mac is the same as your FHCRC user name. You can see your Mac username by typing `whoami` in the terminal window.
- if user names are the same, you can do this: `ssh rhino`.  If not, you'll need to do `ssh myFHCRCuserName@rhino`.   If you're doing this from outside the Hutch, you might need to give the complete address of rhino (rhino.fhcrc.org), and you should be running [webVPN]((https://centernet.fredhutch.org/cn/u/center-it/help-desk/vpn.html)).
- enter your Hutch password. If you've successfully logged on, you'll see a prompt that looks something like this: `rhino02:/home/jayoung`.  This tells us we're on a computer called `rhino02`. 

2. Maybe we'll go through a [simple unix tutorial](unix_intro.md) now, working on rhino for the moment.  

3. When we're ready to do some "real" work that might drain compute resources, we'll ask the rhino computer to give you permission to work on one of the nodes:
- decide how many CPUs you want to ask for, and how many days you'll want access to the node.  1 CPU is enough for very simple stuff. 4 can be nice to run several things simultaneously and all of the available computers have at least 4 CPUs.  Don't run more simultaneous jobs/threads than the number of CPUs you choose!
- from your prompt line on the rhino machine, issue this command:
`grabnode`. You will be asked how many CPUs you want, how much memory and how long you want it for.
- after a few seconds you'll see a bunch of other messages come up, and a prompt that tells you the name of the computer you've been allocated (e.g. the prompt looks like this: "jayoung@gizmok27:~$" when the computer is called gizmok27).
- make sure you keep this window open while you work - when you exit your original `grabnode` window, it ends your permission to work on that gizmo node.  If your computer goes to sleep, it'll also disconnect - there are ways around this, but the simplest one is not to allow the computer to sleep.  I often keep the "htop" command running in this window, to remind me that I should not exit this window.   I open up additional Terminal windows, and use `ssh` from those to actually run commands (see step 4).

4. To log in to a node you’ve ‘grabbed’ from another Terminal window, issue a command that looks something like this:
"ssh gizmoXXX"  (or "ssh jayoung@gizmoXXX" if your Mac and FHCRC user names are different).   You should replace the XXX with whichever computer you've been allocated, e.g. "ssh gizmok27".
The Terminal-Shell-Edit Title menu option can be handy to help keep track of which terminal window is doing what.


Advanced users might want to use the [`screen` or `tmux`](https://sciwiki.fredhutch.org/scicomputing/access_methods/#screen-and-tmux) utilities to manage your login session.

# File storage and access:

You'll probably store smaller files in your own home directory, but any larger files should go in the "fast" storage area: `/fh/fast/malik_h`.  I store pretty much everything in the `fast` disk. 

You probably want to create a folder for yourself in here `/fh/fast/malik_h/user`. Please use your Hutch user name as the folder name (for example, my folder is called `/fh/fast/malik_h/user/jayoung`). For projects/data where >1 lab member is involved, make a folder in `/fh/fast/malik_h/grp`. 

On unix, your home directory will be called something like `/home/jayoung`.  

Large data generated by Hutch Shared Resources (e.g. Illumina sequencing) will appear in `/fh/fast/malik_h/SR`.

Occasionally, we have large files that we want to be visible to people from outside the Hutch. We USED TO use the `/fh/fast/malik_h/pub` for that, via the xfiles server. Example - I have made custom UCSC tracks/hubs and put them in this folder, using URLs that look something like this: http://xfiles.fhcrc.org:7007/malik_h/pub/myHub3/   But that has changed: see [NEW INSTRUCTIONS](https://sciwiki.fredhutch.org/compdemos/ucsc-track-s3/)


## Mounting Hutch storage on a Mac

[This link](https://centernet.fredhutch.org/cn/u/center-it/help-desk/network_drive_paths_mac.html) might help.

Use the Mac's Finder application, and choose `Go` - `Connect to server`. You'll get a menu where you can enter the address of the storage server you want to connect to.

For your home directory, enter this address, replacing jayoung with your own user name:  `smb://FHCRC;jayoung@file.fhcrc.org/home`  

For the `fast` directory, there are two possible different addresses. In Jan 2020, Luna told me that the first way responds quicker. However, links are not visible that way, so for some people/purposes, the second way might be better: scicomp (scientific computing staff) have told me to use the second way.
1.	`smb://center.fhcrc.org/fh/fast/malik_h/`
2.	`smb://FHCRC;jayoung@file.fhcrc.org/fast/malik_h` ()

If you're working on the command line on your Mac (i.e. not working on the cluster, but interacting with your Mac using the Terminal application), the fast directory can be found in `/Volumes/fh/fast/malik_h/` or sometimes `/Volumes/malik_h/`



## Mounting Hutch storage on a PC:

[This link](https://centernet.fredhutch.org/cn/u/center-it/help-desk/network_drive_paths_windows.html) might help.


# Storage on disk (large files)

Our lab gets 5Tb of file storage funded by FHCRC, and above that we'll be charged monthly for disk storage. All of us should keep an eye on any large files we download and generate - it's fine to pay for them if we're actually using them, but we should tidy up files we're not using as much as we can.  (there is another cheaper place we can store "archive" type files that we don't really use, but that system is not so robust yet and requires special access methods)

A command to check file/folder size:  
`du -sm myFolder` = shows total size of everything in myFolder  
`du -sm *`        = shows total size of each item in the current directory  
   (`du` = disk usage, `-s` = summary, `-m` = megabytes)

Some programs or bioinformatic pipelines generate a lot of large intermediate files that can be deleted after the program has run (e.g. de novo assembly programs).



# Using scicomp's installed "environmental modules" (programs)

The scicomp people have installed MANY useful biology tools for us to use on the cluster. See instructions [here](https://sciwiki.fredhutch.org/scicomputing/compute_environments/).

To list all available modules (installed in /app and /app/easybuild): `module avail`

To look for modules whose name contains ‘bowtie’: `module avail bowtie` (it will find partial matches and it’s not case-sensitive).

If there's a program you're interested in trying that doesn't seem to be available, ask me if I have installed it for our lab - I have a few additional tools that scicomp haven't installed as modules. You can also ask them to install things - sometimes they respond very quickly, sometimes it takes a few days or a week.  If you're feeling brave you can install it yourself, for your own use.

To load a module (i.e. make it ready to use). This only applies in your current terminal window, and will reset when you exit that window:  
    `module load Cufflinks`  
That loads the default version of cufflinks, usually the most recent version that's been installed. In the output of `module avail`, the D indicates the default version of each tool, which can change with time. To make your work more reproducible, you probably want to get in the habit of choosing a specific version of each tool for each project, and sticking to it:  
`module load Cufflinks/2.2.1-foss-2018b`

To see all modules currently loaded in your session: `module list`

To unload all modules currently loaded in your session: `module purge`

# Running commands on the cluster in batch mode

More info [here](https://sciwiki.fredhutch.org/scicomputing/compute_jobs/), but basically you do this:

1. figure out the command or set of commands you want to run.  Multiple commands can be run one after the other by concatenating them with ";".  For example, to run `command1` then `command2`, you would do this:  
`command1 ; command2`

2. in the simplest form, you can run the job on a cluster node by wrapping the command by issuing a command like this from a rhino/gizmo command line:  
`sbatch --wrap="command1 ; command2"`


3. but we might also want to specify some options for the sbatch job, for example:
```
sbatch --job-name=testJob -t 0-2 --cpus-per-task=4 --partition=largenode --wrap="command1 ; command2"
```
- job-name:  helps track the job using squeue  
- t:  "walltime" = how much time to allot to the job. 0-2 means 0 days, 2 hours.  Default is 3 days.  
- cpus-per-task: number of CPUs  
- partition: FHCRC has a default queue called `campus` - those nodes have 4 CPUs and 32 Gb memory. The `largenode` queue gets you on a machine with 28 CPUs and 768 Gb memory.

4. some useful commands to see if your job is still going, and to monitor cluster usage are:  
`squeue -u jayoung`  
`squeue`  
`hitparade`  
The jobID (a number, shown when you submit the job and in squeue output) can be useful to diagnose problems if your job fails to start:
`scontrol show job myJobID`  
If the commands you ran resulted in any screen output, or gave errors, that output is sometimes (but maybe not always) captured in a file called `slurm-myJobID.out`

5. sometimes we want to wrap a command that includes using a `module load` function - this is how we do it (it looks  conplex, but it works)
`sbatch  --wrap="/bin/bash -c \"source /app/Lmod/lmod/lmod/init/bash; module load BWA; bwa\""`

6. an alternative, tidier way to run a series of commands within a single sbatch job is to create a "shell script". It’s a plain text file that contains all the commands you want to run – in this very simple example the file is called `myCommands.sh` and here are the contents:
```
#!/bin/bash
echo ‘hello 1’
echo ‘hello 2’
```

If we want to load modules within a shell script, we need to include a ‘source’ line exactly like the second line in this more complicated shell script:
```
#!/bin/bash
source /app/lmod/lmod/init/profile
module load Bowtie2
bowtie2 --help
```
Here’s how you actually run the shell script:  
`sbatch myCommands.sh`  
or, with additional options:  
`sbatch --cpus-per-task=8 -t 1-0 --job-name=myJob1 myCommands.sh`


7. if there's a problem, use scancel to kill a job or jobs. The simplest way is to use the job ID: `scancel myJobID`
You can also do this to cancel all your running jobs - e.g. to cancel all jobs for user `jayoung` you do this:  `scancel -u jayoung `

# Running a bunch of sbatch jobs, on multiple samples

If you want to get fancy, and run the same command on a set of inputs, you might want to set up a **"slurm script"** (a.k.a. **"sbatch script"**) that will look something like this:

```
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
```

A version of that script with detailed comment lines to help you understand it is [here](../example_scripts/simple_slurm_script.sh). 

To run an slurm script: `sbatch myScript.sbatch`

Some key points:
- the lines beginning `#SBATCH` specify resources to request for each sbatch job.  
- the other words in capital letters are *variables* that we define using the = symbol.  Once we've defined variables, we can use their values with `${}` notation.
- notice that `SAMPLE_IDS` has paretheses and specifies several items - it's an array, and we can run some commands on each item in the array, using the `SINGLE_ID` variable.

# Shared databases and installed programs

## databases:

There's a big folder where I keep databases that might be useful for several of us: `/fh/fast/malik_h/grp/public_databases`. The file `database_list.txt` contains a list of many (but not all) of the databases in that folder. 

## installed programs:

It's a bit of a mess, but there are a bunch of programs installed in these locations:
```
/fh/fast/malik_h/grp/malik_lab_shared
/fh/fast/malik_h/grp/malik_lab_shared/bin
/fh/fast/malik_h/grp/malik_lab_shared/linux_gizmo/bin
```

If a program comes with documentation files, I put those in this folder: `/fh/fast/malik_h/grp/malik_lab_shared/help`

Some tools are in subdirectories, e.g. NCBI's edirect tools are here:
`/fh/fast/malik_h/grp/jayoung/malik_lab_shared/linux_gizmo/bin/edirect/`

# Additional notes for PC users

`MobaXterm` (or `Xming`) are a couple of options to allow you to run X-windows programs.

Be careful of line break format - it's different between PC and unix. Can mess things up.
