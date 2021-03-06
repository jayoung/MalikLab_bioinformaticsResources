# BLAST via the command line.
We can blast any NCBI database, contacting NCBI remotely.  We can also format a set of sequences you have and blast those. 

On the Hutch servers, you probably want to use the module that’s been installed for us:

```
module load BLAST+/2.10.1-gompi-2020a
```
(and after you’re done with blast searches for a while, you might want to do `module purge`)

Show command-line help: `blastn -help`

To blast a remote NCBI database (notice the `-remote` option):
```
blastn -remote -query myNuclQuerySeq.fa -db nr -out myNuclQuerySeq.fasta.blastnNR -task blastn
```

To format your own fasta format file as a blastable database:
```
makeblastdb -help
makeblastdb -in myNuclSeqs.fa -dbtype nucl -parse_seqids
makeblastdb -in myProtSeqs.fa -dbtype prot -parse_seqids
```

To blast your own database (after formatting it)
```
blastp -db myProtSeqs.fa -query myProtQuerySeq.fa -out myProtQuerySeq.fa.blastpMyProtSeqs

tblastn -query myProtQuerySeq.fa -db myNuclSeqs.fa -out myProtQuerySeq.fa.tblastn

tblastn -query myProtQuerySeq.fa -db myNuclSeqs.fa -out myProtQuerySeq.fa.tblastnPrimate -entrez_query 'human[Organism]'
```

Blast can give output in other formats - explore the `-outfmt` option. e.g. `-outfmt 6` gives a tabular output that you could sort and filter, perhaps in Excel:
```
tblastn -query myProtQuerySeq.fa -db myNuclSeqs.fa -out myProtQuerySeq.fa.tblastn_tableOutput -outfmt 6
```


To extract a sequence (or partial sequence) from a formatted database
We'll use a program called blastdbcmd. (in the BLAST+ module)
Blastdbcmd only works you used the -parse_seqids option when you formatted the database.

Look at options:
```
blastdbcmd -help
```

Get seq1 from the myNuclSeqs.fa file:  
```
blastdbcmd -db myNuclSeqs.fa -entry seq1 -out seq1.fa
```

Get the reverse-complement of bases 101-200 of seq1 from the myNuclSeqs.fa file:  
```
blastdbcmd -db myNuclSeqs.fa -entry seq1 -out seq1.fa -range 101-200 -strand minus -out seq1-200-101.fa
```

# Shared databases
I've got a bunch of databases in a shared Malik lab folder, and most are formatted for blasting. See this file for a list `/fh/fast/malik_h/grp/malik_lab/public_databases/database_list.txt`.  We can easily download more databases if they'd be useful - talk to me about this.

Some useful databases at NCBI (using `-remote`) include:
```
nr
env_nt
htgs
chromosome
wgs
refseq_genomic
```


# Remote homology searches

## HMMer package
Hutch module: `module load HMMER/3.3.2-gompi-2020b`

hmmbuild creates an HMM from a multiple sequence alignment
```
hmmbuild -h
hmmbuild myAln.fa.hmm myAln.fa
```

hmmsearch takes an HMM and searches protein sequences:
```
hmmsearch -h
hmmsearch -A myAln.vs.seqsToSearch.stockholm myAln.fa.hmm seqsToSearch.fa > myAln.vs.seqsToSearch.stockholm.hmmsearch
```

## Other options

Phyre, HHpred, and I'm sure many more

# Multiple sequence alignments

There are MANY programs to do this. Each has MANY parameters you can set and experiment with.

My current favorites are:
- `mafft` (for DNA or protein alignments, very fast) 
- `MACSE` (for in-frame, translation-based alignments, can handle frameshifts). 

In the past I have also used:
- `PRANK` (in-frame nucleotide alignment
- `clustalw` (DNA/protein)
- `translatorx.pl` (in-frame nucleotide alignment, using Muscle, Clustalw, Prank, or mafft to actually perform the alignment), e.g. `translatorx.pl -i myDNAseqs.fa -o myDNAseqs.translatorx`


## mafft

There are many ways to run MAFFT. Here's how you run it using the defaults:
```
module load MAFFT/7.453-GCC-8.3.0-with-extensions
mafft -h
mafft seqsToAlign.fa > alignment.fa
module purge
```

## MACSE

MACSE requires java, so do this before trying to run it:
```
module load Java/11.0.2 
```
and perhaps this after you've finished:
```
module purge
```

Show all MACSE modes:
```
java -jar /fh/fast/malik_h/grp/malik_lab_shared/bin/macse_v2.06.jar -help
```

Show help for one particular mode:
```
java -jar /fh/fast/malik_h/grp/malik_lab_shared/bin/macse_v2.06.jar -prog alignSequences -help
```

To make an alignment from scratch, use `alignSequences` mode:
```
java -jar /fh/fast/malik_h/grp/malik_lab_shared/bin/macse_v2.06.jar -prog alignSequences -out_NT alignment.fa -seq seqsToAlign.fa
```

To add seqs to an alignment, use `enrichAlignment` mode:
```
module purge
module load Java/11.0.2 
java -jar /fh/fast/malik_h/grp/malik_lab_shared/bin/macse_v2.06.jar -prog enrichAlignment -align existingAlignment.fa  -out_NT newAlignment.fa -seq seqsToAdd.fa
module purge
```

## trimming/masking alignments

`trimal` can trim gappy positions out of an alignment (nucleotide or protein):
```
trimal -in myAln.fa -out myAln.trimal.fa -automated1 -htmlout myAln.trimal.html -colnumbering > myAln.trimal.retainedColumns.txt
```

Also: `SWAMP`

# Phylogenetic trees

## PhyML

Phylogenetic analyis by maximum likelihood. See [manual](https://github.com/stephaneguindon/phyml/blob/master/doc/phyml-manual.pdf).

Input file is a multiple sequence alignment in phylip format.

The simplest way to run phyml is single-threaded version (`phyml`). Entering `phyml` on the command line brings up a menu-driven interface to run phyml. It totally works, but I don't like doing it this way - I prefer to figure out the full command-line version to run it, so that I have a record of how I made each trees, letting me reproduce it in future if necessary, or run the same method on another alignment.

Show help for command-line phyml:  `phyml --help`  
Example: amino acid tree, 100 bootstraps, JTT evolutionary model:
```
phyml -i myaln.phyml -d aa --sequential -m JTT --pinv e --alpha e -f e -b 100
```

But this can be slow, especially if we're running a lot of bootstraps, so we can also run the multi-threaded version (`phyml-mpi`). To see all the options:
```
phyml-mpi --help 
```

Example - to run 1000 bootstraps, using the JTT+I+G+F model, in an interactive session (you would do this in a grabnode session where we asked for 4 CPUs):
```
module load OpenMPI 
mpirun -np 4 phyml-mpi -i bigAln10.phyml -d aa --sequential -m JTT --pinv e --alpha e -f e -b 1000
module purge
```
`-np` = num CPUs (1-12, must have requested that number of CPUs when you logged in to gizmo node)  
`-b` = num bootstraps (must be a whole multiple of num CPUs, so I might do 4 CPUs and 100 or 1000 bootstraps)  
`-m` = evolutionary model  
`-f e` = estimate frequencies by counting alignment  
`--pinv e` = estimate proportion of invariable sites  
`--alpha e` = estimate the gamma distribution shape parameter

We would probably want to run that as an sbatch job instead, requesting a node with 10 CPUs:
```
sbatch --job-name=phyml --cpus-per-task=10 --wrap="/bin/bash -c \"source /app/lmod/lmod/init/profile; module load OpenMPI; mpirun -n 10 --oversubscribe phyml-mpi -i bigAln10.phyml -d aa --sequential -m JTT --pinv e --alpha e -f e -b 1000 > bigAln10.phyml_screenOutput.txt\""
```

## Choose best evolutionary model for nucleotide alignments
Use `jmodeltest`. The help file is here:
`/fh/fast/malik_h/grp/malik_lab_shared/jmodeltest-2.1.10/README`

examples are here:
`/fh/fast/malik_h/grp/malik_lab_shared/jmodeltest-2.1.10/example-data`

Input format can be fasta or phylip

First, change directory, because we have to run prottest from a particular directory (pushd is similar to cd, but remembers which directory we were in before, and lets us use "popd" to return to that directory later):
```
pushd /fh/fast/malik_h/grp/malik_lab_shared/jmodeltest-2.1.7
```

We also load the java module
```
module load Java/1.8.0_181
```

get more help, including details on all options (need to do this from the jmodeltest folder):
```
java -jar jModelTest.jar -help
```

Then run jmodeltest. Here's an example on a real dataset:
```
java -jar jModelTest.jar -d /home/jayoung/temp/testJmodel/masterAln22forFig.fa -o /home/jayoung/temp/testJmodel/masterAln22forFig.fa.jmodeltest.out -f -i -g 4 -s 11 -AIC -a -tr 2
```
I specify the full path to my input file and output file (usually put the output file in the same folder as the input file).  
       option  "-tr 2" = use 2 threads.     
                 BE CAREFUL - default setting is 12 threads.
Number of threads should be between 1 and 12, but you must have requested that number of CPUs from the gizmo node you logged in to. (it may tell you on the screen that it is "proceeding without MPI" but it will still use multiple threads).  If you're testing lots of models, it can be a little slow, so multiple threads are useful.

Then return to the original directory (wherever we were when we did pushd)
popd

take a look at the output to see what the best model was: 
more outputFile.txt
e.g.
Model selected: 
   Model = SYM

## choose best evolutionary model for amino acid alignments
Use `prottest` The help file is here: `/fh/fast/malik_h/grp/malik_lab_shared/prottest-3.4.2/README`

examples are here: `/fh/fast/malik_h/grp/malik_lab_shared/prottest-3.4.2/examples`

Input format can be fasta or phylip.  Should not contain * for stop codons. Long sequence names will be truncated, but that's OK.

First, change directory, because we have to run prottest from a particular directory (pushd is similar to cd, but remembers which directory we were in before, and lets us use "popd" to return to that directory later):
```
pushd /fh/fast/malik_h/grp/malik_lab_shared/prottest-3.4-20140123
```

We also load the java module
```
module load Java/1.8.0_181
```

To get more help, including details on all options (need to do this from the prottest folder):
```
java -jar prottest-3.4.jar  -help
```
Then run prottest. The general form of the command I use is as follows. 
```
java -jar prottest-3.4.jar -i alignmentFile.phyml -all-distributions -F -AIC -BIC -tc 0.5 -threads xxxNumberOfThreads > outputFile.txt
```
I specify the full path to my input file and output file (usually put the output file in the same folder as the input file).  Number of threads should be between 1 and 12, but you must have requested that number of CPUs from the gizmo node you logged in to.

Here's an example on a real dataset:
```
java -jar prottest-3.4.jar -i /home/jayoung/PARPs/alignments/bigAln10.phyml -all-distributions -F -AIC -BIC -tc 0.5 -threads 12 > /home/jayoung/PARPs/alignments/bigAln10.phyml.prottest
```
Then return to the original directory (wherever we were when we did pushd)
```
popd
```

take a look at the output to see what the best model was: 
```
more outputFile.txt
grep 'Best' outputFile.txt
```
e.g.
Best model according to AIC: JTT+I+G+F



## File format conversion (alignments, trees)
Many other tools can convert formats, but I have some scripts to do it on the command line. 

### Convert fasta format to phyml format 

(phyml format is very similar to phylip format, I think just without any * or X or ? characters):  `fasta2phyml.pl` - you can copy this script from `/fh/fast/malik_h/user/jayoung/bin`.
Simply call my script on alignment files in fasta format.
```
fasta2phyml.pl alignment1.fasta
```
The script has a switch at the top - can choose to leave sequence names as they are, or can replace names with sequential names (seq0001, seq0002, etc). Name conversions will be recorded in a file called alignment1.fasta.alias

### Restore original sequence names in a phylip format tree file
```
changenamesinphyliptreefileguessaliasfilename.pl myTree.tre
```

### Convert fasta format to nexus format 
```
fasta2nexus.bioperl alignment1.fasta alignment2.fasta 
```
creates .nex file (the alignment) and .alias file (records the name conversions)

This script replaces original sequence names with sequential names (seq0001, seq0002, etc)

### Restore original sequence names in a nexus format tree file
```
changenamesinnexustreefileguessaliasfilename.pl myTree.tre
```

## mrbayes, to make bayesian phylogenies

It's worth looking at the help files, especially the manual, in this folder:
~/malik_lab_shared/help/mrbayes
e.g ~/malik_lab_shared/help/mrbayes/documentation/Manual_MrBayes_v3.2.0_draft.pdf

input sequence alignment files should be in nexus format

To start up an interactive version of mrbayes (single CPU version):
`mb`

We'll type various commands within this program. First set up parameters, and then run the MCMC chains.

To get detailed help within interactive mrbayes:
help             (general help, lists all available commands)
help sumt   (help on sumt command)
help lset      (help on lset command)
etc

Make a tree using mrbayes - interactive mode
Here's an example of running the commands we used to run a full analysis within interactive mrbayes. These are just some example commands - should think about what parameters to use for your data.

1. start up mrbayes
mb

2. read the data from a nexus format file: 
execute myAlignment.nex    

3.  set some parameters:
lset = set the parameters of the likelihood model
prset = set the priors for the phylogenetic model
lset nst=6    (for DNA alignments, sets evolutionary model to GTR)
lset rates=invgamma   (allows a proportion of invariable sites)
prset aamodelpr=mixed (for amino acid alignments, allows the Markov chain to sample different amino acid models and combine results)

4. look at current settings, check they're OK before we run the analysis:
showparams 
showmodel

5.   Actually run the analysis.:
mcmc Ngen=5000000 Stoprule=yes Stopval=0.01
With these parameters we allow it to run for up to 5million generations (a lot! it might take a while). Stoprule and Stopval mean it will stop if it reaches convergence in fewer generations (convergence threshold is a "split" value of 0.01)

6.    summarize parameters of the run - look at these and make sure they look OK:
 sump 

7.   summarize and output the trees, showing all compatible paritions (default is 50% majority rules - I'd rather see all of them, even if they're more poorly supported): 
sumt relburnin=yes burninfrac=0.25 contype=allcompat

8.    quit

There will be many output files. The consensus tree is in this file:
myAlignment.nex.con.tre

Make a tree using mrbayes - putting all commands in a file 
Make a text file containing all the commands you want to run, and read it into mrbayes like this, capturing any screen output (errors, or good output) to a log file:

mb < myMBrun.txt > myMBrun.log.txt
you can monitor how it's getting on using this command: 
tail myMBrun.log.txt

Make a tree using mrbayes on the cluster using sbatch (commands in a file)
sbatch --job-name=MrBayes --wrap="mb < myMBrun.txt > myMBrun.log.txt"

you can monitor how it's getting on using this command: 
tail myMBrun.log.txt

Use the parallel version of mrbayes, putting all commands in a file
Be warned - an older version of mrbayes had bugs in the parallel version, and bugs were not immediately obvious. http://sourceforge.net/p/mrbayes/bugs/1599/
The program ran and gave an output tree, but instead of giving anything like the expected tree, we got a mostly unresolved tree (very very low bootstraps at all nodes).   This bug is said to be fixed in v3.2.3 but I haven't tested it.

substitute mb with mb_multiCPU, and preface the command with mpirun -np 4 (for 4 CPUs).   Number of CPUs must be evenly divisible into the number of chains (a parameter you can set: the default is 8, so 4 CPUs makes sense) 

mpirun -np 4 mb_multiCPU < myMBrun.txt > myMBrun.log.txt 

you can monitor how it's getting on using this command: 
tail myMBrun.log.txt

Similarly, can use sbatch/wrap to run this command on the cluster (make sure you specify --cpus-per-task)

# dotter, to compare two sequences visually

(need to install X11 on your Mac for this to work, if it's not already installed.  Also, in X11 preferences, under Input, make sure "Emulate three button mouse" is on)

For dotter, human-versus-mouse is quite a distant comparison (mouse evolves fast).   You can usually see conservation of exons, but not much in introns or intergenic regions.   Try to choose more closely related species pairs, or use slowly evolving species as reference (I've had OK luck with dog, or elephant).

Look at options: `dotter --help`.  
Notice in the help output that many options have a short form with one dash, e.g. `-b` and long form with two dashes e.g. `--batch-save`.  This form is common.  

Also, notice that after each option, the help page tells you whether it expects:
- a file name (<file>), 
- a whole number (<int>), 
- a number that might not be a whole number (<float>), 
- or nothing (e.g. the -w option - just tells the program to do something or not)

There's also a manual file here: `~/malik_lab/unix_program_documentation/dotter/Dotter_manual.pdf`

Dotter can take time to run for larger files, so we'll first run it in "batch" mode (-b), where it saves output to a file but does not display it:
`dotter -b seq1vsSeq2.dot seq1.fa seq2.fa`

Then we'll display the output ("-l" means "load previously saved dotter output":
`dotter -l seq1vsSeq2.dot seq1.fa seq2.fa &`
Another window should pop up. The "&" at the end of the command simply lets the command run "in the background", giving us the command like back to do other things.

Features files help us understand dotter output (they annotate where known genes are). They're in gff3 format (see the pdf manual). That's a bit annoying to generate by hand, but here's one way to get a gff3 file for gene annotations available through the Ensembl browser.  As an example, let's get annotations for the region of the human genome containing Mx1 and Mx2:
1. Show your region of interest of the human genome (GRCh37 version) in the Ensembl browser: 21:42730351-42834838
2. Choose "Export Data" 
3. Configure the export (this is how I did it, and it seems to work) 
output= Generic Feature Format Version 3
select region:  make sure this is the region you're interested in (strand = 1, probably) and set some options:  
gene information=yes  
transcripts=no  
exons=yes  
introns=no  
coding sequences=yes
4. Click next, choose Text, and save the output to a file. Put that file on fred in a useful folder.

Make sure the annotated sequence you're comparing has the same name as the sequence in the gff3 file (in this case, it should be called "21").

Also, editing that gff3 file to find-replace "gene" with "mRNA" seems to work a bit better.

Now run dotter loading that features file, using the -q option to tell it that the first sequence we looked at starts at position 42730351 relative to the annotations (those give positions on the whole chromosome):
dotter -f seq1geneAnnotations.gff3 -q 42730351 -l seq1vsSeq2.dot seq1.fa seq2.fa &

We'll look at how to use the Greyramp tool to alter the comparison stringency, and how to zoom in and out to reveal more detail.

Another, even more annoying way to generate a features file:
- use UCSC Table Browser to show genes in the region of interest and download as gtf file
- convertUCSCgtfToDotterFeature.pl humanIFITlocusRefGeneTrack.gtf humanIFITlocus.fa


