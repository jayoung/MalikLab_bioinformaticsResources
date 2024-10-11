# Bioinformatics tools (core)

This document contains notes on tools I'm more familiar with. For notes on other tools I haven't used much, see bioinformatics_tools_not_curated.md

# BLAST via the command line.
We can blast any NCBI database, contacting NCBI remotely.  We can also format a set of sequences you have and blast those. 

On the Hutch servers, you probably want to use the module that’s been installed for us:

```
module load BLAST+/2.10.1-gompi-2020a
```
(and after you’re done with blast searches for a while, you might want to do `module purge`)

## Show command-line help: 

`blastn -help`

## To blast a remote NCBI database 
(notice the `-remote` option):
```
blastn -remote -query myNuclQuerySeq.fa -db nr -out myNuclQuerySeq.fa.blastnNR -task blastn
```

## To format your own fasta format file as a blastable database:
```
makeblastdb -help
makeblastdb -in myNuclSeqs.fa -dbtype nucl -parse_seqids
makeblastdb -in myProtSeqs.fa -dbtype prot -parse_seqids
```

## To blast your own database (after formatting it)
```
blastp -db myProtSeqs.fa -query myProtQuerySeq.fa -out myProtQuerySeq.fa.blastpMyProtSeqs

tblastn -query myProtQuerySeq.fa -db myNuclSeqs.fa -out myProtQuerySeq.fa.tblastn
```

## To blast only some species (remote)

Use an entrez-style query via the `-entrez_query` option: e.g. you can use common or latin name of a species or large taxonomic grouping. 

You can use [NCBI's taxonomy database](https://ncbi.nlm.nih.gov/taxonomy) to help figure out the taxonomic terms you're looking for.  

You can also test whether your entrez query is a good one by pasting the entire query into the [nucleotide](https://ncbi.nlm.nih.gov/nucleotide) or [protein](https://ncbi.nlm.nih.gov/protein) database search bars and seeing whether there are a reasonable number of matches, and whether the matches seem like what you wanted.
```
# human using common name
tblastn -remote -query myProtQuerySeq.fa -db nr -out myProtQuerySeq.fa.tblastnNRhuman -entrez_query 'human[Organism]'

# human using latin name
tblastn -remote -query myProtQuerySeq.fa -db nr -out myProtQuerySeq.fa.tblastnNRhuman -entrez_query 'homo sapiens[Organism]'

# rodents
tblastn -remote -query myProtQuerySeq.fa -db nr -out myProtQuerySeq.fa.tblastnNRhuman -entrez_query 'rodentia[Organism]'

# non-human simian primates:
tblastn -remote -query myProtQuerySeq.fa -db nr -out myProtQuerySeq.fa.tblastnNRhuman -entrez_query 'Simiiformes[Organism] NOT homo sapiens[Organism]'
```

## Blast output formats
Blast can give output in other formats - explore the `-outfmt` option. e.g. `-outfmt 6` gives a tabular output that you could sort and filter, perhaps in Excel:
```
tblastn -query myProtQuerySeq.fa -db myNuclSeqs.fa -out myProtQuerySeq.fa.tblastn_tableOutput -outfmt 6
```
If you already have a blast result you obtained using `-remote`, but you want to see that result in a different format, you first want to find the "result ID" you'll find near the top of the blast output, after `RID:`. For example (using grep to find that ID): 
```
grep 'RID: ' myProtQuerySeq.fa.tblastnPrimate
    RID: 1H0R3J40016
```
Then you use that result ID in the blast_formatter command:
```
blast_formatter -rid 1H0R3J40016 -outfmt 6 -out myProtQuerySeq.fa.tblastnPrimate.outfmt6
```
I think results are only retained on the NCBI server for a week or so: if your blast is older than that, you may need to redo it.

## Extracting the sequences of blast matches

### Extracting one sequence
To extract a sequence (or partial sequence) from a formatted database, we'll use a program called blastdbcmd. (also part of the BLAST+ suite of tools)

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

### Extracting multiple sequences

You probably some sort of tab-delimited table of sequence IDs, start/end coordinates, strand, etc. There are a few different ways to get the seqs you want

#### Extracting multiple sequences from a local database

If your database is local, there are many tools you can use to extract those sequences from the database. One of those is `bedtools getfasta` (use `module load BEDTools`) - see instructions [here](https://www.biostars.org/p/56/) - first you make sure your tab-delimited file is in bed format, then you use `bedtools getfasta`.  Note that [bed format](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) is a bit weird - for start coordinates, the first base of the chromosome is 0, so you need to subtract 1 from the start positions you actually want. The end position is coded in the usual way (so you do NOT subtract 1).

#### Extracting multiple sequences from a remote database
If you want sequences from a remote database, it's a little more complex. One solution is given [here](https://www.biostars.org/p/301274/) and looks like this.

First, make a tab-delimited text file of the sequence regions you want to get. Mine looks like this, where columns are:
  - 1. accession 
  - 2. start position
  - 3. end position 
  - 4. strand (1=forward, 2=reverse)  
  
(start and end coordinates use the 'normal' way of counting, not the weird 0-based start that bed files use)
```
NM_002354       10      30      1
NM_004360       10      30      1
NM_001144663    10      30      1
NM_004063       10      30      2
NM_001145024    10      30      2
```

Then you run a command that looks like this, except that you replace `smallAccList_withCoords.txt` and `smallAccList_withCoords.txt.fa` with the names of your input and output files, respectively. It's a long command, so make sure you copy the entire thing before editing:
```
cat seqsToGet.txt | while read -a F ; do wget -q -O - "https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${F[0]}&seq_start=${F[1]}&seq_stop=${F[2]}&strand=${F[3]}&rettype=fasta" ; done > seqsToGet.txt.fa
```


# Shared local databases

I've got a bunch of databases in a shared Malik lab folder, and most are formatted for blasting. See this file for a list `/fh/fast/malik_h/grp/malik_lab/public_databases/database_list.txt`.  It's not a very well-organized file: sorry! Sometimes I list files under the source where I downloaded them from (e.g. UCSC, NCBI, Ensembl), sometimes under their species or lineage.

We can easily download more databases if they'd be useful - talk to me about this.

Some useful databases at NCBI (using `-remote`) include:
```
nr
env_nt
htgs
chromosome
wgs
refseq_genomic
```


# Searching for remote homologs

## HMMer package

### Running it online:
You can run it online via [EMBL-EBI](https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch) against various databases, and you can restrict your search to certain taxonomic databases.

### Running it on the Hutch server, against databases stored locally:

Hutch module: `module load HMMER/3.3.2-gompi-2020b`

`hmmbuild` creates an HMM from a multiple sequence alignment

show hmmbuild usage:
```
hmmbuild -h
```
run hmmbuild:
```
hmmbuild myAln.fa.hmm myAln.fa
```

`hmmsearch` takes an HMM and searches protein sequences:

show hmmsearch usage:
```
hmmsearch -h
```
run hmmsearch (searches a protein database, in fasta format)
```
hmmsearch -A myAln.vs.seqsToSearch.stockholm myAln.fa.hmm seqsToSearch.fa > myAln.vs.seqsToSearch.stockholm.hmmsearch
```

`jackhmmer` does iterative searches, like psi-blast.  Start with a single sequence query, search a database, build an HMM from the results, search the database again, repeat. Might be good for remote homologs but I haven't tried it.

Protein databases to search: this works on any multiple sequence fasta file. Sometimes I search a file I made myself, sometimes I download all predicted proteins for a given species, sometimes we search bigger databases (uniprot etc).

Check out database files we have stored on the server, listed in `/fh/fast/malik_h/grp/public_databases/database_list.txt`. It's not a very well-organized file: sorry! Sometimes I list files under the source where I downloaded them from (e.g. UCSC, NCBI, Ensembl), sometimes under their species or lineage.

You can look through the whole file, or you might use these search terms in a 'find':
- `uniref`
- `Proteome`
- `protein.faa`
- `.pep.`

## Other options

Phyre, HHpred, and I'm sure many more

# Repetitive sequences

## RepeatMasker

Searches sequences against database of known repeats, detects simple repeats. Gives coordinates of repeats, as well as providing a 'masked' version of the sequence where repeats are replaced by NNNs.

There is an [online server](http://www.repeatmasker.org/cgi-bin/WEBRepeatMasker).

You can also run it on the Hutch server:

By default, RepeatMasker assumes your sequences are from human. It's simple to run:
```
RepeatMasker myseqs.fasta
```

Use the `-species` option if your sequences are not from human
```
RepeatMasker -species mus myseqs.fasta
```

## Tandem repeat finder

There is an [online server](https://tandem.bu.edu/trf/trf.html)

Or you can run it on the server. Detailed instructions are [here](https://tandem.bu.edu/trf/trf.unix.help.html). To run using their recommended parameters:
```
module load TRF/4.09-linux64
trf yoursequence.fa 2 7 7 80 10 50 500 -f -d -m
module purge
```
You'll get four output files, whose names begin `yoursequence.fa.2.7.7.80.10.50.500`




# File format conversion (alignments, trees)

Many other tools can convert formats, but I have some scripts to do it on the command line. 

## Convert fasta format to phyml format 

(phyml format is very similar to phylip format, I think just without any * or X or ? characters):  `fasta2phyml.pl` - you can copy this script from `/fh/fast/malik_h/user/jayoung/bin`.
Simply call my script on alignment files in fasta format.
```
fasta2phyml.pl alignment1.fasta
```
The script has a switch at the top - can choose to leave sequence names as they are, or can replace names with sequential names (seq0001, seq0002, etc). Name conversions will be recorded in a file called alignment1.fasta.alias

## Restore original sequence names in a phylip format tree file
```
changenamesinphyliptreefileguessaliasfilename.pl myTree.tre
```

## Convert fasta format to nexus format 
```
fasta2nexus.bioperl alignment1.fasta alignment2.fasta 
```
creates .nex file (the alignment) and .alias file (records the name conversions)

This script replaces original sequence names with sequential names (seq0001, seq0002, etc)

## Restore original sequence names in a nexus format tree file
```
changenamesinnexustreefileguessaliasfilename.pl myTree.tre
```


# Multiple sequence alignments

There are MANY programs to do this. Each has MANY parameters you can set and experiment with.

My current favorites are:
- `mafft` (for DNA or protein alignments, very fast) 
- `MACSE` (for in-frame, translation-based alignments, can handle frameshifts). 

In the past I have also used:
- `PRANK` (in-frame nucleotide alignment). Seems to do well in many published tests.
- `clustalw` (DNA/protein)
- `translatorx.pl` (in-frame nucleotide alignment, using Muscle, Clustalw, Prank, or mafft to actually perform the alignment), e.g. `translatorx.pl -i myDNAseqs.fa -o myDNAseqs.translatorx`

For VIEWING/editing alignments, here are some options:
- Seaview (mac)
- wasabi
- geneious

Could also check out:

- OPCAT - online pipeline to collect orthologs and make codon-aware alignment, also trims poorly aligned start/end regions.  Looks like it's not so useful for local use (?)
- DGINN pipeline (Lucie Etienne)


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

There are block-filter methods (remove entire columns from alignment) and segment-filter methods (remove part of whole of individual sequences).

I think also some hybrid methods?

Figure 1 of [the OrthoMam paper](https://academic.oup.com/nar/article/52/D1/D529/7318103) shows the pipeline they used for alignments - they do a lot of filtering. 15,868 alignments from 190 mammalian genomes. Website isn't working well as of Oct 11 2024.

### Filter out entire sequences
- [PhylteR](https://academic.oup.com/mbe/article/40/11/msad234/7330000)

### Block methods:

- trimal - can trim gappy positions out of an alignment (nucleotide or protein). Example usage:
```
trimal -in myAln.fa -out myAln.trimal.fa -automated1 -htmlout myAln.trimal.html -colnumbering > myAln.trimal.retainedColumns.txt
```
Some papers say trimal isn't as good as GUIDANCE and/or ZORRO

Also:
- BMGE
- GBLOCKs - some papers say this isn't as good as GUIDANCE

### Segment methods
- [HMMCleaner](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-019-1350-2)

### Other methods - maybe segment, maybe hybrid, or I just haven't looked

I think 
- `SWAMP`
- GUIDANCE (I used this in the kinetochore pipeline)
- ZORRO (Wu, Chatterji, Eisen, PLoS ONE, 2012)
- PREQUAL
- MACSEv2 "trimNonHomologousFragments" module


# Phylogenetic trees

[Phylogenetic Biology](http://dunnlab.org/phylogenetic_biology/) online book, but Casey Dunn

## some jargon

"ultrametric"  (right-aligned)

"tanglegram" = a plot to compare two trees, drawing lines between the same taxa in the two trees. 


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

First, change directory, because we have to run jmodeltest from a particular directory (pushd is similar to cd, but remembers which directory we were in before, and lets us use "popd" to return to that directory later):
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

## Choose best evolutionary model for amino acid alignments

Use `prottest` The help file is here: `/fh/fast/malik_h/grp/malik_lab_shared/prottest-3.4.2/README`

examples are here: `/fh/fast/malik_h/grp/malik_lab_shared/prottest-3.4.2/examples`

Input format can be fasta or phylip.  Should not contain * for stop codons. Long sequence names will be truncated, but that's OK.

First, change directory, because we have to run prottest from a particular directory (pushd is similar to cd, but remembers which directory we were in before, and lets us use "popd" to return to that directory later):
```
pushd /fh/fast/malik_h/grp/malik_lab_shared/prottest-3.4.2/
```

We also load the java module
```
module load Java/1.8.0_181
```

To get more help, including details on all options (need to do this from the prottest folder):
```
java -jar prottest-3.4.2.jar  -help
```
Then run prottest. The general form of the command I use is as follows. 
```
java -jar prottest-3.4.2.jar -i alignmentFile.phyml -all-distributions -F -AIC -BIC -tc 0.5 -threads xxxNumberOfThreads > outputFile.txt
```
I specify the full path to my input file and output file (usually put the output file in the same folder as the input file).  Number of threads should be between 1 and 12, but you must have requested that number of CPUs from the gizmo node you logged in to.

Here's an example on a real dataset:
```
java -jar prottest-3.4.2.jar -i /home/jayoung/PARPs/alignments/bigAln10.phyml -all-distributions -F -AIC -BIC -tc 0.5 -threads 12 > /home/jayoung/PARPs/alignments/bigAln10.phyml.prottest
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

## Miscellaneous other phylogeny-related things

There's the BALTIC python suite for drawing/manipulating trees, from Trevor Bedford's lab.

Stephanie Spielman wrote pyvolve - simulation tools


# Analysis of selective pressures

PAML

DataMonkey / Hyphy

[SWAKK](http://bioinformatics.mdanderson.org/main/SWAKK:Overview) website for sliding window dN/dS

# Ancestral reconstruction

For newer notes, see homing endonuclease project. Older notes:

There is a way to do it with PAML’s codeml on your local computer (run M0 and change the RateAncestor setting to 1, output will be in rst file). I tried this in cd ~/pamlPipeline/geneListsDone/list9_publishedPositivelySelectedGenes.txt_output/manualAnalysis/MX2_tryAncestralReconstruction – see NOTES file there.

But it might be easier to use one of these web servers:
1. http://fastml.tau.ac.il/ 
2.  https://academic.oup.com/bib/article/22/4/bbaa337/6042664   (it looks like you can access the linked website only from outside the Hutch firewall – it’s a Czech site, and I think they’re on some kind of list the Hutch is restricting - so you’ll need to be off campus and off the VPN)
3. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4966924/ (the linked server is not responding for me!)

This paper might be useful:
[Ancestral Sequence Reconstruction as a Tool to Detect and Study De Novo Gene Emergence](https://academic.oup.com/gbe/article/16/8/evae151/7713604?login=false&utm_source=etoc&utm_campaign=gbe&utm_medium=email)  Vakirlis, et al, 2024 Genome Biology and Evolution


# Detecting recombination

GARD (DataMonkey / Hyphy)

LDhat (McVean et al 2002).  Ziheng Yang used this in a 2004 paper.  He also used PLATO and PIST but commented that PLATO is susceptible to interpreting different selective pressures on different parts of the alignment as recombination, when it's not really. Can maybe get around that by analyzing only the third codon position.


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





# Ensembl database

It's confusing – first issue is that there are several Ensembl sites for different sets of species, with a little overlap:
‘main site’ = mostly vertebrates, but some info for selected others, e.g. fly, worm: http://uswest.ensembl.org/index.html
fungi: http://fungi.ensembl.org/index.html
other metazoa (mostly invertebrates flies): http://metazoa.ensembl.org/index.html
https://ensemblgenomes.org/
Non-vertebrates:  https://ensemblgenomes.org/
The multiple sites will eventually get merged in a future Ensembl release.  
Ensembl databases house a TON of information, but for now we are just looking at their phylogenetic treees. These have been generated in an automated way, so will give you a good first pass view of what’s going on with orthology and paralogy, but they might not be totally accurate.

The phylogenetic trees are at two levels:
a.	fairly closely-related species/paralogs (the same species sets found in whichever version of the Ensembl website you’re looking at, like vertebrates). 
b.	‘pan-taxonomic Compara’ - spans a broader range of species: if it was possible to align the genes reasonably well, you may see a tree that includes vertebrates and invertebrates, fungi, plants, protists, etc.  Fewer species in each genus/order.

How to look at Ensembl’s trees:
1.	Open Ensembl website and navigate to your gene of interest
2.	Find a link that says ‘gene tree’ (likely on the left).  For some species you will see two ‘gene tree’ links – one is pan-taxonomic compara, the other for only more related species. 
3.	Scroll to the bottom and click 'view fully expanded tree'.  Spend some time figuring out the display. Your starting gene will be highlighted in red, and parlogs from the same species highlighted in blue. At phylogeny branch points, blue nodes=speciation, red nodes=duplication. The green bars on the right depict the alignment, and can help you see that some sequences in the tree were badly aligned and maybe don’t belong in that location on the tree.  Don’t trust everything you see – duplications in just one species are often genome assembly errors. Hovering over a gene name can get you a link to more info on that gene in a particular species.

Note: unfortunately there is no direct link to pan-taxonomic Compara for human genes (or any vertebrate). This might get fixed in future releases.  It’s annoying.  If your starting point is a human gene, you can use the tree you see on the main Ensembl site to (hopefully) pick out a fly or worm ortholog, get the gene ID, then look it up on the Enseml Metazoa site where you should be able to accesss the pan-taxonomic tree. 

Note: at some point Ensembl changed the algorithms/parameters that group genes before making trees, and now they are more conservative. For genes that are not super well conserved, distant orthologs/paralogs are not always in the same pan-taxonomic compara tree. Sometimes it helps to go to an OLDER version of the Ensembl database to see if there was a larger gene grouping. See http://www.ensembl.info/2018/10/24/changes-to-paralogy-in-release-94/ for details, and maybe use this site: http://jul2018.archive.ensembl.org/index.html 
Note: in April 2021 there is a bug in pan-taxonomic trees where many sequences appear listed as ‘ancestral sequence’ instead of the species they’re actually from. Not useful. This should be fixed when release 51 comes out in early May.


# Genomicus
June 2014: http://www.genomicus.biologie.ens.fr/genomicus-75.02/cgi-bin/search.pl
Based on Ensembl 75 (released Feb 2014, still current as of June 2014)
