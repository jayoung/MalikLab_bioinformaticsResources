# Phylogenetics advice

This doc is for general advice, partly to help me remember the topics I like to introduce people to when we start learning how to do phylogenetics/phylogenomics.

More detailed notes on running some individual algorithms are [here](https://github.com/jayoung/MalikLab_bioinformaticsResources/blob/main/notes/bioinformatics_tools.md)

these notes are chaotic right now

# Evolution and more

You have a new gene and not much is known. What do you do? Look at resources that help you learn about

- conservation/evolution
- expression
- protein structure / domains / motifs
- protein-protein interactions
- genetic information (e.g. knockouts)


# Evolution things you'll want to learn

Using public resources

Obtaining sequences

Homology searches
- blast
- maybe also hmmer etc

Making alignments

Making trees (phylogenetics)

Synteny

Assessing conservation / evolutionary rates / selective pressures

## Evolutionary resources to check out

You can get a lot of information from pre-computed public resources. Usually they're correct, but not always - to be confident about something you'll probably want to do your own analysis. They're especially useful to get a first idea about what's going on with a gene or gene family.



### [UCSC genome browser](https://genome.ucsc.edu/)

For some species, e.g. human, there are some very useful conservation tracks.

#### How to use the browser

If you've never used the browser before, maybe you want to check out some of their [training videos](https://genome.ucsc.edu/training/index.html). 

My quick advice is just to explore. Choose a reference genome for whichever species you are interested in, and search for a gene you care about. Once you're on the browser view:
- notice that if you scroll down there's a list of tracks you can turn on and off (or expand/contract)
- some reference genomes (e.g. human hg38) have tons of tracks, some have very few
- your view will start with a default set of tracks but you'll want to customize it
- in the track-choosing area, each track has a clickable description above the pull-down. That'll give you (a) an explanation of how that track was generated (b) more fine-grained control over what it's showing

#### Useful tracks for evolutionary analysis

"Comparative Genomics" track group. Individual tracks will be called something different for each species/reference genome - this example is for the [human ARC gene](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A142610191%2D142615337&hgsid=4161575949_FsSKUfATEVjjyghhhYGmwxFUAqAv) in the hg38 version of the genome assembly.
- "UCSC 100 vertebrates" - Click on the blue link with that name, and use the "+" box next to "species selection" to turn on ALL the species, then submit.  This represents 100 genome assemblies all aligned back to the human genome

xxx do I want full description?  maybe

### Ensembl browser and trees

New version of the Ensembl website in summer 2026, but trees are not yet available there, so we need to use the archive version of the website. 

Archive versions group species into different browser websites (e.g. human and other vertebrates you'd use the 'main' browser, and fly/worm I think you use the 'Ensembl Metazoa' browser). Not sure what they're doing for the new version.

Select your gene of interest, try to find a link that says 'gene tree'. You'll probably want to click 'fully expand tree' once you get there. Your query gene's name will be highlighted in red. Any paralogs from the same species will be highlighted in blue.

At some point a few years ago they changed the criteria for grouping genes together before making trees - now they tend to split gene families (or diverged orthologs) apart more than they used to.  If you prefer the bigger trees, try looking in older archives.  See note on front page of Genomicus browser about which Ensembl versions to use.

### Genomicus

Uses Ensembl tree data and genomic positions to create viewable synteny maps

### OrthoDB

### NCBI gene
