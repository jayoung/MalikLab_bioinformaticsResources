# how to run alphafold

## Melody's request (Apt 30, 2024)

Hi Janet,

I wrote a small program for the cluster to run basic alphafold in a really straightforward way. I wanted to ask if you would be willing to test it and let me know if you find any mistakes with the documentation or the program itself (although the program has been tested).

https://research.fredhutch.org/campbell/en/etc/alphafold.html

Jeremy recommended you as someone who might find it useful. If there are others in the Malik lab who at least have cluster access (and possibly less terminal experience) and would be willing to test out the documentation, that would be great, too. I wanted to make it as accessible as possible and I anticipate the most difficult part will be getting the sequence file on the cluster. I also wanted to make sure there aren’t big mistakes before I send it to a larger crowd.

Thanks so much!
Melody



## Trying it

```
cd ~/FH_fast_storage/git_more_repos/MalikLab_bioinformaticsResources/janets_NOTES_forMyself/test_MelodyCampbell_alphafoldScript

# copy script
cp /fh/fast/campbell_m/share_scripts/share_run_alphafold_cluster_v1.0 .

# get a query fasta file
cp ~/FH_fast_storage/miscMalikLabGenesOfInterest/cenH3/getSeqs/human_CENPA_ORF.peps.fa .

# I manually made a version without the stop codon at the end
human_CENPA_ORF.peps.noStop.fa
```

it's an interactive script:
```
./share_run_alphafold_cluster_v1.0
```

I ran it on human_CENPA_ORF.peps.noStop.fa with the monomer model


It is possible to run it on several seqs, in two ways:
- a. asking alphafold to test for multimers (put all seqs in one file)
- b. running alphafold on each separately (specify all files, comma separated)

If you are running mulitmer, put all monomers in the same fasta file.
Example: ITGAM_ITGB2.fasta (file contains sequences for ITGAM and ITGB2 each with a >title and sequence)
 
If you want to run multiple sequentially, list all files separating via a comma.
Example: ICAM1.fasta,SINV.fasta,talin.fasta 

## my feedback

Hi Melody,
 
Sorry I didn’t get to this until today – I just ran my first job.  Making alphafold run is easy, thanks to your script and good instructions, but knowing what to do next is less obvious.
 
All my feedback is biased because of my own pre-existing (lack of) knowledge: I’m decent with computer stuff and almost totally ignorant about protein structure.  I hope to get less ignorant before long – I’ve been thinking about taking that class that Barry and you have been teaching 😊   But if there’s any sort of basic intro materials you recommend for understanding structures, I’d love some pointers (and maybe it’d be good to have a link on this page).
 
A big picture thing (for sure biased by my own skillset) - instead of including so much stuff in this particular page about file transfer, basic unix and cluster-related work, you could point people towards external resources.  I find it a bit distracting here – it might be nice to focus your page on how to run alphafold and how to understand the results.    For example, you could point people towards these resources from Hutch DASL and shorten your notes section – DASL have been trying hard to put together good resources:
Basic unix - https://hutchdatascience.org/training/#introduction-to-command-line
Cluster - https://hutchdatascience.org/training/#cluster-101
 
I’d like to see more of the following sorts of material. Some of these might be bigger questions than you want to address here, and I just need to read and learn more (embarrassingly, it’s been 30 years since I took a biochemistry class). But I suspect you’ll get some users a bit like me, whose PIs have said “why don’t you run alphafold on protein X” and they have no idea what to do next. 
Roughly how long will my cluster job take? (I’m sure it’s a range, but is it a few minutes, a few hours, a few days?)
What do I do with the output files?   I can see from my test run (human CENP-A) there are various pdb files, and various pkl files (as well as the msas dir).  Which file(s) should I look at?    A quick google search suggests I can probably ignore the pkl files, and I know that pdb files are likely useful. I’m assuming it’s the “ranked” ones I should use?
do you have a favorite viewer for pdb files, or a shortlist of viewers?  I found this Mol-star server online– it seems to be intuitive: https://www.rcsb.org/3d-view   but maybe there are better choices
for my output, I loaded three files into that viewer (ranked_0.pdb and _1 and _2). To my eye they look identical. Is there some sort of reporting summary of how different the top-ranked models are from each other? Any sort of scores to help me see whether the top one is way better than the rest or they’re very close in score?  After a bit of googling I figured out how to superimpose them in Molstar and I can see they’re not identical.
Maybe it would be good to provide an example output folder that a user can copy to practice looking at predicted structures (integrin, or whatever)?  If we were to get nonsensical output with our own input seqs, it might help us assess that.
Do you have a favorite reference or two on what alphafold is doing?  And/or any on what we can learn from protein structures?
I am aware that there are collection(s) of pre-computed alphafold outputs that most of us might look at before we run a prediction locally,  e.g. https://alphafold.ebi.ac.uk/ looks super useful. What situations would make me want to run my own prediction rather than use pre-computed versions?  Might be good to link to those resources, so that the cluster doesn’t get bogged down with people running things that exist in the world already.
It looks like there are also web-based ways to run alphafold on your own seqs, (and anyone who’s scared of the cluster/unix will gravitate towards those). E.g. Ching-Ho sent this to our lab today for v3 - https://golgi.sandbox.google.com/about   and it looks like there’s a v2 server too.  Is it useful to have a few words on why you’d use the cluster rather than one of those websites? (I imagine compute limitations, and the annoyance of web-based things for people like me)
Is alphafold robust to weird characters in the protein seq, or do we need to be careful to edit those out? (e.g. mine sometimes contain a * for the stop codon).  
More generally, is alphafold robust (does it always run without crashing)?  If not, maybe tell us where we should look for error messages.   If it gives error messages are they understandable or cryptic? (if cryptic, give us a few clues on how to troubleshoot).  Is it useful to warn people not to run it on hundreds of proteins that are 1000s of amino acids long?
Maybe include a link to alphafold documentation, in case we want to read about the differences between the models (monomer versus multimer is obvious, but the three monomer versions are less obvious)
 
And you probably didn’t want the following level of pickiness (sorry! I can’t help myself) but here you go:
Some typos etc in the webpage:
“You can import a fast file”  should be fasta
“(Please excuse my vocal fry; it's just my voice.)”   you have nothing to apologize for!  (I know it’s easy to be self-conscious about that, but it’s useful to set an example to the trainees not to be apologetic about ourselves, especially as women)
“it is unlikly to be your first name.” - unlikly
A typo in the echo commands from the script:
“If you are running mulitmer,”   - misspelled mulitmer
 
That’s it.   More than you wanted, I think!   Sorry.
 
It was useful for me to start thinking about this. I’ve been curious about the whole alphafold bandwagon but had not jumped on until now – I’m way overdue, so thanks for the nudge!
 
Janet
 