# questions 2022, Oct 26

## NGmerge

**Question:**   

"(1) NGmerge is suddenly giving me an error message that “J” (a Q score of 41) is not accepted - I haven’t yet found an easy workaround, but I’m loathe to switch programs because I really liked it." 

**Ideas**

Other people also had this problem - this [github link](https://github.com/jsh58/NGmerge/issues/4) might give a solution, although I haven't read the details.   Let me know how it goes!


## fasta2counts

**Question:**

"(2) After all my pre-processing (trimming + read merging + quality filtering), I used FastX-collapser to collapse all my reads down to a fasta file with the # of counts for unique sequence (in an effort to avoid loading huge files into R). This gave me an awkward fasta file that looks like this:
```
>1_counts
seq
>2_counts
seq 
```
And I’m not quite sure how to extract that into a table with a seq and counts column. Perhaps there’s a better way to do this in Terminal (with sort uniq?); fastX-collapser was just easy and quick."

**Solution:**

This should be very doable in R - see [fasta2counts_andMore.R](codeBits/forJeannette/forJeannette_2022_Oct26/fasta2counts_andMore.R)

**Alternatives**

It should work fine on moderately sized files, but it's possible that will run too slowly on very big files in R. Try it and see!  If that's true, we could instead make a script to run from the command line, or there are ways to optimize that in R.


## match nuc seqs to pre-defined list

**Question:**

"(3) After this collapsing / counting unique step, I need to diverge my workflow:

(a) for one set, I need to match the sequences to a list of pre-defined sequences I ordered"

**Ideas**

We talked a bit before about the `match()` function in R.  There's another tidyverse function `left_join()` that may also do well with this.   Did you get anywhere with `match()` ?

## translate nuc seqs and get distance from wild-type

**Question:**

"(b) for the other (DMS set), I need to do the following (which the biostatisticians did for me previously and would never respond when I asked for the code):
translate my sequences
define the distance from WT in terms of nt changes and amino acid changes, as well as the position of the first amino acid change"

**Solution:**

Translation, distances, first change positions: also very doable in R. See [fasta2counts_andMore.R](codeBits/forJeannette/forJeannette_2022_Oct26/fasta2counts_andMore.R)


## compare read counts across conditions

**Question:**

"4) Then I need to compare the number of reads for each sequence in each of 4 different conditions"

**Ideas**

R's `left_join()` may work well here, too. 