# Introduction to unix ("the command line")

## General notes

File names/directory names ARE case-sensitive.

File names/directory names should NOT contain unusual characters, for example:  `spaces * , /  \ + ! ( )`. These characters are OK: `. - _`

## Basic unix commands:

`cd` - change directory (directory = folder)  
`ls` - list files in current directory  
`ls -l ` - list files, including their size and date stamp  
`..` - the directory above the current one  
`.` - the current directory  
`*` - one or more wild-card characters  
`?` - a single wild-card character  
`/` - use / to string together longer directory paths (or, if it is at the beginning of a file path, or on it's own, means the very top level of the directory structure).  
`mkdir` - make new directory - specify the name of the new directory (e.g. `mkdir newDirectoryName`)  
`rm` - remove a file or files  
`rmdir` remove a directory (only works if it's empty)  
`rm -r directoryName` - removes a directory and ALL of its contents (-r = recursive)  
`pwd` - show working directory  
`more file.txt` view a text file page-by-page. Use the space bar to see the next page, the return key to see just one more line, or the `q` key to quit looking at the file. (`cat` and `less` also let you look through files)  
`head file.txt` or `tail file.txt` – show the first few or last few lines in a file. The default is 10 lines: to see 20 lines, for example, do this: `head -20 file.txt`


## Examples: 
list all files with names ending in `.bam` in the current directory: `ls *.bam`  
list all files with names starting with `Sample` and ending in `.bam` in the current directory: `ls Sample*.bam`  
list all files in a subdirectory of current directory: `ls tetrahymena/annotations`  
list all bam files, in annotations folder which is a subdirectory of the tetrahymena folder: `ls tetrahymena/annotations/*.bam`  
or, if tetrahymena is a directory above the current directory: `ls ../tetrahymena/annotations/*.bam`  

change working directory to the home directory: `cd`  
change working directory to the folder above the current one: `cd ..`  
change working directory to the tetrahymena/annotations dir: `cd tetrahymena/annotations`  

use the `resize` command, if you resize the window after logging in and things look odd

## Capturing screen output from a unix command

There are two types of output that get put up on the screen:   
- "standard output" = most of the output of any command (also known as stdout)  
- "standard error" = some commands produce warnings/errors that go into this output stream (also known as stderr)

You can save those outputs to files using `>` and similar characters. Examples: 

`ls > myListFiles.txt`  - the `>` lets you capture the standard output of a command and put it in a file. It will overwrite the file if it already exists

`ls >> myListFiles.txt`  - the `>>` lets you capture the standard output of a command and add it to a file. If the file already exists, standard output is added to the end of the file, but if not, it creates the file

`ls > myListFiles.txt 2> myStdErr.txt` - puts stdout into myListFiles.txt and stderr into myStdErr.txt.  (`ls` is a bad example for this because it does not produce any stderr)

`ls > allScreenOutput.txt 2>&1`  - puts both stderr and stdout into allScreenOutput.txt  (`ls` is a bad example for this because it does not produce any stderr)


## Some keyboard shortcuts

(these work on many computers but not all)

The tab key can help you complete file names on the command line, avoiding typing. Might need to hit it twice.

Up/down arrow keys scroll through previously used commands.

escape-then-f  / escape-then-b go forward/backwards through the command line you're typing one word at a time

escape-then-delete - delete one word at a time on the command line


## Next level unix commands, if you’re interested in more:

`grep`: search text/files for a string or pattern, e.g. show me the header lines in a fasta file:  
     `grep ‘>’ mySeqs.fa`  
To learn useful options for grep, google ‘grep man page’ and look through one of the hits.  I often use these options: `-i -v -l -L`

`wc`:  ‘word count’ – counts lines and words and bytes in a file. Example: `wc myDat.txt`

## Linking commands together 
We can use the `|` symbol (the ‘pipe’) to link commands together. Example: show me the header lines in a fasta file and count them:
     `grep ‘>’  mySeqs.fa | wc`

`zcat`:  looks inside a gzip compressed file without changing the file.  
Example: show me the first few lines of a fastq.gz file:  
    `zcat mySeqs.fastq.gz | head`  
Example: count lines of a fastq.gz file:  
    `zcat mySeqs.fastq.gz | wc`

`gunzip / gzip` : compress/uncompress a file (does change the original file)

`history`: shows recent commands, e.g.   
    `history`  
    `history | grep ‘blastn’`  



