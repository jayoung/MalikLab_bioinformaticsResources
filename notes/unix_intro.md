# Introduction to unix ("the command line")

## General notes

File names/directory names ARE case-sensitive.

File names/directory names should NOT contain unusual characters, for example:  `spaces * , /  \ + ! ( ) ~ $`. These characters are OK: `. - _`

## Basic unix commands:

`cd` - change directory (directory = folder)  
`ls` - list files in current directory  
`ls -l` - list files, including their size and date stamp  
`ls -lrt` - list files, including their size and date stamp, and sort in reverse temporal order  
`..` - the directory above the current one  
`.` - the current directory  
`*` - one or more wild-card characters  
`?` - a single wild-card character  
`/` - use / to string together longer directory paths    
`mkdir` - make new directory - specify the name of the new directory (e.g. `mkdir newDirectoryName`)  
`rm` - remove a file or files  
`rmdir` remove a directory (only works if it's empty)  
`rm -r directoryName` - removes a directory and ALL of its contents (-r = recursive)  
`pwd` - show working directory  
`more file.txt` view a text file page-by-page. Use the space bar to see the next page, the return key to see just one more line, or the `q` key to quit looking at the file. (`cat` and `less` also let you look through files)  
`head file.txt` or `tail file.txt` – show the first few or last few lines in a file. The default is 10 lines: to see 20 lines, for example, do this: `head -20 file.txt`  
`ctrl-c` - interrupt/cancel a command that's running

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

## Slightly less basic unix commands:

### `mv` to move or rename files

Examples:
- `mv oldfilename.txt newfilename.txt` - rename a file  
- `mv file1.txt file2.txt my_folder_to_move my_destination_folder` - move one or more files into a folder  

`mv` can be a bit confusing, as it can move files OR rename them, depending on how you use it:
- If the last part of the command is the name of an existing folder (e.g. `my_destination_folder` in the second example above), `mv` will move any other files/folders specified in the command into that folder.
- If the last part of the command doesn't yet exist (e.g. `newfilename.txt` in the first example above) `mv` will rename the file

### `cp` to copy files
Examples:
- `cp oldfilename.txt newfilename.txt` - make a new copy of the file, with a different name
- `cp oldfilename.txt my_destination_folder` - make a new copy of the file, within a different folder (if you don't specify a file name within that folder, the new file will have the same file name as the original)


## Adding options to commands

Many commands have options that change their behavior. A google search can tell you more, especially if you include the word `man` (for manual) and `unix` or `linux` in your search terms  (e.g. search with `man ls unix`).

A typical way to add an option is with the `-` sign, e.g. `ls -l`. 

You can specify multiple options, e.g. `ls -l -r -t`. 

For options specified using a single letter, you can also combine them, e.g. `ls -l -r -t` can also be written `ls -lrt`

Some options have a long form and a short form, and the long form is often specified with two `-` signs, e.g. `ls -t` is the same as `ls --time`.

## Relative versus absolute file paths

The "path" to a file means the file name plus its location on the file server.

A "relative path" means the path relative to your current working directory (e.g. `../tetrahymena/annotations` would mean something different when you're working in different directories).

An "absolute path" starts with the `/` character, and specifies a file's location beginning from the very top level of the computer's directory structure.  Absolute paths are useful as you can use the same path no matter what your working directory is.

## Making links to files

Sometimes, especially if we are working with large files, we don't want to **copy** a file to a second location, but we do want to make a link (=alias) so that it's easy to access. 

We pay for each Gb of data we store, so please try not to make unnecessary copies. Making copies also increases the risk that you'll edit one copy of the file but not another: using links to a single 'master' copy avoids that risk.

Example: Illumina sequencing files from the Hutch facility get put in directories with complex names, inside `/fh/fast/malik_h/SR/ngs/illumina/`. We might want to use those large files from the directory where we're doing a bunch of analysis, so we make a link using the `ln -s` command:
```
cd myAnalysisDir
ln -s /fh/fast/malik_h/SR/ngs/illumina/mhays/181029_D00300_0638_AHNC55BCX2/Unaligned/Project_mhays/* .
```
You can delete links without deleting the original data file.

The `-s` means soft link - it's the only kind of link I ever use. If you want to know more, google 'soft hard link unix'.


## Capturing screen output from a unix command

There are two types of output that get put up on the screen:   
- "standard output" = most of the output of any command (also known as stdout)  
- "standard error" = some commands produce warnings/errors that go into this output stream (also known as stderr)

You can save those outputs to files using `>` and similar characters. Examples: 

`ls > myListFiles.txt`  - the `>` lets you capture the standard output of a command and put it in a file. It will overwrite the file if it already exists

`ls >> myListFiles.txt`  - the `>>` lets you capture the standard output of a command and add it to a file. If the file already exists, standard output is added to the end of the file, but if not, it creates the file

`ls > myListFiles.txt 2> myStdErr.txt` - puts stdout into myListFiles.txt and stderr into myStdErr.txt.  (`ls` is a bad example for this because it does not produce any stderr)

`ls > allScreenOutput.txt 2>&1`  - puts both stderr and stdout into allScreenOutput.txt  (`ls` is a bad example for this because it does not produce any stderr).  

More detail: The `2>&1` (means "Redirect stderr to where stdout is currently going") must come AFTER the initial `> file.txt`. "In other words, the `&1` reuses the file descriptor which stdout currently uses."


## Some keyboard shortcuts

(these work on many computers but not all)

The `tab key` can help you complete file names on the command line, avoiding typing. Might need to hit it twice.

`Up/down arrow keys` scroll through previously used commands.

`escape-then-f` / `escape-then-b` : move cursor forward/backwards through the command line you're typing one word at a time

`ctrl-e` / `ctrl-a` : move cursor to end/start of line

`escape-then-delete` : delete one word at a time on the command line


## Next level unix commands, if you’re interested in more:

`grep`: search text/files for a string or pattern, e.g. show me the header lines in a fasta file:  
     `grep ‘>’ mySeqs.fa`  
To learn useful options for grep, google ‘grep man page’ and look through one of the hits.  I often use these options: `-i -v -l -L`

`wc`:  ‘word count’ – counts lines and words and bytes in a file. Example: `wc myDat.txt`

`diff`: detects and displays any differences between two files.  Example `diff file1.txt file2.txt`

`cmp` - a quicker alternative to `diff` that stops when it finds the first difference: `cmp --silent file1 file2 || echo "files are different"`


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



