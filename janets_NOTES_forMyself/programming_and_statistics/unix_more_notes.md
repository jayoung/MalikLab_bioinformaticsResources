# more unix notes

intro to unix notes are in notes/unix_intro.md 

this contains more stuff that will be too much for the average user 

control-s suspends screen (don't do it)
control-q resumes screen (if you accidentally did control-s)

# capturing output/error

to capture standard error, switch to bash shell (exit to get out again) and 
use 2> to redirect.

>& redirects both STDERR and STDOUT

>> appends instead of overwriting

` backticks to interpret a command 
e.g. grep 'lib' `find /traskdata -name *dgb`

# miscellaneous useful commands

## `ls`

`ls -lwc -l` returns the number of files in a directory

`ls -R` list all subdirectories too

`ls -d` provides additional path information


## `wc` options (counting)

wc -l file.txt
wc -w file.txt
count lines or words (or use -m for characters)

## `mklinks`

`mklinks` "Samples" "Path"
is a script in `/home/jayoung/malik_lab_shared/bin` which takes all filenames in "Samples" list, in
directory "Path" and makes soft links to them

## `du` disk usage

`du -sk filename`
gives filesize in kb (s = summarize, k= kb)

`df -h` (human readable)


## python modules

installing python modules:
  e.g. to install module called pysam:
      pip install --user pysam
(and I have a useful environmental variable set:     PYTHONUSERBASE=/home/jayoung/malik_lab_shared )

## linux versions

uname -a

config.guess to determine system architecture

which linux version am I using?
cat /proc/version

# Libraries

`locate libxml2` (on mercury nodes)

the locate tool is sometimes mlocate, slocate and there may be others.

version numbers: sometimes this works:
/lib/libc.so.6 
Sometimes not...

# shell scripts

`-n` tests to see if the argument is non empty
`-z` tests whether operand has zero length

might be a good idea to add this line to every script:
    `set -eu`
it ensures that if something fails in the script, the WHOLE script will fail and an error message will be printed

or a version that also prints every executed line:
    `set -eux`



In shell scripts, the `$status` variable is useful – it’s a shell variable (in tcsh and csh) indicating whether the last command was successful.

A for loop:
```
foreach dog (1 2 3) 
   echo "$dog"
end
```

Using curly brackets: `$dog` is the same as `${$dog}`, but using the parentheses sometimes help avoid confusion.  

Double quotes and single quotes act differently in terms of how variables are evaluated.



# vi editor

command mode:
x is backspace
:w is save
:wq is quit with save
:q! is quit without save
esc once or twice to get into compose mode


# xemacs editor

to turn on line number
esc-x line-number-mode

# which

e.g. which grep
tells the path of the default 'grep' in use. useful in conjunction with man
pages. to use another version the whole path can be specified

# printing direct from the server

`psprint` sends file to printer

`lpstat -t` lists all print queues

`lpr -Pqueue_name` to print

`lpq -Pqueue_name` shows queue

`lprm -Pprinter_name job_number` removes job from queue


# kill

`ps -ef | grep janet`
to see what the jobs are, then
kill (using first number of row)

jobs -l includes process id


# tar

`tar -xzvf 'filename'` to extract
`tar -cf 'newfilename' 'oldfilename'` to compress

x - extract
z - if needed to unzip too
v - verbose (lists files)
f - from file

with gtar:
`gtar -cf - aceGENOME | gzip -c - > aceGENOME.tar.gz`

`tar -cf - aceGENOME | gzip -c - > aceGENOME.tar.gz`

`gunzip < foo.tar.gz | tar xf -`

`gunzip < foo.tar.gz | tar xfi -`
the i option on tar might help get past directory checksum errors

`tar cvf multiplefilesnotindir.tar *.txt *.scf`

compress with gtar:   `gtar cvzf file.tar.gz filenames`
uncompress with gtar: `gtar xvzf file.tar.gz`

`bunzip2 bwa-0.5.9.tar.bz2 `

Uncompresssing multipart gz file:
cat APMT01.fsa*gz > MesAur1.0.fa.gz.gz
gunzip MesAur1.0.fa.gz.gz
gunzip MesAur1.0.fa.gz
It's a little odd that it's twice gzipped, but perhaps that's just how NCBI does it


# find

find /traskdata -name '*pl' 
finds within directory tree /traskdata *pl files
find /traskdata -name '*pl' -exec grep janet \{\} \;
-type d (directories) or f (files), etc.

find directories within the current dir, and then change their positions
find /path/to/base/dir -type d -exec chmod 755 {} +

e.g. to allow Michelle to look in folders on the fast drive using her Mac's finder, they seem to need to have group sticky bit set, e.g. 
   drwxrws--- 78 jayoung malik_h_grp   2719 Mar  2 11:44 UCSC
but I don't want to set that stick bit on files as it seems to screw them up
so here are some commands that finds the directories (recursively) and changes the permissions
    cd /fh/fast/malik_h/grp/public_databases
    find . -type d -exec chmod g+s {} +

# grep

grep -l option means it prints only filenames containing the matching string, not the lines with the string
grep -h option means it prints only the line containing the matching string, not the filenames

grep -v (excludes those lines containing string)
grep -n (shows line numbers)

# laj and blastz
to view blastz alignments using laj:
java -jar ~/laj.jar -noseq
(then get a menu where you can load in the blastz file)


# editing files

## sed

stream editor - could be really useful. see http://www.grymoire.com/Unix/Sed.html

show a chunk of a file using line numbers:

sed -n 5,8p file   (shows lines 5-8)

https://unix.stackexchange.com/questions/288521/with-the-linux-cat-command-how-do-i-show-only-certain-lines-by-number/288523

## awk 

awk is also useful for streaming through files. that's what I've been using for bam files.


# moving files around

## wget

wget http://www.bcgsc.bc.ca/downloads/fpc_mousemap.tar.gz

wget -r -l1 --no-parent http://bio.cs.washington.edu/assessment/parameters/

## rsync

rsync is a quicker way to transfer files than ftp (and maybe more robust?) this is how I transferred files from fred onto a local Mac.
rsync --verbose  --progress --stats --compress --rsh=/usr/bin/ssh \
      --times \
      jayoung@rhino.fhcrc.org:/home/jayoung/tetrahymena/bwa_MIC_bothBatches/* .
I had to edit my .tcshrc file for this to work (can't have it output any text)

some important options to think about using or not:
--recursive
--copy-links
rsync --verbose --progress --stats --compress --rsh=/usr/bin/ssh --times --recursive --copy-links jayoung@rhino.fhcrc.org:/home/jayoung/tetrahymena/FASTQ_files .

## scp 

copying across computers:
```
scp fred:~/filename.txt .
scp jayoung@fred.fhcrc:~/filename.txt .
```
see http://kb.iu.edu/data/agye.html
When using wildcards (e.g., * and ?) to copy multiple files from a remote system, be sure to enclose the filenames in quotes. This is because the Unix shell, not the scp command, expands unquoted wildcards.


# setfacl and permissions

Not sure if this still works

Show permissions in a more sophisticated way
```
getfacl bin
```

Change permissions for multiple groups at once
```
setfacl -m group:Malik_H_lab:r-x bin
setfacl -R -m group:Malik_H_lab:r-x help
```
-R is recursive


To copy the parent dir's current permissions into default permissions for things within it:
```
getfacl --access bin | setfacl -d -M- bin
```
so, should always do this:

```
setfacl -R -m group:Malik_H_lab:r-x bin
getfacl --access bin | setfacl -d -M- bin
```


# email using `mailx`

to email a file from the command line:
```
cat fourRodents.fa | mailx -s "testSubject" jayoung@fhcrc.org
```
it comes as text in the body of the message - seems there's no easy way to attach it as a file

# latency issues:

For a while the gizmo/rhino cluster had some weird issues with slowness accessing files, running `ls`, etc. To troubleshoot, scicomp were asking me to do this:
```
pwd; date; time ls -l > /dev/null
```
