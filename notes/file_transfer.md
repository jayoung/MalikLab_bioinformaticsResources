# transferring files from one computer to another

There are many ways to do it. It's easy to use the Mac finder to copy small files from your Mac to/from the Hutch servers. But for larger files, we usually choose other methods, because they're faster and more tolerant of being interrupted.

Some advice from the Hutch IT people [here](https://sciwiki.fredhutch.org/scicomputing/store_collaboration/)

## Aspera

Works well. See instructions [here](https://aspera.fhcrc.org/index.html)


## rsync

A command-line tool to transfer files between storage systems.  

rsync is probably unnecessary if you can get Aspera working, but here are some old notes I put together when I was using it.

Rsync works if you have a login on both computers, and if rsync is installed on both systems.  

First figure out if rsync is installed on both systems: login to each, and try this command: `which rsync`.  Hopefully you’ll see output something similar to this (could be different, you just want to see a path that ends in rsync):  
`/usr/bin/rsync`  
If nothing comes up, then rsync is not installed, and you need to use a different option.

The general form of rsync command has four parts. In this example we are working on another computer (i.e. NOT on rhino) and we want to transfer some files from rhino onto the local machine. Here are the four parts of the command:
- a. `rsync`  
- b. the options for rsync (e.g. `--verbose  --progress --stats --compress --rsh=/usr/bin/ssh --times`).   Include `--recursive` if you want to go into subdirectories and get those files as well.  
- c. the files to transfer (in this case, on a remote machine called rhino: `jayoung@rhino.fhcrc.org:/home/jayoung/janet_analysis/annotations/455-Mseq-FASTA2*`)  
- d. where to put them (in this case, the current directory on the local machine, i.e. `.`)

### using rsync to transfer files from fred/rhino system to local computer. 
three example commands are below:
- (i) transfer a few files 
- (ii) transfer a whole folder, including subfolders 
- (iii) transfer a folder that also includes links (aliases)

(i) transfer all files whose name begins with "455-Mseq-FASTA2" from the janet_analysis/annotations directory on fred to the current directory (".") on the local server:
```
rsync --verbose  --progress --stats --compress --rsh=/usr/bin/ssh --times jayoung@rhino.fhcrc.org:/home/jayoung/janet_analysis/annotations/455-Mseq-FASTA2* .
```
(then enter your password)

(ii) transfer a folder, with everything that's inside it. We simply add the "--recursive" option:
```
rsync --verbose  --progress --stats --compress --rsh=/usr/bin/ssh --times --recursive jayoung@rhino.fhcrc.org:/home/jayoung/janet_analysis/tables .
```

(iii) transfer a folder with symbolic links inside (link = alias), we add the "--copy-links" option:
```
rsync --verbose --progress --stats --compress --rsh=/usr/bin/ssh --times --recursive --copy-links jayoung@rhino.fhcrc.org:/home/jayoung/janet_analysis/FASTQ_files .
```

### using rsync to transfer files from local computer to fred/rhino system. 
```
rsync --verbose  --progress --stats --compress --rsh=/usr/bin/ssh --times filesOnMyMac*  jayoung@rhino.fhcrc.org:/home/jayoung
```

