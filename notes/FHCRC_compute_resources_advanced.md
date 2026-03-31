# More sbatch tips and tricks

`scancel` examples: 
- cancel a single job using its numerical ID: `scancel myJobID`
- cancel all your own running jobs - e.g. to cancel all jobs for user `jayoung` you do this:  `scancel -u jayoung`
- cancel a numerical range of job IDs: `scancel 50088736-50088740` 
- or a discontinuous set: `scancel 50088736,50088740,50088780` 
- cancel all pending jobs for a particular user: `scancel -u jayoung -t PD`


Reasons a job might not be running: 
- `Resources` = There aren't resources available to run the job
- `Priority` = Higher priority jobs are queued ahead of this job
- `AssociationResourceLimit` = The job would violate a per-user or per-account limit
- `QOSResourceLimit` = The job would exceed a limit for the partition
- `ReqNodeNotAvail` = A requested node configuration is temporarily unavailable

Info on how [priorities](https://sciwiki.fredhutch.org/scicomputing/compute_job_scheduling/) are calculated.

Show information for a running or recent job:  `scontrol show job <job_id>`

To show info on jobs INCLUDING completed jobs
- sacct -u jayoung
- sacct -A malik_h | more

To show current per user / per account settings:
sacctmgr show qos public format=name,maxtresperuser,maxtrespa

A script to show who has what in the queue:
hitparade


Tools to look at priority and "fairshare" (based on recent usage):
- sshare | head
- sshare | grep 'jayoung'
- sprio -u jayoung

Show summary information on each queue (=partition) that's set up on our cluster: `sinfo -o "%P %C"`. Example output (March 2026):

```
PARTITION CPUS(A/I/O/T)    #### allocated/idle/other/total
campus-new* 5342/1038/304/6684
short 5342/1038/304/6684
restart-new 5342/1038/304/6684
interactive 5342/1038/304/6684
chorus 107/85/64/256
canto 0/108/0/108
```


Run a job in a specific partition: `srun -p queue_name`

To run a job on a specific node (I think?): `srun --nodelist hyraxD84 --pty bash`. There's also a way to exclude certain nodes, I think.



# Restoring files from backup

If you accidentally delete something, or want a previous version, go to the directory where those files lived, and try this command:
```
ls .snapshot
```
You'll see backups from various times that you can copy as you like. 

## Advice from Agus (March 2018) on restoring from backup

If we deleted a directory's entire contents:

Find out what time and day that the directory was last modified.
```
ls -al /fh/fast/malik_h/grp/michelle_yeast/fromMichelle
total 85
drwxrws---  2 jayoung malik_h_grp   0 Mar 15 17:49 .
drwxrws--- 17 jayoung malik_h_grp 650 Mar 15 17:13 ..
```
Remove the existing empty directory.
```
rmdir /fh/fast/malik_h/grp/michelle_yeast/fromMichelle
```
Copy the directory from the snapshot:
```
cp -pr /fh/fast/malik_h/grp/michelle_yeast/.snapshot/fast_hourly_2018-03-15-_17-00/fromMichelle /fh/fast/malik_h/grp/michelle_yeast
```

# Compiling code on the cluster, Aug 2020

There are a mix of architectures now in the cluster -make sure to compile on the ‘lowest common denominator’ nodes, the F class nodes. Use:
```
grabnode --constraint=gizmof 
sbatch --constraint=gizmog 
```

# A fix for slow SMB connections

https://osxdaily.com/2015/04/17/fix-slow-folder-populating-cloudkit-macosx/

on mac:  
```
rm ~/Library/Caches/CloudKit/CloudKitMetadata*;killall cloudd
```
   (may need to do it through Finder – permissions may be odd)

# Understanding load on a node
`htop`   
`uptime`

