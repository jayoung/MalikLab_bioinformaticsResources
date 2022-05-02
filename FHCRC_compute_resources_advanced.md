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

