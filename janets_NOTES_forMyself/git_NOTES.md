xxxx there are more git notes somewhere else!

## Permissions

The Mac's smb mount of /fh/fast shows permissions very oddly.  For example, permissions on the same file look like this:

From gizmo:
```
jayoung@gizmok29:~/FH_fast_storage/cromwell-home/FH_WDL101_Cromwell$ ls -l README.md 
-rw-rw---- 1 jayoung malik_h_grp 1486 Nov 30 12:21 README.md
```
From the Mac:
```
bsjanetmacmalik:FH_WDL101_Cromwell jayoung$ ls -l README.md 
-rwx------  1 jayoung  staff  1486 Nov 30 12:21 README.md
```

This can sometimes make VScode (running on the Mac) think that there have been a lot of changes to a file, when it's simply that it understands the permissions differently than it should.

To turn off permissions checking globally:
```
git config --global core.filemode false
```

and/or for an individual repo:
```
git config core.filemode false
```

There will be situations where I DO want to track file permissions. A [Stack Overflow post](https://stackoverflow.com/questions/1257592/how-do-i-remove-files-saying-old-mode-100755-new-mode-100644-from-unstaged-cha/1257613#1257613) suggests a couple of ways to do that:
- toggle the config back and forth for particular commits (`git config core.filemode false` / `git config core.filemode true`)
- directly modify permissions in the repo using a command something like this: `git update-index --chmod=(+|-)x <path>`

## Git config

To show config:
```
git config --list

# global settings stored here:
~/.gitconfig 

# local settings stored here: 
.git/config
```
It is possible to edit those files, but we can also use `git config` to do it. Examples of modifiying a particular setting (permissions checking), globally or locally:
```
git config --global core.filemode false
git config core.filemode false
```







## Forking a repo

There's a 'fork' button on the github website that lets me add a fork of any repo to my own account (https://github.com/jayoung).  I can fork just one branch, or all of them.

Once I've forked the repo, I can either clone it to make a local copy to work on, or for an existing clone, I can change the location it will push/pull from:

Example: 
```
cd ~/FH_fast_storage/cromwell-home/FH_WDL101_Cromwell

# I originally cloned the repo from the master version:
git clone https://github.com/fhdsl/FH_WDL101_Cromwell.git

# Later, I decided I wanted to sync to my fork rather than the master, so we set the remote URL:
git remote set-url origin https://github.com/jayoung/FH_WDL101_Cromwell.git

# show the remote URL:
git remote -v
   # origin  https://github.com/jayoung/FH_WDL101_Cromwell.git (fetch)
   # origin  https://github.com/jayoung/FH_WDL101_Cromwell.git (push)
```

Assuming the fork exists, push/pull/status etc will now look at my fork not the master repo.

## Pull requests

I can use a button on the github website from my own fork ([example](https://github.com/jayoung/FH_WDL101_Cromwell.git)) to submit a pull request. It's under the 'Contribute' dropdown menu.

