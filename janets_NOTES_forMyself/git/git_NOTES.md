xxxx there are more git notes somewhere else!

## Big picture


### git status versus fetch versus pull [advice](https://www.freecodecamp.org/news/git-fetch-vs-pull/)

There generally are at least three copies of a project on your workstation:
- one copy is your own repository with your own commit history (the already saved one, so to say).
- the second copy is your working copy where you are editing and building (not committed yet to your repo).
- the third copy is your local “cached” copy of a remote repository (probably the original from where you cloned yours).

`git status` only compares the local working copy to the local repo (i.e. looks for things you haven't committed yet)

`git fetch` is the command that tells your local git to retrieve the latest meta-data info from the remote, but it does not copy those changes to the current working directory. If you did `fetch` then `status` you would see any changes that were in the remote that aren't yet in your working dir.

`git pull` brings metadata from the remote AND merges those changes into the working directory


## Something I wrote for Grant on showing code in different branches

OK, this is a confusing thing with git to do with branches. I'm still learning myself! I had a quick chat with Maria yesterday and she said that Rasi's lab no longer even uses branches at all. I'm thinking about whether to just abandon that and have us both work in the "main" branch.

might be easiest to clear this up with a quick chat.... but I'm also kind of explaining here: it's basically because the new code is only on the janet-mapability branch - it only becomes part of the main branch when you approve my pull request and merge the branches (note that a "pull request" is different from "git pull" - that's confusing naming)

whether you're viewing through the website, or using local files, you are always viewing the files that are part of a particular branch (and you can switch between which branch you're viewing)

The easiest way to look at a chosen branch (e.g. janet-mapability) is on the web (use the dropdown menu under 'main'). I'd do that for now, if I were you. You can see how it all looks, decide if you want to merge into the main branch or not. After merging, you can do a local pull and you'll have your own copy of the file.

(only read below if you want to be confused!)

However, if you want to switch your local copy switch between branches here's how. Even in writing this, I learned some stuff and had to think really hard. It's a horribly confusing thing. But here goes:

If you switch branches locally, you do have to be careful not to get confused about which branch you're working on. The git status command always clarifies that without actually performing any actions. Then you can

before you do any actions, make sure you're on the branch you think you are (git status) then commit any changes you've made (using the usual three step protocol -  git add --all . then git commit -m "my message" then git push)
next you can sync to make sure your local copy has everything in the remote repo: git pull - if I understand right this does two things: A. get the files from remote for the current working branch and B. get metadata for all branches (file names/modification dates/sizes)
then you can switch your working copy ("checkout") so that files from the branch of interest are the ones you see in your working directory:  git checkout my_branch_name. Yes, this is a bit scary! Files may appear/disappear in your rhino folder when you do this . Don't panic - the files that disappeared will come back if you switch back to the other branch you were on before. (if you don't know the branch names you can find it either on the website, or you can use the git branch -av command to list all local and remote branches - the remote ones begin with "origin").



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

## Git and Rstudio

There's a git/SVN panel in preferences. I think you need to be working in an "Rproject" to use git (?).  Can start a new Rproject by letting Rstudio clone a repo.  Can also start a new repo for an Rproject from Rstudio.  I did NOT see a way to associate an existing Rproject with an existing repo.  I already had git installed on the mac and linked to my github account, and Rstudio somehow worked with my stored credentials (at laest on the work computer).

See docs [here](https://support.posit.co/hc/en-us/articles/200532077)

