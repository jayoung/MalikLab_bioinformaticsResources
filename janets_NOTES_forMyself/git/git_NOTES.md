# Basic notes

Check whether the local repository matches the remote one: `git status`

`git reflog` – shows a list of commits

To list git configuration: `git config –list`

To change preferred editor (for all projects):  `git config --global core.editor "vim"`

To show which files were added during a particular commit:
`git show --pretty="" --name-only d09fb87` (where d09fb87 is the commit ID)

# removing a cached file:

`git rm --cached otherDataFiles/kyleFowler/GEO/suppData/GSE49977_Rec12_multimap.txt`

`git rm --cached otherDataFiles/kyleFowler/GEO/suppData/GSE49977_Rec12_unique.txt`

`git filter-branch -f --index-filter 'git rm -r --cached --ignore-unmatch otherDataFiles/kyleFowler/GEO/suppData' HEAD`

A more general way to clear the cache and get only the right files on the github server:
1.	Make sure all changes are committed
2.	git rm -rf --cached .
3.	git add --all .
4.	git commit
5.	git push

# removing a big file from commit

Common situation:  I did `git add`, I did `git commit` and I got halfway through `git push` before I got a message that I have file(s) that are too big to store on github.  So I need to add that file to `.gitignore`, undo the commit and the add, and redo the commit.

1. add file to .gitignore
2. git reset HEAD~   
3. git add --all .
4. git commit -m "my commit message"
5. git push

# branches

often a software development approach is to use one branch per issue, and merge that branch into main once the issue is complete/fixed.

‘origin/main’ is the main branch, the version stored at github.com
‘main’ is the main branch, stored on the local computer.

We can make new branches, for example if the code is working well, but we want to experiment a bit, we might make a branch.  We can ‘checkout’ specific branches for use.

Later we might then merge our branch with the main branch
Or, we might rebase  - somewhat different from merge.
Merge should take the commits from main and the commits from the branch, and intersperse them.  Rebase might

`git reset hard` can help with deleted things.

# restoring a deleted file:

To restore a deleted file, if I haven’t done a commit afterwards
git checkout HEAD <filename>

tools to examine gitignore
To check why git is ignoring a folder called knownDomesticatedGenes:
git check-ignore -v -- knownDomesticatedGenes

To see every file git is ignoring:
git ls-files -o -i --exclude-standard

# Large files
Github doesn’t want to accept large files. Ideally I include large files/dirs in my .gitignore file before I try to push changes to github. A couple of times I’ve forgotten to do that before I do ‘git push’ and I get an error. Not only do I get an error on that ‘git push’ attempt, but even if I add the file to .gitignore and then try again, it persists in trying to sync the large file. Some combination of these commands seemed to fix it (not sure which ones were key):


# file mode Nov 19 2021
I had a weird issue using VScode where it thought all of my files had uncommitted changes, but it was something to do with permissions and the SMB mount.  I fixed it by issuing two commands from a gizmo command line that prevent git from caring about syncing permissions:

This command changed global options:
git config --global core.fileMode false
it added this line to ~/.gitconfig in the [core] section
		filemode = false

This command changed options within a particular repo (in this case HSV1_Mx_evolution):
git config core.fileMode false
it added this line to .git/config in the [core] section
		fileMode = false

Seemed to fix my weird VScode problem

I also did this on my 16” mac laptop

# OTHER RESOURCES

http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004668
https://swcarpentry.github.io/git-novice/
https://jules32.github.io/2016-07-12-Oxford/git/

# Jan 21 2016, presentation by Paul Litwin of Collaborative Data Services, FHCRC

for help, try: http://git-scm.com/doc

centralized version control (one master copy on a server) versus distributed (various versions of the repository are all peers of each other, and the software tracks most recent changes)

git = the software that manages repositories
github = a website that hosts git-style repositories. Github is free for public repositories, a charge if you want to have private repositories, but FHCRC has an account (ask Dirk). There are also other hosts available, e.g. bitbucket (get 5 private repositories for free)

Once git is installed whichever machine, it can be configured to remember things like your user name, password, preferred text editor (default might be vim - command to exit vim is :q), etc.  Can also install "credential helper" to store name/password for remote servers.

We keep each PROJECT in its own folder. To start a git repository for that project, from the top level of that dir, we do:
   git init
This creates a hidden folder called .git where backups, change logs, etc will be stored

Any files that are in the folder or its subdirs will be backed up and changes will be tracked

Files get backed up in two steps (a) files are staged (or "tracked") (b) changes are committed. Files that have not yet been staged might be known as "untracked" Changes that have been committed to a local repository can then also be backed up to a remote repository copy, or can be "checked out" to copy them back to the local copy. Only committed changes are tracked, etc. When you "commit" that captures the project as a whole.

Some useful commands:
git status

git add (= stage files to temporary staging area - this does NOT commit the changes)
git add . (= stage all files)  (git add --all is the same thing)

git commit  (= commit changes - it will ask for a message to describe that update, probably in an awkward-to-use text editor. If it's vi, can save and quit by doing escape-:w then escape-:q)
git commit -m "my message describing this update"

git log
git log --oneline
git log -p (profile)

git reset

git diff
git diff --cached
git diff --staged

Each commit will have an ID - a checksum number. The most recent commit is also called HEAD. The short version of this ID can be used in various ways, e.g. in git diff or in git reset.

Can make a file called .gitignore if there are files within the project you do NOT want to back up.

Remote repositories: Method A
1. start up a local repository
2. use the github website to make a new repo - it will provide a URL and some example commands to commit to that repository.
3. from the local machine: 
git remote add origin URL
("origin" is convention for the current project)
git push origin master -u

can also edit files directly on the github website

Remote repositories: Method B
Start the repository remotely, and then clone it locally
git clone URL

To interact with the remote repository:
git push
git pull
(git add and git commit work only in the local session - they do not compare the local repo with the remote version)

When you do git pull/push, git makes an attempt to resolve conflicting changes that have been made to different versions of the repository (or at least tries to report them)

Version control: branching

e.g. if you want to maintain a "release" version of the software, and will need to do bug fixes on that, but also want to work on a "version2" that is quite different, and is not ready for public release for a while.  "release" and "version2" are different "branches" of the repository. Branches can be merged/deleted, etc.

"stashing" - a way to temporarily save changes without committing them.

"forking" - creates a personal copy of a public repository that you can work on, and then later if you have bug fixes to contribute  you can "create a pull request" so that the authors can try to pull in your fixes to their current version.

# Credentials (passwords, access tokens)

I was initially using a password to access my github repos but they now want me to use an ‘access token’ (similar to a password).

I first used the github website to set up an access token.

Then I store that access token on each computer I want to access github from.  
https://docs.github.com/en/github/getting-started-with-github/getting-started-with-git/caching-your-github-credentials-in-git

I run this command:
git config credential.helper store
then I do a 
git push
I give it my user name (jayoung) and the access token, and the credentials will be remembered for next time, as it creates a file called ~/.git-credentials

Or, a method from Jenny Smith on Slack to save PAT on the Gizmos:
git credential-store --file ~/.git-credentials store 
This will appear to be "hanging" and won't prompt you for any information. But on the command line, you can enter this:
protocol=https
host=github.com
username=mygithub_username
password=the_personal_access_token_string
This should create a file called ~/.git-credentials stored in your home drive. There is a bit more detail on the git documentation but I didn't think it really described it very well (I had to ask SciComp for help). https://git-scm.com/book/en/v2/Git-Tools-Credential-Storage (edited)

Turned out it was using the stored credentials for one repo but not another. The repo where it DOES use stored credentials has a couple of lines in the .git/config file:  
[credential]
	   helper = store
Try this in the repo where it is not using stored credentials:
git config --global credential.helper store
then I do a 
git push

## Personal access token (new credential method, 2021)

As of Aug 2021 now I have to use a personal access token instead of a password. It’s too long to remember so I have to store it. I tried something suggested at this site but it didn’t seem to work:
https://docs.github.com/en/get-started/getting-started-with-git/caching-your-github-credentials-in-git
https://cli.github.com/manual/gh_auth_login 

1.	from the command line, “gh auth login”
2.	choose GitHub.com
3.	choose HTTPS
4.	say Y to ‘authenticate git with my credentials’
5.	paste my authentication token


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

