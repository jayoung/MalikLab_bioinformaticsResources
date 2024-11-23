# Collaborative Git and GitHub

Chris Lo, DASL, Feb 14 2024 (in person class)

## replit compute environment

we used this site as a compute environment: https://replit.com/@jayoung2/CollaborativeGitGitHubDaSl

Once logged in to replit, used this script to log in to github
```
sh setup.sh 
```
The setup.sh script looks like this:
```
git config --global user.name "Replit_User"
git config --global user.email "address@email.com"
echo "Link for authentication: https://github.com/login/device"
gh auth login
```

Then can do 
```
git clone https://github.com/fhdsl/Collaborative_Git_GitHub_Student_Practice.git
```

## class material

See [slides](2024_DASLclass_Collaborative_Git_and_GitHub.pdf)

We learned about branching and merging, and pull requests

Create a new branch:
- can create using web interface
- convention is to use the goal of the branch to name it (rather than a person's name)
- it is often easiest to only have one person working on each branch
- can have sub-branches of branches

Switch local copy to a different branch:
- `git branch` shows local branches
- `git branch -rv` shows remote branches (v is for verbose, to include most recent commit message)
- `git branch -av` shows local and remote branches (remote branch names begin with `origin/`)
- `git checkout my_branch_name` switches to the named branch
- `git checkout main` switches back to the main branch
- the branch we're currently working on is known as the `HEAD` branch

Seeing branch structure/commit history:
- here's a way to show a crude graph of the branching structure as viewed from our current working branch/local directory: `git log --all --decorate --oneline --graph`
- on the website, we can see a branch structure diagram under 'Insights-Network'.  In some views the arrow points to the parent commit.  With a regular used account this works only for public repos. With 'github pro' this is also available for private repos.

Updating local copy to reflect the remote:
- [`git pull`](https://www.freecodecamp.org/news/git-pull-explained/) does two things:
  - updates the current local working branch.  It fetches (`git fetch`) the new commits and merges (`git merge`) these into your local branch.
  - updates the remote tracking branches for all other branches.


Adding/committing
- I have been doing `git add --all .`
- `git add myFile.txt` is another option for specific files
- `git commit -m "my message for this commit"` would shortcut the vim message adding step

Pull requests (merging branches)
- use the github website. there are various features to choose, including selecting assignees who can approve the pull request
- can (and should) make additional comments during the pull request, to explain what this update does
- when ready, somebody approves the pull request and merges the branches
- can choose to delete the branch or not
- it's good to keep pull requests small and modular

Even after we delete a branch on the remote server, it still exists locally. To delete it locally, we do `git branch --delete my_branch_name`

`git pull` sort of combines two commands: `git fetch` then `git merge`. Doing them separately is a bit safer - if there are changes in the remote repo that conflict with/overwrite changes that have occurred locally, you might lose the local versions.

`git fetch origin` is similar to `git pull` but more subtle.  Not sure I totally understand it, but it looks like `git fetch` (I think by default it does `git fetch origin`) gets the metadata, i.e. tells us what's changed on the remote repo versus the local copy, but it doesn't actually change the files

A **fork** is NOT the same as a **branch** - to fork is to create a new copy of the entire repo in a different user's github account. Sort of a much more extreme version of a branch (think of a fork in the road)


Ways to interact with git:
- command line
- VScode (command line provides finer-grained control)
- Rstudio
- github desktop tool (shows local copy of the repo)

Other things to learn about
- github [wikis](https://docs.github.com/en/communities/documenting-your-project-with-wikis/about-wikis) for hosting documentation

