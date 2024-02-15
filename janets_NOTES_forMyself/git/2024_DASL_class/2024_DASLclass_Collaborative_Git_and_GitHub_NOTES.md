# Collaborative Git and GitHub

Chris Lo, DASL, Feb 14 2024 (in person class)

used this site: https://replit.com/@jayoung2/CollaborativeGitGitHubDaSl

Once logged in to replit, used this script to log in to github
```
sh setup.sh 
```
It looks like this:
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