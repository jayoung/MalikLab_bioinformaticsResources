# set up password-less ssh access to rhino

On the Mac, Terminal, in my home dir, do this:
```
cd
ssh-keygen -t ed25519 -b 4096
    # that creates id_ed25519 and id_ed25519.pub in ~/.ssh
    # do not require a passcode
chmod 400 ~/.ssh/id_ed25519
export USER_AT_HOST="jayoung@rhino03.fhcrc.org"
export PUBKEYPATH="$HOME/.ssh/id_ed25519.pub"
ssh-copy-id -i "$PUBKEYPATH" "$USER_AT_HOST"
     # and entered my password
```

# configure git for automatic sign-in

```
cd
git config --global user.email "jayoung@fredhutch.com"
git config --global user.name "Janet Young"

git config --global core.fileMode false

git credential-store --file ~/.git-credentials store 
    # appears to hang, but enter this:
protocol=https
host=github.com
username=jayoung
password=my_personal_access_token_string

git config --global credential.helper store
```



# installed these applications, and I think I set them up:

```
seaview
dendroscope
geneious - haven't activated free trial yet
igv
slack
zoom
treeviewer
chrome
xquartz
dropbox
xcode
docker
jalview
set up druva backups
vscode
endnote
pymol (no license!)
```

# changed default shell to bash:
- system preferences - users and groups - control-click over my user name, advanced options, change login shell to /bin/bash


# java
see https://osxdaily.com/2024/06/03/how-install-java-mac-m3-m2-m1-apple-silicon/
I've added this to .profile:
export JAVA_HOME=/usr/libexec/java_home


# command-line tools
Terminal:
xcode-select --install


# homebrew (requires command-line tools)
used pkg installer from https://brew.sh
added this to .profile:
export PATH="/opt/homebrew/bin:$PATH"

accept the Xcode license:
sudo xcodebuild -license accept


# wget
from Terminal: 
brew install wget


# figtree - had a problem running it at first, but I think I fixed it. Notes:

Problem:
ERROR launching 'FigTree v1.4.4'.
No suitable Java version found on your system!
This program requires Java  6 or later.
Make sure you install the required Java version

Solution posted here - https://github.com/tofi86/universalJavaApplicationStub/releases/
cd /Applications/FigTree/FigTree\ v1.4.4.app/Contents/MacOS 
mkdir old
mv universalJavaApplicationStub old

wget https://github.com/tofi86/universalJavaApplicationStub/releases/download/v3.3.0/universalJavaApplicationStub-v3.3.0-binary-macos-10.15.zip 
unzip universalJavaApplicationStub-v3.3.0-binary-macos-10.15.zip 
rm universalJavaApplicationStub-v3.3.0-binary-macos-10.15.zip 

## r / rstudio

# regular CRAN packages:
install.packages("tidyverse")
install.packages("janitor")
install.packages("here")
install.packages("ape")
install.packages("cowplot")
install.packages("devtools")
install.packages("openxlsx")
install.packages("kableExtra")
install.packages("palmerpenguins")
install.packages("snakecase")

# bioconductor:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.20")

# core bioC packages (it installs many dependencies of these two):
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

# more packages
BiocManager::install("DiffLogo")
BiocManager::install("EnhancedVolcano")
BiocManager::install("ggmsa")
BiocManager::install("ggseqlogo")
BiocManager::install("ggtreeExtra")
BiocManager::install("MotifDb")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("plyranges")
BiocManager::install("seqLogo")
BiocManager::install("taxize")



# not installing (at least for now)
linebreak? (not sure it exists any more)

netskope client? (Hutch thing to let me get through firewalls for websites in some countries. needed it to look at ggtree website but I don't know if that's true any more)

inkscape?

