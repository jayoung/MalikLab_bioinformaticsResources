# other setup

prevent mac going to sleep?

# sidebar links

Sidebar links periodically disappear - [this](https://discussions.apple.com/thread/253196430) is what I'm trying, to prevent that:

I [mount malik_h at login](https://superuser.com/questions/125189/how-do-i-have-a-network-share-mount-in-mac-os-x-persist-between-reboots) (System Preferences - Login Items & Extensions - drag `malik_h` from Finder into the "Open at Login" list)

Also, in Mac Terminal
```
cd 
mkdir links_to_server_folders
cd links_to_server_folders

ln -s /Volumes/malik_h/ .
ln -s /Volumes/malik_h/user/jayoung/ .
ln -s /Volumes/malik_h/grp/h2av_RNAseq/  .
ln -s /Volumes/malik_h/grp/abo_janet_analysis/ .
ln -s /Volumes/malik_h/grp/worm_mitoNuc/  .
ln -s /Volumes/malik_h/user/jayoung/paml_screen/ .
```

Use Finder to make sidebar shortcuts using those links



# [Safari profiles](https://support.apple.com/en-us/105100)

most bookmarks are shared between profiles, but the 'favorites bar' can be set to be shared or not shared

can customize so that certain web domains are opened in particular profiles

# haven't installed (at least for now)

netskope client? (Hutch thing to let me get through firewalls for websites in some countries. needed it to look at ggtree website but I don't know if that's true any more)

inkscape?

MEGA

linebreak? (not sure it exists any more)

python? (python3 is installed)

haven't signed in to eDynamic dropbox


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

I think that worked. I can sync to github from VScode now, at least for this repository


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
```
export PATH="/opt/homebrew/bin:$PATH"
```

accept the Xcode license:
```
sudo xcodebuild -license accept
```

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
```
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

```

# homebrew version of perl, cpanm and perl modules 

April 2025

## install cpanm 

seems to have worked

```
cpan App::cpanminus
    # choose sudo option
which cpanm
    # /usr/local/bin/cpanm
```

## install perlbrew

to help me manage different versions of perl on the mac
```
sudo cpan App::perlbrew
perlbrew init
```

As instructed after running `perlbrew init`, I add this to my `~/.profile` file: `source ~/perl5/perlbrew/etc/bashrc` using these commands:

```
echo '##### use perlbrew rather than system perl:' >> ~/.profile
echo 'source ~/perl5/perlbrew/etc/bashrc'  >> ~/.profile
echo ''  >> ~/.profile
```
seems to have worked

## install different version


now actually install a different version of perl. The current stable version of perl is v5.40.2

### Try installing perl using perlbrew - failed

```
perlbrew install perl-5.40.2
    # and in another shell:
tail -f ~/perl5/perlbrew/build.perl-5.40.2.log
    # failed: 

    ####
Installation process failed. To spot any issues, check

  /Users/jayoung/perl5/perlbrew/build.perl-5.40.2.log

If some perl tests failed and you still want to install this distribution anyway,
do:

  (cd /Users/jayoung/perl5/perlbrew/build/perl-5.40.2/perl-5.40.2; make install)

You might also want to try upgrading patchperl before trying again:

  perlbrew install-patchperl

Generally, if you need to install a perl distribution known to have minor test
failures, do one of these commands to avoid seeing this message:

  perlbrew --notest install perl-5.40.2
  perlbrew --force install perl-5.40.2

  #### end of   /Users/jayoung/perl5/perlbrew/build.perl-5.40.2.log :

Test Summary Report
-------------------
../lib/locale.t                                                    (Wstat: 6 (Signal: ABRT) Tests: 382 Failed: 0)
  Non-zero wait status: 6
  Parse errors: No plan found in TAP output
Files=2892, Tests=1185301, 728 wallclock secs (21.95 usr  5.04 sys + 234.06 cusr 60.63 csys = 321.68 CPU)
Result: FAIL
Finished test run at Mon Apr 28 13:32:46 2025.
make: *** [test_harness] Error 1
##### Brew Failed #####

```
I didn't try to troubleshoot.



### Try installing perl using tarball

```
## try from a tarball from github/Perl - didn't work
cd /Users/jayoung/perl5/perlbrew/dists/from_github
wget https://github.com/Perl/perl5/archive/refs/tags/v5.40.2.tar.gz
perlbrew install /Users/jayoung/perl5/perlbrew/dists/from_github/v5.40.2.tar.gz
    Unable to determine perl version from archive filename.
    The archive name should look like perl-5.x.y.tar.gz or perl-5.x.y.tar.bz2 or perl-5.x.y.tar.xz

## try from a tarball from cpan
cd /Users/jayoung/perl5/perlbrew/dists/
wget https://www.cpan.org/src/5.0/perl-5.40.2.tar.gz
perlbrew install /Users/jayoung/perl5/perlbrew/dists/from_cpan/perl-5.40.2.tar.gz
    # same error as above
```

### Try installing perl using tarball without tests


```
perlbrew --notest install /Users/jayoung/perl5/perlbrew/dists/from_cpan/perl-5.40.2.tar.gz
    # perl-5.40.2 is successfully installed.
```

Actually switch which perl we run, and install matching cpanm:
```
which perl 
    # /usr/bin/perl
perlbrew switch perl-5.40.2
which perl 
    # /Users/jayoung/perl5/perlbrew/perls/perl-5.40.2/bin/perl
which cpan
    # /Users/jayoung/perl5/perlbrew/perls/perl-5.40.2/bin/cpan
cpan App::cpanminus
which cpanm
    # /Users/jayoung/perl5/perlbrew/perls/perl-5.40.2/bin/cpanm

printenv | grep 'PERL'
PERLBREW_SHELLRC_VERSION=1.01
PERLBREW_VERSION=1.01
PERLBREW_PERL=perl-5.40.2
PERLBREW_ROOT=/Users/jayoung/perl5/perlbrew
PERLBREW_HOME=/Users/jayoung/.perlbrew
PERLBREW_MANPATH=/Users/jayoung/perl5/perlbrew/perls/perl-5.40.2/man
PERLBREW_PATH=/Users/jayoung/perl5/perlbrew/bin:/Users/jayoung/perl5/perlbrew/perls/perl-5.40.2/bin
```


For VScode, make a script called `perl.perlCmd.forVScode.pl` in `/Users/jayoung`, with contents:
```
#!/bin/sh
/Users/jayoung/perl5/perlbrew/perls/perl-5.40.2/bin/perl -I/Users/jayoung/perl5/perlbrew/perls/perl-5.40.2 "$@"
## this is a script I'm adding to try to make VScode able to use perl debugger. See https://github.com/richterger/Perl-LanguageServer/issues/1
```

Also for VScode, I want `Perl::LanguageServer`

First tried this - it tries to install some dependencies too - many fail
```
cpanm Perl::LanguageServer
```

I install a couple of modules (dependencies, I think) on their own, then I retry using sudo

```
sudo cpanm --force inc::latest
sudo cpanm --force Software::License
```

retry using sudo - seems like it worked
```
sudo cpanm Perl::LanguageServer
```

I think VS code is now showing perl scripts OK.

Install a module that cpan wants:
```
sudo cpanm Log::Log4perl
```

# Bioperl


Find the name of the most recent Bioperl version:

```
cpan

d /bioperl/

Reading '/Users/jayoung/.cpan/sources/modules/03modlist.data.gz'
DONE
Writing /Users/jayoung/.cpan/Metadata
Distribution    BOZO/Fry-Lib-BioPerl-0.15.tar.gz
Distribution    CDRAUG/BioPerl-1.7.4.tar.gz
Distribution    CDRAUG/Dist-Zilla-PluginBundle-BioPerl-0.27.tar.gz
Distribution    CJFIELDS/BioPerl-1.007002.tar.gz
Distribution    CJFIELDS/BioPerl-1.6.924.tar.gz
Distribution    CJFIELDS/BioPerl-1.7.7.tar.gz
Distribution    CJFIELDS/BioPerl-1.7.8.tar.gz
Distribution    CJFIELDS/BioPerl-DB-1.006900.tar.gz
Distribution    CJFIELDS/BioPerl-Network-1.006902.tar.gz
Distribution    CJFIELDS/BioPerl-Run-1.007002.tar.gz
Distribution    CJFIELDS/BioPerl-Run-1.007003.tar.gz
Distribution    CJFIELDS/Bundle-BioPerl-2.1.9.tar.gz
12 items found

```

Actually install

First try - it got some of the way but not all:

```
sudo cpanm install CJFIELDS/BioPerl-1.7.8.tar.gz
....

Configuring XML-Writer-0.900 ... OK
Building and testing XML-Writer-0.900 ... OK
Successfully installed XML-Writer-0.900
! Installing the dependencies failed: Module 'XML::LibXML::Reader' is not installed, Module 'XML::LibXML' is not installed
! Bailing out the installation for BioPerl-1.7.8.
68 distributions installed
```

install XML::LibXML and then try Bioperl again
```
sudo cpanm XML::LibXML

Building and testing XML-LibXML-2.0210 ... FAIL
! Installing XML::LibXML failed. See /Users/jayoung/.cpanm/work/1746212689.74612/build.log for details. Retry with --force to force install it.
    # 3/18 tests failed. Try again with --force:
sudo cpanm --force XML::LibXML
    # worked
```


```
sudo cpanm install CJFIELDS/BioPerl-1.7.8.tar.gz
    # looks like it worked
```

test it's installed. It is. Perl scripts can find Bio::SeqIO on my when the first line of the script is `#!/usr/bin/env perl` but not when I use `#!/usr/bin/perl`. But I think I had the opposite problem on gizmo.
```
cd /Volumes/malik_h/user/jayoung/temp/testBioperl 
./testBioperl.pl 
    # first line is    #!/usr/bin/env perl
```
 

However, most of my scripts have `#!/usr/bin/perl` as the first line, in which case it doesn't find Bio::SeqIO
```
/Volumes/malik_h/user/jayoung/bin/seqlength.bioperl D.fasta 
    # first line is   #!/usr/bin/perl

Can't locate Bio/SeqIO.pm in @INC (you may need to install the Bio::SeqIO module) (@INC contains: /Library/Perl/5.34/darwin-thread-multi-2level /Library/Perl/5.34 /Network/Library/Perl/5.34/darwin-thread-multi-2level /Network/Library/Perl/5.34 /Library/Perl/Updates/5.34.1/darwin-thread-multi-2level /Library/Perl/Updates/5.34.1 /System/Library/Perl/5.34/darwin-thread-multi-2level /System/Library/Perl/5.34 /System/Library/Perl/Extras/5.34/darwin-thread-multi-2level /System/Library/Perl/Extras/5.34) at /Volumes/malik_h/user/jayoung/bin/seqlength.bioperl line 5.
```

On the Mac:
```
which perl 
/Users/jayoung/perl5/perlbrew/perls/perl-5.40.2/bin/perl

printenv | grep 'PERL'
PERLBREW_SHELLRC_VERSION=1.01
PERLBREW_VERSION=1.01
PERLBREW_PERL=perl-5.40.2
PERLBREW_ROOT=/Users/jayoung/perl5/perlbrew
PERLBREW_HOME=/Users/jayoung/.perlbrew
PERLBREW_MANPATH=/Users/jayoung/perl5/perlbrew/perls/perl-5.40.2/man
PERLBREW_PATH=/Users/jayoung/perl5/perlbrew/bin:/Users/jayoung/perl5/perlbrew/perls/perl-5.40.2/bin
```

On gizmo:
```
which perl
/home/jayoung/malik_lab_shared/linux_gizmo/bin/perl

PERL_MB_OPT=--install_base "/home/jayoung/malik_lab_shared/perl"
PERL_MM_OPT=INSTALL_BASE=/home/jayoung/malik_lab_shared/perl
PERL_LOCAL_LIB_ROOT=/home/jayoung/malik_lab_shared/perl
PERL5LIB=/home/jayoung/malik_lab_shared/perl:/home/jayoung/malik_lab_shared/perl/lib:/home/jayoung/malik_lab_shared/perl/lib/perl5:/home/jayoung/malik_lab_shared/perl/lib/x86_64-linux-gnu:/home/jayoung/malik_lab_shared/perl/lib/x86_64-linux-gnu/perl:/home/jayoung/malik_lab_shared/perl/lib/x86_64-linux-gnu/perl/5.26.1:/home/jayoung/malik_lab_shared/perl/lib/perl5/x86_64-linux-gnu-thread-multi:/home/jayoung/malik_lab_shared/perl/share/perl/5.26.1:/home/jayoung/malik_lab_shared/perl/share/perl:/home/jayoung/malik_lab_shared/perl/bioperl-live/lib:/home/jayoung/malik_lab_shared/perl/Bio-MolEvol/lib:/home/jayoung/malik_lab_shared/perl/Bio-PopGen/lib:/home/jayoung/malik_lab_shared/perl/Bio-Root/lib:/home/jayoung/malik_lab_shared/perl/bioperl-db/lib:/home/jayoung/malik_lab_shared/perl/Bio-EUtilities/lib:/home/jayoung/malik_lab_shared/perl/bioperl-run/lib:/home/jayoung/malik_lab_shared/perl/ensembl/modules:/home/jayoung/malik_lab_shared/perl/ensembl-io/modules:/home/jayoung/malik_lab_shared/perl/ensembl-compara/modules:/home/jayoung/malik_lab_shared/perl/ensembl-external/modules:/home/jayoung/malik_lab_shared/perl/ensembl-funcgen/modules:/home/jayoung/malik_lab_shared/perl/ensembl-functgenomics/modules:/home/jayoung/malik_lab_shared/perl/ensembl-variation/modules:/home/jayoung/malik_lab_shared/perl/ensemblgenomes-api/modules:/home/jayoung/malik_lab_shared/linux_gizmo/lib/RepeatMasker2025_Apr10/RepeatMasker:
```


xxx so in order to get the scripts I store on the server to work right I need to think about the shebang line. I think on the server `#!/usr/bin/env perl` as a first line does not work because it points to `/home/jayoung/malik_lab_shared/linux_gizmo/bin/perl` and then it can't fund my perl modules, and somehow when I use `#!/usr/bin/perl` it DOES find my modules.  Perhaps all I need to do is make `#!/usr/bin/perl` be higher up in my PATH than `/home/jayoung/malik_lab_shared/linux_gizmo/bin/perl` on the server, then I could use `env perl` on both server and mac?
