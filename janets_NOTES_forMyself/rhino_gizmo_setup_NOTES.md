
# my compute environment

## May 2024

### places I'm installing stuff

`~/bin/` contains a lot of scripts I've written

`/fh/fast/malik_h/grp/malik_lab_shared` has various subdirs. 

`linux_gizmo` has bin, lib, etc etc and is a good place to put stuff I install now.

Specific subdirs of linux_gizmo/bin are in my PATH, as is the directory itself. Where possible I will put groups of tools in a named subdir.

Places that are no longer in my PATH:
`~/malik_lab_shared/bin` has a ton of older installations.  I'm going to gradually get stuff out of there, into linux_gizmo (and subdirs of it) or re-install it.


### login files

Various files are used to set environmental variables and aliases:

### .profile

The `.profile` file will be read on login, if the user's shell is bash, unless `~/.bash_profile` or `~/.bash_login` exist

in my `.profile` file, I do this:
- source `gizmo_login_items.bash` 
- source `.bashrc` 
- add `$HOME/bin` to `PATH`
- define the `NCBI_API_KEY`

```
cd ~
cp /fh/fast/malik_h/grp/malik_lab_shared/help/loginFileTemplates/.profile .
```

### .bashrc

My `.bashrc` file is something I copied from a scicomp file, and does things like add colors to the output of `ls`. 

It also sets up aliases `ll` (`ls -l`) and `la` (`ls -a`)


### gizmo_login_items.bash

The `.profile` file is what directs my login shell to read this `gizmo_login_items.bash` file.

I put most stuff in the file called `gizmo_login_items`, in the `/fh/fast/malik_h/grp/malik_lab_shared/help/loginFileTemplates/` folder. 

It's in a shared folder, so that other people can source it too, if they want to replicate my environment.  

Before May 28, 2024, malik_h_grp had write access to that file, but after that date they only have read access. Also I used to call that file `.gizmocombinedMalikLabrc_bashVersion`.

in my home dir I have a LINK to the master copy of `gizmo_login_items.bash` - to make that link: 
```
cd ~
ln -s /fh/fast/malik_h/malik_lab_shared/help/loginFileTemplates/gizmo_login_items.bash .
```

### .bash_logout

`.bash_logout` simply clears console when I log out for privacy



## spring/summer 2020

old notes, when scicomp [switched rhino/gizmo OS over to bionic flavor of linux](https://sciwiki.fredhutch.org/compdemos/gizmo-bionic-upgrade/)

I needed to rebuild any R libraries

many module names have changed - I will need to change my alias for R, for example

I might need to reinstall other programs I have built


I think I want to have module load BioPerl as part of my login setup

### tools I tested, and they work:

blastn (etc)
paml
phyml
SRA toolkit (updated, and fixed runFastqDump.pl script for new version. seems to work faster now?)
bioperl (I am using my modules)
got rid of all old perl modules, reinstalled Bioperl and various other modules using CPAN
dotter (module load seqtools)
busco and augustus (much more recent version)

I changed perl scripts that launch sbatch jobs that load modules, because now we source a different file to get the module command working:

old way to source a module:
```
$command = "/bin/bash -c \\\"source /app/Lmod/lmod/lmod/init/bash; module load BWA; $command\\\"";
grep -l 'lmod/init/bash' ~/bin/*l
```

new way to source a module:
```
$command = "/bin/bash -c \\\"source /app/lmod/lmod/init/profile; module load BWA; $command\\\"";
```

### perl

/usr/bin/perl versus /usr/bin/env perl :

The main difference is that /usr/bin/perl is a fixed path to one perl (which is usually the one that ships with your operating system), while /usr/bin/env perl will take the first perl in $PATH

### slurm

interactive (grabnode) sessions use salloc. On the old nodes, they do not yield environmental variables that include SLURM in the name. On the new nodes they yield one: SLURM_PARTITION=campus-new


sbatch commands yield many environmental variables that include SLURM in the name, e.g. SLURM_JOB_USER=jayoung : 

cd ~/FH_fast_storage/temp
sbatch --wrap "printenv > env_sbatch_oldNodes.txt" 
cd ~/FH_fast_storage/temp
sbatch --wrap "printenv > env_sbatch_newNodes.txt" 

perhaps I want to set up a login file that does extra things when I'm in an interactive session (e.g. module load BioPerl, module load BLAST+, module load HMMER


## old notes, from when I switched to using bash

### June 4th, 2014

look at env just after I changed over to using to bash, before I have a .bashrc in my account:

using all default paths, etc.

```
jayoung@rhino01:~$ env
MODULE_VERSION_STACK=3.2.10
TERM=xterm-256color
SHELL=/bin/bash
XDG_SESSION_COOKIE=407fc484443fff6d52fe47b400000109-1401910508.382272-1440569804
SSH_CLIENT=140.107.19.110 51747 22
SSH_TTY=/dev/pts/92
USER=jayoung
LD_LIBRARY_PATH=/app/lib:/app/lib64:/usr/local/lib
MODULE_VERSION=3.2.10
MAIL=/var/mail/jayoung
PATH=/opt/moab/bin:/app/Modules/3.2.10/bin:/app/bin:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
PWD=/home/jayoung
LANG=en_US.UTF-8
MODULEPATH=/app/Modules/versions:/app/Modules/$MODULE_VERSION/modulefiles:/app/Modules/modulefiles
LOADEDMODULES=
SCRATCH_LOCAL=/loc/scratch
KRB5CCNAME=FILE:/tmp/krb5cc_31490_mMYiBS8633
SHLVL=1
HOME=/home/jayoung
LOGNAME=jayoung
SSH_CONNECTION=140.107.19.110 51747 140.107.218.55 22
MODULESHOME=/app/Modules/3.2.10
PROMPT_COMMAND=k5start -H 360
DISPLAY=localhost:28.0
SCRATCH=/mnt/rhodium/scratch
module=() {  eval `/app/Modules/$MODULE_VERSION/bin/modulecmd bash $*`
}
_=/usr/bin/env
jayoung@rhino01:~$ groups
groups: cannot find name for group ID 197870
197870 haldaemon g_jayoung g_cprodrig Tapscott_S Tapscott_S_CompBio Koelle_D_ngs Malik_H_lab tapscott_s_grp
```

### June 4th, 2014

look at env just after I changed to bash, and after I copied these three files to my account:

```
cp /etc/skel/.bash_logout .
cp /etc/skel/.profile .
cp /etc/skel/.bashrc .
```
     (.profile was necessary to have .bashrc run upon login)
```
jayoung@rhino03:~$ env
MODULE_VERSION_STACK=3.2.10
TERM=xterm-256color
SHELL=/bin/bash
XDG_SESSION_COOKIE=7def8f4c305901c3f820d9c7000002cf-1401912456.593570-338380126
SSH_CLIENT=140.107.19.110 51907 22
SSH_TTY=/dev/pts/16
USER=jayoung
LS_COLORS=rs=0:di=01;34:ln=01;36:mh=00:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=40;31;01:su=37;41:sg=30;43:ca=30;41:tw=30;42:ow=34;42:st=37;44:ex=01;32:*.tar=01;31:*.tgz=01;31:*.arj=01;31:*.taz=01;31:*.lzh=01;31:*.lzma=01;31:*.tlz=01;31:*.txz=01;31:*.zip=01;31:*.z=01;31:*.Z=01;31:*.dz=01;31:*.gz=01;31:*.lz=01;31:*.xz=01;31:*.bz2=01;31:*.bz=01;31:*.tbz=01;31:*.tbz2=01;31:*.tz=01;31:*.deb=01;31:*.rpm=01;31:*.jar=01;31:*.war=01;31:*.ear=01;31:*.sar=01;31:*.rar=01;31:*.ace=01;31:*.zoo=01;31:*.cpio=01;31:*.7z=01;31:*.rz=01;31:*.jpg=01;35:*.jpeg=01;35:*.gif=01;35:*.bmp=01;35:*.pbm=01;35:*.pgm=01;35:*.ppm=01;35:*.tga=01;35:*.xbm=01;35:*.xpm=01;35:*.tif=01;35:*.tiff=01;35:*.png=01;35:*.svg=01;35:*.svgz=01;35:*.mng=01;35:*.pcx=01;35:*.mov=01;35:*.mpg=01;35:*.mpeg=01;35:*.m2v=01;35:*.mkv=01;35:*.webm=01;35:*.ogm=01;35:*.mp4=01;35:*.m4v=01;35:*.mp4v=01;35:*.vob=01;35:*.qt=01;35:*.nuv=01;35:*.wmv=01;35:*.asf=01;35:*.rm=01;35:*.rmvb=01;35:*.flc=01;35:*.avi=01;35:*.fli=01;35:*.flv=01;35:*.gl=01;35:*.dl=01;35:*.xcf=01;35:*.xwd=01;35:*.yuv=01;35:*.cgm=01;35:*.emf=01;35:*.axv=01;35:*.anx=01;35:*.ogv=01;35:*.ogx=01;35:*.aac=00;36:*.au=00;36:*.flac=00;36:*.mid=00;36:*.midi=00;36:*.mka=00;36:*.mp3=00;36:*.mpc=00;36:*.ogg=00;36:*.ra=00;36:*.wav=00;36:*.axa=00;36:*.oga=00;36:*.spx=00;36:*.xspf=00;36:
LD_LIBRARY_PATH=/app/lib:/app/lib64:/usr/local/lib
MODULE_VERSION=3.2.10
MAIL=/var/mail/jayoung
PATH=/home/jayoung/bin:/opt/moab/bin:/app/Modules/3.2.10/bin:/app/bin:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
PWD=/home/jayoung
LANG=en_US.UTF-8
MODULEPATH=/app/Modules/versions:/app/Modules/$MODULE_VERSION/modulefiles:/app/Modules/modulefiles
LOADEDMODULES=
SCRATCH_LOCAL=/loc/scratch
KRB5CCNAME=FILE:/tmp/krb5cc_31490_cMzpK13269
SHLVL=1
HOME=/home/jayoung
LOGNAME=jayoung
SSH_CONNECTION=140.107.19.110 51907 140.107.218.57 22
MODULESHOME=/app/Modules/3.2.10
LESSOPEN=| /usr/bin/lesspipe %s
PROMPT_COMMAND=k5start -H 360
DISPLAY=localhost:46.0
LESSCLOSE=/usr/bin/lesspipe %s %s
SCRATCH=/mnt/rhodium/scratch
module=() {  eval `/app/Modules/$MODULE_VERSION/bin/modulecmd bash $*`
}
_=/usr/bin/env
```

### June 4th, 2014

look at env after running `.gizmocombinedMalikLabrc_bashVersion` on top of the default login files:

```
MODULE_VERSION_STACK=3.2.10
MANPATH=/home/jayoung/malik_lab_shared/linux_gizmo/man:/home/jayoung/malik_lab_shared/linux_gizmo/bin/qt/doc/man
TERM=xterm-256color
SHELL=/bin/bash
WISECONFIGDIR=/home/jayoung/malik_lab_shared/lib/wisecfg
XDG_SESSION_COOKIE=7def8f4c305901c3f820d9c7000002cf-1401912456.593570-338380126
SSH_CLIENT=140.107.19.110 51907 22
PERL5LIB=/home/jayoung/malik_lab_shared/perl/bioperl-live:/home/jayoung/malik_lab_shared/perl/Bio-Root/lib:/home/jayoung/malik_lab_shared/perl/Bio-EUtilities/lib:/home/jayoung/malik_lab_shared/perl/ensembl/modules:/home/jayoung/malik_lab_shared/perl/ensembl-compara/modules:/home/jayoung/malik_lab_shared/perl/ensembl-external/modules:/home/jayoung/malik_lab_shared/perl/ensembl-functgenomics/modules:/home/jayoung/malik_lab_shared/perl/ensembl-variation/modules:/home/jayoung/malik_lab_shared/perl/bioperl-run/lib:/home/jayoung/malik_lab_shared/perl/lib:/home/jayoung/malik_lab_shared/lib/perl5:/home/jayoung/malik_lab_shared/perl/lib/perl5/x86_64-linux-gnu-thread-multi:/home/jayoung/malik_lab/public_databases/Ensembl/VariantEffectPredictorCache/plugins:/home/jayoung/malik_lab_shared/lib/RepeatMasker2014_Feb5_linux_64/RepeatMasker
QTDIR=/home/jayoung/malik_lab_shared/linux_gizmo/bin/qt
LAGAN_DIR=/home/jayoung/malik_lab_shared/linux_gizmo/bin/lagan
SSH_TTY=/dev/pts/16
CONSED_PARAMETERS=/home/jayoung/malik_lab_shared/lib/.consedrc.linux64
PYTHONUSERBASE=/home/jayoung/malik_lab_shared
USER=jayoung
PHRED_PARAMETER_FILE=/home/jayoung/malik_lab_shared/lib/phredpar.dat
LS_COLORS=rs=0:di=01;34:ln=01;36:mh=00:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=40;31;01:su=37;41:sg=30;43:ca=30;41:tw=30;42:ow=34;42:st=37;44:ex=01;32:*.tar=01;31:*.tgz=01;31:*.arj=01;31:*.taz=01;31:*.lzh=01;31:*.lzma=01;31:*.tlz=01;31:*.txz=01;31:*.zip=01;31:*.z=01;31:*.Z=01;31:*.dz=01;31:*.gz=01;31:*.lz=01;31:*.xz=01;31:*.bz2=01;31:*.bz=01;31:*.tbz=01;31:*.tbz2=01;31:*.tz=01;31:*.deb=01;31:*.rpm=01;31:*.jar=01;31:*.war=01;31:*.ear=01;31:*.sar=01;31:*.rar=01;31:*.ace=01;31:*.zoo=01;31:*.cpio=01;31:*.7z=01;31:*.rz=01;31:*.jpg=01;35:*.jpeg=01;35:*.gif=01;35:*.bmp=01;35:*.pbm=01;35:*.pgm=01;35:*.ppm=01;35:*.tga=01;35:*.xbm=01;35:*.xpm=01;35:*.tif=01;35:*.tiff=01;35:*.png=01;35:*.svg=01;35:*.svgz=01;35:*.mng=01;35:*.pcx=01;35:*.mov=01;35:*.mpg=01;35:*.mpeg=01;35:*.m2v=01;35:*.mkv=01;35:*.webm=01;35:*.ogm=01;35:*.mp4=01;35:*.m4v=01;35:*.mp4v=01;35:*.vob=01;35:*.qt=01;35:*.nuv=01;35:*.wmv=01;35:*.asf=01;35:*.rm=01;35:*.rmvb=01;35:*.flc=01;35:*.avi=01;35:*.fli=01;35:*.flv=01;35:*.gl=01;35:*.dl=01;35:*.xcf=01;35:*.xwd=01;35:*.yuv=01;35:*.cgm=01;35:*.emf=01;35:*.axv=01;35:*.anx=01;35:*.ogv=01;35:*.ogx=01;35:*.aac=00;36:*.au=00;36:*.flac=00;36:*.mid=00;36:*.midi=00;36:*.mka=00;36:*.mp3=00;36:*.mpc=00;36:*.ogg=00;36:*.ra=00;36:*.wav=00;36:*.axa=00;36:*.oga=00;36:*.spx=00;36:*.xspf=00;36:
LD_LIBRARY_PATH=/app/lib:/app/lib64:/usr/local/lib:/lib:/home/jayoung/malik_lab_shared/linux_gizmo/lib:/home/jayoung/malik_lab_shared/lib:/home/jayoung/malik_lab_shared/lib/root:/home/jayoung/malik_lab_shared/lib_linux:/usr/lib64:/home/jayoung/malik_lab_shared/lib_linux_64:/home/jayoung/malik_lab_shared/linux_gizmo/bin/FuzzyClustering/Matlab_Compiler_Runtime/v79/runtime/glnxa64:/home/jayoung/malik_lab_shared/linux_gizmo/bin/FuzzyClustering/Matlab_Compiler_Runtime/v79/extern/include:/home/jayoung/malik_lab_shared/linux_gizmo/local/boost:/home/jayoung/malik_lab_shared/linux_gizmo/local/boost/bin.v2/libs:/home/jayoung/malik_lab_shared/linux_gizmo/bin/qt/lib
WUBLASTMAT=/home/jayoung/malik_lab_shared/linux_gizmo/share/blast/wu_blast/matrix
MODULE_VERSION=3.2.10
MAIL=/var/mail/jayoung
PATH=.:/usr/sbin:/home/jayoung/bin:/home/jayoung/malik_lab_shared/linux_gizmo/bin:/home/jayoung/malik_lab_shared/bin:/home/jayoung/malik_lab_shared/scripts:/home/jayoung/malik_lab_shared/linux_gizmo/bin/VelvetOptimiser-2.2.5:/home/jayoung/malik_lab_shared/linux_gizmo/bin/consed/bin:/home/jayoung/malik_lab_shared/linux_gizmo/bin/SRA_toolkit:/home/jayoung/malik_lab_shared/linux_gizmo/bin/FuzzyClustering/bin:/home/jayoung/malik_lab_shared/linux_gizmo/bin/meme/bin:/home/jayoung/malik_lab_shared/linux_gizmo/bin/weeder:/home/jayoung/malik_lab_shared/linux_gizmo/bin/lagan:/home/jayoung/malik_lab_shared/linux_gizmo/bin/lagan/utils:/home/jayoung/malik_lab_shared/linux_gizmo/bin/phyme/code/bin:/home/jayoung/malik_lab_shared/linux_gizmo/bin/MUMmer3.23:/home/jayoung/malik_lab_shared/linux_gizmo/bin/picard-tools:/home/jayoung/malik_lab_shared/linux_gizmo/bin/qt/bin:/home/jayoung/bin:/opt/moab/bin:/app/Modules/3.2.10/bin:/app/bin:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
LD_RUN_PATH=/home/jayoung/malik_lab_shared/lib
PWD=/home/jayoung
LANG=en_US.UTF-8
MODULEPATH=/app/Modules/versions:/app/Modules/$MODULE_VERSION/modulefiles:/app/Modules/modulefiles
LOADEDMODULES=
SCRATCH_LOCAL=/loc/scratch
ROOTSYS=/home/jayoung/malik_lab_shared
KRB5CCNAME=FILE:/tmp/krb5cc_31490_cMzpK13269
WUBLASTFILTER=/home/jayoung/malik_lab_shared/linux_gizmo/share/blast/wu_blast/filter
EMBOSS_ACDROOT=/home/jayoung/malik_lab_shared/linux_gizmo/share/EMBOSS/acd
SHLVL=1
HOME=/home/jayoung
CONSED_HOME=/home/jayoung/malik_lab_shared/linux_gizmo/bin/consed
EMBOSS_DATA=/home/jayoung/malik_lab_shared/share/linux_gizmo/EMBOSS/data
BOOST_LIBRARYDIR=/home/jayoung/malik_lab_shared/linux_gizmo/lib
LOGNAME=jayoung
BOOST_ROOT=/home/jayoung/malik_lab_shared/linux_gizmo/include/boost
SSH_CONNECTION=140.107.19.110 51907 140.107.218.57 22
MODULESHOME=/app/Modules/3.2.10
PKG_CONFIG_PATH=/home/jayoung/malik_lab_shared/linux_gizmo/lib/pkgconfig
LESSOPEN=| /usr/bin/lesspipe %s
R_HOME=/home/jayoung/malik_lab_shared/linux_gizmo/lib/R
PROMPT_COMMAND=k5start -H 360
DISPLAY=localhost:46.0
R_LIBS_USER=/home/jayoung/malik_lab_shared/linux_gizmo/lib/R/library
LESSCLOSE=/usr/bin/lesspipe %s %s
SCRATCH=/mnt/rhodium/scratch
module=() {  eval `/app/Modules/$MODULE_VERSION/bin/modulecmd bash $*`
}
_=/usr/bin/env
```
