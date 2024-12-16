# VScode extension "Perl" by Gerald Richter. 

Should allow debugging.

There's a lot of information on the extensions panel, to help installation

## Setting up rhino access for VScode

I ran these commands on both my desktop and laptop macs:
```
ssh-keygen -t ed25519 -b 4096
    # that creates id_ed25519 and id_ed25519.pub in ~/.ssh
chmod 400 ~/.ssh/id_ed25519
export USER_AT_HOST="jayoung@rhino03.fhcrc.org"
export PUBKEYPATH="$HOME/.ssh/id_ed25519.pub"
ssh-copy-id -i "$PUBKEYPATH" "$USER_AT_HOST"
     # and entered my password
```
that means I no longer have to enter my password when I do ssh 'jayoung@rhino03.fhcrc.org' (or my alias rhino3), and when I connect from VScode.

Now I can run VScode from within a rhino ssh session. From VScode: command-shift-P brings up the command palette, and I do "Remote-SSH: Connect to Host" and add rhino to the list (ssh jayoung@rhino03.fhcrc.org). Now I can connect directly for each VScode window.

VScode extensions run on the mac or on rhino, so I have to install them in both places.

## VScode settings

sync is ON on both my laptop and my desktop mac (and from a window running through ssh-rhino).  I think that means the json file etc gets synced, e.g. for specifying location of python

in a rhino-based VScode session on the desktop computer (no folder open) I messed around and created a 'profile'. Don't think it is syncing to anything? Not sure what that did

Within each git repo, I do NOT want to sync the .vscode folder that gets put there (because it is being synced in a different way, through VScode application itself).  Make sure that `vscode` is in the `.gitignore` file, and if I accidentally sync it, here are the steps to remove it from the repo:
```
git add --all .
git rm -r .vscode
git commit -m "remove .vscode"
git push
```

## Python
on my mac (desktop and laptop) I have Python 3.9.6

on rhino - maybe I need to uninstall this, as it's not working:
which python3
/home/jayoung/malik_lab_shared/linux_gizmo/bin/python3
python3
python3: error while loading shared libraries: libpython3.7m.so.1.0: cannot open shared object file: No such file or directory


on rhino there are many modules, e.g. fhPython/3.9.6-foss-2021b 

but running module load fhPython/3.9.6-foss-2021b 


### to run python within vscode, running on my mac:

I find code-runner's settings.json file and change 

change this
"python": "python -u",
    or something similar so that it looks like this
"python": "python3",

change this

"python": "python -u",
or something similar so that it looks like this

"python": "python3",


## testing a bash chunk

I get "Code language not supported or defined." when I try to run this with the little 'play' button
```bash
ls
```
Same here. But putting this line within a python script does work (in a Mac-based VScode window)
```python
print("Hello World")
```


## Perl (2022 notes)

I have the perl debugger working (I think?) on my work mac, at least in the bat_reproduction VScode workspace, but probably more generally. 

It is currently working using perl modules installed on the mac.  There IS also a way to get it going using ssh so that I could test on rhino using the same perl modules installed there, but I'm not messing with that.

Notes:
- I installed Perl::LanguageServer (`sudo cpanm --force Perl::LanguageServer`, on the work mac)
- I edited some stuff in the global settings.json file (`/Users/jayoung/Library/Application\ Support/Code/User/settings.json`):
    - this is in the top level of the json: 
```
    "perl.perlCmd": "/Users/jayoung/perl.perlCmd.forVScode.pl",
    "perl.perlInc": [ 
        "/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/lib/5.34.0", 
        "/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/lib/site_perl/5.34.0", 
        "/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/lib/site_perl/5.34.0/darwin-2level", 
        "/Users/jayoung/perl5/lib/perl5" ]
```
   - and so is this:
```
    "launch": {
        "version": "0.2.0",
        "configurations": [
            {
                "type": "perl",
                "request": "launch",
                "name": "Perl-Debug local",
                "program": "${workspaceFolder}/${relativeFile}",
                "stopOnEntry": true,
                "reloadModules": true,
                "args": [],
                "env": {"PERL5LIB":[
                    "/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/lib/5.34.0", 
                    "/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/lib/site_perl/5.34.0", 
                    "/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/lib/site_perl/5.34.0/darwin-2level", 
                    "/Users/jayoung/perl5/lib/perl5"]},
            }
        ]
    },
```

- I made a small script (`/Users/jayoung/perl.perlCmd.forVScode.pl` that looks like this:
```
#!/bin/sh
/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/bin/perl -I/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0 "$@"
## this is a script I'm adding to try to make VScode able to use perl debugger. See https://github.com/richterger/Perl-LanguageServer/issues/1
```

### general perl stuff on my work mac:
```
which perl
/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/bin/perl

which cpanm
/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/bin/cpanm

printenv | grep 'PERL'
PERLBREW_SHELLRC_VERSION=0.94
PERLBREW_VERSION=0.94
PERLBREW_PERL=perl-5.34.0
PERLBREW_ROOT=/Users/jayoung/perl5/perlbrew
PERLBREW_HOME=/Users/jayoung/.perlbrew
PERLBREW_MANPATH=/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/man
PERLBREW_PATH=/Users/jayoung/perl5/perlbrew/bin:/Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/bin

```

cpanm installs stuff to /Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/lib

Jan 30, 2023:
```
sudo cpanm --force Perl::LanguageServer
   # this installed in /Users/jayoung/perl5/perlbrew/perls/perl-5.34.0/lib/site_perl/5.34.0/Perl/
```

For all `cpanm` installs on my work mac, I think I should use `sudo cpanm --force`
