# VScode extension "Perl" by Gerald Richter. 

Should allow debugging.

There's a lot of information on the extensions panel, to help installation

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



## general perl stuff on my work mac:
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
