# Installations

## rhino/gizmo server

see /fh/fast/malik_h/user/jayoung/source_codes/pythonModules/pythonModules_INSTALLATION_NOTES_JY.txt

## Mac desktop (iMac 2017) (installations done Dec 2023)

First I install the latest version of python, using `brew`. This installs a SECOND copy of python3 in a non-default location, along with `pip` and `setuptools`:
```
brew install python
```
Didn't work the first time. Got a message advising me to install Apple Command Line Tools  - I did it, and ittook quite a while
```
xcode-select --install
```

Then try this agin: 
```
brew install python
```
Got much further, but I got this error:
```
==> Pouring python@3.11--3.11.6_1.ventura.bottle.tar.gz
Error: An unexpected error occurred during the `brew link` step
The formula built, but is not symlinked into /usr/local
Permission denied @ dir_s_mkdir - /usr/local/Frameworks
Error: Permission denied @ dir_s_mkdir - /usr/local/Frameworks
```
Sort that out as follows:
```
sudo mkdir /usr/local/Frameworks
sudo chown $(whoami):admin /usr/local/Frameworks
brew link python@3.11
```


Python locations and links aren't quite as I'd like - try this:
```
brew uninstall python
brew install python
# messages include

Unversioned symlinks `python`, `python-config`, `pip` etc. pointing to
`python3`, `python3-config`, `pip3` etc., respectively, have been installed into
  /usr/local/opt/python@3.11/libexec/bin

You can install Python packages with
  pip3 install <package>
They will install into the site-package directory
  /usr/local/lib/python3.11/site-packages

==> `brew cleanup` has not been run in the last 30 days, running now...
Disable this behaviour by setting HOMEBREW_NO_INSTALL_CLEANUP.
Hide these hints with HOMEBREW_NO_ENV_HINTS (see `man brew`).
Removing: /usr/local/Cellar/ca-certificates/2021-10-26... (3 files, 208.5KB)
Removing: /usr/local/Cellar/openssl@3/3.0.1... (6,420 files, 28.1MB)
Error: Permission denied @ apply2files - /usr/local/share/ghostscript/9.16/Resource/CIDFont/HiraKakuPro-W3
```

Try 
```
brew doctor
    # got many messages. But I'll do this to try to fix the brew cleanup error above
sudo rm /usr/local/share/ghostscript/9.16/Resource/CIDFont/Hira*
brew cleanup
    # then
brew uninstall python
brew install python
    # now no errors
```


I want to add `/usr/local/opt/python@3.11/libexec/bin` to my PATH so that I can use `pip` to mean `pip3` and `python` to mean `python3` (in each case, this new version I just installed). So I modify `.profile` to add the following line:
```
export PATH="/usr/local/opt/python@3.11/libexec/bin:$PATH"
```

Now pip and python are both in my PATH
```
which python
    # /usr/local/opt/python@3.11/libexec/bin/python
which pip
    # /usr/local/opt/python@3.11/libexec/bin/pip

which python3
    # /usr/local/bin/python3
ls -l /usr/local/bin/python3
lrwxr-xr-x  1 jayoung  admin  42 Dec 29 14:26 /usr/local/bin/python3 -> ../Cellar/python@3.11/3.11.6_1/bin/python3
```

I think this is what I'd like the header of my phython scripts to have, on my Mac laptop at least:
`/usr/local/bin/python3`

Now I can install Biopython:
```
pip install biopython
   # Installing collected packages: numpy, biopython
   # Successfully installed biopython-1.82 numpy-1.26.2
```

where did it get installed?
```
# show which packages are installed and where they are:
python3 -m pip list -v
```

packages are here: `/usr/local/lib/python3.11/site-packages`

So I want to add that to PYTHONPATH - I add this to my .profile file:
```
#### JY adding to (or creating) PYTHONPATH:
if [ -z ${PYTHONPATH+x} ] 
then
    echo "        making new PYTHONPATH"
    export PYTHONPATH=/usr/local/lib/python3.11/site-packages
else
    echo "        adding to PYTHONPATH"
    export PYTHONPATH=/usr/local/lib/python3.11/site-packages:${PYTHONPATH}
fi
```

More useful packages:
```
pip install ipykernel # for python notebooks
pip install notebook # for python notebooks
pip install pandas # to data tables
pip install matplotlib # for graphs
```




## Mac laptop (newer OS) (installations done Dec 2023)

My mac laptop did not have pip on it, so I do this;

First I install the latest version of python, using `brew`. This installs a SECOND copy of python3 in a non-default location, along with `pip` and `setuptools`:
```
brew install python
```
Messages include:
```
==> python@3.11
Python has been installed as
  /usr/local/bin/python3

Unversioned symlinks `python`, `python-config`, `pip` etc. pointing to
`python3`, `python3-config`, `pip3` etc., respectively, have been installed into
  /usr/local/opt/python@3.11/libexec/bin

You can install Python packages with
  pip3 install <package>
They will install into the site-package directory
  /usr/local/lib/python3.11/site-packages
```

I did this (not sure if it actually helped) to make symlinks again:
```
brew unlink python && brew link python
```

I want to add `/usr/local/opt/python@3.11/libexec/bin` to my PATH so that I can use `pip` to mean `pip3` and `python` to mean `python3` (in each case, this new version I just installed). So I modify `.profile` to add the following line:
```
export PATH="/usr/local/opt/python@3.11/libexec/bin:$PATH"
```

Now pip and python are both in my PATH
```
which python
    # /usr/local/opt/python@3.11/libexec/bin/python
which pip
    # /usr/local/opt/python@3.11/libexec/bin/pip
```

```
which python3
    # /usr/local/bin/python3
ls -l /usr/local/bin/python3
lrwxr-xr-x  1 jayoung  admin  42 Dec 28 15:31 /usr/local/bin/python3 -> ../Cellar/python@3.11/3.11.6_1/bin/python3
```

I think this is what I'd like the header of my phython scripts to have, on my Mac laptop at least:
`/usr/local/bin/python3`

Now I can install Biopython:
```
pip install biopython
   # Installing collected packages: numpy, biopython
   # Successfully installed biopython-1.82 numpy-1.26.2
```

where did it get installed?
```
# show which packages are installed and where they are:
python3 -m pip list -v
```

on my mac laptop, packages are here: `/usr/local/lib/python3.11/site-packages`

So I want to add that to PYTHONPATH - I add this to my .profile file:
```
#### JY adding to (or creating) PYTHONPATH:
if [ -z ${PYTHONPATH+x} ] 
then
    echo "        making new PYTHONPATH"
    export PYTHONPATH=/usr/local/lib/python3.11/site-packages
else
    echo "        adding to PYTHONPATH"
    export PYTHONPATH=/usr/local/lib/python3.11/site-packages:${PYTHONPATH}
fi
```


