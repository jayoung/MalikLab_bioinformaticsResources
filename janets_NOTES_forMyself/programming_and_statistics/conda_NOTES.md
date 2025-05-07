# using conda 

May 2025:  Dan Tenenbaum recommends using `uv` instead of conda. Dan's advice (with a note or two):

----
IMO Conda is bad and we should not be using it. I notice they provide an alternate set of installation instructions using pip. Can you try that?

You should probably create a virtual environment to install this in.

uv is a new package manager that can work just like pip but it’s faster, so something like this:

```
cd where-you-want-your-virtual-env
module purge

# I needed to add this line:
module load fhPython/3.11.3-foss-2023a

module load uv/0.6.4

uv venv # creates .venv dir within the current folder
source .venv/bin/activate
git clone https://github.com/maxplanck-ie/snakepipes
cd snakepipes

## Dan's install command didn't work for me:
# uv pip install .

## but this did:
uv pip install 'snakepipes @ .'

```

Then you should have the snakePipes command.

In subsequent sessions just run the “source” line to activate the venv.
-----



"Conda is an open-source package and environment management system." 

So it's for packages AND virtual envs.

If it's just python packages, I can use pip as an alternative, I think. But if something uses a virtual environment, I think I need to use conda.

As of 2024 the licensing agreements changed to do with the anaconda repository (where it installs from).  The Hutch has a workaround - https://conda-forge.fredhutch.org/ 



My `~/.condarc` file looked like this before May 2 2025:
```
channels:
  - conda-forge
  - bioconda
channel_priority: strict
```

On May 2 2025 I edited it - now it looks like this:
```
channels:
  - conda-forge
  - bioconda
channel_priority: strict
channel_alias: https://conda-forge.fredhutch.org
```


On May 2 2025, here are the contents of `~/.conda`:

```
tree ~/.conda
.conda
├── envs
└── pkgs
    ├── cache
    │   ├── 09cdf8bf.json
    │   ├── 09cdf8bf.q
    │   ├── 258f2e84.json
    │   ├── 258f2e84.q
    │   ├── 2a957770.json
    │   ├── 2a957770.q
    │   ├── 497deca9.json
    │   ├── 497deca9.q
    │   ├── ffeee55f.json
    │   └── ffeee55f.q
    ├── urls
    └── urls.txt
3 directories, 12 files
```

