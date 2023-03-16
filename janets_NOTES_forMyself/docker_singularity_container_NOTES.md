# Resources

An [intro](https://www.freecodecamp.org/news/a-beginner-friendly-introduction-to-containers-vms-and-docker-79a9e3e119b)

Using containers [at the Hutch](https://sciwiki.fredhutch.org/compdemos/Singularity/#access-to-storage-on-fhfast-and-fhscratch-inside-singularity-containers)

Singularity [intro](https://sylabs.io/guides/3.5/user-guide/quick_start.html)

Erick Matsen's [intro](http://erick.matsen.org/2018/04/19/docker.html)

# my first docker 

Trying this for the [pamlWrapper project](https://github.com/jayoung/pamlWrapper/blob/main/buildContainer/docker_container_NOTES.md)

# Docker getting-started example
On my Mac, I first installed the Docker desktop client (that also installed command-line docker).

Then, in Terminal, I start up a container that has their tutorial materials
```
docker run -d -p 80:80 docker/getting-started
```
A docker container appears in the Docker-desktop app's dashboard. I can use the app to start a shell within that container, and to inspect logs, or to delete containers.

From the app's control panel I also clicked, 'run in browser' - that opened up the tutorial material. That let me download a folder of stuff called 'app': it's an example to-do-list app they wrote as a demo.

They have me create a simple Dockerfile within that folder. 

## Build the container
then from the same folder I run this:
```
cd /Users/jayoung/dockerStuff/getting-started/app
docker build -t getting-started .
```
The `-t` lets us provide a tag (= a human-readable name for the container).  
The `.` means look for a `Dockerfile` in the current dir.  

## Run the container
```
docker run -dp 3000:3000 getting-started
```
`-d` means run the new container in "detached" mode (in the background)  
`-p 3000:3000` creates a mapping between the host's port 3000 to the container's port 3000. Without the port mapping, we wouldn't be able to access the application.

Then I can open web browser to http://localhost:3000 and the app is running.

## Changing a container

We would:
- make edits to the code
- stop and remove the container
- rebuild it and restart it

# Managing my docker images and containers (or use Mac gui):
```
# containers
docker ps -a
docker stop <the-container-id>
docker rm <the-container-id>
# Or stop and remove a container in a single command by adding the "force" flag to the docker rm command: 
docker rm -f <the-container-id>

# images:
docker images
docker system df
docker rmi <the-image-id>
```

# Using containers on rhino/gizmo

Rhino/gizmo don't allow docker containers. We first have to convert the docker image to a singularity/apptainer image.

Apptainer is the same thing as singularity, just a newer version. In Feb 2023 Dan told me "I did experience some problems using the latest version (and reported them) but 1.0.1 seems to work fine"

Example:
```
cd ~/FH_fast_storage/paml_screen/pamlWrapper/buildContainer

module purge
module load Apptainer/1.1.6

apptainer build paml_wrapper-v1.3.5.sif docker://jayoungfhcrc/paml_wrapper:version1.3.5
```
A file appears called paml_wrapper-v1.3.5.sif - this file can be used just like a docker image. Example, to run a command within the container (`--cleanenv` ensures the environment outside the container isn't used within the container):
```
singularity exec --cleanenv paml_wrapper-v1.3.5.sif /pamlWrapper/scripts/pw_testScript.bioperl
```

That can go in a shell script for use in sbatch.

Or to get a shell in the container:
```
apptainer run --cleanenv paml_wrapper-v1.3.5.sif
```
my rhino home dir (`/home/jayoung`) is available within the container but `/fh/fast` is not.


To get a shell in the container, mounting PWD and starting the shell within it:
```
apptainer shell --cleanenv --bind $(pwd):/mnt -H /mnt /fh/fast/malik_h/grp/malik_lab_shared/singularityImages/paml_wrapper-v1.3.6.sif
```

To run a command within the container, mounting PWD:
```
apptainer exec --cleanenv --bind $(pwd):/mnt -H /mnt /fh/fast/malik_h/grp/malik_lab_shared/singularityImages/paml_wrapper-v1.3.6.sif /pamlWrapper/scripts/pw_makeTreeAndRunPAML.pl CENPA_primates_aln2a_only5seqs.fa 
```

```
module purge
```

# Understanding containers, kernels, operating systems

The container still runs on the kernel of whichever machine it's running on: that's what `uname -a` will show from within the container.  

If we're running on a Mac or Windows computer, the computer needs a virtual machine in order to make the container work, rather than running directly on the Mac kernel, so `uname -a` shows something like this:
```
root@74ceeece7d11:/# uname -a
Linux 74ceeece7d11 5.15.49-linuxkit #1 SMP Tue Sep 13 07:51:46 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux
```  

On top of the kernel, there's an operating system (OS).  To see what OS we're using, do this: `cat /etc/*-release`.  That will be whatever was the base OS used when building the docker image, and is often different from the OS used outside the container. 

Example:  inside paml_wrapper-v1.3.5 container, `cat /etc/*-release` includes the info `DISTRIB_DESCRIPTION="Ubuntu 14.04.4 LTS"` (whether I'm running it via Docker on my Mac, or via apptainer/singularity on gizmo/rhino)/

Outside the container, on gizmo/rhino, `cat /etc/*-release` includes `DISTRIB_DESCRIPTION="Ubuntu 18.04.6 LTS"` (and I can't even run the cat command on the mac, because those release files don't exist)



# public images

```
docker search continuumio
```

# Notes from my paml_wrapper work

see docker notes halfway down [this doc](https://github.com/jayoung/pamlWrapper/blob/main/README.md)

Build an image from a `Dockerfile`:
```
cd /Users/jayoung/gitProjects/pamlWrapper/buildContainer
docker build -t paml_wrapper .
```

Test a docker image for vulnerabilities: not sure I really need to do this. With paml_wrapper it DOES find various problems, I think all with the base debian system. None are labeled 'critical'. See [here](https://docs.docker.com/get-started/09_image_best/). I think I'll just proceed.
```
docker scan paml_wrapper
```

Show how the build layers were added:
```
docker image history paml_wrapper
```

Open a shell in that core, mounting the entire current dir into workingDir:
```
docker run -v `pwd`:/workingDir -it paml_wrapper
```
Can mount any local directory, but it must be an ABSOLUTE path for this to work right (at least on my mac)

Or, to start the container in detached mode so I can come in and out of it more easily, we add the -d switch
```
docker run -v `pwd`:/workingDir -itd paml_wrapper
docker ps     # to get the container-id
docker exec -it <container-id> /bin/bash
   # then do the work
exit
```

To later remove the container: 
```
docker rm -f 163df768287c
```


# Rasilab example

Bottorff et al:  
[git repo](https://github.com/rasilab/bottorff_2022)

There's a Dockerfile at the top level of that repo that will create our Docker container.

It starts with a miniconda base.

Downloads and upacks some pre-compiled binaries and adds their location to PATH

It sets up the conda enviroment

It uses conda to install some stuff, including R and some packages. It uses some .yml files in the /install dir to do that.

When someone wants to run the code, they first install the software: `docker build -t bottorff_2022 .` This creates a docker image (??).

Then you can run code, e.g. on a mac:
```
docker run --rm -it -v $(pwd):/workspace bottorff_2022 /bin/bash
cd workspace 
sh run_everything.sh
```
or, more complex, can run on the FHCRC cluster via sigularity - see [here](https://github.com/rasilab/bottorff_2022)