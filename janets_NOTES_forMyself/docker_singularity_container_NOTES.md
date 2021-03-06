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