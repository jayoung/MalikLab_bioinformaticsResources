# Second try, Nov 2022

Now I use Amy's [newer instructions](https://hutchdatascience.org/FH_WDL101_Cromwell/introduction.html) from Nov 2022 

First, clone Amy's repo and make my own branch, so I can suggest changes via pull requests
```
cd ~/FH_fast_storage/cromwell-home
git clone https://github.com/fhdsl/FH_WDL101_Cromwell.git
```

# git forking/pull requests

Made some edits to the FH_WDL101_Cromwell instructions.  First I do `add/commit` for my my changes in preparation for a later pull request:
```
git add --all --verbose .
git commit
    # [main e549788] typos and suggested edits
```

Can't (or don't want to) push those changes to the main/master repo (https://github.com/fhdsl/FH_WDL101_Cromwell). Instead I want my own fork of the repo that I mess with, and submit pull requests from.

Within the local copy of the repo, I change the public location for push/pulls/etc:
```
git remote set-url origin https://github.com/jayoung/FH_WDL101_Cromwell.git

# check that worked
git remote -v
   # origin  https://github.com/jayoung/FH_WDL101_Cromwell.git (fetch)
   # origin  https://github.com/jayoung/FH_WDL101_Cromwell.git (push)

# try a push.  Doesn't work YET, because I didn't actually FORK the repo, I had just CLONED it
git push
   # remote: Repository not found.
   # fatal: repository 'https://github.com/jayoung/FH_WDL101_Cromwell.git/' not found
```

I used the 'fork' button on the github website to create a fork of the repo in https://github.com/jayoung 

Now I can do `git push` successfully, and the commit is pushed to https://github.com/jayoung/FH_WDL101_Cromwell.git.   

Then, from https://github.com/jayoung/FH_WDL101_Cromwell.git I can use a button to submit a pull request.




# notes working through the instructions

My AWS credentials are stored in `.aws/credentials`

I previously set up a database container and a database for my cromwell jobs and checked it's still running using the [Hutch myDB dashboard](https://mydb.fredhutch.org/list_containers/).  See [notes](tryWDLworkflows_v1.md).

I get a fresh copy of the `diy-cromwell-server` repo:
```
cd ~/FH_fast_storage/cromwell-home
rm -rf diy-cromwell-server/
git clone --branch main https://github.com/FredHutch/diy-cromwell-server.git

# make a template config file:
mv cromUserConfig.txt cromUserConfig.txt.backup 
cp ./diy-cromwell-server/cromUserConfig.txt .
    # and edit
```

Start up cromwell server:

Running this shell script looks like it actually wipes out the diy-cromwell-server subdir and clones the git repo all over again.  Maybe the user needs a warning not to use that dir for any personal notes / example files

```
cp ./diy-cromwell-server/cromwell.sh .
chmod +x cromwell.sh
./cromwell.sh cromUserConfig.txt
```
The output (Nov30) looked like this:
```
Your configuration details have been found...
Getting an updated copy of Cromwell configs from GitHub...
Setting up all required directories...
Detecting existence of AWS credentials...
Credentials found, setting appropriate configuration...
Requesting resources from SLURM for your server...
Submitted batch job 5041293
Your Cromwell server is attempting to start up on node/port gizmok169:54243.  It can take up to 2 minutes prior to the port being open for use by the shiny app at https://cromwellapp.fredhutch.org or via the R package fh.wdlR. If you encounter errors, you may want to check your server logs at /fh/fast/malik_h/user/jayoung/cromwell-home/server-logs to see if Cromwell was unable to start up.
Go have fun now.
```
If I want to stop the server, I use `scancel 5041293`

It crashed an hour or so after I started it (hadn't even used it yet), with a memory issue:
```
more server-logs/cromwell_5041293.out 
To execute cromwell, run: java -jar $EBROOTCROMWELL/cromwell.jar

To execute womtool, run: java -jar $EBROOTCROMWELL/wamtool.jar
[2022-11-30 12:55:07,34] [warn] modifyDataType will lose primary key/autoincrement/not null settings for mysql.  Use <sql> and re-speci
fy all configuration if this is the case
[2022-11-30 12:55:07,35] [warn] modifyDataType will lose primary key/autoincrement/not null settings for mysql.  Use <sql> and re-speci
fy all configuration if this is the case
[2022-11-30 12:55:07,35] [warn] modifyDataType will lose primary key/autoincrement/not null settings for mysql.  Use <sql> and re-speci
fy all configuration if this is the case
[2022-11-30 12:55:07,35] [warn] modifyDataType will lose primary key/autoincrement/not null settings for mysql.  Use <sql> and re-speci
fy all configuration if this is the case
OpenJDK 64-Bit Server VM warning: INFO: os::commit_memory(0x000014b73b9fa000, 16384, 0) failed; error='Not enough space' (errno=12)
#
# There is insufficient memory for the Java Runtime Environment to continue.
# Native memory allocation (mmap) failed to map 16384 bytes for committing reserved memory.
# An error report file with more information is saved as:
# /fh/fast/malik_h/user/jayoung/cromwell-home/hs_err_pid6812.log
[thread 21195 also had an error]
```

Start it again (Nov 30, 3.09pm):
```
./cromwell.sh cromUserConfig.txt
```
Output:
```
Your configuration details have been found...
Getting an updated copy of Cromwell configs from GitHub...
Setting up all required directories...
Detecting existence of AWS credentials...
Credentials found, setting appropriate configuration...
Requesting resources from SLURM for your server...
Submitted batch job 5065329
Your Cromwell server is attempting to start up on node/port gizmok25:39433.  It can take up to 2 minutes prior to the port being open for use by the shiny app at https://cromwellapp.fredhutch.org or via the R package fh.wdlR. If you encounter errors, you may want to check your server logs at /fh/fast/malik_h/user/jayoung/cromwell-home/server-logs to see if Cromwell was unable to start up.
Go have fun now.
```

Then I use a web browser to connect to the [shiny app](https://cromwellapp.fredhutch.org) that interacts with my cromwell database.


Note: every time we run `./cromwell.sh` it pulls a fresh copy of the `diy-cromwell-server` repo.  

It also starts up a server using the `./diy-cromwell-server/cromwellServer.sh` script, and a config file that will be either `./diy-cromwell-server/fh-S3-cromwell.conf` (if AWS credentials are found) or `./diy-cromwell-server/fh-cromwell.conf`.   The `cromwellServer.sh` script does some work to set up the environment for cromwell jobs, including `module purge`, `module --ignore-cache load cromwell/84`.



# Troubleshooting tg-wdl-VariantCaller


I add this line to the `task annovar {` block: `unset PERL5LIB` to get past an error


works!

all the example workflows work.

next step is try to make my own?

