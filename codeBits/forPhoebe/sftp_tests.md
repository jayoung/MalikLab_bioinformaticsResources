# sftp 

## Motuz:
Genewiz/Azenta provide Illumina sequencing services, and they provide the data via sftp (rather than regular ftp).  This makes it hard to use wget or sbatch approaches to parallelizing downloads. 

It turns out the answer is to use the Hutch's [Motuz](https://motuz.fredhutch.org) file transfer service.  Simple website that allows transfers via many protocols, including sftp, directly from external site onto `/fh/fast`.   It is working OK-ish, although it seems like there are sporadic authentication failures.  Might be due to trying to have too many connections open at once. It can handle some level of simultaneous downloads, but perhaps not as many as we asked for.  


## command-line sftp

before we realized that Motuz would work: interactive sftp worked OK but was very slow. In order to parallelize, we were trying to get the command-line password-less access set up. Here are my notes from when we were trying that.   I suspect our failure may have been because the remote host isn't set up for us to have ssh-like access, only sftp-like access.

Here's the [documentation](https://linux.die.net/man/1/sftp) for `sftp`.  There are various ways to use it:

### 1. basic interactive usage

`sftp remote.host.name` or `sftp user@remote.host.name`

It should ask for a password, and once you're in, youcan  use commands like `cd`, `mget` and `get` interactively.

This should work fine. However, because it requires you to enter user name and password, it'll be harder to parallelize using sbatch scripts.

### 2. set up password-less access

In order to run `sftp` in batch mode, we first need to set up password-free access to the remote machine. This is a one-time thing for each remote machine. I think we do that using the [`ssh-keygen`](https://linux.die.net/man/1/ssh-keygen) command, using two steps:
- generate an ssh key pair - see these [instructions](https://linuxiac.com/generate-ssh-key-pair/). This creates two files, a private key file (`~/.ssh/id_rsa`) and a public key file (`~/.ssh/id_rsa.pub`) 
- put the public key on the remote machine - see these [instructions](https://linuxiac.com/ssh-login-without-password/)

Perhaps we test whether that worked, using a command like `sftp remote.host.name`.

Once that's done, we should be able to connect to the remote server without using a username/password. 


#### 2a. use a batchfile

We might put the sftp commands we want to run in one or more `batchfiles`. I think the `batchfile` would be a plain text file that looks something like this:
```
cd xx_dir_xx
get xx_file1_xx
get xx_file2_xx
```

Once we have the batch file(s) we run sftp like this:

```
sftp -b batchfile user@host 
```

We could make one batchfile for each large file we want to get, and put this sftp command into an sbatch array script.

#### 2b. use echo to mimic a batchfile

With this method it might be even easier to use an sbatch array script, so let's try to make this work.

We would supply the commands to run using `stdin`, rather than using a `batchfile` ("A batchfile of '-' may be used to indicate standard input."). 
```
echo "cd xx_dir_xx ; get xx_file1_xx  ; get xx_file2_xx" | sftp -b - user@host 
```

(for password-free sftp, it's possible we need to mess with this `sftp` option: `-F ssh_config`, but I don't think so)
