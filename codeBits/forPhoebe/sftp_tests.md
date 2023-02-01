# sftp 

`sftp` [manpage](https://linux.die.net/man/1/sftp)

There are various ways to use sftp

## 1. basic interactive usage

`sftp user@host`.  

It should ask for a password, then you use commands like `cd`, `mget` and `get` interactively.

This should work fine, and will be easy, but will require you to enter user name and password, making it harder to parallelize using sbatch.

## 2. set up password-less access

In order to run `sftp` in batch mode, we first need to set up password-free access to the remote machine. This is a one-time thing for each remote machine. I think we do that using the [`ssh-keygen`](https://linux.die.net/man/1/ssh-keygen) command? There are two steps:
- generate an ssh key pair - see these [instructions](https://linuxiac.com/generate-ssh-key-pair/) 
- put the public key on the remote machine - see these [instructions](https://linuxiac.com/ssh-login-without-password/). 

Once that's done, we should be able to connect without using a username/password.

### 2a. use a batchfile

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

### 2b. use echo to mimic a batchfile

With this method it might be even easier to use an sbatch array script, so let's try to make this work.

We would supply the commands to run using `stdin`, rather than using a `batchfile` ("A batchfile of '-' may be used to indicate standard input."). 
```
echo "cd xx_dir_xx ; get xx_file1_xx  ; get xx_file2_xx" | sftp -b - user@host 
```

(for password-free sftp, it's possible we need to mess with this `sftp` option: `-F ssh_config`, but I don't think so)
