# documentation

https://sciwiki.fredhutch.org/scicomputing/access_credentials/#gui-instructions

# Ways to access AWS storage

## command-line client
```
cd ~/FH_fast_storage/temp/testAWS
module load awscli

echo hello | aws s3 cp - s3://fh-pi-malik-h-eco/hello.txt

aws s3 ls s3://fh-pi-malik-h-eco/

# note that you NEED the / after the dir name to list contents - these two things are different:
aws s3 ls s3://fh-pi-malik-h-eco/user
aws s3 ls s3://fh-pi-malik-h-eco/user/

# neither of these behave as I would like, to list >1 dir:
aws s3 ls s3://fh-pi-malik-h-eco/*
aws s3 ls s3://fh-pi-malik-h-eco/*/

module purge
```

## setting AWS credentials for command-line use (a one-time thing)

they're stored in `.aws/credentials`
```
cd 
module purge
module load awscli
aws configure
module purge
```

## Motuz tool 

The Hutch's [Motuz tool](https://motuz.fredhutch.org) helps us see what's stored in AWS, transfer files between `/fh/fast` and AWS

# costs

Hutch rates shown [here](https://centernet.fredhutch.org/cn/u/center-it/scicomp/data-resource.html)

Hutch [billing page](https://grafana.fredhutch.org/d/dy5I3SIMk/data-core-storage-usage?orgId=1&var-storage_type=All&var-identifier=malik_h&var-division=All) for overall summary of our use

A [granular view](https://storage-hotspots.fhcrc.org) of what we're storing. Can download a csv file and look.


Rates July 18, 2022:
Fast storage is free for up to 1 TiB and is then $30/TiB/month
Economy storage is free for up to 100 TiB and is then $3/TiB/month


# AWS web dashboard

Web dashboard is [here](https://us-west-2.console.aws.amazon.com/console/home?nc2=h_ct&region=us-west-2&src=header-signin#).  I set up two factor authentication - open Google authenticator on my phone to get the code.
