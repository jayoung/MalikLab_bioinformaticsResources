# documentation


# Nov 2025 notes
[AWS access credentials](https://sciwiki.fredhutch.org/scicomputing/access_credentials/#amazon-web-services-aws)
- login to AWS on browser via Hutch MyApps, AWS portal
- can get access keys for `fh-pi-malik-h` there (click the triangle next to `fh-pi-malik-h`)

Most of these notes are from July 2022. Updated some in Nov 2025

Also done Nov 2025:
- SSO credentials last for 12 hours
- updated `~/.aws/credentials` (using info found through portal login)

```
which aws
/app/bin/aws

module purge
module load awscli/2.17.54-GCCcore-12.3.0 

which aws
/app/software/awscli/2.17.54-GCCcore-12.3.0/bin/aws

aws configure sso

SSO session name (Recommended): my-sso
SSO start URL [None]: https://d-92674cb6d7.awsapps.com/start/#
SSO region [None]: us-west-2
SSO registration scopes [sso:account:access]:
Attempting to automatically open the SSO authorization page in your default browser.
If the browser does not open or you wish to use a different device to authorize this request, open the following URL:

https://d-92674cb6d7.awsapps.com/start/#/device

Then enter the code:

FZLB-BVFP
The only AWS account available to you is: 891846948394
Using the account ID 891846948394
The only role available to you is: PowerUserWRoleAccess
Using the role name "PowerUserWRoleAccess"
CLI default client Region [us-west-2]:
CLI default output format [None]:
CLI profile name [PowerUserWRoleAccess-891846948394]: default

To use this profile, specify the profile name using --profile, as shown:

aws s3 ls --profile default
```

Motuz works for transfer.

Some command line stuff:
```
module purge
module load awscli
# the / on the end is important otherwise we don't see the contents
aws s3 ls s3://fh-pi-malik-h-eco/user/

# to remove stuff:
aws s3 rm 
aws s3 rm --recursive
```

# Earlier notes

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

## setting AWS credentials for command-line use 

they're stored in `.aws/credentials`. They'll expire about once a day, so update that file when you get new ones.

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


[Understanding folder sizes better](https://bobbyhadz.com/blog/aws-s3-get-folder-size) - can do it through AWS web dashboard or through command-line

Rates July 18, 2022:
Fast storage is free for up to 1 TiB and is then $30/TiB/month
Economy storage is free for up to 100 TiB and is then $3/TiB/month


# AWS web dashboard

Web dashboard is [here](https://us-west-2.console.aws.amazon.com/console/home?nc2=h_ct&region=us-west-2&src=header-signin#).  I set up two factor authentication - open Google authenticator on my phone to get the code.

The costs we see on the AWS web dashboard are NOT the costs our lab actually gets charged (the Hutch subsidizes costs.)

From Center IT (see email from Brikti Kahsay, July 20, 2022):  
The data in cost explorer isn’t the rate you get charged, that’s the rate AWS charges to FH (which gets paid from CIT’s budget and then we recoup a portion of that through chargebacks).  
CIT does chargebacks at a subsidized rate for economy storage ($3/TB/month). While Cost Explorer is useful to broadly understand where costs are incurring, it’s not useful for understanding how much you will pay at the end of the day.  
You are correct that Simple Storage Service (called S3) is economy storage. There are two buckets with the name -eco in the title in your account – those are the buckets FH subsidizes.  
Certain services, like Virtual Private Cloud (VPC) are security services installed by my team. We don’t charge for those, but they’re important for maintaining the security of your account. You may see CloudWatch in there as well – same situation there (it’s just a logging and monitoring service). Feel free to disregard.  
If you log into your AWS account, go to the top of the homepage and type “S3”. You should get a list of all your S3 buckets (there should be 2 of them). Just above that list you should see a button that says “View Storage Lens Dashboard”. If you click on that, that will give you breakdowns and reports of your usage, including size of the data in each bucket.
 