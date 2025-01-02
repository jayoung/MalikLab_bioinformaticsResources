alias rhino="ssh -Y rhino.fhcrc.org"
alias rhino1="ssh -Y rhino01.fhcrc.org"
alias rhino2="ssh -Y rhino02.fhcrc.org"
alias rhino3="ssh -Y rhino03.fhcrc.org"
alias testJY="echo testingTesting"


# use perlbrew rather than system perl:
source ~/perl5/perlbrew/etc/bashrc

### set PAML_WRAPPER_HOME
export PAML_WRAPPER_HOME=/Users/jayoung/gitProjects/pamlWrapper

### add to PATH
export PATH="/usr/local/opt/python@3.11/libexec/bin:$PATH"
export PATH=$PATH:$HOME/bin
export PATH=$PATH:$HOME/bin/ensembl-vep-release-104.3
# needed for bioconductor to see vep:
export PATH=$PATH:$HOME/bin/ensembl-vep-release-104.3/vep
export PATH=$PATH:/Library/Frameworks/R.framework/Versions/Current/Resources
export PATH=$PATH:/Applications/Apollo-2.0.5/bin


######## JY adding to (or creating) PERL5LIB for general perl modules
# first path element
export TEMP_PL=$HOME/perl5/perlbrew/perls/perl-5.34.0/lib/5.34.0
# adding more path elements
export TEMP_PL=$HOME/perl5/perlbrew/perls/perl-5.34.0/lib/site_perl/5.34.0:$TEMP_PL
export TEMP_PL=$HOME/perl5/perlbrew/perls/perl-5.34.0/lib/site_perl/5.34.0/darwin-2level:$TEMP_PL

if [ -z ${PERL5LIB+x} ] 
then
    echo "        making new PERL5LIB"
    export PERL5LIB=$TEMP_PL
else
    echo "        adding to PERL5LIB"
    export PERL5LIB=$TEMP_PL:${PERL5LIB}
fi
unset TEMP_PL



####### JY adding to (or creating) PYTHONPATH:
if [ -z ${PYTHONPATH+x} ] 
then
    echo "        making new PYTHONPATH"
    export PYTHONPATH=/usr/local/lib/python3.11/site-packages
else
    echo "        adding to PYTHONPATH"
    export PYTHONPATH=/usr/local/lib/python3.11/site-packages:${PYTHONPATH}
fi




##### application specific things:
export ENTREZ_KEY=e9e901b80c0eedeb4712c5f0ca0646793e08



###### set DOCKER_BUILDKIT: allows RUN --mount=type=bind" commands in Dockerfiles. And much more: https://github.com/moby/buildkit/blob/master/frontend/dockerfile/docs/syntax.md
export DOCKER_BUILDKIT=1

######## other

# Your previous /Users/jayoung/.profile file was backed up as /Users/jayoung/.profile.macports-saved_2017-07-20_at_18:24:45
##
# MacPorts Installer addition on 2017-07-20_at_18:24:45: adding an appropriate PATH variable for use with MacPorts.
export PATH="/opt/local/bin:/opt/local/sbin:$PATH"
# Finished adapting your PATH environment variable for use with MacPorts.

#THIS MUST BE AT THE END OF THE FILE FOR SDKMAN TO WORK!!!
export SDKMAN_DIR="/Users/jayoung/.sdkman"
[[ -s "/Users/jayoung/.sdkman/bin/sdkman-init.sh" ]] && source "/Users/jayoung/.sdkman/bin/sdkman-init.sh"

# MacPorts Installer addition on 2022-04-19_at_16:50:14: adding an appropriate MANPATH variable for use with MacPorts.
export MANPATH="/opt/local/share/man:$MANPATH"
# Finished adapting your MANPATH environment variable for use with MacPorts.

