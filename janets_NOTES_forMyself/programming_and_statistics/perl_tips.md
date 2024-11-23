# Perl tips

www.perl.com

man pages (xman)

in regular expressions:
^ means the beginning of the line
$ means the end of the line
[^t] means match anything except t
\s is a whitespace

$string =~ s/\s.+$// will take only the first word of the line

\G match starting from where last m//g left off.

e.g. 
my @array = $string =~ m/\G.{80}/g;
my $newstring = join("\n",@array);

puts a newline in every 80 characters

+ match any number of the preceding character
? match until the first occurrence of the next character
putting the match expression in brackets means you can get it out using $1, $2 
(numbers refer to order of open bracket - brackets can be nested)

$/ is a special variable - use it to define what a file gets split on in a while loop.
the default is /n. to split on one or more empty lines, $/ = ""
undef $/ will get the entire contents of the source file.

print UC ($bacname);
gives uppercase version of bacname (if I get command UC not found then 
need to do as 'use' command at top of script - look on perl web page to find it)

to find out what class an object is in:

```
my $class1 = ref $thisaln;
```

or:
```
use Scalar::Util qw(blessed);
my $class2 = blessed ($thisaln);
```

You can use $#fred to get the index value of the last element of @fred.

decimal places:
```
$a = 0.255;
$b = sprintf("%.2f", $a);
print "Unrounded: $a\nRounded: $b\n";
printf "Unrounded: $a\nRounded: %.2f\n", $a;
```
the %.2f formats to 2 decimal places 

byte order, endianness:
```
use Config;
print "byteorder " . $Config{ byteorder } . "\n";
print "unpack " . unpack 'I', "\x01\x02\x03\x04" . "\n";
print "\ndone\n\n";
```

# perl modules

## module installation

### new method - cpanm

see notes in ~/FH_fast_storage/source_codes/CPAN/aa_JY_CPAN_NOTES.txt

### old method, on the mac

use "sudo" to run command as super-user

```
perl Makefile.PL 
make
make test
sudo make install
```

### old method 1 rhino/gizmo, the Makefile.PL method

installing perl modules:
download from CPAN, uncompress
 
perl Makefile.PL LIB=/home/jayoung/malik_lab_shared/perl/lib INSTALLMAN1DIR=/home/jayoung/malik_lab_shared/man/man1 INSTALLMAN3DIR=/home/jayoung/malik_lab_shared/man/man3 INSTALLBIN=/home/jayoung/malik_lab_shared/bin INSTALLSCRIPT=/home/jayoung/malik_lab_shared/scripts INSTALLSITEMAN3DIR=/home/jayoung/malik_lab_shared/man/man3

or with some more things specified:

perl Makefile.PL LIB=/home/jayoung/malik_lab_shared/perl/lib PREFIX=/home/jayoung/malik_lab_shared/perl/lib INSTALLSITELIB=/home/jayoung/malik_lab_shared/perl/lib INSTALLMAN1DIR=/home/jayoung/malik_lab_shared/man/man1 INSTALLMAN3DIR=/home/jayoung/malik_lab_shared/man/man3 INSTALLBIN=/home/jayoung/malik_lab_shared/perl/bin INSTALLSCRIPT=/home/jayoung/malik_lab_shared/scripts 

make 
make test
make install

### old method 2 rhino/gizmo, the Build.PL method

or, the build method:
perl Build.PL LIB=/home/jayoung/malik_lab_shared/perl/lib PREFIX=/home/jayoung/malik_lab_shared/perl/lib INSTALLSITELIB=/home/jayoung/malik_lab_shared/perl/lib INSTALLMAN1DIR=/home/jayoung/malik_lab_shared/man/man1 INSTALLMAN3DIR=/home/jayoung/malik_lab_shared/man/man3 INSTALLBIN=/home/jayoung/malik_lab_shared/perl/bin INSTALLSCRIPT=/home/jayoung/malik_lab_shared/perl/scripts
./Build
./Build test
./Build install


## module versions

find out version of a perl module:
> perl -MModule -e 'print "$Module::VERSION\n"'
For example, with Net::SSLeay:
> perl -MNet::SSLeay -e 'print "$Net::SSLeay::VERSION\n"'
1.25

to check versions in bioperl:
print "root version ", $Bio::Root::Version::VERSION, "\n";
print "PAML version ", $Bio::Tools::Phylo::PAML::VERSION, "\n";

print "Codeml version ";
printf "%vd\n", $Bio::Tools::Run::Phylo::PAML::Codeml::VERSION, "\n\n";

## module locations

to find out [where a module is being loaded from](http://www.unix.org.ua/orelly/perl/prog3/ch30_02.htm
):
```
perldoc -l warnings
perldoc -l R
perldoc -l RReferences
```
or, within a script:
```
use Bio::Factory::EMBOSS;
print "\n" . $INC{"Bio/Factory/EMBOSS.pm"} . "\n\n";
```

if it says library -lgdbm is needed, the library to look for is libgdbm
(setenv LD_LIBRARY_PATH may be useful here)

# installing perl itself on a mac

## notes from June 2019

1. figure out whether perl is installed?  (if it is, this command report where it is installed, if it is not, will report nothing):
    which perl
Sounds like it really should be installed on newer macs

2. check that perl is working (this also happens to report the installed version):
perl -v
   
3. here's one way to check whether various modules are installed:
     perldoc -l moduleName
If the module is installed, the output will report where it installed, if not it will say "No documentation found"

e.g. "warnings" module should be installed
perldoc -l warnings
    /System/Library/Perl/5.18/warnings.pm
    
e.g. Bio::SeqIO is a basic bioperl module and is probably not installed:
perldoc -l Bio::SeqIO
    No documentation found for "Bio::SeqIO".

(and from here on I'm following instructions from this website - might be useful if we need to troubleshooot:
   http://www.bioperl.org/wiki/Installing_Bioperl_for_Unix#INSTALLING_BIOPERL_THE_EASY_WAY
I'm skipping the "Install expat" stage from those instructions, as I think this is already installed on OS X)


4. upgrade/install the "cpan" program - this is prgram that makes it easy to install perl modules

run this command:
   perl -MCPAN -e shell

I was then asked three questions, and gave these answers:

a. Q=Would you like to configure as much as possible automatically? [yes] 
   A=yes

b. Q=What approach do you want?  (Choose 'local::lib', 'sudo' or 'manual') [local::lib] 
   A=sudo

c. Q=Would you like me to automatically choose some CPAN mirror sites for you? (This means connecting to the Internet) [yes] 
   A=yes

It did some stuff, and then gave me a cpan shell prompt:
cpan[1]> 
at which I did this:

install Bundle::CPAN
   (had to give it my mac password to complete the install)
   (got a weird prompt "Enter arithmetic or Perl expression: exit" and pressed return)
   (looks like it worked, so next command is to quit cpan shell:
q

5. Use the "cpan" program Install/upgrade Module::Build, and make it your preferred installer:

cpan
install Module::Build
o conf prefer_installer MB
o conf commit
q


6. install bioperl using cpan
cpan
   Then find the name of the most recent Bioperl version:
d /bioperl/
   I got this output, from which I choose the most recent version of BioPerl (BioPerl-1.6.924.tar.gz):
Reading '/Users/jayoung/.cpan/Metadata'
  Database was generated on Thu, 18 Jun 2015 22:17:02 GMT
Distribution    BOZO/Fry-Lib-BioPerl-0.15.tar.gz
Distribution    CDRAUG/Dist-Zilla-PluginBundle-BioPerl-0.20.tar.gz
Distribution    CJFIELDS/BioPerl-1.6.901.tar.gz
Distribution    CJFIELDS/BioPerl-1.6.923.tar.gz
Distribution    CJFIELDS/BioPerl-1.6.924.tar.gz
Distribution    CJFIELDS/BioPerl-DB-1.006900.tar.gz
Distribution    CJFIELDS/BioPerl-Network-1.006902.tar.gz
Distribution    CJFIELDS/BioPerl-Run-1.006900.tar.gz
Distribution    CJFIELDS/Bundle-BioPerl-2.1.9.tar.gz
Distribution    CJFIELDS/Dist-Zilla-PluginBundle-BioPerl-0.23.tar.gz
Distribution    RBUELS/Dist-Zilla-PluginBundle-Bioperl-0.01.tar.gz
11 items found

Then I install that:
install CJFIELDS/BioPerl-1.6.924.tar.gz

I got some warnings, but I proceeded (I'll deal with those later):

  requires:
    !  Data::Stag is not installed
  build_requires:
    !  Test::Most is not installed
  recommends:
    *  Algorithm::Munkres is not installed
    *  Array::Compare is not installed
    *  Bio::Phylo is not installed
    *  Convert::Binary::C is not installed
    *  GD is not installed
    *  Graph is not installed
    *  GraphViz is not installed
    *  HTML::TableExtract is not installed
    *  PostScript::TextBlock is not installed
    *  SOAP::Lite is not installed
    *  SVG is not installed
    *  SVG::Graph is not installed
    *  Set::Scalar is not installed
    *  Sort::Naturally is not installed
    *  Spreadsheet::ParseExcel is not installed
    *  XML::DOM is not installed
    *  XML::DOM::XPath is not installed
    *  XML::Parser::PerlSAX is not installed
    *  XML::SAX::Writer is not installed
    *  XML::Twig is not installed

Checking optional features...
EntrezGene............disabled
  requires:
    ! Bio::ASN1::EntrezGene is not installed
MySQL Tests...........disabled
  requires:
    ! DBD::mysql is not installed
Pg Tests..............disabled
  requires:
    ! DBD::Pg is not installed


Q=Do you want to run the Bio::DB::GFF or Bio::DB::SeqFeature::Store live database tests? y/n [n] 
A=n

Q=Install [a]ll BioPerl scripts, [n]one, or choose groups [i]nteractively? [a] 
A=a

Q=Do you want to run tests that require connection to servers across the internet
(likely to cause some failures)? y/n [n] 
A=n

gave it my Mac password

seemed like it worked

then I install the modules listed above as being not installed:
install Test::Most Array::Compare Bio::Phylo Convert::Binary::C GD Graph GraphViz HTML::TableExtract SOAP::Lite SVG SVG::Graph Set::Scalar Sort::Naturally Spreadsheet::ParseExcel XML::DOM XML::DOM::XPath XML::Parser::PerlSAX XML::SAX::Writer XML::Twig Bio::ASN1::EntrezGene DBD::mysql DBD::Pg PostScript::TextBlock Algorithm::Munkres

there may be some errors - don't think they're important: let's try running some Bioperl scripts and troubleshoort errors later.  I'm guessing those errors just affect a very tiny part of the Bioperl functionality that we're not likely to use.

q


I found that perldoc didn't work after I did that installation, so I fixed it like this:
sudo chmod uga+x /usr/bin/perldoc

then I tested whether Bioperl is installed, and it is:
perldoc -l Bio::SeqIO
   output should be: /Library/Perl/5.18/Bio/SeqIO.pm

in future, if you try to run a perl script and you which looks like this at the beginning:
    Can't locate Algorithm/Munkres.pm in @INC (you may need to install the Algorithm::Munkres module) 

then you can do install that module using the "cpan" program": 
cpan 
install Algorithm::Munkres
   (and you may need to give it your password
q (to quite cpan)



7. then install Bio-EUtilities (it is not part of the main Bioperl release and we need it for the getgenbankJY_eutils.bioperl script) 

unfortunately it does not seem to install well using the cpan method, so we do it manually instead:

a. download a file for the module into some kind of temporaray folder:
http://search.cpan.org/CPAN/authors/id/C/CJ/CJFIELDS/Bio-EUtilities-1.73.tar.gz
   (or click the download link near the top: http://search.cpan.org/dist/Bio-EUtilities/)

b. cd to that temp folder, uncompress the tar.gz file
tar -xzf Bio-EUtilities-1.73.tar.gz
cd Bio-EUtilities-1.73

c. run these commands
perl Makefile.PL  
make
make test
sudo make install
   (and supply your password)

this probably worked - try running the getgenbankJY_eutils.bioperl script now. can now delete the temporary folder

