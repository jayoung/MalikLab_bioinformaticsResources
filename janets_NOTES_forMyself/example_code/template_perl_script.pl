#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long 'HelpMessage';


#### set up defaults for the options
my $infile = "default_value";
my $option2 = 99;
my $verbose;


## GetOptions syntax:  https://perldoc.perl.org/Getopt/Long.html

## = for mandatory options 
## : for optional options (defaults of "" or 0 will be supplied if default is not specified)
## s for string, i for numeric
## or, with nothing after the option name, it's a flag
GetOptions ("in=s"       => \$infile, 
            "option2=i"  => \$option2, 
            "verbose",
            "help"       => sub { HelpMessage(0) } ) or HelpMessage(1);


######## Usage in POD format (used by HelpMessage):

=head1 NAME

template_perl_script.pl - description of what it does

=head1 SYNOPSIS

  --option1    description of options, including
  --option2    including what the defaults are and whether they're mandatory or optional
  --help,-h       Print this help

=head1 VERSION

0.01

=cut

################


### check args
if($infile eq "") {
    print "\n\nERROR - you must specify the input file\n";
    HelpMessage(1);
}
if(!-e $infile) {
    print "\n\nERROR - the input file you specified does not exist\n";
    HelpMessage(1);
}

## figure out output file name
my $outfile = $infile;
$outfile =~ s/\.fa$//;
$outfile .= ".processed.fa";


## to process a regular text file
# open (IN, "< $infile") || die "can't open $infile\n";
# open (OUT, "> $outfile");
# while (<IN>) {
#     my $line = $_; chomp $line;
#     ## do something!
# }
# close IN;
# close OUT;

## to process a fasta file
my $seqIN = Bio::SeqIO->new(-file => "< $file", '-format' => 'fasta');
my $seqOUT = Bio::SeqIO->new(-file => "> $out", '-format' => 'fasta');
while (my $seq = $seqIN ->next_seq) {
    my $id = $seq->display_id();
    my $letters = $seq->seq(); 
    my $desc = $seq->desc();
    my $newseqobj = Bio::Seq->new(-seq => $letters, -id => $desc);
    $seqOUT->write_seq($newseqobj);
}




