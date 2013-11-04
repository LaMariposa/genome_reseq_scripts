#!/usr/bin/perl

#readpairxml.pl by Megan Supple
#created 1 November 2013
#last modified 1 November 2013

#script to generate trace archive xml entires for read pairs from a fasta file for one of the pairs

#usage readpairxml.pl <in.fasta>

use lib $ENV{PERL5LIB};
 
use strict;
use warnings;
use Data::Dumper;

my $usage = "usage: readpairxml.pl <in.fasta>
arguments (required):
	<in.fasta> is the fasta file for one set of paired sequencing reads
";

die "$usage" unless (@ARGV == 1);

#read in command line arguments
my ($infasta)=@ARGV;

#open input fasta and output xml file
open(INFASTA, $infasta)||die "can't open input fasta file. $!\n";
my @prefix=split(/\./, $infasta);

#my $outxml=join('.',$prefix[0],"readpairs.xml");
my $outxml="temp.pairs.xml";
open(OUTXML, ">$outxml");

#print header
print OUTXML "\t<LIBRARY ID=\"small\" NAME=\"small\" MIN=\"200\" MAX=\"500\">\n";

#make an entry for each read pair
#loop through the input fasta file
my $r=1;
my $line;
while($line=<INFASTA>)
	{
	  #get seq id
	  chomp $line;
	  my @entry=split(">",$line);
	  my @entry2=split("/",$entry[1]);
	  my $read=$entry2[0];
	  #print entry based on header
	  print OUTXML "\t\t<INSERT ID=\"seq_$r\" NAME=\"SEQ$r\">
			\t\t\t<SEQUENCE ID=\"seq1\" NAME=\"\@$read/1\"/>
			\t\t\t<SEQUENCE ID=\"seq2\" NAME=\"\@$read/2\"/>
			</INSERT>\n";
	  #read in and ignore sequence
	  $line=<INFASTA>;
	  $r++;
	}

#print footer
print OUTXML "\t</LIBRARY>";




close INFASTA;
close OUTXML;
