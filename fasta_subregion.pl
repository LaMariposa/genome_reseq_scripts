#!/usr/bin/perl


#fasta_subregion.pl by Megan Supple 03 April 2012
#script to extract a specified region from entries in an aligned fasta file
#usage fasta_subregion.pl input.fasta start_coord stop_coord


use lib $ENV{PERL5LIB};
use strict;
use warnings;
use Getopt::Long;


#read in command line arguments
my ($infasta, $start1, $stop)=@ARGV;


my $usage = "Usage: fasta_subregion.pl <infile.fasta> <start> <stop>
options:
        <infile.fasta>          an input aligned fasta file
        <start>			first position to extract
        <stop>			last position to extract
";

die "$usage" unless (@ARGV == 3);
my ($infile) = @ARGV;


#open input file and output file
open(IN, $infasta)||die "can't open input fasta file. $!\n";
open(OUT, ">$infasta$start1-$stop.fasta");

#need to shift the frame
my $start=$start1-1;

#declare variables
my $line;
my $seq;
my $header;
my $width=50;

#read in each line of the file
while($line=<IN>)
	{
	  chomp($line);
	  if ($line =~ m/^>.*/)
		{
		  #if there is sequence
		  if ($seq)
			{
		  	  #print header to fasta file
			  print OUT "$header:$start1-$stop\n";
			  #print sequence of interest in lines of $width
			  my $length=$width; 
			  for (my $pos=$start; $pos<$stop;$pos+=$width)
				{
				  if ($pos+$width>$stop){$length=$stop-$pos;}
				  print OUT substr($seq, $pos, $length), "\n";
				}
			}
		  #set header
		  $header=$line;
		  #clear sequence
		  $seq="";
		}
		else
		{
		  $seq .=$line
		}
	}

			  #print final entry
                          #print header to fasta file
                          print OUT "$header:$start1-$stop\n";
                          #print sequence of interest in lines of $width
                          my $end=$width;
                          for (my $pos=$start; $pos<$stop;$pos+=$width)
                                {
                                  if ($pos+$width>$stop){$end=$stop-$pos;}
                                  print OUT substr($seq, $pos, $end), "\n";
                                }

close IN;
close OUT;
