#!/usr/bin/perl

use lib $ENV{PERL5LIB};
use strict;
use warnings;

my ($infasta)=@ARGV;
my $line;


my @M=("A","C");
my @R=("A","G");
my @W=("A","T");
my @S=("C","G");
my @Y=("C","T");
my @K=("G","T");
my @B=("C","G","T");
my @D=("A","G","T");
my @H=("A","C","T");
my @V=("A","C","G");



#parse input fasta file, looking for ambiguous codes
open(INFASTA, $infasta)||die "can't open input fasta file. $!\n";
open(OUT,">$infasta.randiupac.fasta");
while($line=<INFASTA>)
	{
	  #test if header line
	  if ($line=~/^>./){print OUT $line;}
	     else
		{
		  if ($line=~/[MRWSYKBDHV]/)
			{
			  #line contains an ambigous base, so go character by character
			  my @chars=split("",$line);
			  foreach(@chars)
				{
				  if ($_=~/[ACGT]/){print OUT "$_";}
				    elsif ($_=~"M"){print OUT $M[rand @M];}
                                    elsif ($_=~"R"){print OUT $R[rand @R];}
                                    elsif ($_=~"W"){print OUT $W[rand @W];}
                                    elsif ($_=~"S"){print OUT $S[rand @S];}
                                    elsif ($_=~"Y"){print OUT $Y[rand @Y];}
                                    elsif ($_=~"K"){print OUT $K[rand @K];}
                                    elsif ($_=~"B"){print OUT $B[rand @B];}
                                    elsif ($_=~"D"){print OUT $D[rand @D];}
                                    elsif ($_=~"H"){print OUT $H[rand @H];}
                                    elsif ($_=~"V"){print OUT $V[rand @V];}
				    else {print OUT "$_";}
				}
			}
		     else {print OUT $line;}
		}
	}
