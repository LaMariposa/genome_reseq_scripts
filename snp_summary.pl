#!/usr/bin/perl

#snp_summary.pl by Megan Supple, 2 Mar 2012
#script to summarize snp calls from the matrix genotype file

#usage  snp_summary.pl <in.matrix> <ref.fasta>
	#<in.matrix> is the matrix genotype format file generated from vcf2format.pl
	#<ref.fasta> is a fasta reference file

#output is
	#STDOUT 		the number of positions with variation in any of the samples 
	#genotyped_summary.txt	number of positions genotyped per sample
	#het_summary.txt	number of heterozygote snps per sample
	#snp_summary.txt	number of snps relative to the reference per sample

use lib $ENV{PERL5LIB}; 
use strict; 
use warnings;

my $usage = "Usage:  snp_summary.pl <in.matrix> <ref.fasta>
options (required):
	<in.matrix>	a matrix of genotypes (formatted as in vcf2format.pl)
	<ref.fasta>	a reference sequence in fasta format
";

die "$usage" unless (@ARGV == 2);

#read in command line arguments
my ($ingenos, $inref)=@ARGV;

#open input file and output file for snp summar
open(GENOS, $ingenos)||die "can't open input fst file. $!\n";
open(OUTGENO, ">genotyped_summary.txt");
open(OUTSNP, ">snp_summary.txt");
open(OUTHET, ">het_summary.txt");

#read in reference sequence
open(REF, $inref)||die "can't open input fst file. $!\n";
my %ref;
my $contig;
my $seq="";

while (my $fasta=<REF>)
	{
	  chomp($fasta);
	  #check if new contig, if so add to hash
	  if ($fasta=~ m/^>.*/)
		{
		  #record contig name
		  $contig=substr($fasta, 1);
		  #reset seq to null
		  $seq="";
		}
		else
		{
		  #add seq
		  $seq.=$fasta;
		}
	  #add contig and seq to hash
	  $ref{$contig}=$seq;	    
	}
#my @list=keys(%ref);
#print @list;

#print input parameters to output file
print OUTGENO "input genotypes file=$ingenos\nnumber of positions genotyped\n";
print OUTSNP "input genotypes file=$ingenos\nnumber of snps relative to the referece\n";
print OUTHET "input genotypes file=$ingenos\nnumber of heterozygote snps\n";

#read header info
my $line=<GENOS>;
my @samples=split(" ", $line);
#remove non sample information
shift(@samples);
shift(@samples);
my $n_samples=@samples;
print "number of samples=$n_samples\n";
#print headers
print "contig\tnum_snps_among_samples\n";
my $samples=join("\t",@samples);
print OUTGENO "contig\t$samples\n";
print OUTSNP "contig\t$samples\n";
print OUTHET "contig\t$samples\n";

my @genoed=(0)x$n_samples;
my $genoed;
my @hets=(0)x$n_samples;
my $hets;
my @ref_snps=(0)x$n_samples;
my $ref_snps;
my @temp_genos=();
my $pop_snps=0;
my $current_contig="";
my @contig_seq;

#process input file until EOF
while($line=<GENOS>)
	{
	  #break up line into component parts
	  my @entry=split(" ", $line);
	  
	  #if new contig
	  if ($entry[0] ne $current_contig)
		{
		  #print out last contig results, unless first contig
		  if ($current_contig ne "") 
			{
			  $genoed=join("\t",@genoed);
			  $ref_snps=join("\t",@ref_snps);
			  $hets=join("\t",@hets);
			  print OUTGENO "$current_contig\t$genoed\n";
			  print OUTSNP "$current_contig\t$ref_snps\n";
			  print OUTHET "$current_contig\t$hets\n";
			  print "$current_contig\t$pop_snps\n";
			}
		  #zero values and start with new contig information
		  $current_contig=$entry[0];
		  $pop_snps=0;
		  @contig_seq=split("",$ref{$current_contig});

		  #loop over individuals and start new counts
                  for (my $i=0; $i<$n_samples; $i++)
                        {
			  @genoed=(0) x $n_samples;
			  @ref_snps=(0) x $n_samples;
			  @hets=(0)x$n_samples;

                          if ($entry[$i+2] ne "N") 
				{
				  #increment genoed
				  $genoed[$i]++;
				  #determine if SNP relative to reference
				  if ($entry[$i+2] ne @contig_seq[$entry[1]+1]) 
					{
					  $ref_snps[$i]++;
					}
                                  #make list of genotypes, excluding Ns
                                  push (@temp_genos, $entry[$i+2]);
				}
			  if ($entry[$i+2] ne "N" && $entry[$i+2] ne "A" && $entry[$i+2] ne "T" && $entry[$i+2] ne "C" && $entry[$i+2] ne "G") {$hets[$i]++}
                        }
                  #determine if there is variation at this position within the population
                  my %string = map { $_, 1 } @temp_genos;
                        if (keys %string != 1 && @temp_genos) {$pop_snps++;}
		}
	  #if same contig
		else
		{
		  my @temp_genos=();
		  #loop over individuals and increment counts
          	  for (my $i=0; $i<$n_samples; $i++)
			{
			  if ($entry[$i+2] ne "N") 
				{
				  #increment genoed
				  $genoed[$i]+=1;
                                  #determine if SNP relative to reference
                                  if ($entry[$i+2] ne @contig_seq[$entry[1]-1])
                                        {
					  $ref_snps[$i]++;
                                        }

				  #make list of genotypes, excluding Ns
				  push (@temp_genos, $entry[$i+2]);
				}
			  if ($entry[$i+2] ne "N" && $entry[$i+2] ne "A" && $entry[$i+2] ne "T" && $entry[$i+2] ne "C" && $entry[$i+2] ne "G") {$hets[$i]++}
			  
			}	

		  #determine if there is variation at this position in any sample
		  my %string = map { $_, 1 } @temp_genos;
			if (keys %string != 1 && @temp_genos) {$pop_snps++;}
		}

	  
	}

#print out final results
$genoed=join("\t",@genoed);
$hets=join("\t",@hets);
$ref_snps=join("\t",@ref_snps);
print OUTGENO "$current_contig\t$genoed\n";
print OUTSNP "$current_contig\t$ref_snps\n";
print OUTHET "$current_contig\t$hets\n";
print "$current_contig\t$pop_snps\n";

close GENOS;
close OUTGENO;
close OUTSNP;
close OUTHET;
print "done!\n";










 












