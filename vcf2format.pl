#!/usr/bin/perl

#vcf2format.pl by Megan Supple 27 May 2013
#last updated 8 Nov 2013
#script to convert vcf file(s) into various formats
	#option to specifiy a single contig
	#use IUPAC ambiguity codes
	#possible output formats
		#fasta
		#matric

#usage:    vcf2format.pl [options] <infile>
	#<infile> is a text file that lists each race and vcf files, each file on a single line
	#race1 race1.vcf
        #race2 race2.vcf
	
use strict;
use warnings;

use lib $ENV{PERL5LIB};
use FormatGenos;
use Bio::PopGen::Population;
use Getopt::Long;
use Data::Dumper;
	

#read in command line arguments
my $contig;	#contig of interest
my $fasta;	#generate a fasta file
my $matrix;	#generates a matrix file

GetOptions (	"contig=s" => \$contig,
		"fasta" => \$fasta,
		"matrix" => \$matrix
  	);
		
my $usage = "Usage: vcf2format.pl [options] <infile>
options:
	<infile>		an infile that lists each race and vcf file, each file on a single line 
	-contig <string>	specify a single contig
	-fasta			generate fasta file
	-matrix			generate a matrix file
";

die "$usage" unless (@ARGV == 1);
my ($infile) = @ARGV;


#declare some variables
my @pop_names;                  #array with names of all pops
my @pops=();                    #data structure with genotype information


#open and read in input file
open(INFILE, $infile)||die "can't open input file. $!\n";
my @input=<INFILE>;
close INFILE;
	
#check and see if any output format is specified
($fasta || $matrix)||die "no output format specified\n";

#create a hash of contig names and sizes from the header of the first vcf
print "gathering contig information\n";
#open vcf file
my @info=split(" ",$input[0]);
open(VCF, $info[1])||die "can't open vcf file. $!\n";
#ignore header lines until contig information
my $line;
do {$line=<VCF>;} until ($line =~ m/^##contig=<ID=.*/);
#while contig information, record lengths
my %contigs;
while ($line =~ m/^##contig=<ID=.*/)
        {
         my @contig_info=split("=",$line);
         my @contig_id=split(",",$contig_info[2]);
         my @contig_size=split(">",$contig_info[3]);
         #add to hash
         $contigs{$contig_id[0]} = $contig_size[0];
         $line = <VCF>;
        }
close VCF;


#if an input contig is specified, remove other contigs from the contig hash
if ($contig)
	{
	  my $size=$contigs{$contig};
       		($size)||die "specified contig does not exist in first vcf header\n";
	  for (keys %contigs) { delete $contigs{$_}; }
	  $contigs{$contig}=$size;
	}
 
#get list of contigs
my @contig_list=keys(%contigs);


#get genos and generate output for each contig in turn (to reduce memory usage)
while (my ($contig, $size)=each(%contigs))
	{
	  print "processing contig=$contig ($size bp)\n";

	  #process each vcf
	  for (my $i=0; $i<@input; $i++)
        	{
         	  #get race and file name from input
         	  my @info=split(" ",$input[$i]);

         	  #call subroutine to process the vcf file
         	  print "processing $info[1]\n";
         	  my $indivs_p=FormatGenos::vcf2geno($info[1],$contig);

         	  #create the population (if it is new) and push into @pops
         	  print "creating population\n";
         	  my $pop_flag=@pop_names; #set flag for new population
         	  
		  #check if new population, if old population--reset pop flag
         	  for (my $j=0; $j<@pop_names; $j++)
                	{
                 	  if ($info[0] eq $pop_names[$j]) {$pop_flag=$j}
                	}

         	  #if new population, add population object and add name to pop_names array
         	  if ($pop_flag==@pop_names)
                	{
                	  push(@pops, Bio::PopGen::Population->new(-name => $info[0]));
                 	  push(@pop_names, $info[0]);
                	}

         	  #add each individual to the population (get number of individuals processed from size of array pointed to, pop from pop_flag)
	 	  print "adding individuals to the population\n";
         	  for (my $k=0; $k<@$indivs_p; $k++)
                	{
		 	  #add individual
                 	  $pops[$pop_flag]->add_Individual($$indivs_p[$k]);
                	}
        	}	


	  #get a list of individuals
	  my @individuals=();
	  #for each population
	  for (my $i=0; $i<@pop_names; $i++)
		{
		  my @inds=$pops[$i]->get_Individuals();
		  #for each individual
		  my $count=$pops[$i]->get_number_individuals;
		  for (my $j=0; $j<$count; $j++)
			{
			  my $ind=$inds[$j]->unique_id;
			  push (@individuals, $ind);
			}	
		}


	  #generate fasta files
	  if ($fasta)
		{
		  print "generating fasta filess\n";
		  #call module to generate the fasta file for the contig
        	  FormatGenos::createFasta($contig,$size,\@pop_names,\@pops);
		}
          #generate matrix files
          if ($matrix)
                {
                  print "generating matrix files\n";
                  #call module to generate the matrix file for the contig
		  FormatGenos::createMatrix($contig,$size,\@pop_names,\@pops);
                }
	  #clear data from the contig
	  @pops=(); @pop_names=();	
	}


print "DONE!!!\n";

exit;

