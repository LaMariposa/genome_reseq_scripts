#!/usr/bin/perl

#sample_consensus.pl by Megan Supple 07May2014
#script to generate a consensus for specified samples from an aligned fasta file
#usage sample_consensus.pl input.fasta sample.list
	#sample.list is a space delimited file with one sample per row and 
	#col1=population, col2=phenotype, col3=sampleID
	#population(col1) will be ignored, but is used for consistency with other scripts (eg fasta2popgen.pl)

use lib $ENV{PERL5LIB};
use strict;
use warnings;

use Getopt::Long; 
use List::MoreUtils qw(uniq);
use Bio::AlignIO;
use Bio::SeqIO;

#read in command line arguments
my ($infasta, $insamples)=@ARGV;

my $usage = "Usage: sample_consensus.pl <infile.fasta> <sample.list>
options:
<infile.fasta> an input aligned fasta file
<sample.list> is a space delimited file with one sample per row and
        	col1=population(ignored), col2=phenotype, col3=sampleID
";

die "$usage" unless (@ARGV == 2);

my @phenos;  	#list of phenos
my %samples;	#has of sample and associated phenotype

#read in sample and phenotype information
	open(PHENO, $insamples)||die "can't open input sample file. $!\n";
	#parse to read in phenotypes
	while(my $line=<PHENO>)
		{
		  chomp($line);
		  my $sample=(split /\ /, $line)[2];
		  my $pheno=(split /\ /, $line)[1];
		  $samples{$sample}=$pheno;
		}

#generate a fasta for each phenotype
#convert the Ns in the fasta file to ?, so they will be ignored when generating a consensus
        my %handles;
        my $current_pheno;

	#open fasta file and new out fasta file
	@phenos=values %samples;	
	my @pheno_list=uniq @phenos;

	#for each pheno, make an output file
	for (@pheno_list)
		{
		  open (my $fh, ">", "$infasta.$_") || die $!;
		  $handles{$_}=$fh;
		}

	#parse infasta
	open(IN, $infasta)||die "can't open input fasta file. $!\n";
	while(my $line=<IN>)
		{
		  #if id line, determine pheno and write to outfile
		  if ($line=~m/^>.*/) 
			{
			  my $id=substr $line,1,-1;
			  $current_pheno=$samples{$id};
			  print {$handles{$current_pheno}} $line;
			}
		    #if line with no Ns, print line
		    elsif ($line!~/N/) {print {$handles{$current_pheno}} $line;}
		    else
			{
			  $line=~ s/N/?/g;
			  print {$handles{$current_pheno}} $line;
			}
		}
	
	close IN;
	#close pheno fastas
	for (@pheno_list)
		{
		  my $fh=$handles{$_};
		  close $fh;
		}


#open the output file
#determine base file name to use for output file
my $basename=(split /\./, $infasta)[-2];
my $outfile=$basename.".consensus.fasta";
my $out=Bio::SeqIO->new(-file => ">$outfile",
                        -format => 'fasta');


#for each phenotype, generate a consensus
	foreach my $ph (@pheno_list)
		{
		  my $infasta_pheno=$infasta . "." . $ph;
		  
		  #read in the aligned fasta
		  my $str=Bio::AlignIO->new(-file => "$infasta_pheno", 
			  		    -format => 'fasta');
		  my $aln=$str->next_aln();

		  #generate the consensus
		  my $con_aln=$aln->consensus_iupac();

		  #print to the output file
		  my $header=$basename . "_consensus_" . $ph;
		  my $seq=Bio::Seq->new(-seq => "$con_aln",
		   		        -display_id => $header);
		  $out->write_seq($seq);
}

