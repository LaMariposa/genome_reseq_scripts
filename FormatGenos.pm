package FormatGenos;

use strict;
use warnings;

use Data::Dumper;
use Bio::PopGen::Individual;
use Bio::PopGen::Population;
use Bio::PopGen::Genotype;

#make hash for IUPAC codes
my(%iupac)=(
'AA' => 'A',
'CC' => 'C',
'GG' => 'G',
'TT' => 'T',
'AC' => 'M',
'CA' => 'M',
'AG' => 'R',
'GA' => 'R',
'AT' => 'W',
'TA' => 'W',
'CG' => 'S',
'GC' => 'S',
'CT' => 'Y',
'TC' => 'Y',
'GT' => 'K',
'TG' => 'K'
);

############################################################################################################
#vcf2geno
#parses vcf data into individual genotype objects
#input is a single vcf
#returns a pointer to an array of pointers to individual genotype objects for each individual in the vcf file
sub vcf2geno
	{
	 my($vcf,$contig)=@_;
	 my $line;
	 my $ft_flag;

	 #open vcf file
	 open(VCF, $vcf)||die "can't open vcf file $vcf. $!\n";	
	 
	 #ignore header lines (marked by ##) of vcf file	
	 do {$line=<VCF>;} until ($line !~ m/^##.*/);
	 	
	 #process individual information	
	 my @entry=split(" ", $line);
	 #determine number of individuals sampled
	 my $n_samples = @entry - 9;
	 #create an array of pointers, with each pointer pointing to an object for an individual
	 my @indivs=();
	 for (my $i=0; $i<$n_samples; $i++)
	 	{
		 #create object for the individual and push onto indivs array 	
		  push (@indivs, Bio::PopGen::Individual->new(-unique_id => $entry[$i+9],
								 -genotypes => []));
	 	}	 

	 #process each line until EOF
	 while($line=<VCF>) 
	 	{
		 #read each line and enter genotypes into individuals	 
		 my @entry=split(" ", $line);
		 #if not contig of interest, ignore
		 if ($entry[0] eq $contig)
			{
		 	 #ignore locus if didn't pass
			 if ($entry[6] eq "PASS")
				{
	  		 	  my $marker="$entry[0]_$entry[1]";
			 	  #enter genotype for each sample
		 		  for (my $i=0; $i<$n_samples; $i++)
		 			{
				 	  #check for ft tag (indicating some individual genotypes failed)
				 	  if ($entry[8] =~ m/FT/){$ft_flag=1}
				  	  else {$ft_flag=0}
				 	  #only make entry if genotype present and pass ft
				 	  #only make entry if genotyped (skip ./. entries)
				 	  if ($entry[$i+9] !~ m/^\..*/)#do if entry does not begin with a "."
				 		{
					 	  #only make entry if no ft or pass ft
					 	  if ($entry[8] !~ m/FT/ || $entry[$i+9] =~ m/PASS/)
							{
							 #get alleles
							 my $allele0=$entry[3];
							 my $allele1=$entry[4];
				 
							 #get indiv geno
							 my @geno=split("",$entry[$i+9]);
						 
							 #first allele
							 my $allele_a;
							 if ($geno[0]==0) {$allele_a=$allele0}
							   elsif ($geno[0]==1) {$allele_a=$allele1}
							   else {print "uh oh, something went horribly wrong with the first allele"}
				   				   
							 #second allele
							 my $allele_b;  
							 if ($geno[2]==0) {$allele_b=$allele0}
							   elsif ($geno[2]==1) {$allele_b=$allele1}
							   else {print "uh oh, something went horribly wrong with the second allele"}
					 
							 $indivs[$i]->add_Genotype(Bio::PopGen::Genotype->new(-alleles => [$allele_a, $allele_b],
											                         -marker_name => "$entry[0]_$entry[1]" ));
							}
						}					                             
			 		}	
		 		}
			}
		}
	 close VCF;
	 return \@indivs;
	}	 


######################################################################################
#createFasta
#prints out a fasta file for genotypes
#input contig name, contig size, pointer to population names, pointer to population information
#outputs a fasta file
#returns nothing
sub createFasta
        {
	  my($contig,$contig_size,$pop_names_p,$pops_p)=@_;
	
          #open output files for the contig
          open(FASTA,">$contig.samples.fasta");
          
	  #loop over the populations
          for (my $k=0; $k<@$pop_names_p; $k++)    #k tracks the population
                {
                  #loop over the individuals in the population
                  #get a list of individual in the population
                  my @inds=@$pops_p[$k]->get_Individuals();
                  my $count=@$pops_p[$k]->get_number_individuals;
                  #for each individual
                  for (my $l=0; $l<$count; $l++)        #l tracks individual within a population
                       {
                          #get information for the individual
                          my $race=@$pops_p[$k]->name;
                          my $ind_id=$inds[$l]->unique_id;
                          my @ind=@$pops_p[$k]->get_Individuals(-unique_id=>$ind_id);
                          #print header for individual
                          print FASTA ">$ind_id\_$race\n";
                          #loop over position
                          for (my $j=1; $j<=$contig_size; $j++) #j tacks position
                                {
                                  #get genotypes
                                  my @genotypes=eval{$ind[0]->get_Genotypes(-marker=>"$contig\_$j")};
                                  #if genotyped get alleles, if not genotyped print N
                                  if ($genotypes[0])
                                        {
                                          #genotyped so get alleles
                                          my @alleles=$genotypes[0]->get_Alleles();
                                          #joing the two alleles into a genotype string
                                          my $geno=join("",@alleles);
                                          print FASTA "$iupac{$geno}";
                                        }
                                        else
                                        {
                                          #not genotyped, so print N
                                          print FASTA "N";
                                        }
                                  #format width of fasta
                                  if ($j % 50 ==0 && $contig_size>$j) {print FASTA "\n";}

                                }
                          print FASTA "\n";
                        }
                }

          close FASTA;

	}
1;
