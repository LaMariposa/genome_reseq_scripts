#!/bin/bash

#gatk_genos.sh
#script to call and filter genotypes using GATK's pipeline for multi sample snp calling
#Megan Supple
#25 April 2013

#usage gatk_genos.sh <genos.input> 

#<genos.input> is a text file with entries as in example input file (see example.genos.input or example at end of this script)

#requires GATK, mutt

#produces output in current directory:
	#<race>.snps.filtered.vcf	-***filtered genotype call for downstream analysis***
	#<race>.snps.filtered.vcf.idx	-index file
	#<race>.snps.raw.metrics	-metrics file
	#<race>.snps.raw.vcf		-unfiltered genotype calls
	#<race>.snps.raw.vcf.idx	-index file
	 
#based on Best Practice Variant Detection with the GATK v3 
	#http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3
	#Better: sample-level realignment: dedup and realign
#snps called at race level
	#one script per race


###Let's begin###

#read in input file
source $1
#set output file for stdout
exec &>$race.genos.log

#print some header information
echo -e "running $0 from working directory `pwd`\n"
echo "using input parameters in $1:"
cat $1	
echo -e "\n\nstarting snp calling pipeline for $race from ${#id[@]} samples\n"

#make variable for all bams combined
snpbams=""
for ((a=0; a<${#id[@]}; a++))
 	do
		snpbams="$snpbams -I $pathbams/${id[a]}$bamsuffix"
  	done
echo "snps called from:"  	
echo $snpbams	

#hypercoverage for all samples combined (number of samples * max cov per sample)
hyper=`expr ${#id[@]} \\* $maxcov`

#call SNPs on all samples from the race
echo -e "\ncall snps"
java -jar $gatk -T UnifiedGenotyper -R $reference $snpbams -o $race.snps.raw.vcf -metrics $race.snps.raw.metrics -nt $threads -out_mode EMIT_ALL_CONFIDENT_SITES -stand_call_conf $call_conf -stand_emit_conf $emit_conf --heterozygosity $het
		#out_mode EMIT_ALL_CONFIDENT_SITES prints all confident genotype calls, not just variants
		#-Deep (> 10x coverage per sample) data: we recommend a minimum confidence score threshold of Q30 with an emission threshold of Q10. These Q10-Q30 calls will be emitted filtered out as LowQual.  (gatk website)
		#-stand_emit_conf,--standard_min_confidence_threshold_for_emitting=The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be emitted (and filtered if less than the calling threshold)	
		#-dcov defaults to 250 per sample, When running on projects with many samples at low coverage (e.g. 1000 Genomes with 4x coverage per sample) we usually lower this value to about 10 times the average coverage: '-dcov 40'. 		
		#--heterozygosity 0.025, estimated from the median value of the number of het snps/genotyped site in each sample of previous run						  

#filter SNPs	
echo -e "\nfilter snps"		  		
java -jar $gatk -T VariantFiltration -R $reference -V $race.snps.raw.vcf -o $race.snps.filtered.vcf --filterExpression "QD < 5.0" --filterName qd --filterExpression "FS > 200.0" --filterName strandbias --filterExpression "HRun > 5" --filterName homopolymer --filterExpression "DP > $hyper" --filterName hypercov --genotypeFilterExpression "GQ < $gqual" --genotypeFilterName xqual --genotypeFilterExpression "DP < $mincov" --genotypeFilterName xlowcov --genotypeFilterExpression "DP > $maxcov" --genotypeFilterName xhypercov 1> $race.filter.out
			#--clusterSize <clusterSize> --clusterWindowSize 10 
			#--filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterName HARD_TO_VALIDATE 
			#"QD < 5.0 || HRun > 5 || FS > 200.0"  for snps from gatk website
			#'QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10' (brian)
			#-et NO_ET
	       
echo | mutt -s "snp calling complete for $race" $email
echo -e "\nDONE!!!"	
	

##example input file for genotype calling (25 April 2013)
#
##race name
#race=racename
#
##path to realigned bams
#pathbams=/PATH/TO/BAMS
#
##sample ids (prefix for input bam files)
#id=( sampleA sampleB )
##suffix for bam files
#bamsuffix=.markdup.realign.fixed.reorder.bam
#
##reference that was used for alignments
#reference=/PATH/TO/REF/reference.fasta
#
##paths to programs
#picard=/PATH/TO/PICARD/DIR/picard-tools-1.53
#gatk=/PATH/TO/GATK/DIR/GenomeAnalysisTK-1.2-4-gd9ea764/GenomeAnalysisTK.jar
#
##number of threads to use
#threads=10
#
##variables
#call_conf=30
#emit_conf=10
#mincov=10
#maxcov=100
#het=0.025
#gqual=30.0
#
##email adress for notification when script completes
#email=name@domain.com

	
