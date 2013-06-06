#!/bin/bash

#bam2blastid.sh
#script to call genotypes from bam files individually and blast a specified region to a database
#Megan Supple
#1 June 2013

#usage bam2blastid.sh <blastid.input>
#<blastid.input> is a text file with input bams and blast database information (see example.blastid.input or example at the end of this script)
#requires gatk_genos.sh, vcf2format.pl, fasta_subregion.pl, blast, mutt, and all requirements within these scripts/programs
#output
	#bam2blastid.log				stdout/stderr
	#blastreport.out				blast report (tabular format)
	#<ref_contig>.samples.fasta 			fasta sequences for all samples for the entire reference contig
	#<ref_contig>.samples.fastaStart-stop.fasta	fasta sequences for all samples for the entire reference contig
	#<sample>.snps.filtered.vcf			genotypes relative to the entire reference sequence
	#topblasthits.out				list of the single best blast hit for each sample

#read in input file
source $1
#set output file for stdout
exec &>bam2blastid.log

#create an empty file to list the vcfs
>vcf.list

#create reference variable from ref_path and ref_file
db=$db_path/$db_file

#for each sample
for ((a=0; a<${#ids[@]}; a++))
	do
		echo "analyzing sample ${ids[a]}"
		
		#prep input file for SNP calling by appending new info onto original file
		echo "id=( ${ids[a]} )" > temp.txt
		echo "race=( ${ids[a]} )" >> temp.txt
		cat $1 temp.txt > temp.in
		
		#call genotypes
		gatk_genos.sh temp.in
		
		#generate infile for script to convert vcf to fasta
                echo "${pre_taxa_id[a]} ${ids[a]}.snps.filtered.vcf" >> vcf.list

		#remove a bunch of useless files
		rm ${ids[a]}.filter.out ${ids[a]}.snps.filtered.vcf.idx ${ids[a]}.snps.raw.vcf.idx ${ids[a]}.snps.raw.metrics ${ids[a]}.snps.raw.vcf
		rm temp.txt temp.in
	done
	
	#generate fasta
	#generate contig fasta
	vcf2format.pl -contig $contig -fasta vcf.list
	rm vcf.list
	#pull region of interest out of fasta
	fasta_subregion.pl $contig.samples.fasta $start $end	

	#blast
	#format blast db
	#check if already formated
	if [ $db_type = "nucleotide" ] 
	  then
		if (test -e $db.nhr)
			then echo "skipping formatting $db_type database..already done"
			else echo "formatting $db_type database"
				formatdb -i $db -p F
		fi
	elif [ $db_type = "protein" ] 
	  then
		if (test -e $db.phr)
                       then echo "skipping formatting $db_type database..already done"
                       else echo "formatting $db_type database"
                                formatdb -i $db -p T
		fi
	else 
	  #then 
	  echo "$db_type is not a valid database type.  Specify a "nucleotide" or "protein" type"
	fi

	#blast
	blastall -p $blast_type -d $db -i $contig.samples.fasta$start-$end.fasta -o blastreport.out -m8 -b1 -e $eval
 
	#filter blast output for top hit only for each sample
	sort -u -k1,1 blastreport.out > topblasthits.out

echo | mutt -s "bam2blastid script from $1 complete" $email2
echo DONE!!!


##example input file for bam2blastid.sh (1 June 2013)
#
##input bams
#pathbams=/PATH/TO/BAMS                          #path to realigned bams
#bamsuffix=.markdup.realign.fixed.reorder.bam    #suffix for bam files
#
##sample ids (prefix for input bam files)
#ids=( sampleA sampleB sampleC )
#
##initial species id for each sample
#pre_taxa_id=( speciesA speciesB speciesC )
#
##minimum e-value to report blast hits
#eval=0.001
#
##blast db
#db_path=/PATH/TO/BLAST/DB               #path to blast db
#db_file=blast.db.fa                     #blast db fasta file
#db_type=nucleotide                      #type of blast database (protein or nucleotide)
#
##type of blast search (eg blastn, blastp, blastx, tblastn, tblastx)
#blast_type=blastn
#
##genomic region of interest
#contig="contig_name"            #contig in reference with blast piece
#start=100                       #start position of blast piece
#end=200                         #end position of blast piece
#
##reference that was used for alignments
#reference=/PATH/TO/REF/reference.fasta
#
##paths to programs
#picard=/PATH/TO/PICARD/picard-tools-1.53
#gatk=/PATH/TO/GATK/GenomeAnalysisTK-1.2-4-gd9ea764/GenomeAnalysisTK.jar
#
##number of threads to use
#threads=10
#
##GATK genootyping variables
#call_conf=30
#emit_conf=10
#mincov=10
#maxcov=251      #GATK defaults to downsample every site to <251, so setting to 251 keeps everything
#het=0.025
#gqual=30.0
#
##email adress for notification when script completes
#email2=name@domain.com

