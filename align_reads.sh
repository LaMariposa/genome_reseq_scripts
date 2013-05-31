#!/bin/bash

#align_reads.sh
#script to align paired end reads to a reference using bwa, prepare alignments for SNP calling, and qc alignment.
#Megan Supple
#15 April 2013

#usage align_reads.sh <align.input>

#<align.input> is a text file with entries as in example input file (see example.align.input or example at end of this script)

#requires BWA, samtools, picard, GATK, mutt

#produces output in current directory:
	#align.out	-- output from stdout and stderr
	#gatk_qc 	-- directory containing alignment assessment for each sample
	#for each sample
		#<id>.mapped.sort.bam			--  bwa mapped reads
		#<id>.mapped.sort.bam.bai
		#<id>.markdup.intervals
		#<id>.markdup.realign.fixed.reorder.bai
		#<id>.markdup.realign.fixed.reorder.bam --  ***alignments prepped for SNP calling***
		#<id>.markdups.metrics

#preparing files for SNP calling is based on: 
#Best Practice Variant Detection with the GATK v3 
	#http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3
	#Better: sample-level realignment: dedup and realign
	#no dbSNP, so cannot recalibrate
	#each sample consists of a single library.  If not, additional steps need to be added



###Let's begin###

#read in input file
source $1
#set output file for stdout
exec &>align.log


#create reference variable from ref_path and ref_file
reference=$ref_path/$ref_file


#print some header information
echo "Reference Alignment Pipeline"
echo -e "running $0 from working directory `pwd`\n"
echo "using input parameters in $1:"
cat $1
echo -e "\nBeginning Pipeline..."
echo -e "aligning reads from ${#id[@]} samples to the reference genome: $reference \n"


#need a BWA index for the reference
#check to see if one already exists (test if one of the files is there, and if so, assume the other are as well)
if (test -e $reference.bwt) 
	then echo "skipping indexing...index already exists"
	else echo "creating index"
		 bwa index $reference
fi


#align reads and prep for SNP calling
for ((a=0; a<${#id[@]}; a++))
	do
		#align reads to reference with bwa
		echo -e "\n\naligning sample ${id[a]} with bwa"
		echo "bwa aln $q -t $t -n $n -o $o -e $e -l $l -k $k $reference $files/${file1[a]} > ${id[a]}_1.sai"
			  bwa aln $q -t $t -n $n -o $o -e $e -l $l -k $k $reference $files/${file1[a]} > ${id[a]}_1.sai
		echo "bwa aln $q -t $t -n $n -o $o -e $e -l $l -k $k $reference $files/${file2[a]} > ${id[a]}_2.sai"
			  bwa aln $q -t $t -n $n -o $o -e $e -l $l -k $k $reference $files/${file2[a]} > ${id[a]}_2.sai
  		echo "bwa sampe -f ${id[a]}.sam -r "@RG\tID:${id[a]}\tSM:${sample[a]}\tPL:$platform" $reference ${id[a]}_1.sai ${id[a]}_2.sai $files/${file1[a]} $files/${file2[a]}"
			  bwa sampe -f ${id[a]}.sam -r "@RG\tID:${id[a]}\tSM:${sample[a]}\tPL:$platform" $reference ${id[a]}_1.sai ${id[a]}_2.sai $files/${file1[a]} $files/${file2[a]}
		
		#create mapped only, index, sorted bam
		echo -e "\nfiltering alignments"
		echo "samtools view -Sb -F 4 -o ${id[a]}.mapped.bam ${id[a]}.sam"
			  samtools view -Sb -F 4 -o ${id[a]}.mapped.bam ${id[a]}.sam
		echo "samtools sort ${id[a]}.mapped.bam ${id[a]}.mapped.sort"
			  samtools sort ${id[a]}.mapped.bam ${id[a]}.mapped.sort
		echo "samtools index ${id[a]}.mapped.sort.bam"
			  samtools index ${id[a]}.mapped.sort.bam

		#remove intermediate files
		echo "remove intermediate bwa files"
		rm ${id[a]}_1.sai ${id[a]}_2.sai ${id[a]}.sam ${id[a]}.mapped.bam
		
		#prep alignments with picard and gatk
		echo -e "\nprep alignment files for SNP calling"
		
		#mark duplicates for each library/sample
                echo "mark duplicates for sample ${id[a]}"
                java -jar $picard/MarkDuplicates.jar INPUT=${id[a]}.mapped.sort.bam OUTPUT=${id[a]}.markdup.bam METRICS_FILE=${id[a]}.markdups.metrics VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true

                #realign around indels
                echo -e "\nfind realign region"s
                java -jar $gatk -T RealignerTargetCreator -R $reference -I ${id[a]}.markdup.bam -o ${id[a]}.markdup.intervals
                echo -e "\nrealign"
                java -jar $gatk -T IndelRealigner -R $reference -I ${id[a]}.markdup.bam -targetIntervals ${id[a]}.markdup.intervals -o ${id[a]}.markdup.realign.bam

                #fix pair info
                echo -e "\nfix pair info"
                java -jar $picard/FixMateInformation.jar INPUT=${id[a]}.markdup.realign.bam OUTPUT=${id[a]}.markdup.realign.fixed.bam CREATE_INDEX=true SO=coordinate VALIDATION_STRINGENCY=SILENT

                #reorder
                echo -e "\nreorder"
                java -jar $picard/ReorderSam.jar INPUT=${id[a]}.markdup.realign.fixed.bam OUTPUT=${id[a]}.markdup.realign.fixed.reorder.bam REFERENCE=$reference VALIDATION_STRINGENCY=SILENT

                #index
                echo -e "\nindex"
                java -jar $picard/BuildBamIndex.jar I=${id[a]}.markdup.realign.fixed.reorder.bam VALIDATION_STRINGENCY=SILENT

		#remove intermediate files
		echo -e "\nremove intermediate picard/gatk files"
		rm ${id[a]}.markdup.bam ${id[a]}.markdup.bai ${id[a]}.markdup.realign.bam ${id[a]}.markdup.realign.bai ${id[a]}.markdup.realign.fixed.bam ${id[a]}.markdup.realign.fixed.bai

	done			  

echo | mutt -s "alignments from $1 complete" $email
echo -e "\ndone aligning and preparing for SNP calling"
echo -e "beginning quality assesment of alignments\n" 

#qc alignments			  

#create directory for qc output
mkdir gatk_qc

for ((a=0; a<${#id[@]}; a++))
        do
		echo -e "\nassessing alignments for sample ${id[a]}"
		java -jar $gatk -T FlagStat -R $reference -I ${id[a]}.markdup.realign.fixed.reorder.bam -o gatk_qc/FlagStat_${id[a]}.out
		java -jar $gatk -T CountLoci -R $reference -I ${id[a]}.markdup.realign.fixed.reorder.bam -o gatk_qc/CountLoci_${id[a]}.out
		java -jar $gatk -T CallableLoci -R $reference -I ${id[a]}.markdup.realign.fixed.reorder.bam -o gatk_qc/Callable_${id[a]}.out -summary gatk_qc/Callable_summary_${id[a]}.out -mmq 10 -mbq 20 -minDepth $mincov -maxDepth $maxcov
		java -jar $gatk -T CoarseCoverage -R $reference -I ${id[a]}.markdup.realign.fixed.reorder.bam -o gatk_qc/CoarseCoverage_${id[a]}.out --granularity 1000
	done

echo | mutt -s "alignment and qc script from $1 complete" $email
echo -e "\nDONE!!!"

#example input file to align samples to a reference (April 2013)

##input fastq files (make sure left and right match up in the order)
##path to fastq files
#files=/PATH/TO/FASTQ/DIR/
##left reads
#file1=( sampleA.R1.fastq.gz sampleB.R1.fastq.gz sampleC.R1.fastq.gz )
##right reads
#file2=( sampleA.R2.fastq.gz sampleB.R2.fastq.gz sampleC.R2.fastq.gz )
##quality score encoding
##if your reads are not phred+33/sanger/Illumina 1.8, uncomment this line so bwa uses phred+64 and converts to the standard phred+33 for alignment files
##q="-I"
#
##sam read group tags
#id=( sampleA sampleB sampleC )
#sample=( sampleA sampleB sampleC )
#platform=ILLUMINA
#
##location of reference genome
#ref_path=/PATH/TO/REFERENCE/DIR/
#ref_file=reference.fasta
#
##location of program executables (assumes BWA and samtools are already in your PATH)
#picard=/PATH/TO/PICARD/DIR/picard-tools-1.53
#gatk=/PATH/TO/GATK/DIR/GenomeAnalysisTK-1.2-4-gd9ea764/GenomeAnalysisTK.jar
#
##BWA mapping parameters (set for 100 bp reads and sanger quality scores)
##Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]
#n=8
##Maximum number of gap opens [1]
#o=2
##Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]
#e=3
##Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for .-k 2.. [inf]
#l=35
##Maximum edit distance in the seed [2]
#k=2
##number threads
#t=16
#
##gatk alignment qc parameters
##minimum and maximum coverage per position to be considered callable
#mincov=10
#maxcov=100
#
##your email address to notify you when the alignments are complete and when the script finishes
#email="name@domain.com"
#
