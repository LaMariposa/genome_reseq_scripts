#!/bin/bash

#clean_reseq4assembly.sh
#script to clean whole genome resequencing reads prior to assembly
#Megan Supple
#9 June 2013

#usage: clean_reseq4assembly.sh <clean.reseq.input>
#<clean.reseq.input> is a text file with entries as in example input file (see example.clean.reseq.input or example at the end of this script)

#requires: Trimmomatic, FLASh, JELLYFISH, Quake, FastQC.

#produces output in the current directory: 
	#clean.log			-- output from stdout and stderr
	#resource.log			-- log of computation resources used
	#*stats.txt 			-- stats from error correction
	#<sample>.hist* 		-- kmer distribution
	#*.fastqc*			-- FastQC reports
	#<sample>.pair1.fastq.gz	--**cleaned reads, pair 1**
	#<sample>.pair1.fastq.gz	--**cleaned reads, pair 2**
	#<sample>.single.fastq.gz	--**cleaned reads, unpaired**



#read in input file
source $1
#set output file for stdout
exec &>clean.log


echo "cleaning  data" > $resourcelog

#loop over each sample
for ((a=0; a<${#id[@]}; a++))
	do
		echo -e "\ncleaning sample ${id[a]}"	
		echo -e "\ncleaning sample ${id[a]}" >> $resourcelog
		
		#adapter trim, merge completely overlapping reads, and quality filtering with Trimmomatic
		echo -e "\ntrimmomatic"
		echo -e "\ntrimmomatic" >> resource.log
		if [ $score = "-phred33" ]
		  then
			/usr/bin/time -v -o $resourcelog -a java -classpath $path_trimmomatic org.usadellab.trimmomatic.TrimmomaticPE -threads $threads $score $files/${file1[a]} $files/${file2[a]} ${id[a]}.trim.p1.fastq ${id[a]}.trim.u1.fastq ${id[a]}.trim.p2.fastq ${id[a]}.trim.u2.fastq ILLUMINACLIP:$adapter:$seedMismatches:$palindromeClipThreshold:$simpleClipThreshold HEADCROP:$headCropLen SLIDINGWINDOW:$windowSize:$windowQuality LEADING:$leadQuality TRAILING:$trailQuality 
		elif [ $score = "-phred64" ]
		  then
                        /usr/bin/time -v -o $resourcelog -a java -classpath $path_trimmomatic org.usadellab.trimmomatic.TrimmomaticPE -threads $threads $score $files/${file1[a]} $files/${file2[a]} ${id[a]}.trim.p1.fastq ${id[a]}.trim.u1.fastq ${id[a]}.trim.p2.fastq ${id[a]}.trim.u2.fastq ILLUMINACLIP:$adapter:$seedMismatches:$palindromeClipThreshold:$simpleClipThreshold HEADCROP:$headCropLen SLIDINGWINDOW:$windowSize:$windowQuality LEADING:$leadQuality TRAILING:$trailQuality TOPHRED33
		else
		  echo "$score is not a valid quality score for trimmomatic. Should be -phred33 or -phred64"; exit 1
		fi  

		#merge overlapping reads with FLASH
		echo -e "\nflash"
		echo -e "\nflash" >> $resourcelog
		/usr/bin/time -v -o $resourcelog -a flash -o ${id[a]} -t $threads -m $minOverlap -x $maxMismatchDensity -p 33 -M $maxOverlap ${id[a]}.trim.p1.fastq ${id[a]}.trim.p2.fastq

		#error correct with jellyfish and quake
                #use jellyfish to count kmers
		echo -e "\njellyfish"
		echo -e "\njellyfish" >> $resourcelog
                /usr/bin/time -v -o $resourcelog -a jellyfish count -m $quakeK -q -C -t $threads -s 2000000000 --quality-start 33 --stats=Stats -o ${id[a]}.jelly ${id[a]}.notCombined_1.fastq ${id[a]}.notCombined_2.fastq ${id[a]}.trim.u1.fastq ${id[a]}.trim.u2.fastq ${id[a]}.extendedFrags.fastq

		#check if need to merge
		if ( test -e ${id[a]}.jelly_1 )
			then /usr/bin/time -v -o $resourcelog -a jellyfish qmerge -m $quakeK -s 4000000000 -o ${id[a]}.jf.output ${id[a]}.jelly_*
			else mv ${id[a]}.jelly_0 ${id[a]}.jf.output
		fi
                
		#create counts file
		/usr/bin/time -v -o $resourcelog -a jellyfish qdump -c -t -o ${id[a]}.counts ${id[a]}.jf.output

                #make a file list to input into correct
                touch file.list
                echo "${id[a]}.notCombined_1.fastq ${id[a]}.notCombined_2.fastq" >> file.list
                echo "${id[a]}.trim.u1.fastq" >> file.list
		echo "${id[a]}.trim.u2.fastq" >> file.list
		echo "${id[a]}.extendedFrags.fastq" >> file.list

                #use quake to correct the reads
		echo -e "\nquake"
		echo -e "\nquake" >> $resourcelog
                /usr/bin/time -v -o $resourcelog -a $path_quake/bin/correct -f file.list -k $quakeK -c $cov_cut -l $minlen -m ${id[a]}.counts -p $threads -q 33

		#merge single end fastqs
		echo -e "\ncat"
		echo -e "\ncat" >> $resourcelog
		/usr/bin/time -v -o $resourcelog -a cat ${id[a]}.extendedFrags.cor.fastq ${id[a]}.trim.u2.cor.fastq ${id[a]}.trim.u1.cor.fastq ${id[a]}.notCombined_1.cor_single.fastq ${id[a]}.notCombined_2.cor_single.fastq > ${id[a]}.single.fastq
		
		#rename final paired files
		mv ${id[a]}.notCombined_1.cor.fastq ${id[a]}.pair1.fastq
		mv ${id[a]}.notCombined_2.cor.fastq ${id[a]}.pair2.fastq 	


		#FastQC final files
		echo -e "\nFastQC"
		fastqc -t $threads ${id[a]}.pair1.fastq
		fastqc -t $threads ${id[a]}.pair2.fastq
		fastqc -t $threads ${id[a]}.single.fastq
		
		
		#gzip final files and remove intermediates
		gzip ${id[a]}.pair1.fastq &
		gzip ${id[a]}.pair2.fastq &
		gzip ${id[a]}.single.fastq &

		rm ${id[a]}.trim.u1.fastq ${id[a]}.trim.u2.fastq ${id[a]}.trim.p1.fastq ${id[a]}.trim.p2.fastq ${id[a]}.notCombined_1.fastq ${id[a]}.notCombined_2.fastq ${id[a]}.extendedFrags.fastq ${id[a]}.notCombined_1.cor_single.fastq ${id[a]}.notCombined_2.cor_single.fastq ${id[a]}.trim.u1.cor.fastq ${id[a]}.trim.u2.cor.fastq ${id[a]}.extendedFrags.cor.fastq ${id[a]}.counts error_model.${id[a]}.* ${id[a]}.jf.output file.list Stats

		echo -e "\ndone cleaning sample ${id[a]}\n\n\n" 
 		
	done

echo DONE!!!

echo | mutt -s "reseq_assembly_clean.sh script from $1 complete" $email




##example input file to clean resequencing reads prior to assembly (9 June 2013)
#
##input fastq files (make sure left and right match up in the order)
##path to paired end files
#files=/PATH/TO/FASTQ/FILES
##files
#file1=( sampleA.read1.gz sampleB.read1.gz )
#file2=( sampleB.read2.gz sampleB.read2.gz )
##id
#id=( sampleA sampleB )
#
##general parameter
#minlen=30       #minimum read length to output in final file
#
##trimmomatic parameters
#adapter=/PATH/TO/ADAPTER/FILE/adapter.fasta
#score=-phred33  #-phred33 or -phred64
#headCropLen=0
#seedMismatches=2
#palindromeClipThreshold=30
#simpleClipThreshold=12
#windowSize=4
#windowQuality=20
#leadQuality=10
#trailQuality=10
#path_trimmomatic=/PATH/TO/trimmomatic.jar
#
##flash parameters
#minOverlap=15
#maxOverlap=70 #or calc from readlen, frag len, sd [70 is default for 100bp reads, 180bp frag, 18sd]
#maxMismatchDensity=0.2
#
##Quake/Jellyfish parameters
#quakeK=19
#cov_cut=1
#path_quake=/PATH/TO/QUAKE/DIR/
#
##misc system and log files
#threads=6
#resourcelog=resource.log
#
##email for notification
#email=name@domain.com
