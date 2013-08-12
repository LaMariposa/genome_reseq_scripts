#!/bin/bash

#assemble_params.sh
#script to assemble multiple samples using a variety of parameters 
#Megan Supple
#27 June 2013

#usage: assemble_params.sh
#input information (files, parameters, etc.) is located in the header of the script 

#requires:  ABySS

#produces: a fasta file for each sample for each parameter combination.

#input information
#path to cleaned fastq files (it will look for 3 files for each sample id named: <id>.pair1.fastq.gz, <id>.pair2.fastq.gz, <id>.single.fastq.gz)
files=/PATH/TO/CLEANED/FASTQS
id=( sampleA sampleB sampleC )

#assembly parameters
#kmer
kmer=( 21 31 41 51 61 )
#bubble=bubbleX * kmer
bubbleX=( 3 4 5 )   

email=name@domain.com
resourcelog=resource.log



#begin script
mkdir workingdir
cd workingdir

#loop over samples
for ((a=0; a<${#id[@]}; a++))
	do
		#loop over parameters 
		#kmer
		for ((k=0; k<${#kmer[@]}; k++))
			do

				#assemble with default bubble
				bubble=default
				name=${id[a]}\_k${kmer[k]}_b$bubble
                                echo -e "\n\nassembly=$name sample=${id[a]}  kmer=${kmer[k]}  bubble=$bubble"
				/usr/bin/time -v -o ../$resourcelog -a abyss-pe name=$name k=${kmer[k]} \
					in="$files/${id[a]}.pair1.fastq.gz $files/${id[a]}.pair2.fastq.gz" \
					se=$files/${id[a]}.single.fastq.gz
				cp $name-contigs.fa ..
				cp $name-stats ..
				rm *			

				#vary bubble
				for ((b=0; b<${#bubbleX[@]}; b++))
					do

						#calculate bubble parameter
						bubble=`expr ${kmer[k]} \\* ${bubbleX[b]}`
						name=${id[a]}\_k${kmer[k]}_b$bubble
						echo -e "\n\nassembly=$name sample=${id[a]}  kmer=${kmer[k]}  bubble=$bubble"

						#assemble
                                		 /usr/bin/time -v -o ../$resourcelog -a abyss-pe name=$name k=${kmer[k]} b=$bubble \
							in="$files/${id[a]}.pair1.fastq.gz $files/${id[a]}.pair2.fastq.gz" \
							se=$files/${id[a]}.single.fastq.gz
						cp $name-contigs.fa ..
                                		cp $name-stats ..
                                		rm *

					done
			done
	done

cd ..
rmdir workingdir	
echo | mutt -s "assemblies complete" $email
echo -e "\nDONE!!!"	
