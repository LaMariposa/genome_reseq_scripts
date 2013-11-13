#!/bin/bash

#assemble_samples.sh
#script to assemble a set of samples with ABySS
#Megan Supple
#24 August 2013
#last modified 13 Nov 2013

#usage assemble_samples.sh <assemble.input>
	#<assemble.input is a text file file with entries as in example input file (see example.assemble.input or example at the end of this script)
#requires ABySS
#produces ABySS output for each sample, a lot file, and a file with resource usage

#read in input file
source $1
#set output file for stdout
exec &>assemble.log

#loop over samples
for ((a=0; a<${#id[@]}; a++))
        do
		#make assembly name
		name=${id[a]}\_k$kmer\_b$bubble
                echo -e "\n\n\n\nassembly=$name sample=${id[a]}  kmer=$kmer  bubble=$bubble"

                #assemble
                /usr/bin/time -v -o $resourcelog -a abyss-pe name=$name k=$kmer b=$bubble \
                                                              in="$files/${id[a]}.pair1.fastq.gz $files/${id[a]}.pair2.fastq.gz" \
                                                              se=$files/${id[a]}.single.fastq.gz
	done

echo | mutt -s "assemblies complete" $email
echo -e "\nDONE!!!"


##example input file to assemble a set of samples
#
##path to cleaned fastq files (it will look for 3 files for each id, named <id>.pair1.fastq.gz, <id>.pair2.fastq.gz <id>.single.fastq.gz)
#files=/PATH/TO/FASTQS/
#id=( sampleA sampleB sampleC )
#
##assembly parameters
#kmer=<int>
#bubble=<int>
#
##email adress for notification when script completes
#email=name@domain.com
#
##file to record resource usage
#resourcelog=resource.log

