#!/bin/bash

#scaffold_reseq.part1.sh
#script part 1 to scaffold using Bambus2 with paired reads and a reference genome
#Megan Supple
#25 Oct 2013
#last modified 12 Nov 2013

#usage scaffold_reseq.part1.sh <scaffold.input>
	#<scaffold.input> is a text files with entries as in example in put file (see example.scaffold.input or example at the end of this script)
#requires BWA, samtools, ABySS's abyss-samtoafg.pl, AMOS (dev version), FASTX-Toolkit, MUMmer, other scripts in my genome_reseq_scripts repository
#produces: 
	#<id>.contig		TIGR .contig file
	#<id>.reads.fasta	single file with all sequencing reads in fasta format
	#<id>.evidence.xml	trace arcive xml file of links
	#scaffold.log		log file

#read in input file
source $1
#set output file for stdout
exec &>scaffold.log

#generate tigr .contig file
	#align reads to assembly
	#index assembly as reference
	echo -e "\n\nindexing reference"
	bwa index $assembly
			
	#align paired reads
	echo -e "\n\naligning paired reads"
	bwa aln -t $threads -f reads_1.sai $assembly $read1
	bwa aln -t $threads -f reads_2.sai $assembly $read2
	bwa sampe -f pe.sam -r "@RG\tID:$id\tSM:$id\tPL:Illumina" $assembly reads_1.sai reads_2.sai $read1 $read2
	samtools view -Sb -o pe.bam pe.sam
	samtools sort pe.bam pe.sort

	#align single reads
	echo -e "\n\naligning orphaned reads"
	bwa aln -t $threads -f reads_3.sai $assembly $readsingle
	bwa samse -f se.sam -r "@RG\tID:$id\tSM:$id\tPL:$platform" $assembly reads_3.sai $readsingle
	samtools view -Sb -o se.bam se.sam 
	samtools sort se.bam se.sort 

	#merge single and paired
	echo -e "\n\nmerging paired and orphaned alignments"
	samtools merge $id.bam pe.sort.bam se.sort.bam
	#create sam file
	echo -e "\n\nconverting bam to sam"
	samtools view -h $id.bam > $id.sam

	#create afg file from sam file
	echo -e "\n\ncreating afg file"
	$path_samtoafg/abyss-samtoafg.pl $assembly $id.sam > $id.afg
	
	#create amos bank with afg file
	echo -e "\n\ncreating amos bank"
	$path_amos/bank-transact -cb $id.bnk -m $id.afg
	#copy bank, just in case
	cp -r $id.bnk $id.bnk.COPY

	#create .contig file from bank
	echo -e "\n\ncreating .contig file"
	$path_amos/bank2contig $id.bnk > $id.contig


	#if .contig file is not empty, remove intermediate files
	if [ -s $id.contig ]
	then
		echo "\n\nremoving intermediate files"
		rm $assembly.rsa $assembly.sa $assembly.rbwt $assembly.bwt 
		rm $assembly.rpac $assembly.amb $assembly.ann $assembly.pac
		rm *.sai *.sam *.bam $id.afg 
		rm -r $id.bnk.COPY $id.bnk
	else
		echo "error: .contig file is empty"
	fi


#convert fastqs to a single fasta file
	echo -e "\n\nconvert fastqs to fasta"
	#unzip the reads and convert to fasta and rezip
	gunzip $read1
	fastq_to_fasta -n -Q33 -i ${read1%.gz} -o ${read1%fastq.gz}fasta
	gzip ${read1%.gz}
	
	gunzip $read2
        fastq_to_fasta -n -Q33 -i ${read2%.gz} -o ${read2%fastq.gz}fasta
        gzip ${read2%.gz}

	gunzip $readsingle
	fastq_to_fasta -n -Q33 -i ${readsingle%.gz} -o ${readsingle%fastq.gz}fasta
	gzip ${readsingle%.gz}

	#cat the fasta files
	cat ${read1%fastq.gz}fasta ${read2%fastq.gz}fasta ${readsingle%fastq.gz}fasta > $id.reads.fasta


#generate xml file
	echo -e "\n\ngenerating xml file"
	mydate=$(date)

	#generate file header
	echo -e "<?xml version=\"1.0\" ?>\n" > temp.xml
	echo -e "<EVIDENCE ID=\"$id\"
	          DATE=\"$mydate\"
	          PROJECT=\"$assembly\"
	          PARAMETERS=\"$read1,$read2,$ref\"  >" >> temp.xml

	#make xml pair information
	$path_grs/readpairxml.pl ${read1%fastq.gz}fasta

        #remove the separate fasta files
        rm ${read1%fastq.gz}fasta ${read2%fastq.gz}fasta ${readsingle%fastq.gz}fasta


	#contig information and alignment to reference
	#align to reference
	echo -e "\n\naligning to the reference"
	nucmer $ref $assembly &> nucmer.out
	#convert out.delta to XML for bambus
	show-tiling -x out.delta > mum.xml
	#add mum.xml to current evidence xml
	tail -n+6 mum.xml | head -n-1 > temp.mum
	cat temp.xml temp.pairs.xml temp.mum > $id.evidence.xml
	

	#generate final tag
	echo "</EVIDENCE>" >> $id.evidence.xml

	#remove intermediate files
	rm temp.xml temp.pairs.xml nucmer.out out.delta mum.xml temp.mum

echo | mutt -s "scaffold_reseq.sh script from $1 complete" $email
echo -e "\n\nDONE!!!"


##example input files to scaffold resequencing assemblies
#
##sample id
#id=sampleA
#
##initial assembly
#assembly=/PATH/TO/ASSEMBLY/sampleA-contigs.fa
#
##cleaned, fixed sequencing reads
#read1=/PATH/TO/CLEANED/READS/$id.pair1.fixed.fastq.gz
#read2=/PATH/TO/CLEANED/READS/$id.pair2.fixed.fastq.gz
#readsingle=/PATH/TO/CLEANED/READS/PATH/TO/CLEANED/READS/$id.single.fixed.fastq.gz
#
##reference genome
#ref=/PATH/TO/REFERENCE/DIR/reference.fasta
#
##number of threads
#threads=6
#
##path to non-standard programs
##path to abyss-samtoafg.pl
#path_samtoafg=/PATH/TO/SAMTOAFG/DIR
##path to bin of dev version of amos
#path_amos=/PATH/TO/AMOS/amos/bin
##path to my genome_reseq_scripts repository
#path_grs=/PATH/TO/REPOSITORY/genome_reseq_scripts
#
##email address for mutt notification
#email="name@domain.com"
