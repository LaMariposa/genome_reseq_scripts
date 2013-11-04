#!/bin/bash

#scaffold_reseq.sh
#script to scaffold using Bambus2 with paired reads and a reference genome
#Megan Supple
#25 Oct 2013

#usage
#requires
#produces


#read in input file
source $1
#set output file for stdout
exec &>scaffold.log

#generate tigr .contig file
	#align reads to assembly
	#index assembly as reference
	bwa index $assembly
			
	#align paired reads
	bwa aln -t $threads -f reads_1.sai $assembly $read1
	bwa aln -t $threads -f reads_2.sai $assembly $read2
	bwa sampe -f pe.sam -r "@RG\tID:$id\tSM:$id\tPL:Illumina" $assembly reads_1.sai reads_2.sai $read1 $read2
	samtools view -Sb -o pe.bam pe.sam
	samtools sort pe.bam pe.sort

	#align single reads
	bwa aln -t $threads -f reads_3.sai $assembly $readsingle
	bwa samse -f se.sam -r "@RG\tID:$id\tSM:$id\tPL:$platform" $assembly reads_3.sai $readsingle
	samtools view -Sb -o se.bam se.sam 
	samtools sort se.bam se.sort 

	#merge single and paired
	samtools merge $id.bam pe.sort.bam se.sort.bam
	#create sam file
	samtools view -h $id.bam > $id.sam

	#create afg file from sam file
	$path_samtoafg/abyss-samtoafg.pl $assembly $id.sam > $id.afg
	
	#create amos bank with afg file
	$path_amos/bank-transact -cb $id.bnk -m $id.afg
	#copy bank, just in case
	cp -r $id.bnk $id.bnk.COPY

	#create .contig file from bank
	$path_amos/bank2contig $id.bnk > $id.contig


	#if .contig file is not empty, remove intermediate files
	if [ -s $id.contig ]
	then
		echo "removing intermediate files"
		rm *.rsa *.sa *.rbwt *.bwt *.rpac *.amb *.ann *.pac
		rm *.sai *.sam *.bam $id.afg 
		rm -r $id.bnk.COPY
	else
		echo "error: .contig file is empty"
	fi


#convert fastqs to a single fasta file
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
	mydate=$(date)

	#generate file header
	echo -e "<?xml version=\"1.0\" ?>\n" > temp.xml
	echo -e "<EVIDENCE ID=\"$id\"
	          DATE=\"$mydate\"
	          PROJECT=\"$assembly\"
	          PARAMETERS=\"$read1,$read2,$ref\"  >" >> temp.xml

	#make xml pair information
	$path_grs/readpairxml.pl ${read1%fastq.gz}fasta

        #remove the fasta files
        rm ${read1%fastq.gz}fasta ${read2%fastq.gz}fasta ${readsingle%fastq.gz}fasta


	#contig information and alignment to reference
	#align to reference
	nucmer $ref $assembly &> nucmer.out
	#convert out.delta to XML for bambus
	show-tiling -x out.delta > mum.xml
	#add mum.xml to current evidence xml
	tail -n+6 mum.xml | head -n-1 > temp.mum
	cat temp.xml temp.pairs.xml temp.mum > $id.evidence.xml
	

	#generate final tag
	echo "</EVIDENCE>" >> $id.evidence.xml

#scaffold
	$path_amos/toAmos -o $id.afg -s $id.reads.fasta -c $id.contig -x $id.evidence.xml
	$path_amos/bank-transact -cb $id.full.bnk -m $id.afg
	cp -r $id.full.bnk $id.full.bnk.COPY
	$path_amos/goBambus2 $id.full.bnk $id --all	

echo | mutt -s "scaffold_reseq.sh script from $1 complete" $email
echo DONE!!!

