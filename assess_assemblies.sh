#!/bin/bash

#assess_assemblies.sh
#script to evaluate a set of assemblies with a variety of metrics
#Megan Supple
#29 June 2013

#usage assess_assemblies.sh <assess.input>
#<assess.input> is a text file with entries as in example input file (see example.assess.input or example at the end of this script)

#requires  assemblathon2-analysis git repository, cegma, bwa, samtools, amos 

#produces 
	#<sample>.assessment.summary		**a summary file of key metrics**
	#<sample>*.stats			a variety of summaries
	#<sample>.basic.summary			assemblathon_stats.pl output
	#<assembly>.mapped.sort.bam/bai		an alignment of the assembly to a reference for visual inspection
	

#read in input file
source $1
#set output file for stdout
exec &>assess.log


#create empty files to list the output
>$sample.basic.summary
>general.out
>compass.out
>cegma.complete.out
>cegma.partial.out
>amos.out


echo -e "beginning assessment"
for ((a=0; a<${#fasta[@]}; a++))
	do
		echo -e "\nassessing ${fasta[a]}"
		
		#basic stats with assemblathon script assemblathon_stats.pl 
		echo "assemblathon stats"
		$thon_path/assemblathon_stats.pl -csv -graph -genome_size 400000000 $path2fastas/${fasta[a]} >> $sample.basic.summary 
		#add stats to output file
		tail -n1 ${fasta[a]%fa}csv >> general.out
	
		#assess versus BAC sequences with assemblathon script compass_modified.pl 
		echo "compass"
		$thon_path/compass_modified.pl -s /usr/bin/samtools -n 2000000 $ref $path2fastas/${fasta[a]} > compass.temp
		#add compass results to output file
		grep -A3 "By-row output:" compass.temp | tail -n1 >> compass.out
	
		#cegma
		echo "cegma"
		#fix headers so not just numbers
		awk '/^>/{$0=">contig"++i}1' $path2fastas/${fasta[a]} > temp.fa
		cegma -T $threads -g temp.fa -p $cegma > cegma.out
		#parse output
		echo -n ${fasta[a]} >> cegma.complete.out
                echo -n ${fasta[a]} >> cegma.partial.out
		grep "Complete" output.completeness_report | head -n2 | tail -n1 >> cegma.complete.out
		grep "Partial" output.completeness_report | head -n2 | tail -n1 >> cegma.partial.out

                #remove intermediate files
                rm ${fasta[a]%fa}scaffold.NG50.csv
                chmod 666 compassRun.TEMP.ref
	        rm compassRun.*
                rm temp.fa output.cegma* output.completeness_report cegma.out
	done

echo -e "\nsummarizing output"
#add header to the assemblathon stats
head -n1 ${fasta[0]%fa}csv > head.txt
cat head.txt general.out > $sample.general.stats
rm general.out head.txt
#add header to the compass stats
grep -A2 "By-row output:" compass.temp | tail -n1 > head.txt
cat head.txt compass.out > $sample.compass.stats
rm compass.out head.txt compass.temp
#add header to cegma stats
echo "assembly part/complete #Prots complete%Completeness - complete#Total completeAveragePerKOG complete%Ortho" > head.txt
cat head.txt cegma.complete.out > $sample.cegma.complete.stats
echo "assembly part/complete #Prots part%Completeness - part#Total partAveragePerKOG part%Ortho" > head.txt
cat head.txt cegma.partial.out > $sample.cegma.partial.stats
rm cegma.complete.out cegma.partial.out head.txt



echo -e "\n\nassessing without a reference"
for ((a=0; a<${#fasta[@]}; a++))
        do	
		echo -e "\nassessing ${fasta[a]}"

		#index contigs as reference genome
		bwa index $path2fastas/${fasta[a]}

		#paired
		bwa aln -t $threads $path2fastas/${fasta[a]} $read_path/$sample.pair1.fastq.gz > reads_1.sai
		bwa aln -t $threads $path2fastas/${fasta[a]} $read_path/$sample.pair2.fastq.gz > reads_2.sai
		bwa sampe -f pe.sam -r "@RG\tID:${fasta[a]}_pe\tSM:$sample\tPL:$platform" $path2fastas/${fasta[a]} reads_1.sai reads_2.sai $read_path/$sample.pair1.fastq.gz $read_path/$sample.pair2.fastq.gz
		samtools view -Sb -o pe.bam pe.sam
		samtools sort pe.bam pe.sort				

		#single
		bwa aln -t $threads $path2fastas/${fasta[a]} $read_path/$sample.single.fastq.gz > reads_3.sai
		bwa samse -f se.sam -r "@RG\tID:${fasta[a]}_se\tSM:$sample\tPL:$platform" $path2fastas/${fasta[a]} reads_3.sai $read_path/$sample.single.fastq.gz
		samtools view -Sb -o se.bam se.sam
		samtools sort se.bam se.sort

		#merge
		samtools merge ${fasta[a]}.bam pe.sort.bam se.sort.bam
		#create sam file
		samtools view -h ${fasta[a]}.bam > ${fasta[a]}.sam

		#create afg file
		/storage/data_1/megan/programs/abyss-samtoafg.pl $path2fastas/${fasta[a]} ${fasta[a]}.sam > ${fasta[a]}.afg

		#create amos bank
		/usr/local/bin/amos-3.1.0/bin/bank-transact -cb ${fasta[a]}.bnk -m ${fasta[a]}.afg
		
		#run amos validate
		/usr/local/bin/amos-3.1.0/bin/amosvalidate_fixed ${fasta[a]}.bnk

		#generate file for visualization
		#hawkeye ${fasta[a]}.bnk		

		#summarize results
		#put header
		echo "${fasta[a]}" > temp.amos.out
		#count mate pair unhappiness
		wc -l ${fasta[a]}.ce.feat >> temp.amos.out
		#count bad read coverage
		wc -l ${fasta[a]}.depth.feat >> temp.amos.out
		#count misassemblies with at least 3 of the 7 types and a total of at least 9 features
		awk -F':' '{print $2,$3}' ${fasta[a]}.suspicious.feat | awk '{if ($1>8 && $3>2) print $0}' | wc -l >> temp.amos.out	
                #use alignments to estimate reads mapping to assemblies
                map=$(samtools view -c -F 4 ${fasta[a]}.bam)
                unmap=$(samtools view -c -f 4 ${fasta[a]}.bam)
                map_percent=$(echo $map $unmap | awk '{print $1/($1+$2)}')
                echo $map_percent >> temp.amos.out

		#format output
		cut -d " " -f 1 temp.amos.out | tr $'\n' $' ' >> amos.out
		echo "" >> amos.out
	
		#remove intermediate files
		rm ${fasta[a]}.bnk.runAmos.log ${fasta[a]}.22.n22mers ${fasta[a]}.fasta ${fasta[a]}.afg ${fasta[a]}.*feat ${fasta[a]}.snps
		rm -r ${fasta[a]}.bnk 		
		rm ${fasta[a]}.rsa ${fasta[a]}.sa ${fasta[a]}.rbwt ${fasta[a]}.bwt ${fasta[a]}.rpac ${fasta[a]}.amb ${fasta[a]}.ann ${fasta[a]}.pac
		rm reads_1.sai reads_2.sai pe.sam pe.bam pe.sort.bam reads_3.sai se.sam se.bam se.sort.bam
		rm ${fasta[a]}.bam ${fasta[a]}.sam
		rm temp.amos.out
	done

#add header to the amos stats
echo "assembly amos_mates amos_cov amos_all %map" > head.txt
cat head.txt amos.out > $sample.amos.stats
rm amos.out head.txt



#create summary sheet
echo "summarizing results"
tr ' ' '_' <$sample.general.stats | awk -F"," '{print $1,$4,$37,$55,$57}' > general.temp
tr ' ' '_' <$sample.compass.stats | awk '{print $2, $11, $15, $16, $17}' > compass.temp
awk '{print $1, $4, $7}' $sample.cegma.complete.stats > cegma.complete.temp
awk '{print $1, $4, $7}' $sample.cegma.partial.stats > cegma.partial.temp
pr -m -t -s" "  general.temp compass.temp cegma.complete.temp cegma.partial.temp $sample.amos.stats > $sample.assessment.summary
rm general.temp compass.temp cegma.complete.temp cegma.partial.temp

echo | mutt -s "assessment summary complete" $email -a $sample.assessment.summary 



echo -e "\ngenerate bams for visualization"

#index reference if it doesn't already exist 
if (test -e $ref.bwt)
	then echo "skipping indexing...index already exists"
else echo "creating index"
	/storage/data_1/megan/programs/bwa-0.7.5a/bwa index $ref
fi

#align each assembly to the reference
for ((a=0; a<${#fasta[@]}; a++))
        do

                echo -e "\nassessing ${fasta[a]}"	
		/storage/data_1/megan/programs/bwa-0.7.5a/bwa mem -t $threads $ref $path2fastas/${fasta[a]} > ${fasta[a]}.sam  
		samtools view -Sb -F 4 -o ${fasta[a]}.mapped.bam ${fasta[a]}.sam
		samtools sort ${fasta[a]}.mapped.bam ${fasta[a]}.mapped.sort
		samtools index ${fasta[a]}.mapped.sort.bam
		rm ${fasta[a]}.sam ${fasta[a]}.mapped.bam 

	done


echo | mutt -s "assessment script from $1 complete" $email
echo -e "\nDONE!!!"
