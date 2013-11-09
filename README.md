genome_reseq_scripts
====================
scripts for processing Illumina resequencing data:

align_reads.sh
--------------
Aligns paired end Illumina reads to a reference sequence, prepares alignments for genotyping (marks duplicates and realigns around indels), and examines the quality of the alignments. See example.align.input for example input file. Produces a bam file for each sample to be used for SNP calling. Requires BWA, samtools, picard, GATK.

assemble_params.sh
------------------
Assemble multiple samples using a variety of parameters. Specificy input information (files, parameters, etc) in the header of the script. Produces a fasta file for each sample for each parameter combination. Requires ABySS.

assess_assemblies.sh
--------------------
Evaluates a set of assemblies with a variety of metrics. See example.align.input for example input file. Produces a summary file with key metrics. Requires assemblathon2-analysis git repository, cegma, bwa, samtools, amos.

bam2blastid.sh
--------------
Calls genotypes from bam files individually and blasts a specified genomic region to a database. See example.blastid.input for example input file. Produces a report indicating the top blast hit for each sample. Requires BLAST, plus additional scripts from this repository (gatk_genos.sh, vcf2format.pl, fasta_subregion.pl) and all requirements within these scripts

clean_reseq4assembly.sh
-----------------------
Cleans whole genome resequencing reads prior to assembly. See example.clean.reseq.input for example input file. Produces three cleaned fastq files for each pair of input fastq files (two with pairs, one with orphans). Requires Trimmomatic, FLASh, JELLYFISH, Quake, and FastQC.

fasta_subregion.pl
------------------
Extracts a specified region from entries in an aligned fasta file. Inputs are an aligned fasta file, start position, and stop position. Produces an aligned fasta of the specified region.

gatk_genos.sh
-------------
Calls and filters multi-sample genotypes from a set of alignment (bam) files using GATK's pipeline for multi sample snp calling. See example.genos.input for example input file. Produces a single vcf with filtered genotype calls for each sample. Requires GATK.

snp_summary.pl
--------------
Summarizes genotypes from a matrix genotype file.  Inputs are a matrix of genotypes and a fasta reference.  Outputs the number of SNPs among the samples, the number of genotyped positions per sample, the number of heterozygote SNPs per sample, and the number of SNPs relative to the reference for each sample.

vcf2format.pl
-------------
Convert filtered vcf file(s) (from GATK v1.2) into fasta format. Input is a file with a list of vcf files. Outputs are a fasta file for each contig. Requires FormatGenos.pm from this repository.  **warning: this script is very specific to older versions of GATK and the filtering in gatk_genos.sh
