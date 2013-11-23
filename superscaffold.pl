#!/usr/bin/perl

#superscaffold.pl by Megan Supple 14 Nov 2013
#last updated 14 Nov 2013

#script to scaffold contigs based on mummer alignments

#usage:	superscaffold.pl <in.contigs.fasta> <in.bits.txt>
	#<in.contigs.fasta> is a set of assembled contigs in fasta format
	#<in.bits.txt> is an ordered list of mummer alignments for a single reference contig/scaffold from 
		#nucmer $ref $id.contigs.fasta
		#show-coords -rB out.delta


use strict;
use warnings;


#command line arguments
my $usage = "Usage: superscaffold.pl <in.contigs.fasta> <in.bits.txt>
options [all required]:
	<in.contigs.fasta> 	a set of assembled contigs in fasta format
	<in.bits.txt> 		an ordered list of mummer alignments from nucmer and show-coords -rB
";

die "$usage" unless (@ARGV==2);
my ($contigs, $bits) = @ARGV;


#make a hash of contigs names and sequences
open(REF, $contigs)||die "can't open input contigs file. $!\n";
my %ref;
my $contig;
my $seq="";

while (my $fasta=<REF>)
        {
         chomp($fasta);
         #check if new contig, if so add to hash
         if ($fasta=~ m/^>.*/)
                {
                 #record contig name
		 my @temp=split(" ",$fasta);
                 $contig=substr($temp[0], 1);
                 #reset seq to null
                 $seq="";
                }
                else
                {
                 #add seq
                 $seq.=$fasta;
                }
         #add contig and seq to hash
         $ref{$contig}=$seq;        
        }
close REF;
#while ( my ($key, $value) = each(%ref) ) {
 #       print "$key => $value\n";
  #  }

#open and read in input alignment bits file
open(INBITS, $bits)||die "can't open mummer alignments file. $!\n";


#some variables
#position on reference
my $start;
my $end;
my $last_end;
#previous orientation
my $prev_orient="";
#position on contig
my $cont_start;
my $cont_end;
my $prev_cont_start;
my $prev_cont_end;
#to track contig names
my $current_contig="";
my $previous_contig="";
#string to print to outfile
my $result="";
my $current_seq="";
my $previous_seq="";


#read in first entry in bits file
my $line=<INBITS>;
my @info=split("\t",$line);
#set end so this piece starts as the beginning
if ($info[17] eq "Plus"){$cont_start=$info[6]}
        elsif ($info[17] eq "Minus"){$cont_start=$info[7]}
$last_end=$info[8]-$cont_start;


#process each entry in the bits file
while ($line)
	{
	  #set current contig information
	  @info=split("\t",$line);
          $current_contig=$info[0];
print "contig=$current_contig\n";

	  #get correct sequence orientation
	  if ($info[17] eq "Plus"){$current_seq=$ref{$current_contig}}
		elsif ($info[17] eq "Minus")
			{
			   #on minus strand so reverse complement
                           $current_seq=reverse($ref{$current_contig});
                           $current_seq=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
			}
#print "seq=$temp_seq\n";

	  #see how it relates to previous entry
	  #check to see if they are different contigs
	  if ($current_contig ne $previous_contig)
		{
		  #they are different, so process

		  #calculate extended position of contig relative to reference
                  if ($info[17] eq "Plus"){$cont_start=$info[6];$cont_end=$info[7]}
                        elsif ($info[17] eq "Minus"){$cont_start=$info[7];$cont_end=$info[6]}
                  $start=$info[8]-$cont_start+1;
                  $end=$info[9]+$info[2]-$cont_end;

		  #do they overlap
		  if ($start>$last_end)
			{	
		  	  #no overlap with previous, so print out Ns for the gap
			  my $numNs=$start-$last_end-1;
			  $result.='N' x $numNs;
#$result.="\n";
			  #now print out contig
			  $result.=$current_seq;
#$result.="\n"; 
			}

			else
			{
		  	  #the current and previous overlap, so test various cases
			  #pull out the sequence that is exected to overlap
			  my $overlap=$last_end-$start+1;
print "overlap=$overlap\n";
			  my $bit_previous_seq=substr $previous_seq, -$overlap;
			  my $bit_current_seq=substr $current_seq, 0, $overlap;
#print "prev=$bit_previous_seq\n";
#print "curr=$bit_current_seq\n";
			  			  
			  #check if the sequences are identical
			  if ($bit_current_seq eq $bit_previous_seq)
				{
				  #sequences are identical so remove bit from current and print
print "identical\n";

my $x=substr $current_seq, $overlap;
#print "x=$x\n";
				  $result.=substr $current_seq, $overlap;
#$result.="\n";
				}
				else
				{
print "not identical\n";
				  #are the sequence identical/global match with high similarity
				  #do the edges overlap
				  #is there no match--chop off pieces that should overlap and add double the Ns
				  #chop off the end of the results
#print "result=$result\n";
				  $result=substr($result, 0, -$overlap);
				  #add Ns
#$result.="\n";
				  $result.='N' x $overlap;
				  $result.='N' x $overlap;
#print "contig is=$current_contig\n";
#print "overlap=$overlap\n";
#$result.="\n";
				  #chop off start of current contig and add
				  $result.=substr $current_seq, $overlap;
#$result.="\n";
				}
		  	}
		}
		else 
		{
		  #they are the same
		  #test to see if order and orientation are consistent
		  print "same contig???\n";
		  if ($prev_orient eq $info[17])
			{
			  print "same orientation\n";
print "prevcontstart=$prev_cont_start\n";
print "prevcontend=$prev_cont_end\n";
print "start=$cont_start\n";
print "end=$cont_end\n";
			  if ($info[17] eq "Plus" && $cont_end > $prev_cont_start || $info[17] eq "Minus" && $cont_end < $prev_cont_start)
				{
				  print "order ok\n";
				}
			}
			else {print "same contig, different orientation, check manually\n";}
		}

	  #read in next entry
	  $line=<INBITS>;
	  #set previous contig info
          $prev_cont_start=$cont_start;
          $prev_cont_end=$cont_end;
          $prev_orient=$info[17];
          $previous_contig=$current_contig;
          $last_end=$end;
          $previous_seq=$current_seq;

	}


#open outfile
#make oufile name
my @outfile=split(/\./,$contigs);
open(OUTFASTA, ">$outfile[0].superscaff.fasta");
#print header
print OUTFASTA ">$outfile[0].superscaff.fasta\n";
print OUTFASTA $result;
print OUTFASTA "\n";

close INBITS;
close OUTFASTA;
