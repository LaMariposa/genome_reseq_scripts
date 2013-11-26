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

use List::Util qw( min max );

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
#flag to ignore contig
my $ignore=0;

#read in first entry in bits file
my $line=<INBITS>;
my @info=split("\t",$line);
#set end so this piece starts as the beginning
if ($info[17] eq "Plus"){$last_end=$info[8]-$info[6];}
        elsif ($info[17] eq "Minus"){$last_end=$info[8]-$info[2]+$info[6]-1}


#process each entry in the bits file
while ($line)
	{
	  #set current contig information
	  @info=split("\t",$line);
          $current_contig=$info[0];
	  print "examining contig $previous_contig & $current_contig: ";
	  #set contig position based on orientation
          if ($info[17] eq "Plus"){$cont_start=$info[6];$cont_end=$info[7]}
          	elsif ($info[17] eq "Minus"){$cont_start=$info[7];$cont_end=$info[6]}


	  #get correct sequence orientation
	  if ($info[17] eq "Plus"){$current_seq=$ref{$current_contig}}
		elsif ($info[17] eq "Minus")
			{
			   #on minus strand so reverse complement
                           $current_seq=reverse($ref{$current_contig});
                           $current_seq=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
			}

	  #see how it relates to previous entry
	  #check to see if they are different contigs
	  if ($current_contig ne $previous_contig)
		{
		  #they are different, so process
		  #calculate extended position of contig relative to reference
		  if ($info[17] eq "Plus"){$start=$info[8]-$info[6]+1}
			elsif ($info[17] eq "Minus"){$start=$info[8]-$info[2]+$info[6];}
              	  if ($info[17] eq "Plus"){$end=$info[9]+$info[2]-$info[7];}
			elsif ($info[17] eq "Minus"){$end=$info[9]+$info[7]-1;}	 
		  #do they overlap
		  if ($start>$last_end)
			{
		  	  #no overlap with previous, so print out Ns for the gap
			  print "no overlap\n";
			  my $numNs=$start-$last_end-1;
			  $result.='N' x $numNs;
			  #now print out contig
			  $result.=$current_seq;
			}

			else
			{
		  	  #the current and previous overlap, so test various cases
			  #pull out the sequence that is exected to overlap
			  my $overlap=$last_end-$start+1;
			  print "overlap=$overlap, ";
			  my $bit_previous_seq=substr $previous_seq, -$overlap;
			  my $bit_current_seq=substr $current_seq, 0, $overlap;
			  			  
			  #check if the sequences are identical
			  if ($bit_current_seq eq $bit_previous_seq)
				{
				  #sequences are identical so remove bit from current and print
				  print "identical\n";
				  $result.=substr $current_seq, $overlap;
				}
				else
				{
				  #is there high similarity as pieces are slid apart
				  my $diffs=($bit_current_seq ^ $bit_previous_seq) =~tr/\0//;
				  my $similarity=$diffs/length($bit_current_seq);
				  while (length($bit_current_seq) > 10 && $similarity < 0.9)
					{
					  #aren't similar enough, so slide one
					  $bit_previous_seq=substr $bit_previous_seq, 1;
					  $bit_current_seq=substr $bit_current_seq, 0, -1;
					  #calculate similarity
					  $diffs=($bit_current_seq ^ $bit_previous_seq) =~tr/\0//; 
					  $similarity=$diffs/length($bit_current_seq);
					}

				  #if a good match was found, merge the contigs
				  if (length($bit_current_seq) > 10)
					{
					  #remove matching piece from current and print
					  my $size=length($bit_current_seq);
					  my $percent=$similarity*100;
					  print "a match of length $size and similarity $percent%\n";
					  $result.=substr $current_seq, $size;
					}

					else
					{
					  #no good match found sliding apart, so slide together
					  #reset values
					  my $bit_prev=substr $previous_seq, -$overlap;
                          		  my $bit_curr=substr $current_seq, 0, $overlap;

					  #find size of smallest contig
					  my @lengths=(length($previous_seq), length($current_seq));
					  my $min_len=min(@lengths);
					  
					  #find similarity
					  my $diffs=($bit_curr ^ $bit_prev) =~tr/\0//;
                                  	  my $similarity=$diffs/length($bit_curr);
					  my $slide=1; #tracks how many positions slid

					  #while no match
					  while(length($bit_curr) < $min_len && $similarity < 0.9)
						{
						  #aren't similar enough, so slide one
						  $bit_prev=substr $previous_seq,-$overlap-$slide;
						  $bit_curr=substr $current_seq,0, $overlap+$slide;
						  
						  #calculate similarity
						  $diffs=($bit_curr ^ $bit_prev) =~tr/\0//;
                                          	  $similarity=$diffs/length($bit_curr);
					
						  #update slide
						  $slide++;
						}
					  #if match found, merge by removing matching piece from current 
					  if (length($bit_curr) < $min_len)
						{
                                         	  my $size=length($bit_curr);
                                          	  my $percent=$similarity*100;
						  print "a match of length $size and similarity $percent%\n";
						  $result.=substr $current_seq, $size;
						}
						else
						{
					  	  #there was no match found--chop off pieces that should overlap and add double the Ns
					  	  print "no similarity\n";	

				  	  	  #chop off the end of the results
				  	  	  $result=substr($result, 0, -$overlap);
				  	  	  #add Ns
				  	  	  $result.='N' x $overlap;
				  	  	  $result.='N' x $overlap;
				  	  	  #chop off start of current contig and add
				  	  	  $result.=substr $current_seq, $overlap;
						}
					}
				}
		  	}
		}
		else 
		{
		  #they are the same
		  #test to see if order and orientation are consistent
		  print "same contig, ";
		  if ($prev_orient eq $info[17])
			{
			  if (($info[17] eq "Plus" && $cont_end > $prev_cont_start) || ($info[17] eq "Minus" && $cont_end < $prev_cont_start))
				{
				  print "proper order and orientation\n";
				  #don't do anything but update coordinates
                  		  if ($info[17] eq "Plus"){$start=$info[8]-$info[6]+1}
                  		  	elsif ($info[17] eq "Minus"){$start=$info[8]-$info[2]+$info[6]}
             		          if ($info[17] eq "Plus"){$end=$info[9]+$info[2]-$info[7];}
                	       		elsif ($info[17] eq "Minus"){$end=$info[9]+$info[7]-1;}
				}
				else 
				{
				  print "but different order, check manually\n";
				  #skip contig so don't update coordinates
				  $ignore=1;			  
				}
			}
			else 
			{
			  print "but different orientation, check manually\n";
                          #skip contig so don't update coordinates
                          $ignore=1;
			}
		}

	  #read in next entry
	  $line=<INBITS>;
	  #set previous contig info
	  if ($ignore==0)
		{
          	  $prev_cont_start=$cont_start;
          	  $prev_cont_end=$cont_end;
          	  $prev_orient=$info[17];
          	  $previous_contig=$current_contig;
          	  $last_end=$end;
          	  $previous_seq=$current_seq;
		}
		else
		{
		  #ignore contig and reset flag
		  $ignore=0;
		}

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
