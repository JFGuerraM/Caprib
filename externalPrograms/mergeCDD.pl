#!/usr/bin/perl
#Juan Guerra 18-10-16  
#script to do set operations within species
use strict;
use warnings;
#use Set::Scalar;
#use Data::Dumper;

my $path= pop(@ARGV);
my $hitdata=pop(@ARGV);
my $report=pop(@ARGV);



findCDD();
exit;

sub indices{ 
	my $countHits=0; 	#counting every proteins in the file
	my $REP;
	my ($dictP,$dictM); #dictionary for proteins and mutations
	my ($protein, $refseq);
	open($REP, $report);
	while(<$REP>){
		chomp;
		$a=$_;
		my @values =split(";", $a);
		if (($a !~ /^Protein/)&&(defined $values[4])){
			$countHits++; 			
			$protein= $values[0];
			$refseq=$values[1];			
			$dictP->{$refseq}=$protein;			
		}		
	}
	close($REP);

	return $dictP;
}

sub mutations{ 
	my $countHits=0; 	#counting every proteins in the file
	my $REP;
	my $dictP={};
	my $dictM={}; #dictionary for proteins and mutations
	my ($protein, $refseq);
	open($REP, $report);
	while(<$REP>){
		chomp;
		$a=$_;
		my @values =split(";", $a);
		if (($a !~ /^Protein/)&&(defined $values[4])){
			$countHits++; 			
			$protein= $values[4];
			$refseq=$values[1];			
			$dictP->{$refseq}=$protein;			
		}		
	}
	close($REP);
	
	return $dictP;
}


sub findCDD{
	my $countHits=0; 	#counting every cdd in the file
	my $HIT;
	my $dictP=indices();	
	my ($fho, $output);
	my $dictM=mutations();
	my ($protein, $refseq, $max, $min);
	my @values;
	$output= "$path";
	open($HIT, $hitdata);
	open($fho,">$output");
	print $fho "Protein\tMutation position\tQuery\tHit type\tPSSM-ID\tFrom\tTo\tE-Value\tBitscore\tAccession\tShort name\tIncomplete\tSuperfamily\n";

	while(<$HIT>){
		chomp;
		$a=$_;
		my @line =split("\t", $a);
		if ($a =~ /^Q\#/){
			$countHits++;			
			$refseq=(split(" ", $line[0]))[2]; #on enleve le Q# - 
			$min=$line[3];
			$max=$line[4];
			
			if (exists ($dictM->{$refseq})){
				@values = %$dictM{$refseq};
				my @mutations=split(",",($values[1]));
				my $ref=%$dictP{$refseq};
				foreach my $i(@mutations){
					$i =~ s/[^a-zA-Z0-9_]//g;
					if (($i >= $min) && ($i <= $max)){						
						print $fho "$ref \t$i\t".$a."\n";						
					}					
				}				
			}			
		}		
	}
	close($HIT);
	close $fho;
	print $countHits;
}