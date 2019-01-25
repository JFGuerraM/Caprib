#!/usr/bin/perl
#Juan Guerra 18-10-16  
#script to merge caprib analyses with cdd ncbi features results
use strict;
use warnings;
#use Set::Scalar;
use Data::Dumper;

#usage: mergeCDDfeature.pl ISvsD.csv ISvsDfeatdata.txt merged.tsv

my $path= pop(@ARGV);
my $hitdata=pop(@ARGV);
my $report=pop(@ARGV);

findCDD();
exit;

sub indices{ 
	my $countHits=0; 	#counting every proteins in the file
	my $REP;
	my %dictP; #dictionary for proteins and mutations
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
			$dictP{$refseq}=$protein;			
		}		
	}
	close($REP);

	return %dictP;
}


sub mutations{ 
	my $countHits=0; 	#counting every proteins in the file
	my $REP;	
	my %dictM={}; #dictionary for proteins and mutations
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
			$dictM{$refseq}=$protein;			
		}		
	}
	close($REP);
	
	return %dictM;
}



sub findCDD{
	my $countHits=0; 	#counting every cdd in the file
	my $HIT;
	my %dictP=indices();	
	my ($fho, $output);
	my %dictM=mutations();
	my ($refseq, $max, $min);
	my (@values, @coordinate);
	$output= "$path";#Elongata"mergedfetMtb.txt";#	
	open($HIT, $hitdata);
	open($fho,">$output");
	print $fho "Protein\tMutation position\tQuery\tType\tTitle\tcoordinates\tcomplete size\tmapped size\tsource domain\n";
								
	while(<$HIT>){
		chomp;
		$a=$_;		
		my @line =split("\t", $a);
		if ($a =~ /^Q\#/){
			$countHits++;			
			$refseq=(split(" ", $line[0]))[2]; 			
			#we take coordinates, 3 cases: A34, p45 ..., P100; A34-P45; P100
			if($line[3] =~/,/g){				
				@coordinate=split(",", $line[3]);	
				$min= undef;
			} 
			elsif($line[3] =~/-/g){
				
				$min=(split("-", $line[3]))[0];
				$min =~ s/[a-zA-Z_]//g;				
				$max=(split("-", $line[3]))[1];
				$max =~ s/[a-zA-Z_]//g;
			}
			elsif($line[3] =~/[a-zA-Z_]/){								
				@coordinate=($line[3]);
				$min= undef;				
			}
			
			if (exists ($dictM{$refseq})){
				
				@values = $dictM{$refseq};
				
				my @mut=();
				if(scalar @values > 1){
					@mut=split(",",($values[1]));
				}else{
					@mut =$values[0];
				}				
				my $ref=$dictP{$refseq}; #locustag
				
				foreach my $i(@mut){
					$i =~ s/[^a-zA-Z0-9_]//g;
					
					if (defined $min && ($i >= $min) && ($i <= $max)){						
						print $fho "$ref \t$i\t".$a."\n";						
					}
					else{
						foreach my $j(@coordinate){
							$j =~ s/[a-zA-Z_]//g;							
							if ($i == $j){
								print $fho "$ref \t$i\t".$a."\n";
							
							}
						}
					}
				}				
			}			
		}		
	}
	close($HIT);
	close $fho;
	print "count= ".$countHits;
}