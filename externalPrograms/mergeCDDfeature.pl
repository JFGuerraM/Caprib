#!/usr/bin/perl
#Juan Guerra 18-10-16  
#script to do set operations within species
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
	my ($dictP,$dictM); #dictionary for proteins and mutations
	my ($protein, $refseq);
	open($REP, $report);
	while(<$REP>){
		chomp;
		$a=$_;
		my @values =split(";", $a);
		if (($a !~ /^Protein/)&&(scalar @values >5)){
			$countHits++; 			
			$protein= $values[0];
			$refseq=$values[1];			
			$dictP->{$refseq}=$protein;			
		}
		
	}
	close($REP);
	print "Reporte con $countHits proteinas\n";
	#print $countHits;
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
		
		if (($a !~ /^Protein/)&&(scalar @values >3)){
			$countHits++; 			
			$protein= $values[4];
			$refseq=$values[1];			
			$dictP->{$refseq}=$protein;			
		}
		if ($a =~ /Stats: counting variations/){
			last;
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
	my ($refseq, $max, $min);
	my (@values, @coordinate);
	$output= "$path";#Elongata"mergedfetMtb.txt";#
	print $hitdata;	
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
			print $refseq."\n";
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
			
			if (exists ($dictM->{$refseq})){
				
				@values = %$dictM{$refseq};
				print "values1: @values[1]\n";
				my @mutations=split(",",($values[1]));
				print @mutations;
				my $ref=%$dictP{$refseq}; #locustag
				print $ref."\n";
				foreach my $i(@mutations){
					$i =~ s/[^a-zA-Z0-9_]//g;
					
					if (defined $min && ($i >= $min) && ($i <= $max)){
						print "min: ".$min."\n";
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