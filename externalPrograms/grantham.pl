#!/usr/bin/perl -w 
#Juan Guerra 04-06-2019

use strict;
use warnings;

# Usage: ./grantham.pl results.csv indice seuil
#indice est la colonne de l'organisme que tu veux filtrer, s'il s'agit du premier organisme alors indice = 1 et ainsi de suite...
#seuil est la valeur de grantham a filtrer

my $seuil= pop(@ARGV);;
my $indice= pop(@ARGV);
my $inputf = pop(@ARGV);



my $output="mutationsGrantham.csv";


my $REPORT;
my $fho;
my $line;
open($REPORT, $inputf);

my $skip=1;
my $count=0;
my @var;
my $organisms_count=0;
my $proteins_candidate=0;
my $grantham;
my $ex;
my $variation;

open ($fho, ">$output");
print $fho "Protein;name;mutation;Grantham;EX\n";
while (<$REPORT>) {
	chomp;
	$line =$_;
	#nombre d'organismes:
	if($line=~/^Group B:/){
		@var=split(";",$line);
		$organisms_count=scalar (split(",",$var[1]));
		print "nombre d'organismes: $organisms_count\n";
	}
	#nombre de protéines
	if($line=~/^Proteins candidates:/){
		@var=split(";",$line);
		
		print "nombre de protéines: $var[1]\n";
	}
	#on prend que les protéines 
	if ($skip) {		
		if($line=~/^Protein;/){			
			$skip=0;
		}else{
			next;
		}
	}
	elsif($line=~/^Stats:/){		
		$skip=1;
	}
	#finalement on trouve
	elsif(!$skip && !/^$/ ){		#/^\s*$/ length($line)!=0
		@var=split(";",$line);
		$count++;
		my $protein=$var[0];
		my $name=$var[2];
		my $mutations=$var[4+$indice];

		$mutations =~ s/(\{|\})//g;
		if($mutations=~/,/){
				my @each_mutation=split(",",$mutations);
				foreach my $x (@each_mutation) {
					$grantham=(split(":",$x))[1];
					$ex=(split(":",$x))[2];
					$variation=(split(":",$x))[0]; #AA/AA
					if ($grantham>=$seuil) {
						print $fho "$protein;$name;$variation;$grantham;$ex\n";
					}
					
				}
				
		}else{
			$grantham=(split(":",$mutations))[1];
			$ex=(split(":",$mutations))[2];
			$variation=(split(":",$mutations))[0]; #AA/AA
			if ($grantham>=$seuil) {
						print $fho "$protein;$name;$variation;$grantham;$ex\n";
					}

		}
		
	}
	
}
print $count;
close ($REPORT);
close ($fho);
