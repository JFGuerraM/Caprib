#!/usr/bin/perl -w 
#Juan Guerra 12-08-16  rewrite 02-03-18
use strict;
use Data::Dumper;

my $identityLimit= pop(@ARGV);
my $output=pop(@ARGV); 
my $inputf = pop(@ARGV);


blastProteins($inputf, $identityLimit);

sub blastProteins{

$inputf=$_[0];
$identityLimit=$_[1];
#File handler
my ($BLAST, $fho);


open($BLAST, $inputf);
open ($fho, ">$output");
print $fho "Proteine;length Sbjct;Identity;similarity;e-value;AA_identical;AA_similar;AA_diference;gaps Query;gaps Subject; Stop codon\n"; #CSV head	

#protein name 
my $queryR;
#query, subject, evalue, %identity, protein length, %similar
my($q, $s, $ev, $longueurP, $sim);
#auxilar variables to the split sequence line
my($aux0, $aux1, $aux2, $aux3);
#to extract the sequence
my ($start, $end);
my ($firstQ,$lastQ,$seqQ, $seqS, $sequence);     
my $id=0;
#to know that we're between query ans subject sequence
my $control2=0;



my $line;
#Sometimes we have more than one hit for the same protein
#in order to know if we're in the first we add to control that start at 0. 

my $control;

while(<$BLAST>){
	chomp;
	$line=$_;	
	#if new protein founded	we take the name
	if ($line =~ /^Query=(.*)/){ 
		$queryR=$1;			
		$q=0;
		$s=0;
		$control=0;			
	}
	#if we arrive to score line then we take the evalue "ev"
	elsif($line=~  /^ Score =/){		
		if (!$control){
			$ev = (split)[7]; 
			$ev =~ s/,//;
			
			#print $ev,"\n";
		}		
		elsif ($control >= 1){
			$q=$s=1;
		}
	}
	#To take the %indentity "id",reference protein lenght "longueurP" and similar score "sim" 
	elsif($line=~ /^ Identities/){		
		$control++;		
		if ($control==1){
			$id = (split)[3];
			$id =~ s/(\(|\)|\%|,)//g;
			$longueurP=(split)[2]; #query protein lenght
			$longueurP=~ s/((.*)\/)//g;
			$sim= (split)[7];
			$sim =~ s/(\(|\)|\%|,)//g;
			
		}		
		elsif ($control >= 1){
			$q=$s=1;
		}			
	}
	
	#if no hits found and identityLimit =0
	elsif($line=~ /^\*\*\*\*\* No hits found \*\*\*\*\*/){		
		$id = 0;
		$sim= 0;		
		$longueurP=0;
		$ev=0;
		$q=$s=1;
		if ($id >= $identityLimit){
			print $fho "$queryR;$longueurP;$id;$sim;$ev\n";
		}
		
				
	}
	
	elsif(($line=~ /^Query /)&& (!$q)){
		
		($aux0, $aux1, $aux2, $aux3) = split;  # split separa $_ en 4 escalares separados por espacios en blanco
		if(!$q && $control==1){	$firstQ = $aux1;	}		
		$start=index($line, $aux2);
		$end=length($aux2);
		$lastQ = $aux3;		
		$seqQ .= $aux2;		
		$q++;  # we are reading 'query' Q=1
		$s=0;		
		$control++;	
		
		$control2=1;
		
	}	
	elsif(($line=~ /Sbjct/)&& (!$s))	{			
		$seqS.= substr($_,$start,$end);	
		$s++;
		$q=0;	
		#$control=0;
		$control2=0;
	} 
	#we're reading the comparative line 
	elsif ($control2){
		my $var=$_;				
		$sequence .=substr($var,$start,$end);		
		$control2=1;
		
	}
	if(($id < $identityLimit)&&(/^Effective/ && $s )){
		$control++;
		$s=1;
		$sequence= '';
		$q = 1;
		$seqQ='';
		$seqS='';
		$queryR=""
	}
	#we print only filtered information for the first hit
	if (($id >= $identityLimit)&&(/^Effective/ && $s )){
		#first informations: "Proteine;length Sbjct;Identity;similarity;e-value"		
		#now we need to extract information in the compare line
		#so, if there are a sequence, we'll not take if we don't found hits id=0
		
		if (length($sequence)>0){
			#Counters to follow amino acids positions in the query, subject,comparative line, gaps
			my $counting =$firstQ;
			my $countgapsQ =$firstQ;
			my $countgapsS=$firstQ;
			my $countGi=$firstQ; # nous dit la position du premier acide aminé dans le blast de la proteine  
			my ($begin, $final)=0; #varibles para determinar consecutivos
			
			my $controlgap=0;
			#to store each position
			my @gap=();
			my @gapOut=();
			my @aaStop =();
			my (@aaidentical,@aasimilar,@aadiference)=();
			my(@aasimilar_L,@aadiference_L);
			my (@gapQ, @gapS, @trunqueS)=();
			#to take every gap index in the sequence
			my @indicesGapSeq=(); #pour prendre les indices de chaque gap in the sequence
			my  @aas=(); 
			my @aasOut=();
			#we split our three sequences	
			my @seq=split("",$sequence);
			my @seqQ=split("",$seqQ);
			my @seqS=split("",$seqS);
			
			#we walk through the query sequence i order to take gap in position 
			#to store it in gapQ list 
			foreach my $letter(@seqQ){									
				if($letter=~ /\-/){					
					push((@gapQ,$countgapsQ));						
				}
				$countgapsQ++;
			}
			#we walk through the subject sequence i order to take gap out position 
			#to store it in gapS list or in the stop codon list @trunqueS 
			foreach my $letter(@seqS){					
				if($letter=~ /\-/){
					push(@gapS,$countgapsS);						
				}
				elsif($letter=~ /\*/){
					push(@trunqueS,$countgapsS);
						
				}					
				$countgapsS++;
			}
			#if query don't start at position 1 means the amino acids before are different
			if ($firstQ!=1){
				my @tempList=1..$firstQ-1;
				foreach my$i(@tempList){
					push(@aadiference,$i);
				}																	
			}
			#now we'll work in the comparative sequence
			#
			foreach my $letter(@seq){					
				#print $letter;
				if($letter=~ /\s/){	#if maybe there is a gap...										
					for (my $i = 0; $i < scalar @gapQ; $i++){	#if the gap is in the Query take the gap query list						
						#print $countGi ;
						if ($gapQ[$i] == $countGi){								
							#print $countGi, "\n" ;								
							if (($i < (scalar @gapQ-1) &&(($gapQ[$i]+1) == $gapQ[$i+1]))){				#si son consecutivos ej [12]+1 =[13]
								push(@aas, $seqS[($countGi)-($firstQ)]);																		
							}
							else{
								$begin=$counting-1; 
								push(@aas, $seqS[($countGi)-($firstQ)]);									
								push(@gap,($begin, @aas));									
								@aas=();
							}								
							$controlgap++;							
						}						
					}
					for (my $i = 0; $i < scalar @gapS; $i++){	#if the gap is in the Subject take the gap subject list 
						
						if ($gapS[$i] == $countGi){								
																				
							if (($i < (scalar @gapS-1) &&(($gapS[$i]+1) == $gapS[$i+1]))){			#si son consecutivos ej [12]+1 =[13]
								push(@aasOut, $seqQ[($countGi)-($firstQ)]);	
								push(@indicesGapSeq, ($counting-1));#$dif ($countGi)-($firstQ)
								}
							else{
								if (scalar(@indicesGapSeq)>0){	$begin=shift(@indicesGapSeq)} 
								else {	$begin=$counting-1;	}																		
								push(@aasOut, $seqQ[($countGi)-($firstQ)]);									
								push(@gapOut,($begin, @aasOut));									
								@aasOut=();
								@indicesGapSeq=();
							}								
							$controlgap++;
							$counting++;								
						}						
					}						
					for (my $i = 0; $i < scalar @trunqueS; $i++){
						if ($trunqueS[$i] == $countGi){	
							push(@aaStop, $counting);
							$controlgap++;																
						}							
					}						
						
					if (!$controlgap){
						push(@aadiference_L,"$counting".":"."$seqQ[$countGi-($firstQ)]\/$seqS[$countGi-($firstQ)]");
						$counting ++;
					}						
					$controlgap=0;											
				}
				#if identical amino acid 						
				elsif ($letter=~ /[A-Z]/){
					push(@aaidentical,$counting);
					$counting ++;	
				}
				#if similar amino acid	
				elsif($letter=~ /\+/){
					
					push(@aasimilar_L,"$counting".":"."$seqQ[$countGi-($firstQ)]\/$seqS[$countGi-($firstQ)]");
					$counting ++;
				}
				$countGi++;														
			}
			print $fho "$queryR;$longueurP;$id;$sim;$ev";
			print $fho ";".(join",",@aaidentical).";".(join",",@aasimilar_L).";".
			(join",",@aadiference_L).";".(join",",@gap).";".(join",",@gapOut).";".
			(join",",@aaStop)."\n";			
			$control++;
			$s=1;
			$sequence= '';
			$q = 1;
			$seqQ='';
			$seqS='';
			$queryR="";
		}		
	}	
}
close ($BLAST);
close ($fho);
}