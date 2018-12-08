#!/usr/bin/perl -w

use strict;
use warnings;
use Config;

# INPUT:    a text file with the protein TABLE from NCBI
#           a text file with the protein FASTA from NCBI (the whole genome)
# ToDO:     combine the table with the fasta
# OUTPUT:   a text file with the locus tag and the associated sequence
# HINT:     NCBI - table has the following format:
#           Location | Strand | Length | PID | Gene | Synonym | Code | COG | Product
#           NCBI - fasta has the following format:
#           >gi|GeneID|ref|...|...

# usage-definition
my $usage= << "USE";
  USAGE:   perl "combine_TABLE_FASTA.pl" protein_fasta protein_table path
  options: none
  example: perl "combine_TABLE_FASTA.pl" protein_fasta.txt protein_table.txt projectName/Fasta
USE

# get the information from the program call
my $requiredArgs = 3;
my $argcount     = $#ARGV+1;
$argcount eq $requiredArgs or die "$usage";

my $path= pop(@ARGV);
my $table = pop(@ARGV);
my $fasta = pop(@ARGV);


print $Config{osname}."\n";

my $calcul=calculate_fasta ();
exit;

sub calculate_table {
  my %combination;
  my $TABLE;
	my ($locus,$refseq);
	open($TABLE, $table);

    while (<$TABLE>) {
		my $line = $_;
		if ($line =~ /^CDS/) {
			my @line_array = split("\t", $line);
			$refseq=$line_array[10];
			#print $refseq."\n";
			$locus= $line_array[16];

            $combination{$refseq}= $locus;
        }
    }
    close ($TABLE);
    return %combination;
}

# calculate_fasta - parses a NCBI fasta file
# @return - hash with the gi number and the sequence

sub calculate_fasta {
    my $FASTA;
    my $sequence;
	  my %combination = calculate_table ();
    open($FASTA, $fasta);
    my $output_name;
  	if(($Config{osname} eq "linux")||($Config{osname} eq "darwin")){
  		$output_name = "$path\/combined.fasta";
  	}
  	else{
  		$output_name = "$path\\combined.fasta";
  	}

    my $OUTPUT;
	  open($OUTPUT, ">$output_name");
    my $count=0;
    
    while (<$FASTA>) {
      my $a;
      my $line = $_;

      # regex which looks for the header line in fasta-file:
      if ($line =~ />/ ){
        my $var=(split)[0];
		    my $gi=substr $var, 1;
        if(defined $sequence){
          print $OUTPUT $sequence."\n";
        }
        if (exists $combination{$gi}){
          
          $a = $combination{$gi};
				  $count+=1;
				  print $OUTPUT ">$a\n";
          $sequence = "";
        }       
      }
      elsif ($line =~ /(\w+)/) {
        # concatenate each sequence line
        $sequence = $sequence.$1;
      }

    }
    print "File with $count proteins\n";
    
	print $OUTPUT $sequence."\n";
    close ($FASTA);
    close ($OUTPUT);
}
