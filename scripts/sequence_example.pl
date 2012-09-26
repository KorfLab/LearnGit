use warnings;
use strict;
use Sequences;

#Create a new Fasta type by handing it the filename
my $fasta = new Fasta("../Sequences/sequences.fa.gz");
my $seq = $fasta->getNext;  #imports first sequence of fasta file as Sequence type
my $seqs = $fasta->getAll;  #imports the remaining fasta as Sequences type


#Get and print header and sequence from Sequence Type
my $header = $seq->header;  #Get header string
my $sequence = $seq->sequence;  #Get sequence string
my $seq_ref  = $seq->sequence_ref; #Get reference to sequence string
print "$header\n$sequence\n" . $$seq_ref . "\n";

#Print Sequence header and sequence in Fasta format directly
$seq->print();

#Iterate through Sequences and print them twice
for(my $i=0;$i<$seqs->size;$i++){
    my $sequence = $seqs->[$i];  #Get ith sequence and print it
    $sequence->print;
    
    $sequence = $seqs->pop;  #Pop last sequence off
    $sequence->print; #print the sequence
}
