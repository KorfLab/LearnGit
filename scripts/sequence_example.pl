#!/usr/bin/perl
use warnings;
use strict;
use Sequence;

#Open fasta file/files for import.   
my $fasta = open_fasta(["../Sequences/test_seq1.fa",
                        "../Sequences/test_seq2.fa",
                        "../Sequences/test_seq3.fa"]);

#Import one fasta file at a time and print to stdout
while (my $seq=get_next_fasta($fasta)){
    my $header = $seq->{HEADER};
    my $sequence = $seq->{SEQUENCE};
    print ">$header\n$sequence\n";
}


# how to shuffle a sequence
my $seq = "ACGTACGTACGT";
my $shuffled = shuffle_seq($seq);
print "\n$seq -> shuffle -> $shuffled\n\n";


#Open fasta file/files for import
$fasta = open_fasta("../Sequences/sequences.fa.gz");



#Import all at once
my $seqs = get_all_fastas($fasta);



#Print each sequence in $seqs
foreach my $sq (@$seqs){
    my $header = $sq->{HEADER};
    my $sequence = $sq->{SEQUENCE};
    print ">$header\n$sequence\n";
}

