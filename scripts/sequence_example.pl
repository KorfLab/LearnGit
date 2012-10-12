#!/usr/bin/perl
use warnings;
use strict;
use Sequence;

#Open fasta file/files for import.   
my $fasta = open_fasta("../sequences/test_seq1.fa",
                        "../sequences/test_seq2.fa",
                        "../sequences/test_seq3.fa");

#Import one fasta file at a time and print to stdout
while (my $seq = get_next_fasta($fasta)){
    my $header = $seq->{HEADER};
    my $sequence = $seq->{SEQUENCE};
    print ">$header\n$sequence\n";
}


# how to shuffle a sequence
my $seq = "ACGTACGTACGT";
my $shuffled = shuffle_seq($seq);
print "\n$seq -> shuffle -> $shuffled\n\n";


# reverse sequence
$seq = "AAAACCCCGGGGTTTT";
my $reversed = reverse_seq($seq);
print "\n$seq -> reversed -> $reversed\n\n";

#Open fasta file/files for import
$fasta = open_fasta("../sequences/sequences.fa.gz");



#Import all at once
my $seqs = get_all_fastas($fasta);



#Print each sequence in $seqs
foreach my $sq (@$seqs){
    my $header = $sq->{HEADER};
    my $sequence = $sq->{SEQUENCE};
    print ">$header\n";
}

# Make random sequence of 1000 length
my $rand_dna = rand_seq(1000, "dna");
my $rand_rna = rand_seq(1000, "rna");
my $rand_pro = rand_seq(1000, "protein");

my %cust_ref = ("k" => 0.1, "o" => 0.2, "r" => 0.05, "f" => 0.1, "l" => 0.15, "a" => 0.25, "b" => 0.15);
my $rand_custom = rand_seq(1000, "custom", \%cust_ref);
my $rand_file   = rand_seq(1000, "file"  , "probtable.tsv");

print "dna; $rand_dna\nrna: $rand_rna\npro: $rand_pro\ncustom: $rand_custom\nfile:$rand_file\n";
