#!/usr/bin/perl
use warnings;
use strict;
use Sequence;


#---------------------- Fasta Import ---------------------------#

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


#---------------------- Genbank Import ---------------------------#
my $genbank = open_genbank("../sequences/sequences.gbk.gz");

#Import one fasta file at a time and print to stdout
while (my $seq = get_next_genbank($genbank)){
    my $header = $seq->{HEADER};
    my $sequence = $seq->{SEQUENCE};
    print "$header\n";
}


#---------------------- Shuffle A Sequence ---------------------------#

# how to shuffle a sequence
my $seq = "ACGTACGTACGT";
my $shuffled = shuffle_seq($seq);
print "\n$seq -> shuffle -> $shuffled\n\n";


#---------------------- Reverse A Sequence ---------------------------#

# reverse sequence
$seq = "AAAACCCCGGGGTTTT";
my $reversed = reverse_seq($seq);
print "\n$seq -> reversed -> $reversed\n\n";


#-------------------- Generate Ranom Sequence -_----------------------#

# Make random sequence of 1000 length
my $rand_dna = rand_seq(100, "dna");
my $rand_rna = rand_seq(100, "rna");
my $rand_pro = rand_seq(100, "protein");

<<<<<<< HEAD
my %cust_ref = ("k" => 10, "o" => 0.5, "r" => 0.5, "f" => 0.5, "l" => 1, "a" => 2, "b" => 2);
my $rand_custom = rand_seq(1000, "custom", \%cust_ref);

print "dna; $rand_dna\nrna: $rand_rna\npro: $rand_pro\ncustom: $rand_custom\n";
=======
my %cust_ref = ("k" => 0.1, "o" => 0.2, "r" => 0.05, "f" => 0.1, "l" => 0.15, "a" => 0.25, "b" => 0.15);
my $rand_custom = rand_seq(100, "custom", \%cust_ref);
my $rand_file   = rand_seq(100, "file"  , "probtable.tsv");
print "\n\nFunction rand_seq(\$length, \$type[dna rna protein custom file], [NA or \$custom hash], [optional: File]): Random sequence generator\n";
print "dna: $rand_dna\nrna: $rand_rna\npro: $rand_pro\ncustom: $rand_custom\nfile:$rand_file\n";

# Translate a sequence and reverse translate it
my $dnatoprotein = translate_codon($rand_dna);
my $proteintodna = rev_translate_codon($dnatoprotein);
print "\n\nFunction translate_codon(\$dna): Translate a sequence\nFunction rev_translate_codon(\$protein): Reverse translate protein seq\n";
print "dna:$rand_dna\ndna to protein: $dnatoprotein\nprotein to dna: $proteintodna\n\n";

#----------------------- Check for Whitespace ------------------------#

# reverse sequence
$seq = "AAAACCCCGGGG
TTTT -AAA";
#$seq = "AAAACCCCGGGGTTTT -AAA!";
my $nowhitespace_seq = clean_sequences($seq);
print "\nSequence before cleaning:\n$seq\n\nSequence after cleaning:\n$nowhitespace_seq \n\n";


