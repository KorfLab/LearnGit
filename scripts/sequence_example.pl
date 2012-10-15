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

my %cust_ref = ("k" => 10, "o" => 0.5, "r" => 0.5, "f" => 0.5, "l" => 1, "a" => 2, "b" => 2);
my $rand_custom = rand_seq(1000, "custom", \%cust_ref);

print "\n\nFunction rand_seq(\$length, \$type[dna rna protein custom], [NA or \$custom hash]\n";
print "dna; $rand_dna\nrna: $rand_rna\npro: $rand_pro\ncustom: $rand_custom\n";

#----------------------- Translate DNA seq to Protein ------------------------#
# Translate a sequence and reverse translate it
my $dnatoprotein = translate_codon($rand_dna);
my $proteintodna = rev_translate_codon($dnatoprotein);
print "\n\nFunction translate_codon(\$dna): Translate a sequence\nFunction rev_translate_codon(\$protein): Reverse translate protein seq\n";
print "dna:$rand_dna\ndna to protein: $dnatoprotein\nprotein to dna: $proteintodna\n\n";

#----------------------- Shannon Entropy ------------------------#
my $entropy = entropy_shannon($rand_dna);
my $entropy_case = entropy_shannon($rand_dna, "case");
print "\n\nFunction entropy_shannon(\$dna, [case]): Calculate Shannon's Entropy of a dna sequence\n";
print "Shannon Entropy of sequence rand_dna():\ncase insensitive: $entropy\ncase sensitive: $entropy_case\n";

#----------------------- Check for Whitespace ------------------------#
$seq = "AAAACCCCGGGGTTTT -AAA";
#$seq = "AAAACCCCGGGGTTTT -AAA!";
my $nowhitespace_seq = clean_sequences($seq);
print "\n\nFunction clean_sequence(\$seq): Cleans whitespace from a DNA, RNA, or protein sequence\n";
print "\nSequence before cleaning:\n$seq\n\nSequence after cleaning:\n$nowhitespace_seq \n\n";

#-------------------- Needleman-Wunsch Alignment ---------------------#
my $alignseq1 = "AAATTTCCCGGG";
my $alignseq2 = "AATTCCCG";
my $MATCH = 1;
my $MISMATCH = -4;
my $GAP = -1;
my ($aligned1, $aligned2) = nw_align($alignseq1, $alignseq2, $MATCH, $MISMATCH, $GAP);
print "\n\nFunction nw_align(\$seq): Aligns two sequences using the Needleman-Wunsch Alignment\n";
print "\nBefore Needleman-Wunsch Alignment\nSeq1:\t$alignseq1\nSeq2:\t$alignseq2\n";
print "\nAfter Needleman-Wunsch Alignment\nAligned Seq1:\t$aligned1\nAligned Seq2:\t$aligned2\n";
