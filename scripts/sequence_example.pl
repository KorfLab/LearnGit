#!/usr/bin/perl
use warnings;
use strict;
use Sequence;


#---------------------- Fasta Import ---------------------------#

#Open fasta file/files for import.   
my $fasta = Sequence::open_fasta("../sequences/test_seq1.fa",
                        "../sequences/test_seq2.fa",
                        "../sequences/test_seq3.fa");

#Import one fasta file at a time and print to stdout
while (my $seq = Sequence::get_next_fasta($fasta)){
    my $header = $seq->{HEADER};
    my $sequence = $seq->{SEQUENCE};
    print ">$header\n$sequence\n";
}

#Open fasta file/files for import
$fasta = Sequence::open_fasta("../sequences/sequences.fa.gz");

#Import all at once
my $seqs = Sequence::get_all_fastas($fasta);

#Print each sequence in $seqs
foreach my $sq (@$seqs){
    my $header = $sq->{HEADER};
    my $sequence = $sq->{SEQUENCE};
    print ">$header\n";
}


#---------------------- Genbank Import ---------------------------#
my $genbank = Sequence::open_genbank("../sequences/sequences.gbk.gz");

#Import one fasta file at a time and print to stdout
while (my $seq = Sequence::get_next_genbank($genbank)){
    my $header = $seq->{HEADER};
    my $sequence = $seq->{SEQUENCE};
    print ">$header\n";
}


#---------------------- Shuffle A Sequence ---------------------------#

# how to shuffle a sequence
my $seq = "ACGTACGTACGT";
my $shuffled = Sequence::shuffle_seq($seq);
print "\n$seq -> shuffle -> $shuffled\n\n";


#---------------------- Reverse A Sequence ---------------------------#

# reverse sequence
$seq = "AAAACCCCGGGGTTTT";
my $reversed = Sequence::reverse_seq($seq);
print "\n$seq -> reversed -> $reversed\n\n";


#---------------------- Create Random Seq ---------------------------#

# Make random sequence of 100 length
my $rand_seq_length = 1E2;
my $rand_dna = Sequence::create_rand_seq($rand_seq_length, "dna");
my $rand_rna = Sequence::create_rand_seq($rand_seq_length, "rna");
my $rand_pro = Sequence::create_rand_seq($rand_seq_length, "protein");

my %cust_ref = (
	"A" => 100, "T" => 100, "G" => 100, "C" => 100, "N" => 10,
	"k" => 1, "o" => 1, "r" => 1, "f" => 1, "l" => 1, "a" => 1, "b" => 1
	);
my $rand_custom = Sequence::create_rand_seq($rand_seq_length, "custom", \%cust_ref);

print "\n\nFunction rand_seq(\$length, \$type[dna rna protein custom], [optional: \$custom hash])\n";
print "$rand_seq_length length of create_rand_seq()\ndna: $rand_dna\nrna: $rand_rna\npro: $rand_pro
custom hash: 
my \%cust_ref = (
	\"A\" => 100, \"T\" => 100, \"G\" => 100, \"C\" => 100, \"N\" => 10,
	\"k\" => 1, \"o\" => 1, \"r\" => 1, \"f\" => 1, \"l\" => 1, \"a\" => 1, \"b\" => 1)
	\;
custom dna: $rand_custom\n";

#CHECK FOR RANDOMNESS of 100000bp seq#
$rand_seq_length = 1E5;
$rand_dna = Sequence::create_rand_seq($rand_seq_length, "dna");
my %rand;
my $ldna = length($rand_dna);
for (my $i = 0; $i < $ldna; $i++) {
	my $nuc = substr($rand_dna, $i, 1);
	$rand{$nuc}++;
	$nuc = substr($rand_dna, $i, 2) if $i < $ldna-1;
	$rand{$nuc}++ if $i < $ldna-1;
	$nuc = substr($rand_dna, $i, 3) if $i < $ldna-2;
	$rand{$nuc}++ if $i < $ldna-2;
}
my $DEBUG = 0;
if ($DEBUG == 1) {
	print "DEBUG create_rand_seq(): How random is the seq? (Tested on 1E5 bp DNA for 1-3 mer)
log odd ratio is a function of log(observed/expected), where if it's random enough, obs = exp, and value would be 0 (log(1) = 0)\n";
	foreach my $nuc (sort keys %rand) {
		$rand{$nuc} /= $ldna if length($nuc) == 1;
		$rand{$nuc} /= ($ldna-1) if length($nuc) == 2;
		$rand{$nuc} /= ($ldna-2) if length($nuc) == 3;
		printf "$nuc\t$rand{$nuc} (log odd = %.4f)\n", log($rand{$nuc} / 0.25**(length($nuc)));
	}
}
#----------------------- Create random seq based on kmer freq ------------------------#
print "\n\nFunction create_rand_seq_kmer(\$seq_length, \\\%kmer_table)\nGenerate \$seq_length bp sequence based on hash reference \\\%kmer_table\n";
###### NOT COMPLETE YET VVVVVVVVVVVVVVVVVVV
# How to use create_rand_seq_kmer:
#I have a table of a biased nucleotide frequencies
my %biased_nuc = (A => 0.5, T => 0.4, G => 0.07, C => 0.03);
# Based on the biased nuc freq table, I use create_rand_seq() to generate random sequence of 100k length
my $seq_kmer = Sequence::create_rand_seq(100000, "custom", \%biased_nuc);
# Then use my quick dirty count_kmer function to make a kmer hash reference table based on the sequence above
my $kmer = Sequence::count_kmer(3, $seq_kmer);
# Then based on the kmer hash reference table, I make sequence of 100k length
$seq_kmer = Sequence::create_rand_seq_kmer(100000, $kmer); # Final product #
my $seq_100bp = substr($seq_kmer, 0, 100);
print "create_rand_seq_kmer first 100bp result: $seq_100bp\n";

# DEBUG: As a test of how random it is, create_rand_seq_kmer output kmer table must be similar as input kmer table.
# I use 3 mer so that it's not printing crazy amount of kmer
my $kmer2 = Sequence::count_kmer(3, $seq_kmer);

my %kmer = %{$kmer};
my %kmer2 = %{$kmer2};
$DEBUG = 0;
if ($DEBUG == 1) {
	print "DEBUG create_rand_seq_kmer(): How random is the seq? (Tested on 100000 bp DNA for 3 mer)
log odd ratio is a function of log(observed/expected), where if it's random enough, obs = exp, and value would be 0 (log(1) = 0)
kmer\tIn\tOut\tlogodd\n";
	foreach my $nuc (sort {$kmer{$b} <=> $kmer{$a}} keys %kmer) {
		$kmer2{$nuc} = 0 if not defined($kmer2{$nuc});
		my $logodd = $kmer2{$nuc} == 0 ? "inf" : log($kmer{$nuc} / $kmer2{$nuc});
	        print "$nuc\t$kmer{$nuc}\t$kmer2{$nuc}\t$logodd\n";
	}
}
######## NOT COMPLETE YET ^^^^^^^^^^^^^^^^^^^^

#----------------------- Translate DNA seq to Protein ------------------------#
# Translate a sequence and reverse translate it
print "\n\nFunction translate_codon(\$dna): Translate a sequence\nFunction rev_translate_codon(\$protein): Reverse translate protein seq\n";
my $dnatoprotein = Sequence::translate_codon($rand_custom);
my $proteintodna = Sequence::rev_translate_codon($dnatoprotein);
print "dna:$rand_custom\ndna to protein: $dnatoprotein\nprotein to dna: $proteintodna\n\n";


#----------------------- Shannon Entropy ------------------------#
my $Sh_entropy = Sequence::entropy_shannon($rand_custom);
print "\n\nFunction entropy_shannon(\$dna): Calculate Shannon's Entropy of any sequence (e.g. DNA)\n";
print "DNA = $rand_custom\n";
print "Shannon Entropy of sequence rand_dna(\$rand_custom): $Sh_entropy\n";


#----------------------- Check for Whitespace ------------------------#
$seq = "AAAACCCCGGGGTTTT -AAA";
#$seq = "AAAACCCCGGGGTTTT -AAA!";
my $nowhitespace_seq = Sequence::clean_sequence($seq);
print "\n\nFunction clean_sequence(\$seq): Cleans whitespace from a DNA, RNA, or protein sequence\n";
print "\nSequence before cleaning:\n$seq\n\nSequence after cleaning:\n$nowhitespace_seq \n\n";


#---------------------- Check the Type of a Sequence ---------------------------#

$seq = "ACGTACGTACGT";
my $seq_type = Sequence::check_type($seq);
print "\n$seq -> check type -> $seq_type\n\n";


#-------------------- Needleman-Wunsch Alignment ---------------------#
my $alignseq1 = "AAATTTCCCGGG";
my $alignseq2 = "AATTCCCG";
my $MATCH = 1;
my $MISMATCH = -4;
my $GAP = -1;
my ($aligned1, $aligned2) = Sequence::nw_align($alignseq1, $alignseq2, $MATCH, $MISMATCH, $GAP);
print "\n\nFunction nw_align(\$seq): Aligns two sequences using the Needleman-Wunsch Alignment\n";
print "\nBefore Needleman-Wunsch Alignment\nSeq1:\t$alignseq1\nSeq2:\t$alignseq2\n";
print "\nAfter Needleman-Wunsch Alignment\nAligned Seq1:\t$aligned1\nAligned Seq2:\t$aligned2\n";


#----------------------- Extract subsequence ------------------------#
print "\n\n";
$seq = "TACAATCCGCCGGCGGTTTT";

my ($start, $length) = (2, 4);
print "Sequence: $seq\n";
print "Extracting from position $start, length $length\n";
my $subseq = Sequence::extract_subsequence($seq, $start, $length);
print "Subsequence: $subseq\n\n";

# negative offset
($start, $length) = (-5, 3);
print "Sequence: $seq\n";
print "Extracting from position $start, length $length\n";
$subseq = Sequence::extract_subsequence($seq, $start, $length);
print "Subsequence: $subseq\n\n";

# no length specified
$start  = 1;
print "Sequence: $seq\n";
print "Extracting from position $start\n";
$subseq = Sequence::extract_subsequence($seq, $start);
print "Subsequence: $subseq\n\n";


# Trying to extract too much sequence
($start, $length) = (5, 100);
print "Sequence: $seq\n";
print "Extracting from position $start, length $length\n";
$subseq = Sequence::extract_subsequence($seq, $start, $length);
print "Subsequence: $subseq\n\n";
