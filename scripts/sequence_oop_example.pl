#!/usr/bin/perl
#
# Name:        sequence_oop_example.pl
# Description: A script to demonstrate the usage of Sequence_OOP functions
# Author:      Korf Lab
# License:     <BSD, MIT, GPL, CC, PD, choose one>

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Sequence_OOP;

################
# command line #
################

my %opt;
getopts('hp:', \%opt);

my ($var1, $var2) = @ARGV;

# global variables and other initializations

#############
# main loop #
#############


#----------------------------  Create a new Sequence  ---------------------------------#

my $sequence    = new Sequence();  #empty sequence;
$sequence       = new Sequence("ACGTACGT");
$sequence       = new Sequence("ACGTACGT", TYPE => "DNA");
$sequence       = new Sequence(SEQUENCE => "ACGTACGTACGT");
$sequence       = new Sequence(HEADER   => "Empty Sequence");
$sequence       = new Sequence(HEADER   => "This is my Sequence",
                               SEQUENCE => "ACGTACGTACGT");
$sequence       = new Sequence(SEQUENCE => "CGGAAATTTGGT",
                               TYPE     => "PRO");

#---------------------------- Format Sequence as Fasta ------------------------#
print $sequence->as_fasta_string . "\n";
$sequence->print_fasta;
#----------------------------  Using sequence accessor method -----------------#

print $sequence->sequence . "\n";   #print sequence

$sequence->sequence("ACGTACGT");    #assign sequence
print $sequence->sequence ."\n";    #print ACGTACGT\n

print $sequence . "\n";     #print using overloaded '""'
print length $sequence;     #print length using overloaded '""'
print "\n";

#---------------------------  Using header accessor method --------------------#

$sequence->header("My Header");     #Setting header
print $sequence->header . "\n";     #Getting header

#---------------------------  Using length method  ----------------------------#

print $sequence->length . "\n";     #get length using Sequence::length

#---------------------------  Using Sequence Concatenation --------------------#
my $first_seq   = new Sequence ("ACGT");    #create new Sequence
my $second_seq  = new Sequence ("NNNN");    #create new Sequence
my $string      = "AUUG";  #create string

$first_seq->concatenate("__TT");
print $first_seq->sequence . "\t=\tACGT__TT\n";   #print $first_seq

#Using overloaded operator '.='
$first_seq .= $string;  #Concatenate AUUG to ACGTACGT
print $first_seq->sequence . "\t=\tACGT__TTAUUG\n";   #print $first_seq


#Using concatenate_copy
my $new_seq = $first_seq->concatenate_copy($second_seq);
print $new_seq->sequence . "\t=\tACGT__TTAUUGNNNN\n";   #print $new_seq

#Using overloaded operator '.'
$new_seq = $first_seq . $second_seq;
print $new_seq->sequence . "\t=\tACGT__TTAUUGNNNN\n";   #print $new_seq

#Concatenating with string using overloaded '.'
$new_seq = "ata_" . $first_seq;
print $new_seq->sequence . "\t=\tata_ACGT__TTAUUG\n";   #print $new_seq


$new_seq = $first_seq . "_tat";
print $new_seq->sequence . "\t=\tACGT__TTAUUG_tat\n";   #print $new_seq


#-------------------------- Using Sequence copy -------------------------------#
$new_seq = $first_seq->copy;  #Copy $first_seq with copy subroutine

$new_seq = $first_seq;  #Copy $first_seq using overloaded '='

$new_seq .= "ACGT";

print "Original Sequence: $first_seq\n";
print "Copied Sequence after concat: " . $new_seq . "\n";

#----------------------------Generate random sequence----------------------------------#
my %biased_kmer = ("A" => 0.5, "T" => 0.3, "G" => 0.05, "C" => 0.15);
my $DNA = Sequence_OOP::generate_random_sequence(100000, "custom", "DNA", \%biased_kmer);

#_---------------------------Generate kmer sequence-----------------------------------#
my %kmer = %{Sequence_OOP::Sequence::count_kmer(3, $DNA->{sequence})};
$DNA = Sequence_OOP::generate_kmer_sequence(100000, \%kmer);

#------------------------------Translate----------------------------------------------#
my $Translate = Sequence_OOP::Sequence::translate($DNA);

#----------------------------------Clean----------------------------------------------#
my $dirtyDNA = "ASERTWB D\nT";
$DNA->clean;


############
# clean up #
############

###############
# subroutines #
###############

sub foo {}

BEGIN {} # remove if you don't need
END {} # remove if you don't need

=head1 Documentation

Include extra documentation here if you like

=cut

__END__

free form comments and thoughts
