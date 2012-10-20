#!/usr/bin/perl
#
# Name:        <script_name.pl>
# Description: <one line description>
# Author:      <your name>
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

die "
usage: <program_name> [options] <arguments...>
options:
  -h
  -p <parameter>
"; # adjust as necessary

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
