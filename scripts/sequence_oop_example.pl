#!/usr/bin/perl

# 

# Name:        <script_name.pl>

# Description: <one line description>

# Author:      <your name>

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

" unless @ARGV == 0; # adjust as necessary



my %opt;

getopts('hp:', \%opt);



my ($var1, $var2) = @ARGV;



# global variables and other initializations



#############

# main loop #

#############


#----------------------------Generate random sequence----------------------------------#
my %biased_kmer = ("A" => 0.5, "T" => 0.3, "G" => 0.05, "C" => 0.15);
my $DNA = Sequence_OOP::generate_random_sequence(100000, "custom", "DNA", \%biased_kmer);

#_---------------------------Generate kmer sequence-----------------------------------#
my %kmer = %{Sequence_OOP::Sequence::count_kmer(3, $DNA->{sequence})};
$DNA = Sequence_OOP::generate_kmer_sequence(100000, \%kmer);
my $DEBUG = 0;
#DEBUG#
#How random is random?
if ($DEBUG==1) {
	my %kmer2 = %{Sequence_OOP::Sequence::count_kmer(3, $DNA->{sequence})};
	print "kmer\tIN\tOUT\tLOGODD\n";
	foreach my $kmer (sort {$kmer{$b} <=> $kmer{$a}} keys %kmer) {
	        my $logodd = $kmer2{$kmer} == 0 ? "inf" : log($kmer{$kmer}/$kmer2{$kmer});
		        print "$kmer\t$kmer{$kmer}\t$kmer2{$kmer}\t$logodd\n";
	}
}
#------------------------------Translate----------------------------------------------#
#my $Translate = Sequence_OOP::Sequence::translate($DNA);





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

