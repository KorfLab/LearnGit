#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use getNCBI;

die "Usage: $0 <FASTA Filename> <Genbank Filename> [NCBI accessions]" unless @ARGV>=3;

my $fasta_file = shift @ARGV;
my $genbank_file = shift @ARGV;

open FASTA, ">>$fasta_file" or die "Couldn't open $fasta_file for appending.\n";
open GENBANK, ">>$genbank_file" or die "Couldn't open $genbank_file for appending.\n";

foreach my $acc (@ARGV){
    my $fasta_sequence  = getNCBIFasta($acc);
    my $gbk_sequence    = getNCBIGenbank($acc);
    
    if (defined $fasta_sequence){
        $fasta_sequence=~s/\n*$//;
        print FASTA $fasta_sequence . "\n";
    }
    else{
        print "Couldn't get $acc FASTA record from NCBI\n";
    }
    
    if (defined $gbk_sequence){
        $gbk_sequence=~s/\n*$//;
        print GENBANK $gbk_sequence . "\n";
    }
    else{
        print "Couldn't get $acc GENBANK record from NCBI\n";
    }
}