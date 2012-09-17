#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use LWP::Simple;

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

################ Subroutines ##########################

#Get NCBI Genbank file
#Query NCBI given an accession.
#Then returns a Genbank as scalar
sub getNCBIGenbank{
    my $acc = shift;
    
    my $accn = $acc . "[accn]";
    
    #assemble the esearch URL
    my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    my $url = $base . "esearch.fcgi?db=nucleotide&term=$accn&usehistory=n";
    
    #post the esearch URL
    my $output = get($url);
    
    #parse WebEnv and QueryKey
    my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
    my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
    
    #assemble the efetch URL
    $url = $base . "efetch.fcgi?db=nucleotide&query_key=$key&WebEnv=$web";
    $url .= "&rettype=gb&retmode=text";
    
    #post the efetch URL
    my $genbank = get($url);
    
    if ($genbank=~/\<\?xml/){
        return;
    }
    
    return $genbank;
}


#Get NCBI FASTA file
#Query NCBI given an accession.
#Then returns a FASTA as scalar
sub getNCBIFasta{
    my $acc = shift;
    
    my $accn = $acc . "[accn]";
    
    #assemble the esearch URL
    my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    my $url = $base . "esearch.fcgi?db=nucleotide&term=$accn&usehistory=n";
    
    #post the esearch URL
    my $output = get($url);
    
    #parse WebEnv and QueryKey
    my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
    my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
    
    #assemble the efetch URL
    $url = $base . "efetch.fcgi?db=nucleotide&query_key=$key&WebEnv=$web";
    $url .= "&rettype=fasta&retmode=text";
    
    #post the efetch URL
    my $fasta = get($url);
    
    if ($fasta=~/^</){
        return;
    }
    
    return $fasta;
}