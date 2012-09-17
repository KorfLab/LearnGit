package getNCBI;
use strict;
use warnings;
use LWP::Simple;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(getNCBIGenbank getNCBIFasta);


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



