package Sequence;
use strict;
use warnings;

sub new {
    my ($class,$file) = @_;
    
    my $fh;
    if (ref $file ne "GLOB"){
        if (ref $file ne "SCALAR"){
            if (!-e $file){
                die "$file doesn't exist\n";
            }
            else{
                open $fh , "<" , $file or die "Couldn't open " . $file . " for reading fasta file\n";
            }
        }
        else{
            open $fh , "<" , $$file or die "Couldn't open " . $$file . " for reading fasta file\n";
        }
    }
    else{
        $fh = $file;
    }
    
    my $self = bless {}, $class;
    $self->{FH}=$fh;
    
    #Check first line
    my $first_line = <$fh>;
    
    
    return $self;
}



1;
