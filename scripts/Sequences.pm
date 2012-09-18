use strict;
use warnings;

my $NUCLEOTIDE="ACGTURYSWKMBDHVN._";
my $AMINO_ACID="ABCDEFGHIKLMNPQRSTVWXYZ.*";



###############################################################################
package Fasta;

#Create new Fasta object from File or GLOB
#Todo: implement gzip compatibility
sub new {
    my ($class,$file) = @_;
    
    my $fh;
    if (ref $file eq "GLOB"){
        $fh=$file;
    }
    elsif (ref $file eq "SCALAR"){
        open $fh , "<" , $$file or die "Couldn't open " . $$file . " for reading fasta file\n";
    }
    else{
        if (!-e $file){
            die "$file doesn't exist\n";
        }
        else{
            open $fh , "<" , $file or die "Couldn't open " . $file . " for reading fasta file\n";
        }
    }
    
    my $self = bless {}, $class;
    $self->{FH}=$fh;
    
    #Check first line for valid Fasta mark
    my $first_line = <$fh>;
    chomp $first_line;
    if ($first_line=~/^>/ || $first_line=~/^;/){
        $self->{LAST}=$first_line;
        seek($fh,0,0);
    }
    else{
        die "Not FASTA format file\n";
    }
    
    return $self;
}

#Get the next sequence in the file and return a Sequence object
#Todo: test using FASTA format with > or ;\n;...\n
sub getNext{
    my $self = shift;
    my $fh = $self->{FH};
    my $seq = new Sequence();
    while(<$fh>){
        chomp;
        if (/^;/){
            next;
        }
        elsif (/^>/){
            $seq->set_header($self->{LAST});
            $self->{LAST}=$_;
            return $seq;
        }
        else{
            $seq->{sequence}.=$_;
        }
    }
    
    if (length $seq->{sequence}==0){
        return;
    }
    
    $seq->set_header($self->{LAST});
    return $seq;
}

#Get all the sequences in the file and return a Sequences object
sub getAll{
    my $self = shift;
    my $fh = $self->{FH};
    my $seqs = new Sequences();
    my $seq = new Sequence();
    while(<$fh>){
        chomp;
        if (/^;/){
            next;
        }
        elsif (/^>/){
            $seq->set_header($self->{LAST});
            $self->{LAST}=$_;
            $seqs->push_sequence($seq);
            $seq=new Sequence();
        }
        else{
            $seq->{sequence}.=$_;
        }
    }
    
    if (length $seq->{sequence}==0){
        return $seqs;
    }
    
    $seq->set_header($self->{LAST});
    $seqs->push_sequence($seq);
    return $seqs;
}

###############################################################################
package Genbank;

#Create a new Genbank object
sub new{
    
}

#Get the next sequence in the file and return a Sequence object
sub getNext{
    
}

#Get all the sequences in the file and return a Sequences object
sub getAll{
    
}





###############################################################################
package Sequence;

#Note may think about keeping track of whether sequence is DNA,RNA, or Amino acid

#Create a new sequence object
sub new{
    my ($class,$header,$seq)=@_;
    my $self = bless {sequence=>"",header=>""}, $class;
    if (defined $header){
        $self->set_header($header);
        $self->set_sequence($seq);
    }
    return $self
}

#Return the lenght of the sequence
sub length{
    my $self=shift;
    return length $self->{sequence};
}

#Set the sequence
sub set_sequence{
    my ($self,$sequence)=@_;
    $self->{sequence}=$sequence;
    return;
}

#Set the header of the sequence
sub set_header{
    my ($self,$header)=@_;
    $self->{header}=$header;
    return;
}

#Return the sequence
sub get_sequence{
    my $self=shift;
    return $self->{sequence};
}

#Return a reference to the Sequence
sub get_sequence_ref{
    my $self=shift;
    return \$self->{sequence};
}

#Return the header
sub get_header{
    my $self=shift;
    return $self->{header};
}


#Create a copy of the Sequence object
sub copy{
    my $self=shift;
    my $copy = new Sequence($self->{header},$self->{sequence});
    return $copy;
}

#Clean the sequence of invalid and whitespace characters
sub clean{
    
}

#Check sequence is valid DNA,RNA, or amino acid sequence
sub check{
    
}

#Reverse sequence and return a new Sequence object
sub reverse{

}

#Complement the sequence and return a new Sequence object
sub complement{
    
}

#Reverse complement the sequence and return a new Sequence object
sub reverse_complement{
    
}

#Extract subsequence given a range and return a new Sequence object
sub extract{
    
}

#Translate sequence
sub translate{
    
}

#Count given k-mer size and return hash {AAA=>10,AAC=>2...}
sub count_kmer{
    
}

#Match Regex or string to sequence and return array of matching positions
sub pattern_match{
    
}

#Shuffle sequence
sub shuffle{
    
}

#print sequence (FASTA format)
sub print{
    
}

###############################################################################
#Sequences Package
#Array of sequence
package Sequences;

sub new{
    my $class=shift;
    my $self = bless [], $class;
    return $self
}

sub push_sequence{
    my ($self,$seq)=@_;
    push @$self,$seq;
    return $self;
}

sub pop_sequence{
    my $self=shift;
    return pop @$self;
}

sub size{
    my $self=>shift;
    return scalar @$self;
}

sub clean{
    
}

sub check{
    
}

sub reverse{

}

sub complement{
    
}

sub reverse_complement{
    
}

sub extract{
    
}

sub translate{
    
}

sub count_kmer{
    
}

sub pattern_match{
    
}

sub shuffle{
    
}

sub print{
    
}


###############################################################################
package main;

sub generate_random_sequence{
    
}

sub generate_kmer_sequence{
    
}

sub global_alignment{
    
}

sub local_alignment{
    
}

sub calculate_entropy{
    
}

sub calculate_relative_entropy{
    
}

1;
