use strict;
use warnings;

my $NUCLEOTIDE="ACGTURYSWKMBDHVN._";
my $AMINO_ACID="ABCDEFGHIKLMNPQRSTVWXYZ.*";

###############################################################################
package Sequence;

#Note may think about keeping track of whether sequence is DNA,RNA, or Amino acid

#Create a new sequence object
sub new{
    my ($class,$header,$seq)=@_;
    my $self = bless {sequence=>undef,header=>undef}, $class;
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
sub sequence{
    my $self=shift;
    return $self->{sequence};
}

#Return a reference to the Sequence
sub sequence_ref{
    my $self=shift;
    return \$self->{sequence};
}

#Return the header
sub header{
    my $self=shift;
    return $self->{header};
}


#Create a copy of the Sequence object
sub copy{
    my $self=shift;
    my $copy = new Sequence($self->{header},$self->{sequence});
    return $copy;
}

#Return sequence as Fasta formatted string
sub asString{
    my $self=shift;
    return ">" . $self->{header} . "\n" . $self->{sequence} . "\n";
}

#print sequence (FASTA format)
sub print{
    my $self=shift;
    print $self->asString;
    return;
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

###############################################################################
#Sequences Package
#Array of sequence class
package Sequences;

sub new{
    my $class=shift;
    my $self = bless [], $class;
    return $self
}

sub push{
    my ($self,$seq)=@_;
    CORE::push @$self,$seq;
    return $self;
}

sub pop{
    my $self=CORE::shift;
    return CORE::pop @$self;
}

sub shift{
    my $self=CORE::shift;
    return CORE::shift @$self;
}

sub unshift{
    my ($self,$seq)=@_;
    CORE::unshift @$self,$seq;
    return $self;
}

sub size{
    my $self=CORE::shift;
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
package Fasta;
use IO::Uncompress::Gunzip;

#Import fasta formatted sequences individually or collectively from a filehandle
#, filename, or gzipped fasta file

#Create new Fasta object from File or GLOB
#implements gunzip compatibility using the IO::Uncompress::Gunzip core module
sub new {
    my ($class,$file) = @_;
    
    my $fh;
    if (ref $file eq "GLOB"){
        $fh=$file;
    }
    else{
        my $filename;
        if (ref $file eq "SCALAR"){
            $filename = $$file;
        }
        else{
            $filename=$file;
        }
        
        if (!-e $filename){
            die "$filename doesn't exist\n";
        }
        else{
            if ($filename=~/\.gz$/){
                $fh = new IO::Uncompress::Gunzip $filename
                    or die "Couldn't open " . $filename . " for reading fasta file\n";
            }
            else{
                open $fh , "<" , $filename or die "Couldn't open " . $filename . " for reading fasta file\n";
            }
        }
    }
    
    my $self = bless {}, $class;  #Bless anonymous hash as FASTA 
    $self->{FH}=$fh; #assign filehandle in FASTA
    
    #Check first line for valid Fasta mark
    my $first_line = <$fh>;
    chomp $first_line;
    if ($first_line=~/^>/ || $first_line=~/^;/){ #First line could be > or ;
        $first_line=~s/^>|^;//;
        $self->{LAST}=$first_line;      #Store that line so we don't have to
                                        #re-read or move file position pointer.
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
            s/^>|^;//;
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
    
    while(!eof($fh)){
        my $seq = $self->getNext;
        if ($seq->length==0){
            last;
        }
        else{
            $seqs->push($seq);
        }
    }
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
