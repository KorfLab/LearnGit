use strict;
use warnings;

my $NUCLEOTIDE="ACGTURYSWKMBDHVN._";
my $AMINO_ACID="ABCDEFGHIKLMNPQRSTVWXYZ.*";

#Need to check passed data type in the functions to make sure it is a Sequence
#or a compatible datatype.


###############################################################################
package Sequence;

#Note may think about keeping track of whether sequence is DNA,RNA, or Amino acid

#Create a new sequence object
sub new{
    my ($class,$header,$seq)=@_;
    my $self = bless {sequence=>undef,header=>undef}, $class;
    if (defined $header){
        $self->set_header($header);
    }
    
    if (defined $seq){
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
    
    if (ref $sequence eq "SCALAR"){
        $self->{sequence}=$$sequence;
    }
    elsif (ref $sequence){
        warn "Trying to set_sequence with non-scalar type\n";
    }
    else{
        $self->{sequence}=$sequence;
    }
    return $self;
}

#Set the header of the sequence
sub set_header{
    my ($self,$header)=@_;
    
    if (ref $header eq "SCALAR"){
        $self->{header}=$$header;
    }
    elsif (ref $header){
        warn "Trying to set_header with non-scalar type\n";
    }
    else{
        $self->{header}=$header;
    }
    
    return $self;
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

#Create an new Sequences type
sub new{
    my $class=shift;
    my $self = bless [], $class;
    
    return $self
}


#Push a Sequence into Sequences
sub push_seq{
    my ($self,$seq)=@_;
    if (ref $seq ne "Sequence"){
        warn "Tried to push a non-Sequence type on to Sequences\n";
        return $self;
    }
    push @$self,$seq;
    return $self;
}

#Pop a Sequence out of Sequences
sub pop_seq{
    my $self=shift;
    return pop @$self;
}

#Shift a Sequence out of Sequences
sub shift_seq{
    my $self=shift;
    return shift @$self;
}

#Unshift a Sequence into Sequences
sub unshift_seq{
    my ($self,$seq)=@_;
    if (ref $seq ne "Sequence"){
        warn "Tried to unshift a non-Sequence type on to Sequences\n";
        return $self;
    }
    
    unshift @$self,$seq;
    return $self;
}

#Get the number of sequences in Sequences
sub size{
    my $self=shift;
    return scalar @$self;
}

#Get the Sequence at the ith position
sub at{
    my ($self,$i)=@_;
    if ($i>=@$self){
        warn "Position $i: Uninitialized Sequence within Sequences\n";
        return new Sequence();
    }
    return $self->[$i];
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
# Format specifications:
#       http://en.wikipedia.org/wiki/FASTA_format
#       http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml

################ Constructor ################

#Create new Fasta object from File or GLOB
#implements gunzip compatibility using the IO::Uncompress::Gunzip core module
sub new {
    my ($class,$file) = @_;
    my $self = bless {FH=>undef, LAST=>undef, FILES=>undef, CURRENT_FILE=>undef}, $class;  #Bless anonymous hash as FASTA 
    if (defined $file){
        if (ref $file eq "ARRAY"){
            $self->{FILES}=$file;
            $self->{CURRENT_FILE}=0;
            $self->set_filehandle($self->{FILES}->[0]);
        }
        else{
            $self->set_filehandle($file);
        }
    }
    return $self;
}


################ Methods ################

#Sets the filehandle for Fasta (Opens filehandle if filename is provided)
sub set_filehandle{
    my ($self,$file)=@_;
    
    #Clear old filehandle and last header line
    if (exists $self->{FH} && defined $self->{FH}){
        close ($self->{FH});
        $self->{FH} = undef;
        $self->{LAST} = undef;
    }
    
    if (!defined $file){
        return;
    }
    
    my $fh;
    
    #Make assignment of filehandle or open filehandle
    if (ref $file eq "GLOB"){
        $fh=$file;
    }
    else{
        my $filename;
        if (ref $file eq "SCALAR"){
            $filename = $$file;  #Get filename from scalar_ref
        }
        else{
            $filename=$file;
        }
        
        #Check to make sure file exists and is readable.  Open filehandle
        if (!-e $filename && !-r $filename){
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
    
    if (!defined $fh){
        warn "Fasta::getNext() called on undefined filehandle\n";
        return;
    }
    
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
    
    $seq->set_header($self->{LAST});
    
    if ((length $seq->{sequence}==0) && (length $seq->{header}==0)){
        return;
    }
    
    if (eof($fh)){
        $self->{LAST}=undef;
        if (defined $self->{FILES}){
            $self->{CURRENT_FILE}++;
            if ($self->{CURRENT_FILE}< scalar @{$self->{FILES}}){
                $self->set_filehandle($self->{FILES}->[$self->{CURRENT_FILE}]);    
            }
            else{
                $self->{CURRENT_FILE}=undef;
                $self->{FILES}=undef;
            }
        }
    }
    
    return $seq;
}

#Get all the sequences in the file and return a Sequences object
sub get_all_fastas{
    my $self = shift;
    my $fh = $self->{FH};
    
    if (!defined $fh){
        warn "Fasta::getAll() called on undefined filehandle\n";
        return;
    }
    
    my $seqs = new Sequences();
    
    while(my $seq = $self->getNext){
        $seqs->push_seq($seq);
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
