package Fasta;
use strict;
use warnings;
use IO::Uncompress::Gunzip;
use Sequence_OOP;

#Delete these comments and provide POD for documentation


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


1;
