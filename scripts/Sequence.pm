use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
my $VERSION = 0.2;

#Open fasta file and check first line using _set_filehandle
#Can accept Typeglob, Lexical filehandle, single filename or multiple filenames (array or array_ref);
#Fasta can be either unzipped or gzipped fasta file
#Outputs a data structure that is used by get_next_fasta(...) and get_all_fastas(...)
sub open_fasta {
    my @files = @_;
    my $file_array_size = scalar @files;
    my $file;
    
    my $fasta;
    
    #Initialize Fasta data structure
    $fasta = {FILES => undef, CURRENT_FILE => undef, FH => undef, LAST => undef};
    
    if ($file_array_size==0){
        die "No File or Filehandles provided to function open_fasta().\n";
    }
    elsif ($file_array_size==1){
        $file = $files[0];
    }
    else{
        $fasta->{FILES}=\@files;
        $fasta->{CURRENT_FILE}=0;
        my $file = $files[0];
    }
    
    if (ref $file eq "ARRAY"){
        $fasta->{CURRENT_FILE}=0;
        _set_filehandle($fasta, $file->[0]);
        $fasta->{FILES}=$file;
    }
    else{
        _set_filehandle($fasta, $file);
    }
    
    return $fasta;
}

#Sets the filehandle for Fasta (Opens filehandle if filename is provided)
sub _set_filehandle{
    my ($fasta,$file)=@_;
    
    #Clear old filehandle and last header line
    if (exists $fasta->{FH} && defined $fasta->{FH}){
        close ($fasta->{FH});
        $fasta->{FH} = undef;
        $fasta->{LAST} = undef;
    }
    
    if (!defined $file){
        return 0;
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
    
   
    $fasta->{FH}=$fh; #assign filehandle in FASTA
    
    #Check first line for valid Fasta mark
    my $first_line = <$fh>;
    chomp $first_line;
    if ($first_line=~/^>/ || $first_line=~/^;/){ #First line could be > or ;
        $first_line=~s/^>|^;//;
        $fasta->{LAST}=$first_line;      #Store that line so we don't have to
                                        #re-read or move file position pointer.
    }
    else{
        die "Not FASTA format file\n";
    }
    
    return 1; 
}


#Get the next sequence in the file and return the sequence and header
#Todo: test using FASTA format with > or ;\n;...\n
sub get_next_fasta{
    my $fasta = shift;
    
    #Assign filehandle to temp variable
    my $fh = $fasta->{FH};
    
    #Check that filehandle is defined
    if (!defined $fh){
        warn "get_next_fasta() called on undefined filehandle\nOpen";
        return;
    }
    
    #initialize sequence hash
    my $seq = {HEADER => undef, SEQUENCE => undef};

    while(<$fh>){
        chomp;
        if (/^;/){ #If line starts with ";" ignore the line
            next;
        }
        elsif (/^>/){ #If line starts with ">" assign header to LAST
            s/^>//;  #Remove ">" line
            $seq->{HEADER} = $fasta->{LAST};
            $fasta->{LAST} = $_;
            return $seq;
        }
        else{
            $seq->{SEQUENCE}.=$_;  #Assign line to sequence
        }
    }
    
    #Handle when we reach EOF
    $seq->{HEADER} = $fasta->{LAST};
    
    #If sequence and header have zero length return nothing
    if ((length $seq->{SEQUENCE}==0) && (length $seq->{HEADER}==0)){
        return;
    }
    
    #If EOF move to next file and initialize the header in _set_filehandle
    if (eof($fh)){
        $fasta->{LAST}=undef;
        if (defined $fasta->{FILES}){
            $fasta->{CURRENT_FILE}++;
            if ($fasta->{CURRENT_FILE} < scalar @{$fasta->{FILES}}){
                _set_filehandle($fasta, $fasta->{FILES}->[$fasta->{CURRENT_FILE}]);    
            }
            else{
                $fasta->{CURRENT_FILE}=undef;
                $fasta->{FILES}=undef;
            }
        }
    }
    
    return $seq;
}

#Get all the sequences in the file and return a array_ref to array of hashes{sequences and headers}
sub get_all_fastas{
    my $fasta = shift;
    my $fh = $fasta->{FH};
    
    if (!defined $fh){
        warn "get_all_fastas() called on undefined filehandle\n";
        return;
    }
    
    my $seqs = [];
    
    #Use get_next_fasta to import all fastas
    while(my $seq = get_next_fasta($fasta)){
        push @$seqs,$seq;
    }
    
    return $seqs;
}

sub complement {
	# Returns the complement of the input $seq of type $nuc_type, which can be either "DNA" or "RNA"
    # Given no value for $nuc_type, $nuc_type defaults to "DNA"
	my ($seq, $nuc_type) = @_;
	$nuc_type //= "DNA"; # If $nuc_type is undefined, set to DNA
	if ($nuc_type eq "DNA") {$seq =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/}
	elsif ($nuc_type eq "RNA") {$seq =~ tr/AUGCYRKMBDHVaugcyrkmbdhv/UACGRYMKVHDBuacgrymkvhdb/}
	return $seq;
}

sub rev_comp {
	# Returns the reverse complement of the input $seq of type $nuc_type, which can be either "DNA" or "RNA"
    # Given no value for $nuc_type, $nuc_type defaults to "DNA"
	# Identical to complement(), except that the input sequence order is reversed before being complemented
	my ($seq, $nuc_type) = @_;
	$nuc_type //= "DNA"; # If $nuc_type is undefined, set to DNA
	$seq = reverse($seq);
	$seq = complement($seq, $nuc_type);
	return $seq;
}


1;