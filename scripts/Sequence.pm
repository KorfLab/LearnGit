use strict;
use warnings;
use List::Util qw(sum);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
my $VERSION = 0.1;

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

sub kl_distance {
	# These are two hash references which contain the kmer counts for two different sequences
	my ($P_ref,$Q_ref) = @_;
	
	#If any value does not exist or equals 0
	map { if(!exists $P_ref->{$_}){$P_ref->{$_}=0;}} keys %$Q_ref;
	map { if(!exists $Q_ref->{$_}){$Q_ref->{$_}=0;}} keys %$P_ref;
	
	#If either of the above conditions is met, each distribution will be increased by 1 to normalize the additional or already existing zero value
	foreach my $count (keys %$P_ref) {
		if ($P_ref->{$count}==0 or $Q_ref->{$count}==0) {
			foreach my $position (keys %$P_ref) {
				$P_ref->{$position} += 1;
				$Q_ref->{$position} += 1;
			}
		}
	}
	
	my $sumP = sum(values %$P_ref);
	my $sumQ = sum(values %$Q_ref);
	
	my %P_distn;
	my %Q_distn;
	
	#Within the @P_distn array, each key is assigned a value that is set equal to the positional value ($_) divided by the sum(using the map command)
	@P_distn{keys %$P_ref} = map {$_ /= $sumP} values %$P_ref;
	@Q_distn{keys %$Q_ref} = map {$_ /= $sumQ} values %$Q_ref;
	
	#Create hash for KL distance
	my %KL;
	
	my $kl_pq = 0;
	my $kl_qp = 0;
	foreach my $seq (keys %P_distn){
		$kl_pq += $P_distn{$seq} + log ($P_distn{$seq}/$Q_distn{$seq});
		$kl_qp += $Q_distn{$seq} + log ($Q_distn{$seq}/$P_distn{$seq});
	}
	
	$KL{"D(P||Q)"} = $kl_pq;
	$KL{"D(Q||P)"} = $kl_qp;
	
	return (\%KL);
}


1;