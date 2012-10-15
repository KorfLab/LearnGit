use strict;
use warnings;
use List::Util qw(shuffle);
use IO::Uncompress::Gunzip;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
my %Translation = (
        'AAA' => 'K', 'AAC' => 'N', 'AAG' => 'K', 'AAT' => 'N',
        'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',
        'AGA' => 'R', 'AGC' => 'S', 'AGG' => 'R', 'AGT' => 'S',
        'ATA' => 'I', 'ATC' => 'I', 'ATG' => 'M', 'ATT' => 'I',
        'CAA' => 'Q', 'CAC' => 'H', 'CAG' => 'Q', 'CAT' => 'H',
        'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',
        'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',
        'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L',
        'GAA' => 'E', 'GAC' => 'D', 'GAG' => 'E', 'GAT' => 'D',
        'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',
        'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',
        'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',
        'TAA' => '*', 'TAC' => 'Y', 'TAG' => '*', 'TAT' => 'Y',
        'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S',
        'TGA' => '*', 'TGC' => 'C', 'TGG' => 'W', 'TGT' => 'C',
        'TTA' => 'L', 'TTC' => 'F', 'TTG' => 'L', 'TTT' => 'F'
);
my %Rev_Translation;
foreach my $nuc (keys %Translation) {
	$Rev_Translation{$Translation{$nuc}} = $nuc;
}

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
        $file = $files[0];
    }
    
    if (ref $file eq "ARRAY"){
        $fasta->{CURRENT_FILE}=0;
        _set_filehandle($fasta, $file->[0]);
		_init_fasta_first_line($fasta);
        $fasta->{FILES}=$file;
    }
    else{
        _set_filehandle($fasta, $file);
		_init_fasta_first_line($fasta);
    }
    
    return $fasta;
}

#Sets the filehandle for Fasta or Genbank files (Opens filehandle if filename is provided)
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
    
    #Make assignment of filehandle or open filehandle
    if (ref $file eq "GLOB"){
        $fasta->{FH} = $file;
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
                $fasta->{FH} = new IO::Uncompress::Gunzip $filename
                    or die "Couldn't open " . $filename . " for reading fasta file\n";
            }
            else{
                open my $fh , "<" , $filename or die "Couldn't open " . $filename . " for reading fasta file\n";
				$fasta->{FH} = $fh;
            }
        }
    }
    
    return 1; 
}


#Open fasta file and initialize the import of first line of fasta file
sub _init_fasta_first_line{
	my $fasta = shift;
	
	my $fh = $fasta->{FH}; #assign filehandle in FASTA
    
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
	return;
}


#Get the next sequence in the file and return the sequence and header
#Todo: test using FASTA format with > or ;\n;...\n  (Old format defined in FASTA)
sub get_next_fasta{
    my $fasta = shift;
    
    #Assign filehandle to temp variable
    my $fh = $fasta->{FH};
    
    #Check that filehandle is defined
    if (!defined $fh){
        warn "get_next_fasta() called on undefined filehandle\n";
        return;
    }
    
    #initialize sequence hash
    my $seq = {HEADER => $fasta->{LAST}, SEQUENCE => undef};

    while(<$fh>){
        chomp;
        if (/^;/){ #If line starts with ";" ignore the line
            next;
        }
        elsif (/^>/){ #If line starts with ">" assign header to LAST
            #$seq->{HEADER} = $fasta->{LAST};
			
			s/^>//;  #Remove ">" line
            $fasta->{LAST} = $_;
            return $seq;
        }
        else{
            $seq->{SEQUENCE}.=$_;  #Assign line to sequence
        }
    }
	
    
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
				_init_fasta_first_line($fasta);
            }
            else{
                $fasta->{CURRENT_FILE}=undef;
                $fasta->{FILES}=undef;
            }
        }
    }
    
    return $seq;
}

#Get all the sequences in the file and return a array_ref to array of
#hashes{sequences and headers}
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


#Open fasta file and check first line using _set_filehandle
#Can accept Typeglob, Lexical filehandle, single filename or multiple filenames (array or array_ref);
#Fasta can be either unzipped or gzipped fasta file
#Outputs a data structure that is used by get_next_fasta(...) and get_all_fastas(...)
sub open_genbank{
	my @files = @_;
    my $file_array_size = scalar @files;
    my $file;
    
    my $genbank;
    
    #Initialize Fasta data structure
    $genbank = {FILES => undef, CURRENT_FILE => undef, FH => undef, LAST => undef};
    
    if ($file_array_size==0){
        die "No File or Filehandles provided to function open_fasta().\n";
    }
    elsif ($file_array_size==1){
        $file = $files[0];
    }
    else{
        $genbank->{FILES}=\@files;
        $genbank->{CURRENT_FILE}=0;
        $file = $files[0];
    }
    
    if (ref $file eq "ARRAY"){
        $genbank->{CURRENT_FILE}=0;
        _set_filehandle($genbank, $file->[0]);
		_init_genbank_first_line($genbank);
        $genbank->{FILES}=$file;
    }
    else{
        _set_filehandle($genbank, $file);
		_init_genbank_first_line($genbank);
    }
    
    return $genbank;
}

#Open fasta file and initialize the import of first line of fasta file
sub _init_genbank_first_line{
	my $genbank = shift;
	
	my $fh = $genbank->{FH}; #assign filehandle in FASTA
    
    #Check first line for valid Fasta mark
    my $first_line = <$fh>;
    chomp $first_line;
    if ($first_line=~/^LOCUS/){ #First line could be > or ;
        $first_line=~s/^LOCUS\s+//;
        $genbank->{LAST}=$first_line;      #Store that line so we don't have to
                                        #re-read or move file position pointer.
    }
    else{
        die "Genbank File should start with first line \"LOCUS\".\nThis file format is not supported.\n
		The first line is :\n$first_line\n";
    }
	return;
}


#Get the next sequence in the file and return the sequence and header
#Todo: test using FASTA format with > or ;\n;...\n  (Old format defined in FASTA)
sub get_next_genbank{
    my $genbank = shift;
    
    #Assign filehandle to temp variable
    my $fh = $genbank->{FH};
    
    #Check that filehandle is defined
    if (!defined $fh){
        warn "get_next_fasta() called on undefined filehandle\nOpen";
        return;
    }
    
    #initialize sequence hash
    my $seq = {HEADER => $genbank->{LAST}, SEQUENCE => undef};
	
	my $seen_origin=0;

    while(<$fh>){
        chomp;
        if (/^LOCUS/){ #If line starts with ">" assign header to LAST
            s/^LOCUS\s+//;  #Remove ">" line
            #$seq->{HEADER} = $genbank->{LAST};
            $genbank->{LAST} = $_;
            return $seq;
        }
		elsif (/^ACCESSION/){
			s/^ACCESSION\s+//;
			$seq->{HEADER} = "ACCESSION\t" . $_  . "\tLOCUS\t" . $seq->{HEADER};
		}
		elsif (/^ORIGIN/){
			$seen_origin=1;
		}
		elsif (/^\/\//){
			$seen_origin=0;
		}
        elsif ($seen_origin){
			s/[0-9\s]+//g;
            $seq->{SEQUENCE}.=uc($_);  #Assign line to sequence
        }
		else{
			next;
		}
    }
    
    #If sequence and header have zero length return nothing
    if ((length $seq->{SEQUENCE}==0) && (length $seq->{HEADER}==0)){
        return;
    }
    
    #If EOF move to next file and initialize the header in _set_filehandle
    if (eof($fh)){
        $genbank->{LAST}=undef;
        if (defined $genbank->{FILES}){
            $genbank->{CURRENT_FILE}++;
            if ($genbank->{CURRENT_FILE} < scalar @{$genbank->{FILES}}){
                _set_filehandle($genbank, $genbank->{FILES}->[$genbank->{CURRENT_FILE}]);
				_init_genbank_first_line($genbank);
            }
            else{
                $genbank->{CURRENT_FILE}=undef;
                $genbank->{FILES}=undef;
            }
        }
    }
    
    return $seq;
}

#Get all the sequences in the genbank file and return a array_ref to array of
#hashes{SEQUENCE=>seq and HEADER=>header}
sub get_all_genbanks{
    my $genbank = shift;
    my $fh = $genbank->{FH};
    
    if (!defined $fh){
        warn "get_all_genbanks() called on undefined filehandle\n";
        return;
    }
    
    my $seqs = [];
    
    #Use get_next_fasta to import all fastas
    while(my $seq = get_next_fasta($genbank)){
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
	else {
		warn "An invalid value was provided for \$nuc_type. Please enter either 'DNA' or 'RNA'";
		return "Error";
		}
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

# Shuffle a DNA/RNA/protein sequence
sub shuffle_seq{
	my ($seq) = @_;
	my @seq = split(//, $seq);
	my @shuffled_seq = shuffle(@seq);
	return(join('', @shuffled_seq));
}

# reverse a sequence
sub reverse_seq {
	my ($seq) = @_;
	return (reverse($seq));
}


# Returns a random sequence of type $type and length $length. 
# $type (case insensitive) can be "dna" (default), "rna", "protein", or "custom".
# dna, rna, and protein will give equal weight to each alphabet (DNA - A/T/G/C, Protein - 20 amino acids).
# "custom" will have a third input hash reference $hash, which will have the alphabet as key and weight as value
# Please don't make the weight in crazy big number. Keep them less than 10 if possible.
sub rand_seq {
	my ($length, $type) = @_;
	my $ref = $_[2] if $type =~ /^custom$/i;
	my @possible_type = qw(dna rna protein custom);
	die "usage: rand_seq(length [int], type [dna, rna, protein, custom], [custom: hash of ref])\n" unless defined($length);
	die "Error at subroutine rand_seq: length must be positive integer (your input: $length).\n" unless $length =~ /^\d+$/;
	$type = "dna" if not defined($type);
	die "Error at subroutine rand_seq: type must be dna/rna/protein/custom/file (your input: $type).\n" unless grep(/^$type$/, @possible_type);
	die "Error at subroutine rand_seq: type is custom but you have undefined ref hash or ref hash contain zero keys/value.\n" if ($type =~ /^custom$/i and keys %{$ref} == 0);

	# Get reference based on type #
	my %ref;
	# DNA/RNA
	%ref = ("A" => 1, "a" => 1, "G" => 1, "g" => 1, "T" => 1, "t" => 1, "C" => 1, "c" => 1) if $type =~ /^dna$/i or $type =~ /^rna$/i;
	# Protein
	%ref = (
	"A" => 1, "R" => 1, "N" => 1, "D" => 1, "C" => 1, 
	"Q" => 1, "E" => 1, "G" => 1, "H" => 1, "I" => 1, 
	"L" => 1, "K" => 1, "M" => 1, "F" => 1, "P" => 1, 
	"S" => 1, "T" => 1, "W" => 1, "Y" => 1, "V" => 1
	) if $type =~ /^protein$/i;
	# Custom
	%ref = %{$ref} if $type =~ /^custom$/i;
	
	# Randomize #
	# First make a dummy sequence at $pool that has length $lmax (the bigger lmax, the more accurate)
	my ($lmax, $seq, $pool, $prob) = (99999);
	foreach my $alphabet (sort keys %ref) {
		my $value = int($ref{$alphabet} * $lmax);
		for (my $i = 0; $i < $value; $i++) {
			$pool .= $alphabet;
		}
	}
	
	# Then randomly take alphabets from $pool to make $seq 
	my $lseq = 0;
	while ($lseq  < $length) {
		$seq .= substr($pool, int(rand(length($pool))), 1);
		$lseq = length($seq);
	}
	$seq =~ tr/T/U/ if $type =~ /^rna$/i;
	return($seq);
}


# Take a case-insensitive sequence input $seq and optional start position $start (default 0)
# And there is also a third option of what should undefined be (X, N, 0)
# This was found somewhere at Ian/Keith's code so credit is theirs
sub translate_codon {
        my ($seq, $start, $undef) = @_;
	$start = 0 if not defined($start);
	$undef = 0 if not defined($undef);
	die "start position of translation has to be positive integer\n" if $start !~ /^\d+$/;

	$seq = uc($seq);
	$seq =~ s/U/T/; #RNA
        my $trans = "";
        for (my $i = $start; $i < length($seq); $i+=3) {
                my $codon = substr($seq, $i, 3);
                last if length($codon) < 3;
                if (not exists $Translation{$codon}) {$trans .= $undef}
                else                                 {$trans .= $Translation{$codon}}

        }
        return ($trans);
}

sub rev_translate_codon {
        my ($seq, $undef) = @_;
	$undef = "NNN" if not defined($undef);

	$seq = uc($seq);
	print "SEQ = $seq\n";
        my $trans = "";
        for (my $i = 0; $i < length($seq); $i++) {
                my $amino = substr($seq, $i, 1);
                if (not exists $Rev_Translation{$amino}) {$trans .= $undef}
                else                                     {$trans .= $Rev_Translation{$amino}}
        }
        return ($trans);
}

sub entropy_shannon {
	my ($seq, $case) = @_;
	$seq = uc($seq) if defined($case);
	my %seq;
	for (my $i = 0; $i < length($seq); $i++) {
		my $nuc = substr($seq, $i, 1);
		$seq{'nuc'}{$nuc}++;
		$seq{'tot'}++;
	}
	foreach my $nuc (keys %{$seq{'nuc'}}) {
		$seq{'nuc'}{$nuc} /= $seq{'tot'};
		$seq{'ent'} += -1 * $seq{'nuc'}{$nuc} * (log($seq{'nuc'}{$nuc}) / log (2));
	}
	return($seq{'ent'});
}


# Cleans a sequence(s) of white space or other characters that aren't valid. 
# Takes a sequence (string) and cleans it for unwanted characters
sub clean_sequence{
	my ($sequence) = @_;
	$sequence =~ s/\s//g;
	$sequence = uc($sequence);
	if ($sequence =~ m/[^ACGTURYSWKMBDHVN\*\-EFILPQ]/) {
		die "Nonstandard characters found in this sequence";
	}
	return ($sequence);
}



# Aligns two sequences using the Needleman-Wunsch Alignment algorithm 
sub nw_align{
	# get sequences passed values
	my ($seq1, $seq2, $MATCH, $MISMATCH, $GAP) = @_;

	# Initialize matrix
	my @matrix;
	$matrix[0][0]{score} = 0;
	$matrix[0][0]{pointer} = "none";
	for(my $j = 1; $j <= length($seq1); $j++) {
	    $matrix[0][$j]{score}   = $GAP * $j;
	    $matrix[0][$j]{pointer} = "left";
	}
	for (my $i = 1; $i <= length($seq2); $i++) {
	    $matrix[$i][0]{score}   = $GAP * $i;
	    $matrix[$i][0]{pointer} = "up";
	}

	# fill
	for(my $i = 1; $i <= length($seq2); $i++) {
	    for(my $j = 1; $j <= length($seq1); $j++) {
	        my ($diagonal_score, $left_score, $up_score);
	
	        # calculate match score
	        my $letter1 = substr($seq1, $j-1, 1);
	        my $letter2 = substr($seq2, $i-1, 1);                            
	        if ($letter1 eq $letter2) {
	            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
	        }
	        else {
	            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
	        }

	        # calculate gap scores
	        $up_score   = $matrix[$i-1][$j]{score} + $GAP;
	        $left_score = $matrix[$i][$j-1]{score} + $GAP;

	        # choose best score
	        if ($diagonal_score >= $up_score) {
	            if ($diagonal_score >= $left_score) {
	                $matrix[$i][$j]{score}   = $diagonal_score;
	                $matrix[$i][$j]{pointer} = "diagonal";
	            }
	        else {
	                $matrix[$i][$j]{score}   = $left_score;
	                $matrix[$i][$j]{pointer} = "left";
	            }
	        } else {
	            if ($up_score >= $left_score) {
	                $matrix[$i][$j]{score}   = $up_score;
	                $matrix[$i][$j]{pointer} = "up";
	            }
	            else {
	                $matrix[$i][$j]{score}   = $left_score;
	                $matrix[$i][$j]{pointer} = "left";
	            }
	        }
	    }
	}

	# trace-back
	my $align1 = "";
	my $align2 = "";

	# start at last cell of matrix
	my $j = length($seq1);
	my $i = length($seq2);

	my $loopflag = 1;
	while ($loopflag) {
	    if ($matrix[$i][$j]{pointer} eq "diagonal") {
	        $align1 .= substr($seq1, $j-1, 1);
	        $align2 .= substr($seq2, $i-1, 1);
	        $i--;
	        $j--;
	    }
	    elsif ($matrix[$i][$j]{pointer} eq "left") {
	        $align1 .= substr($seq1, $j-1, 1);
	        $align2 .= "-";
	        $j--;
	    }
	    elsif ($matrix[$i][$j]{pointer} eq "up") {
	        $align1 .= "-";
	        $align2 .= substr($seq2, $i-1, 1);
	        $i--;
	    }    
	    # ends at first cell of matrix
	    if ($matrix[$i][$j]{pointer} eq "none"){
	    	$loopflag = 0; 
	    }
	}

	$align1 = reverse $align1;
	$align2 = reverse $align2;
	return ($align1, $align2);
}


# extract_subsequence, extract and return subsequence from larger sequence
sub extract_subsequence{
	my ($seq, $start, $length) = @_;
	my $subseq;
	
	# need to have two versions, based on whether length is specified
	if (not defined $length){
		$subseq = substr($seq, $start);
		return($subseq);
	}

	# check that length of subseq won't exceed end of $seq, warn if so
	if (($start + $length) > length($seq)){
		warn "WARNING: You are trying to extract beyond range of sequence (\$start = $start, \$length = $length)\n";
	} 
	$subseq = substr($seq, $start, $length);
	
	return($subseq);
}


1;


