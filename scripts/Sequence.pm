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
        $genbank->{LAST}=$first_line;      #Store that line so we don't have to
                                        #re-read or move file position pointer.
    }
    else{
        die "Genbank File should start with first line \"LOCUS\".\nThis file format is not supported.\n
		The first line is :\n$first_line\n";
    }
	return;
}


#Get the next genbank record in the file and return the sequence, header, and
#annotation
#Input: variable from open_genbank(..)
#Output: hash ref with {SEQUENCE->..., HEADER->..., ANNOTATION->...}
sub get_next_genbank{
    my $genbank = shift;
    
    #Assign filehandle to temp variable
    my $fh = $genbank->{FH};
	if (eof($fh)){
		return;
	}
    
    #Check that filehandle is defined
    if (!defined $fh){
        warn "get_next_fasta() called on undefined filehandle\nOpen";
        return;
    }
    
    #initialize sequence assign Header information
	if ($genbank->{LAST} !~ /^LOCUS/){
		warn "Genbank file begins without LOCUS on first line of record\n" .
		"Line:\t" . $genbank->{LAST};
		$genbank->{LAST} = <$fh>;
		return;
	}
    
	#Parse genbank locus Field (first line of record)
	chomp $genbank->{LAST};
	my $seq = {HEADER => $genbank->{LAST}, SEQUENCE => undef};
	_add_genbank_locus($fh,$seq,$genbank->{LAST});
	$genbank->{LAST} = <$fh>;
	
	my $seen_origin=0;	
	
    while(defined $genbank->{LAST}){
		my $line = $genbank->{LAST};
        chomp $line;

		if ($line =~ /^SOURCE/){
			$genbank->{LAST} = _add_genbank_source($fh,$seq,$line);
		}
		elsif ($line =~ /^REFERENCE/){
			$genbank->{LAST} = _add_genbank_reference($fh,$seq,$line);
		}
		elsif ($line =~ /^FEATURES/){
			$genbank->{LAST} = _add_genbank_feature($fh,$seq);
		}
		elsif ($line =~ /^ORIGIN/){
			$seen_origin=1;
			$genbank->{LAST} = <$fh>;
			chomp $genbank->{LAST};
		}
		elsif ($line =~ /^\/\//){  #end of record
			if (!eof($fh)){
				$genbank->{LAST} = <$fh>;
				chomp $genbank->{LAST};
			}
			
			last;
		}
        elsif ($seen_origin){
			$line =~ s/[0-9\s]+//g;
            $seq->{SEQUENCE}.=uc($line);  #Assign line to sequence
			$genbank->{LAST} = <$fh>;
			chomp $genbank->{LAST};
        }
		elsif ($line=~/^(\w+)/){  #Import general field information
			$genbank->{LAST} = _add_genbank_general($fh,$seq,$1,$line);
		}
		else{
			$genbank->{LAST} = <$fh>;
			chomp $genbank->{LAST};
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

sub _add_genbank_locus{
	my ($fh,$seq,$line) = @_;
	$line =~ s/^LOCUS\s+//;  #Remove ">" line
	$line =~ m/^(\S+)\s+(\d+)\s+bp\s+(\S+)\s+(linear|circular)*\s+(PRI|ROD|MAM|VRT|INV|PLN|BCT|VRL|PHG|SYN|UNA|EST|PAT|STS|GSS|HTG|HTC|ENV)\s+(\S+)$/g;
	$seq->{ANNOTATION}->{LOCUS}->{NAME} = $1;
	$seq->{ANNOTATION}->{LOCUS}->{LENGTH} = $2;
	$seq->{ANNOTATION}->{LOCUS}->{MOLECULE_TYPE} = $3;	
	$seq->{ANNOTATION}->{LOCUS}->{MOLECULE_STRUCT} = $4;
	$seq->{ANNOTATION}->{LOCUS}->{GENBANK_DIVISION} = $5;
	$seq->{ANNOTATION}->{LOCUS}->{MODIFICATION_DATE} = $6;
	
	return;
}

#Parse the Record Field and any lines after that don't have a field ID
#Modifies $seq and adds annotation to the Annotation hash
#Returns the last line read from the file
#Must return last line because gzipped files can't use seek.
sub _add_genbank_general{
	my ($fh,$seq,$type,$line,) = @_;
	$line =~ s/^$type\s+//;
	$seq->{ANNOTATION}->{$type} = $line;
	
	$line = <$fh>;
	chomp $line;
	
	while ($line=~/^\s+/){
		$line =~ s/^\s+//;
		$seq->{ANNOTATION}->{$type} .= " " . $line;
		$line = <$fh>;
		chomp $line;
	}
		
	return $line;
}


#Parse the Source Field of the record and the sub-field Organism
#Modifies $seq and adds annotation to the Annotation hash
#Returns the last line read from the file
sub _add_genbank_source{
	my ($fh,$seq,$line) = @_;
	$line =~ s/^SOURCE\s+//;
	$seq->{ANNOTATION}->{SOURCE} = $line;
	
	my $position = tell $fh;
	$line = <$fh>;
	chomp $line;
	
	if ($line=~/^\s+ORGANISM/){
		$line =~ s/^\s+ORGANISM\s+//;
		$line = _add_genbank_general($fh,$seq,"ORGANISM",$line);
	}
	
	return $line;
}

#Parse the Refernce ID field and associated sub-fields
#Modifies $seq and adds annotation to the Annotation hash
#Returns the last line read from the file
sub _add_genbank_reference{
	my ($fh,$seq,$line) = @_;
	$line =~ s/^REFERENCE\s+//;
	$line =~ m/^(\d+)(\s+(\(.*\)))?/;
	my $ref_index = $1;
	my $ref_coord = (defined $3) ? $3 : q();
	
	$seq->{ANNOTATION}->{REFERENCES}->[$ref_index]->{COORD} = $ref_coord;
	
	$line = <$fh>;
	chomp $line;
	
	my $parent = q();
	
	while ($line=~/^s+/){
		if ($line=~/^\s+AUTHORS/){
			$parent = "AUTHORS";
			$line =~ s/^\s+AUTHORS\s+//;
			$seq->{ANNOTATION}->{REFERENCES}->[$ref_index]->{AUTHORS} = $line;
		}
		elsif ($line=~/^\s+CONSRTM/){
			$parent = "TITLE";
			$line =~ s/^\s+CONSRTM\s+//;
			$seq->{ANNOTATION}->{REFERENCES}->[$ref_index]->{CONSRTM} = $line;
		}
		elsif ($line=~/^\s+TITLE/){
			$parent = "TITLE";
			$line =~ s/^\s+TITLE\s+//;
			$seq->{ANNOTATION}->{REFERENCES}->[$ref_index]->{TITLE} = $line;
		}
		elsif ($line=~/^\s+JOURNAL/){
			$parent = "JOURNAL";
			$line =~ s/^\s+JOURNAL\s+//;
			$seq->{ANNOTATION}->{REFERENCES}->[$ref_index]->{JOURNAL} = $line;
		}
		elsif ($line=~/^\s+PUBMED/){
			$parent = "PUBMED";
			$line =~ s/^\s+PUBMED\s+//;
			$seq->{ANNOTATION}->{REFERENCES}->[$ref_index]->{PUBMED} = $line;
		}
		elsif ($line=~/^REFERENCE/){
			last;
		}
		else{
			$line=~s/^\s+//;
			$seq->{ANNOTATION}->{REFERENCES}->[$ref_index]->{$parent} .= " " . $line;
		}
		
		$line = <$fh>;
		chomp $line;
	}
	
	return $line;
}


#Parse the Feature fields and associated sub-fields
#Adds features to array of Features in Annotation hash
#Modifies $seq and adds annotation to the Annotation hash
#Returns the last line read from the file
sub _add_genbank_feature{
	my ($fh,$seq) = @_;
	
	my $line = <$fh>;
	chomp $line;
	
	while ($line=~/^\s{5}(\S+)\s+(\S+)/){
		my $feature;
		$feature->{FEATURE}=$1;
		$feature->{COORD}=$2;
		$line = _parse_genbank_feature($fh,$feature);
		
		push @{$seq->{ANNOTATION}->{FEATURES}},$feature;
		chomp $line;
	}
	
	return $line;
	
}


#Parse individual features sub-fields
#Modifies $feature hash_ref
#Returns the last line read from the file
sub _parse_genbank_feature{
	my ($fh,$feature) = @_;
	my $line = <$fh>;
	chomp $line;
	
	my $parent = "COORD";
	
	while ($line =~ m/^\s{21}(.+)/){
		my $group = $1;
		
		if ($group =~ /^\/(\S+)\=\"(.*)\"$/){
			$parent = uc($1);
			$feature->{$parent} = $2;
		}
		else{
			$line =~ s/^\s+//;
			$feature->{$parent} .= " " . $line;
		}
		
		$line = <$fh>;
		chomp $line;
	}
	
	return $line;
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
# $type (case insensitive) can be "dna" (default), "rna", "protein", "custom", or "file".
# dna, rna, and protein will give equal probability to each alphabet (DNA - A/T/G/C has 25%, Protein - 1/22).
# "custom" will have a third input hash reference $hash, which will have the alphabet as key and probability as value
# "file" will have a third option $fh, which is a tab/comma delimited format file with alphabet at column 1 with equal 
# probability, unless you define the frequency/probability yourself at col 2, e.g:
# A	0.01
# B	0.99
# or
# A,9
# B,1
# C,100
sub rand_seq {
	my ($length, $type, $extra) = @_;
	my @possible_type = qw(dna rna protein custom file);
	die "usage: rand_seq(length [int], type [dna, rna, protein, custom ,file], extra)\n" unless defined($length);
	die "Error at subroutine rand_seq: length must be positive integer (your input: $length).\n" unless $length =~ /^\d+$/;
	$type = "dna" if not defined($type);
	die "Error at subroutine rand_seq: type must be dna/rna/protein/custom/file (your input: $type).\n" unless grep(/^$type$/, @possible_type);
	die "Error at subroutine rand_seq: type is custom but you have undefined ref hash or ref hash contain zero keys/value.\n" if ($type =~ /^custom$/i and keys %{$extra} == 0);
	die "Error at subroutine rand_seq: type is file but you give no file input.\n" if ($type =~ /^file$/i and not defined($extra));

	# Get reference based on type #
	my %ref;
	# DNA/RNA
	%ref = ("A" => 0.25, "G" => 0.25, "T" => 0.25, "C" => 0.25) if $type =~ /^dna$/i or $type =~ /^rna$/i;
	# Protein
	%ref = (
	"A" => 1/20, "R" => 1/20, "N" => 1/20, "D" => 1/20, "C" => 1/20, 
	"Q" => 1/20, "E" => 1/20, "G" => 1/20, "H" => 1/20, "I" => 1/20, 
	"L" => 1/20, "K" => 1/20, "M" => 1/20, "F" => 1/20, "P" => 1/20, 
	"S" => 1/20, "T" => 1/20, "W" => 1/20, "Y" => 1/20, "V" => 1/20
	) if $type =~ /^protein$/i;
	# Custom
	%ref = %{$extra} if $type =~ /^custom$/i;
	# File
	if ($type =~ /^file$/i) {
		open (my $in, "<", $extra) or die "Error at subroutine rand_seq: Cannot read from $extra: $!.\n";
		my @ref = <$in>;
		chomp(@ref);
		my $print_once = 0;
		for (my $i = 0; $i < @ref; $i++) {
			my $ref = $ref[$i];
			$ref = uc($ref);
			my ($alphabet, $value) = split(",", $ref) if $ref =~ /\,/;
			   ($alphabet, $value) = split("\t", $ref) if $ref =~ /\t/;
			   ($alphabet)	       = $ref if $ref !~ /\,/ and $ref !~ /\t/;

			# If the file is delimited by another format (space, underline) we tell user
			if (not defined($value) and $alphabet =~ /.+\d+\.\d+/ and $print_once == 0) {
				print "Possible error at subroutine rand_seq: Looks like you have file format error, file must be tab/comma delimited. Example correct Format:\nALPHABET1\t0.01\nALPHABET2\t0.99\n";
				$print_once = 1;
			}

			elsif (defined($value)) {
				# To be honest, the program will work correctly even if "value" is not probability.
				# It just makes more sense for me (also my OCD) that the value is in probability format.
				# Feel free to comment this die function out if you don't care
				#die "Error at subroutine rand_seq: Probability must be positive fraction between 0-1! (your input: $alphabet $value)\n" unless ($value eq 0 or $value eq 1 or $value =~ /^\d*\.\d+$/ or $value =~ /^\d+\.*\d*\/\d+\.*\d*$/) and $value <= 1 and $value >= 0;
				
				# Nevermind, it's better to use weighted instead of probability
				$ref{$alphabet} = $value;
			}
			else {
				$ref{$alphabet} = 1 / @ref;
			}			
		}
	}

	
	# Randomize #
	
	# Pre-check to make sure probability is not more than 1. Feel free to comment the die function below if you don't care
	#my $probtest = 0;
	#foreach my $alphabet (sort keys %ref) {
	#	$probtest += $ref{$alphabet};
	#}
	#die "Error at subroutine rand_seq: Total probability is more than 1 ($probtest)\n" unless $probtest < 1 or $probtest eq 1;
	# Since this is weighted we don't need probability.

	# FIrst make a dummy sequence at $pool that has length $lmax (the bigger lmax, the more accurate)
	my ($lmax, $seq, $pool, $prob) = (999999);
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
        my ($seq, $start, $undef) = @_;
	$start = 0 if not defined($start);
	$undef = "NNN" if not defined($undef);
	die "start position of translation has to be positive integer\n" if $start !~ /^\d+$/;

	$seq = uc($seq);
	print "SEQ = $seq\n";
        my $trans = "";
        for (my $i = $start; $i < length($seq); $i++) {
                my $amino = substr($seq, $i, 1);
                if (not exists $Rev_Translation{$amino}) {$trans .= $undef}
                else                                     {$trans .= $Rev_Translation{$amino}}
        }
        return ($trans);
}

# Cleans a sequence(s) of white space or other characters that aren't valid. 
# Takes a sequence (string) and cleans it for unwanted characters
sub clean_sequences{
	my ($sequence) = @_;
	$sequence =~ s/\w//;
	$sequence = uc($sequence);
#	$sequence ~= s/[a-z]/[A-Z]/;	
	if($sequence =~ m/[^ACGTURYSWKMBDHVN.-EFILPQ]/){
		die "Nonstandard characters found in this sequence";
	}
	return ($sequence);
}

1;


