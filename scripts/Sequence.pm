package Sequence;
use strict;
use warnings;
use List::Util qw(shuffle);
use IO::Uncompress::Gunzip;
use Inline C => 'DATA';


my %Translation = (
	'AAA' => 'K',	'AAC' => 'N',	'AAG' => 'K',	'AAT' => 'N',
	'AAR' => 'K',	'AAY' => 'N',	'ACA' => 'T',	'ACC' => 'T',
	'ACG' => 'T',	'ACT' => 'T',	'ACR' => 'T',	'ACY' => 'T',
	'ACK' => 'T',	'ACM' => 'T',	'ACW' => 'T',	'ACS' => 'T',
	'ACB' => 'T',	'ACD' => 'T',	'ACH' => 'T',	'ACV' => 'T',
	'ACN' => 'T',	'AGA' => 'R',	'AGC' => 'S',	'AGG' => 'R',
	'AGT' => 'S',	'AGR' => 'R',	'AGY' => 'S',	'ATA' => 'I',
	'ATC' => 'I',	'ATG' => 'M',	'ATT' => 'I',	'ATY' => 'I',
	'ATM' => 'I',	'ATW' => 'I',	'ATH' => 'I',	'CAA' => 'Q',
	'CAC' => 'H',	'CAG' => 'Q',	'CAT' => 'H',	'CAR' => 'Q',
	'CAY' => 'H',	'CCA' => 'P',	'CCC' => 'P',	'CCG' => 'P',
	'CCT' => 'P',	'CCR' => 'P',	'CCY' => 'P',	'CCK' => 'P',
	'CCM' => 'P',	'CCW' => 'P',	'CCS' => 'P',	'CCB' => 'P',
	'CCD' => 'P',	'CCH' => 'P',	'CCV' => 'P',	'CCN' => 'P',
	'CGA' => 'R',	'CGC' => 'R',	'CGG' => 'R',	'CGT' => 'R',
	'CGR' => 'R',	'CGY' => 'R',	'CGK' => 'R',	'CGM' => 'R',
	'CGW' => 'R',	'CGS' => 'R',	'CGB' => 'R',	'CGD' => 'R',
	'CGH' => 'R',	'CGV' => 'R',	'CGN' => 'R',	'CTA' => 'L',
	'CTC' => 'L',	'CTG' => 'L',	'CTT' => 'L',	'CTR' => 'L',
	'CTY' => 'L',	'CTK' => 'L',	'CTM' => 'L',	'CTW' => 'L',
	'CTS' => 'L',	'CTB' => 'L',	'CTD' => 'L',	'CTH' => 'L',
	'CTV' => 'L',	'CTN' => 'L',	'GAA' => 'E',	'GAC' => 'D',
	'GAG' => 'E',	'GAT' => 'D',	'GAR' => 'E',	'GAY' => 'D',
	'GCA' => 'A',	'GCC' => 'A',	'GCG' => 'A',	'GCT' => 'A',
	'GCR' => 'A',	'GCY' => 'A',	'GCK' => 'A',	'GCM' => 'A',
	'GCW' => 'A',	'GCS' => 'A',	'GCB' => 'A',	'GCD' => 'A',
	'GCH' => 'A',	'GCV' => 'A',	'GCN' => 'A',	'GGA' => 'G',
	'GGC' => 'G',	'GGG' => 'G',	'GGT' => 'G',	'GGR' => 'G',
	'GGY' => 'G',	'GGK' => 'G',	'GGM' => 'G',	'GGW' => 'G',
	'GGS' => 'G',	'GGB' => 'G',	'GGD' => 'G',	'GGH' => 'G',
	'GGV' => 'G',	'GGN' => 'G',	'GTA' => 'V',	'GTC' => 'V',
	'GTG' => 'V',	'GTT' => 'V',	'GTR' => 'V',	'GTY' => 'V',
	'GTK' => 'V',	'GTM' => 'V',	'GTW' => 'V',	'GTS' => 'V',
	'GTB' => 'V',	'GTD' => 'V',	'GTH' => 'V',	'GTV' => 'V',
	'GTN' => 'V',	'TAA' => '*',	'TAC' => 'Y',	'TAG' => '*',
	'TAT' => 'Y',	'TAR' => '*',	'TAY' => 'Y',	'TCA' => 'S',
	'TCC' => 'S',	'TCG' => 'S',	'TCT' => 'S',	'TCR' => 'S',
	'TCY' => 'S',	'TCK' => 'S',	'TCM' => 'S',	'TCW' => 'S',
	'TCS' => 'S',	'TCB' => 'S',	'TCD' => 'S',	'TCH' => 'S',
	'TCV' => 'S',	'TCN' => 'S',	'TGA' => '*',	'TGC' => 'C',
	'TGG' => 'W',	'TGT' => 'C',	'TGY' => 'C',	'TTA' => 'L',
	'TTC' => 'F',	'TTG' => 'L',	'TTT' => 'F',	'TTR' => 'L',
	'TTY' => 'F',	'TRA' => '*',	'YTA' => 'L',	'YTG' => 'L',
	'YTR' => 'L',	'MGA' => 'R',	'MGG' => 'R',	'MGR' => 'R',
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
# $type (case insensitive) can be "dna" (default), "rna", "protein", or "custom".
# dna, rna, and protein will give equal weight to each alphabet (DNA - A/T/G/C, Protein - 20 amino acids).
# "custom" will have a third input hash reference $hash, which will have the alphabet as key and weight as value
# Please don't make the weight in crazy big number. Keep them less than 10 if possible.
sub rand_seq {
	my ($length, $type) = @_;
	my $ref = $_[2] if $type =~ /^custom$/i;
	my @possible_type = qw(dna rna protein custom);
	die "usage: rand_seq(length [int], type [dna, rna, protein, custom], [custom: hash of ref])\n" unless defined($length);
	die "Error at subroutine rand_seq: length must be positive integer (your input: $length).\n" unless $length =~ /^\d+E{0,1}\d*$/i;
	$type = "dna" if not defined($type);
	die "Error at subroutine rand_seq: type must be dna/rna/protein/custom/file (your input: $type).\n" unless grep(/^$type$/, @possible_type);
	die "Error at subroutine rand_seq: type is custom but you have undefined ref hash or ref hash contain zero keys/value.\n" if ($type =~ /^custom$/i and keys %{$ref} == 0);

	# Get reference based on type #
	my %ref;
	# DNA/RNA
	%ref = ("A" => 1, "G" => 1, "T" => 1, "C" => 1) if $type =~ /^dna$/i or $type =~ /^rna$/i;
	#%ref = ("A" => 1, "a" => 1, "G" => 1, "g" => 1, "T" => 1, "t" => 1, "C" => 1, "c" => 1) if $type =~ /^dna$/i or $type =~ /^rna$/i; #case insensitive
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
	# Stepped CDF #
	my %cdf;

	# Get total of weights #
	my $total = 0;
	foreach my $char (sort {$ref{$a} <=> $ref{$b}} keys %ref) {
		my $value = $ref{$char};
		$total += $value;
	}

	# Get steps of weights #
	my $curr_val = 0;
	foreach my $char (sort {$ref{$a} <=> $ref{$b}} keys %ref) {
		$cdf{$char}{'min'} = $curr_val;
		$curr_val += $ref{$char}/$total;
		$cdf{$char}{'max'} = $curr_val;
	}
	
	# Make sequence #
	my ($lseq, $seq) = 0;
	while ($lseq < $length) {
		my $rnum = rand();
		foreach my $char (sort {$cdf{$a} <=> $cdf{$b}} keys %cdf) {
			if ($rnum > $cdf{$char}{'min'} and $rnum <= $cdf{$char}{'max'}) {
				$seq .= $char;
				$lseq++;
				last;
			}
		}
	}

	# Correct for RNA #
	$seq =~ tr/Tt/Uu/ if $type =~ /^rna$/i;
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
        for (my $i = $start; $i < length($seq)-3; $i+=3) {
                my $codon = substr($seq, $i, 3);
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

# Calculate Shannon's Entropy
# Default is case insensitive
sub entropy_shannon {
	my ($seq, $case) = @_;
	$seq = uc($seq) if not defined($case);
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

__DATA__
__C__

void align_sw (const char *s1, const char *s2, double m, double n, double g) {
	int i, j, l1, l2, max_i, max_j, min_i, min_j, pos, alen;
	double d, v, h, max_s = 0;
	double ** sm;
	char ** tm;
	char *a1, *a2, *a3;
		
	l1 = strlen(s1), l2 = strlen(s2);
	
	/* allocate matrices */
	sm = malloc(sizeof(double*) * (l1+1));
	tm = malloc(sizeof(double*) * (l1+1));
	for (i = 0; i <= l1; i++)  sm[i] = malloc(sizeof(double) * (l2+1));
	for (i = 0; i <= l1; i++)  tm[i] = malloc(sizeof(double) * (l2+1));
		
	/* init first column and row */
	for (i = 0; i <= l1; i++) sm[i][0] = 0, tm[i][0] = '.';
	for (j = 1; j <= l2; j++) sm[0][j] = 0, tm[0][j] = '.';
		
	/* fill matrix */
	max_i = max_j = 0;
	for (i = 1; i <= l1; i++) {
		for (j = 1; j <= l2; j++) {
			if (s1[i-1] == s2[j-1]) d = sm[i-1][j-1] + m;
			else                    d = sm[i-1][j-1] + n;
			v = sm[i-1][ j ] + g;
			h = sm[ i ][j-1] + g;
						
			if (d > h && d > v && d > 0) sm[i][j] = d, tm[i][j] = 'd';
			else if (h > v && h > 0)     sm[i][j] = h, tm[i][j] = 'h';
			else if (v > 0)              sm[i][j] = v, tm[i][j] = 'v';
			else                         sm[i][j] = 0, tm[i][j] = '.';
			
			if (d > max_s) max_s = d, max_i = i, max_j = j;
		}
	}
	
	/* allocate alignment strings */
	a1 = malloc((l1 + l2 + 1) * sizeof(char));
	a2 = malloc((l1 + l2 + 1) * sizeof(char));
	a3 = malloc((l1 + l2 + 1) * sizeof(char));
	
	/* traceback */
	i = max_i, j = max_j, pos = 0;
	while (sm[i][j] > 0) {
		if (i == 0 || j == 0) {
			fprintf(stderr, "found 1\n");
		}
		if (pos >= (l1 + l2)) {
			fprintf(stderr, "found 2\n");
		}
		a1[pos] = s1[i-1]; 
		a2[pos] = s2[j-1];
		min_i = i;
		min_j = j;
		pos++;
		if      (tm[i][j] == 'd') i--, j--;
		else if (tm[i][j] == 'h') j--;
		else if (tm[i][j] == 'v') i--;
	}
	a1[pos] = '\0', a2[pos] = '\0', a3[pos] = '\0';
			
	/* reverse strings */
	for (i = 0; i < strlen(a1); i++) a3[pos -i -1] = a1[i];
	for (i = 0; i < strlen(a2); i++) a1[pos -i -1] = a2[i];
	alen = strlen(a1);
		
	/* create return values */
	Inline_Stack_Vars;
	Inline_Stack_Reset;
	Inline_Stack_Push(sv_2mortal(newSVnv(max_s)));
	Inline_Stack_Push(sv_2mortal(newSVpv(a3, alen)));
	Inline_Stack_Push(sv_2mortal(newSVpv(a1, alen)));
	Inline_Stack_Push(sv_2mortal(newSViv(min_i)));
	Inline_Stack_Push(sv_2mortal(newSViv(max_i)));
	Inline_Stack_Push(sv_2mortal(newSViv(min_j)));
	Inline_Stack_Push(sv_2mortal(newSViv(max_j)));
	Inline_Stack_Done;
		
	/* clean up */
	free(a1); free(a2); free(a3);
	for (i = 0; i <= l1; i++) free(sm[i]), free(tm[i]);
	free(sm); free(tm);
	
}

