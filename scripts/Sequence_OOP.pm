use strict;
use warnings;

my $NUCLEOTIDE="ACGTURYSWKMBDHVN._";
my $AMINO_ACID="ABCDEFGHIKLMNPQRSTVWXYZ.*";
my %Translation = (
        'AAA' => 'K',   'AAC' => 'N',   'AAG' => 'K',   'AAT' => 'N',
        'AAR' => 'K',   'AAY' => 'N',   'ACA' => 'T',   'ACC' => 'T',
        'ACG' => 'T',   'ACT' => 'T',   'ACR' => 'T',   'ACY' => 'T',
        'ACK' => 'T',   'ACM' => 'T',   'ACW' => 'T',   'ACS' => 'T',
        'ACB' => 'T',   'ACD' => 'T',   'ACH' => 'T',   'ACV' => 'T',
        'ACN' => 'T',   'AGA' => 'R',   'AGC' => 'S',   'AGG' => 'R',
        'AGT' => 'S',   'AGR' => 'R',   'AGY' => 'S',   'ATA' => 'I',
        'ATC' => 'I',   'ATG' => 'M',   'ATT' => 'I',   'ATY' => 'I',
        'ATM' => 'I',   'ATW' => 'I',   'ATH' => 'I',   'CAA' => 'Q',
        'CAC' => 'H',   'CAG' => 'Q',   'CAT' => 'H',   'CAR' => 'Q',
        'CAY' => 'H',   'CCA' => 'P',   'CCC' => 'P',   'CCG' => 'P',
        'CCT' => 'P',   'CCR' => 'P',   'CCY' => 'P',   'CCK' => 'P',
        'CCM' => 'P',   'CCW' => 'P',   'CCS' => 'P',   'CCB' => 'P',
        'CCD' => 'P',   'CCH' => 'P',   'CCV' => 'P',   'CCN' => 'P',
        'CGA' => 'R',   'CGC' => 'R',   'CGG' => 'R',   'CGT' => 'R',
        'CGR' => 'R',   'CGY' => 'R',   'CGK' => 'R',   'CGM' => 'R',
        'CGW' => 'R',   'CGS' => 'R',   'CGB' => 'R',   'CGD' => 'R',
        'CGH' => 'R',   'CGV' => 'R',   'CGN' => 'R',   'CTA' => 'L',
        'CTC' => 'L',   'CTG' => 'L',   'CTT' => 'L',   'CTR' => 'L',
        'CTY' => 'L',   'CTK' => 'L',   'CTM' => 'L',   'CTW' => 'L',
        'CTS' => 'L',   'CTB' => 'L',   'CTD' => 'L',   'CTH' => 'L',
        'CTV' => 'L',   'CTN' => 'L',   'GAA' => 'E',   'GAC' => 'D',
        'GAG' => 'E',   'GAT' => 'D',   'GAR' => 'E',   'GAY' => 'D',
        'GCA' => 'A',   'GCC' => 'A',   'GCG' => 'A',   'GCT' => 'A',
        'GCR' => 'A',   'GCY' => 'A',   'GCK' => 'A',   'GCM' => 'A',
        'GCW' => 'A',   'GCS' => 'A',   'GCB' => 'A',   'GCD' => 'A',
        'GCH' => 'A',   'GCV' => 'A',   'GCN' => 'A',   'GGA' => 'G',
        'GGC' => 'G',   'GGG' => 'G',   'GGT' => 'G',   'GGR' => 'G',
        'GGY' => 'G',   'GGK' => 'G',   'GGM' => 'G',   'GGW' => 'G',
        'GGS' => 'G',   'GGB' => 'G',   'GGD' => 'G',   'GGH' => 'G',
        'GGV' => 'G',   'GGN' => 'G',   'GTA' => 'V',   'GTC' => 'V',
        'GTG' => 'V',   'GTT' => 'V',   'GTR' => 'V',   'GTY' => 'V',
        'GTK' => 'V',   'GTM' => 'V',   'GTW' => 'V',   'GTS' => 'V',
        'GTB' => 'V',   'GTD' => 'V',   'GTH' => 'V',   'GTV' => 'V',
        'GTN' => 'V',   'TAA' => '*',   'TAC' => 'Y',   'TAG' => '*',
        'TAT' => 'Y',   'TAR' => '*',   'TAY' => 'Y',   'TCA' => 'S',
        'TCC' => 'S',   'TCG' => 'S',   'TCT' => 'S',   'TCR' => 'S',
        'TCY' => 'S',   'TCK' => 'S',   'TCM' => 'S',   'TCW' => 'S',
        'TCS' => 'S',   'TCB' => 'S',   'TCD' => 'S',   'TCH' => 'S',
        'TCV' => 'S',   'TCN' => 'S',   'TGA' => '*',   'TGC' => 'C',
        'TGG' => 'W',   'TGT' => 'C',   'TGY' => 'C',   'TTA' => 'L',
        'TTC' => 'F',   'TTG' => 'L',   'TTT' => 'F',   'TTR' => 'L',
        'TTY' => 'F',   'TRA' => '*',   'YTA' => 'L',   'YTG' => 'L',
        'YTR' => 'L',   'MGA' => 'R',   'MGG' => 'R',   'MGR' => 'R',
);
my %Rev_Translation;
foreach my $nuc (keys %Translation) {
        $Rev_Translation{$Translation{$nuc}} = $nuc;
}
#Need to check passed data type in the functions to make sure it is a Sequence
#or a compatible datatype.

package Sequence;

#Create a new sequence object
sub new{
    my ($class,$header,$seq,$seq_type)=@_;
    my $self = bless {sequence=>undef,header=>undef,seq_type=>undef}, $class;
    if (defined $header){
        $self->set_header($header);
    }
    
    if (defined $seq){
        $self->set_sequence($seq);
    }
    
    if (defined $seq_type){
        # $seq_type can be "DNA", "RNA", or "PRO"
        $self->set_seq_type($seq_type);
    }
    return $self
}

#Return the lenght of the sequence
sub length{
    my $self=shift;
    return length $self->{sequence};
}

sub complement {
	# Returns the complement of $seq as a new Sequence object
	my $self = shift;
    if (undef $self->{seq_type}) {
        warn "Trying to use comp without a value assigned for \$seq_type\n";
        return $self;
    }
    if ($self->{seq_type} eq "PRO") {
        warn "Cannot find the complement for a protein";
        return $self;
    }
    
    my $tmp_seq = $self->{seq};
	if    ($self->{seq_type} eq "DNA") {$tmp_seq =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/}
	elsif ($self->{seq_type} eq "RNA") {$tmp_seq =~ tr/AUGCYRKMBDHVaugcyrkmbdhv/UACGRYMKVHDBuacgrymkvhdb/}
    my $new_obj = new Sequence($self->{header}, $tmp_seq, $self->{seq_type});
	return $self, $new_obj;
}    
sub rev_comp {
	# Returns the reverse complement of $seq as a new Sequence object
	my $self = shift;
	if (undef $self->{seq_type}) {
        warn "Trying to use rev_comp without a value assigned for \$seq_type\n";
        return $self;
    }
	elsif ($self->{seq_type} eq "PRO") {
        warn "Cannot find the reverse complement for a protein";
        return $self;
    }
    
	my $tmp_seq = reverse($self->{seq});
	if    ($self->{seq_type} eq "DNA") {$tmp_seq =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/}
	elsif ($self->{seq_type} eq "RNA") {$tmp_seq =~ tr/AUGCYRKMBDHVaugcyrkmbdhv/UACGRYMKVHDBuacgrymkvhdb/}
    my $new_obj = new Sequence($self->{header}, $tmp_seq, $self->{seq_type});
	return $self, $new_obj;
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

sub set_seq_type{
    my ($self, $seq_type)=@_;
    if (ref $seq_type eq "SCALAR") {
        if    ($seq_type eq "DNA" or $seq_type eq "D")  {$self->{seq_type} = "DNA"}
        elsif ($seq_type eq "RNA" or $seq_type eq "R")  {$self->{seq_type} = "RNA"}
        elsif ($seq_type eq "PRO" or $seq_type eq "P")  {$self->{seq_type} = "PRO"}
        else  {warn "Invalid \$seq_type provided to set_seq_type\n"}
    }
    else {warn "Trying to use set_seq_type with non-scalar type for \$seq_type\n"}
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

#Extract subsequence given a range and return a new Sequence object
sub extract{
    
}

#Translate sequence
sub translate {
	my ($self) = @_;
	if (undef $self->{seq_type} or $self->{seq_type} ne "DNA" or $self->{seq_type} ne "RNA") {
        	warn "Translate error: sequence type must be either DNA/RNA\n";
		return $self;
	}
	else {
		my $seq = $self->{sequence};
		$seq =~ tr/Uu/Tt/;
		my $trans;
		for (my $i = 0; $i < CORE::length($seq)-3; $i+=3) {
			my $codon = substr($seq, $i, 3);
			if (not exists $Translation{$codon}) {$trans .= "0"}
			else                                 {$trans .= $Translation{$codon}}
		}
		my $new_obj = new Sequence($self->{header}, $trans, "PROTEIN");
		return ($self, $new_obj);
	}
   
}

#Count given k-mer size and return hash {AAA=>10,AAC=>2...}
###############################################
# NON OOP ABBY'S COUNT_KMER DELETE THIS LATER #
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#
sub count_kmer {
        my ($k_length, $seq) = @_;
        my $seq_length = length($seq);
                        
        # Make sure that the sequence length is not shorter than the k-mer length
        if ($seq_length < $k_length) {
                print STDERR "Sequence length ($seq_length) is shorter than k-mer length ($k_length)";
                return;
        } else {
                # Hash to store kmers and associated counts
                my %count;
                        
                # Read through sequence with a window size of $length
                for (my $i = 0; $i <= $seq_length - $k_length; $i++) {
                        my $kmer = substr($seq, $i, $k_length);
                                
                        # Tally counts for kmer   
                        $count{$kmer}++;
                }

                # Return reference to count hash
                return(\%count);
        }
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
# NON OOP ABBY'S COUNT_KMER DELETE THIS LATER #
###############################################

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

sub translate {
        my ($self) = @_;
        if (undef $self->{seq_type} or $self->{seq_type} ne "DNA" or $self->{seq_type} ne "RNA") {
                warn "Translate error: sequence type must be either DNA/RNA\n";
                return $self;
        }
        else {
                my $seq = $self->{sequence};
                $seq =~ tr/Uu/Tt/;
                my $trans;
                for (my $i = 0; $i < CORE::length($seq)-3; $i+=3) {
                        my $codon = substr($seq, $i, 3);
                        if (not exists $Translation{$codon}) {$trans .= "0"}
                        else                                 {$trans .= $Translation{$codon}}
                }
                my $new_obj = new Sequence($self->{header}, $trans, "PRO");
                return ($self, $new_obj); 
        }
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

        my ($target_seq_length, $target_seq_type, $seq_type, $custom_char_weight) = @_;
        die "Error at subroutine rand_seq: length must be positive integer (your input: $target_seq_length).\n" unless $target_seq_length =~ /^\d+E{0,1}\d*$/i;
	
        $target_seq_type = "dna" if not defined($target_seq_type);

        # Define char weight reference based on type #
        my %char_weight;
        # DNA
        if (lc($target_seq_type) eq "dna") {%char_weight = ("A" => 1, "G" => 1, "T" => 1, "C" => 1)}
        # RNA
        elsif (lc($target_seq_type) eq "rna") {%char_weight = ("A" => 1, "G" => 1, "U" => 1, "C" => 1)}
        # Protein
        elsif (lc($target_seq_type) eq "protein") {%char_weight = (
        "A" => 1, "R" => 1, "N" => 1, "D" => 1, "C" => 1,
        "Q" => 1, "E" => 1, "G" => 1, "H" => 1, "I" => 1,
        "L" => 1, "K" => 1, "M" => 1, "F" => 1, "P" => 1,
        "S" => 1, "T" => 1, "W" => 1, "Y" => 1, "V" => 1)}
        # Custom
        elsif (lc($target_seq_type) eq "custom") {
                die "Your type input is \"custom\" but you don't define custom hash\n" unless defined($custom_char_weight);
                %char_weight = %{$custom_char_weight};
        }
        # Else, please go die
        else {die "Please input a valid sequence type [dna|rna|protein|custom]\n"}
        # Randomize #
        # Create a hash to store Stepped Cumulative Probability Distribution (SCPD)#
        my %scpd;

        # Define steps of the probabilities #
        # E.g. if weights for A, C, G, T are 10, 10, 10, 30, we will make a distribution of boundaries
        # starting from 0 + initial weight (10), ending at total sum of weights
        # A = 10 (0+10)
        # C = 20 (10+10)
        # G = 30 (20+10)
        # T = 60 (30+30)
        my $step_boundary = 0;
        foreach my $char (sort {$char_weight{$a} <=> $char_weight{$b}} keys %char_weight) {
                $step_boundary += $char_weight{$char};
                $scpd{$char} = $step_boundary;
        }

        # Create random sequence #
        my $random_seq;
        POSITION: for(my $i = 0; $i < $target_seq_length; $i++){
                my $random_num = rand() * $step_boundary ; # $step_boundary is the total weight

                # Loop over each characer and ask if random number is less than the weight boundary
                foreach my $char (sort {$scpd{$a} <=> $scpd{$b}} keys %scpd) {
                        if ($random_num <= $scpd{$char}) {
                                $random_seq .= $char;
                                next POSITION;
                        }
                }
        }
	my $new_obj = new Sequence("Randomized_seq", $random_seq, $seq_type);
        return($new_obj);
    
}

sub generate_kmer_sequence{
        my ($seq_length, $kmer_table) = @_;
        die "Error at create_rand_seq_kmer: sequence length must be positive integer\n" unless defined($seq_length) and $seq_length =~ /^\d+E{0,1}\d*$/;
        die "Error at create_rand_seq_kmer: please input a hash refernce of kmer table\n" unless defined($kmer_table) and (keys %{$kmer_table} > 0);
        my %kmer_table = %{$kmer_table};

        # Create %weight hash that is used to store total weight of each nucleotide given nothing and
        # total weight of each nucleotide given previous nucleotide based on %kmer_table.
        # The weight is going to be used later to create random sequence using Markov chain dinucleotide.
        my %weight;
        foreach my $kmer (sort keys %kmer_table) {
                my $weight = $kmer_table{$kmer};
                for (my $i = 0; $i < length($kmer); $i++) {
                        my $nuc1 = substr($kmer, $i, 1);
                        my $nuc2 = substr($kmer, $i+1, 1) if $i < length($kmer) - 1;
                        $weight{weight}{$nuc1} = $weight / length($kmer);
                        $weight{$nuc1}{weight}{$nuc2} = $weight / (length($kmer) - 1) if $i < length($kmer) - 1;
                }
        }

        # Now use Stepped Cumulative Distribution Probability (weight)
        # Create stepped boundaries based on total weight of each nucleotide ($nuc2), given another nucleotide ($nuc)
        # There are two types of hash which stores boundary steps for each $nuc:
        # 1. %init_kmer hash for first nucleotide of seq, which doesn't prev nuc on the sequence (only contain $nuc)
        # 2. %kmer hash for the rest of seq, since weight of each nucleotide($nuc2) depend on prev nucleotide ($nuc)
        my %init_kmer;
        my %kmer;
        foreach my $nuc (sort keys %{$weight{weight}}) {
                $init_kmer{total} += $weight{weight}{$nuc};
                $init_kmer{weight}{$nuc} = $init_kmer{total};
                foreach my $nuc2 (sort keys %{$weight{$nuc}{weight}}) {
                        $kmer{$nuc}{total} += $weight{$nuc}{weight}{$nuc2};
                        $kmer{$nuc}{weight}{$nuc2} = $kmer{$nuc}{total};
                }
        }

       # Sequence is stored at $seq
        my $seq;
        POSITION: for (my $i = 0; $i < $seq_length; $i++) {
                my $rand_number = rand();

                # First nuc of sequence, use %init_kmer
                if ($i == 0) {
                        foreach my $nuc (sort {$init_kmer{weight}{$a} <=> $init_kmer{weight}{$b}} keys %{$init_kmer{weight}}) {
                                $seq .= $nuc and next POSITION if $rand_number * $init_kmer{total} < $init_kmer{weight}{$nuc};
                        }
                }

                # Rest of sequence, use %kmer
                else {
                        # Get prev nucleotide
                        my $prev_nuc = substr($seq, $i-1, 1);
                        # Get each nucleotide boundary given previous nucleotide, sorted smallest => biggest 
                        foreach my $nuc (sort {$kmer{$prev_nuc}{weight}{$a} <=> $kmer{$prev_nuc}{weight}{$b}} keys %{$kmer{$prev_nuc}{weight}}) {
                                $seq .= $nuc and next POSITION if $rand_number * $kmer{$prev_nuc}{total} < $kmer{$prev_nuc}{weight}{$nuc};
                        }
                }

        }
	
	my $new_obj = new Sequence("Randomized_seq_kmer", $seq, "DNA");
        return($new_obj);
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

