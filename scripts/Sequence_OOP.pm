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

package Sequence;
use overload '.' => \&concatenate_copy;
use overload '.=' => \&concatenate;
use overload '=' => \&copy;
use overload '""' => \&sequence;

=head2 Sequence Constructor

Sequence::new constructor subroutine

Examples:

    my $seq = new Sequence(); #empty sequence
    
    my $seq = new Sequence("Bare Sequence");

    my $seq = new Sequence( -header   => "Fasta Header or Sequence Name,
                            -sequence => "Sequence ACTDDA");
                            
    my $seq = new Sequence( -header   => "Fasta Header or Sequence Name,
                            -sequence => "Sequence ACTDDA",
                            -type     => "PROT");


Accepted types are "UNDEFINED"  , "UNDEF" (DEFAULT)
                   "UNKNOWN"    , "UNK"
                   "RNA"        , "R"
                   "DNA"        , "D"
                   "PROT"       , "P" , "AA"
                   
(Types initially created by Matt Porter.  Paul Lott added UNDEFINED and UNKNOWN)

UNDEFINED = sequence that hasn't been check and type hasn't been defined

UNKNOWN = Has been checked and the sequence isn't RNA, DNA, or PROT  


=cut


sub new{
    my %default=("TYPE"     =>  "UNDEFINED",
                 "HEADER"   =>  q(),
                 "SEQUENCE" =>  q());
    my %arg;
    my $class = shift;
    
    if (scalar @_ % 2 == 0){
        %arg=(%default,@_);
    }
    else{        
        my $seq = shift;
        %arg = (%default,@_);
        $arg{SEQUENCE} = $seq;
     }
    my $self = bless {sequence  =>  $arg{SEQUENCE},
                      header    =>  $arg{HEADER},
                      seq_type  =>  $arg{TYPE}} , $class;
    #TODO:  once the set_type function is fixed.  I'll set the type using set_type();
    
    return $self;
}


=head2 length

Return the length of the sequence

Example:
my $length = $sequence->length;

=cut

sub length{
    my $self=shift;
    return CORE::length $self->{sequence};
}


=head2 concatenate

concatenate allows user to concatenate two Sequences together or a string to the
Sequence instance.  Provides functionality for '.=' operator

Example:  
$seq->concatenate("ACGTA"); #concatenate a string to end of Sequence instance
$seq->concatenate($alt_seq);  #concatenate another Sequence instance

Example through overloaded operator '.=':

$seq .= "ACGTA";
$seq .= $alt_seq;


=cut

sub concatenate{
    my ($self, $r_self)= @_;
    
    if (ref $self eq "Sequence" && ref $r_self eq "Sequence"){
        $self->{sequence} .=  $r_self->{sequence};
    }
    elsif (ref $self eq "Sequence" && !(ref $r_self)){
        $self->{sequence} .= $r_self;
    }
    
    return $self
}

=head2 concatenate_copy

concatenate_copy subroutine allows user to concatenate two Sequences, or
Sequence and string together.   Returns a new Sequence instance.  Provides
functionality for '.' perl concatenate operator.

Example:
my $new_seq = $seq->concatenate_copy($alt_seq);
my $new_seq = $seq->concatenate_copy("ACGT");

Example through overloaded operator '.':

my $new_seq = $seq . "ACGT";    #Concatenate a Sequence and string
my $new_seq = "ACGT" . $seq;    #Concatenate a Sequence and string
my $new_seq = $seq . $alt_seq;  #Concatenate two Sequence instances

=cut

sub concatenate_copy{
    my ($self, $concat_self, $order) = @_;
    
    my $new_seq = new Sequence();
    
    
     if ($order){
        if (ref $concat_self eq "Sequence"){
            $new_seq->{sequence} = $concat_self->{sequence} . $self->{sequence};
        }
        else{
            $new_seq->{sequence} = $concat_self . $self->{sequence};
        }
    }
    else{
         if (ref $concat_self eq "Sequence"){
            $new_seq->{sequence} = $self->{sequence} . $concat_self->{sequence} ;
        }
        else{
            $new_seq->{sequence} = $self->{sequence} . $concat_self;
        }
    }
    
    return $new_seq;
}

=head2 copy

Copy a Sequence instance and return the new Sequence instance

Example:
my $new_seq = $seq->copy();

Example using overloaded operator '=':
my $new_seq = $seq;

=cut

sub copy{
    my $self = shift;
    my $new_seq = new Sequence( SEQUENCE=>$self->sequence, HEADER=>$self->{header}, TYPE=>$self->{seq_type});
    return $new_seq;
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


=head2 sequence accessor

Accessor for sequence

If used without passing parameters, the function will return the sequence string.
If passed a string, it will assign the string to the Sequence;

Example:

    print $seq->sequence . "\n";   #passes sequence to print

    $seq->sequence("ACGTACGT")  # assigns sequence the passed value;

=cut

sub sequence{
    my ($self,$sequence)=@_;
    
    if (!defined $sequence){
        return $self->{sequence};
    }
    
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

=head2 sequence_ref

Return a reference to sequence string
For long sequenes passing a reference to function is much faster than passing
by copy.   This allows quick access to a reference to sequence

Example:

example_function($seq->sequence_ref);  #Passing a reference to sequence to function

=cut

sub sequence_ref{
    my $self=shift;
    return \$self->{sequence};
}

=head2 header

Access or Returns the header of the Sequence object
If parameter provided, header will be changed to provided parameter.
If no parameters provided, function will return the header of Sequence object.

Example:

$seq->header("My sequence definition"); #Set the header
my $string = $seq->header;  #Assign $string the header string from Sequence object

=cut

sub header{
    my ($self,$header)=@_;
    
    if (!defined $header){
        return $self->{header};
    }
    
    if (ref $header eq "SCALAR"){
        $self->{header}=$$header;
    }
    elsif (ref $header){
        warn "Trying to set header with non-scalar type\n";
    }
    else{
        $self->{header}=$header;
    }
    
    return $self;
}

=head2 as_fasta_string

Return sequence as a Fasta formatted string

=cut

sub as_fasta_string{
    my $self=shift;
    return ">" . $self->{header} . "\n" . $self->{sequence} . "\n";
}

=head2 print_fasta

Print the Sequence object as Fasta
If passed a GLOB, it will print the Fasta formatted Sequence to the GLOB

Example:
    $seq->print_fasta;
    $seq->print_fasta(\*STDOUT);

=cut
sub print_fasta{
    my ($self,$fh)=@_;
    if (defined $fh && ref $fh eq "GLOB"){
        print $fh $self->as_fasta_string;
    }
    else{
        print $self->as_fasta_string;
    }
    
    return;
}

#Clean the sequence of invalid and whitespace characters
sub clean{
 	#my ($sequence) = @_;
    my $self=shift;
 	my $sequence = $self->sequence;
	$sequence =~ s/\s//g;
	$sequence = uc($sequence);
	if ($sequence =~ m/[^ACGTURYSWKMBDHVN\*\-EFILPQ]/) {
		die "Nonstandard characters found in this sequence";
	}
	$self->set_sequence($sequence);
	return;   
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
	if (not defined($self->{seq_type}) and $self->{seq_type} ne "DNA" and $self->{seq_type} ne "RNA") {
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
	my $new_obj = new Sequence($random_seq, "random_seq", $seq_type);
	my $new_obj = new Sequence(	SEQUENCE 	=> $random_seq, 
					HEADER 		=> "random_seq",
					TYPE		=> $seq_type);
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
	my $new_obj = new Sequence(	SEQUENCE 	=> $seq, 
					HEADER 		=> "random_kmer_seq",
					TYPE		=> "DNA");
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

