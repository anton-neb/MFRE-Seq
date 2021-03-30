#!/usr/bin/perl -w

use strict;

# Given a FASTA file of reads (output from process_reads.pl) and a set of lengths to check, search
# each of those lengths separately for conserved motifs.  This is done by examining an ungapped
# alignment of reads column by column and comparing the nucleotide distribution at that column with
# expected background distributions of each IUPAC base symbol using KL-divergence.

# Input is FASTA file of reads (output from process_reads.pl); name of the overall output file
# from the pipeline; the file name of a temporary FASTA file in which all reads derived from motifs
# found here have been removed (designated by the pipeline as "temp_reads.fasta"); a list of lengths
# to check (output from lengths_to_check.pl); a minimum redundancy (output from read_stats.pl);
# and the iteration number from the pipeline (for annotation of the output file).

# General procedure:
#   1. Get reads from FASTA input file and keep only those of designated lengths that are also
#         pCCMD reads and of sufficient copy number.
#   2. For each designated length separately, attempt to derive a conserved motif from the reads.
#   3. Score motifs based on specificity and identify the motif with the highest score.
#   4. Remove all reads created by cutting at that motif and pass the reduced set as a FASTA
#         file back to the pipeline.  All other motifs found in step 2, if they are real, should
#         be found again in this reduced set in the next iteration of the pipeline.


### DEFINE VARIABLES

# Input FASTA file, output filenames, a space-separated list of lengths to examine, and a minimum
# copy number from which to derive lengths.
my( $read_file, $outfile, $fasta_output, $length_string, $min_copyno, $iteration ) = @ARGV;

# Parameters for deriving motifs
$min_copyno = 2 if ($min_copyno == 1);  # ensure a minimum redundancy of 2
my $MIN_READS = 20;                     # and minimum number of unique reads

# Parse the length string into a hash
my @lengths = split ' ', $length_string;
my %length_hash;
foreach( @lengths ) {
    $length_hash{$_} = 0;
}

# Hashes for read sequence data
my %nonred_sequences;  # all non-redundant sequences from the input file; sequence => length
my %seqs_for_output;   # for output FASTA file; sequence => FASTA_header
my %nt_fraction;       # stores fractions of A, C, G, and T in all reads examined


### READ DATA FROM FILE

# Read each reference-matching entry in the FASTA file
my( $seqs_all, $seqs_rightlength, $seqs_mincopy, $seqs_used );  # keep track of numbers of reads
my( $thislen, $thiscopies, $header );      # length, copy number, and header of present read
open( IN, $read_file ) || die "Failed to open $read_file: $!\n";
while (<IN>) {
    chomp;

    # Get read length, number of copies, and complete line from ID line
    if ( /^>/ ) {
	last unless ( /reference.*length (\d+)\|copies (\d+)$/ );
	$thislen = $1;
        $thiscopies  = $2;
	$header = $_;
	next;
    }

    # Ignore all reads not of one of the designated lengths, and all reads of low copy number
    $seqs_all++;            # all sequences examined
    next unless (exists ($length_hash{ $thislen }));   # keep only pre-defined lengths
    $seqs_rightlength++;    # all sequences that are of the lengths in question
    next unless ($thiscopies >= $min_copyno);          # keep only reads with sufficient copy number
    $seqs_mincopy++;        # all sequences of right length and of sufficient copy number

    # Ignore all reads that are not pCCMD
    my $seq = $_;
    my $top_c_pos = $thislen - 17;    # 0-based location of C in a pCCMD read
    my $bot_c_pos = 16;               # 0-based location of G in a pCCMD read
    next unless (substr( $seq, $top_c_pos, 1 ) eq 'C' && substr( $seq, $bot_c_pos, 1 ) eq 'G');
    $seqs_used++;           # all sequences of right length and copy number and also pCCMD-reads

    # Store data for all reads passing the above tests
    $nonred_sequences{ $seq } = $thislen;         # store full sequence and its length
    $seqs_for_output{ $seq } = $header;           # store full sequence and header for output file
    
    # For all reads actually used to derive motifs, add number of each nt to running totals
    for my $nt qw(A C G T) {
	$nt_fraction{$nt} += () = $seq =~ /$nt/gi;
    }
}
close IN;

# Convert total nt numbers to fractions, to be used in KL-divergence for background
my $totalnt;
for my $nt qw(A C G T) {
    $totalnt += $nt_fraction{$nt};
}
for my $nt qw(A C G T) {
    $nt_fraction{$nt} = $nt_fraction{$nt}/$totalnt;
}


### APPEND LENGTH INFO TO PIPELINE OUTPUT FILE

open( OUT, ">>$outfile" ) || die "Failed to append to $outfile: $!\n";
print OUT "--------------------\n";
print OUT "READ STATISTICS (find_motifs.pl)\n";
print OUT "--------------------\n";
print OUT "Total no. of non-redundant sequences:            $seqs_all\n";
print OUT "No. of these also with lengths in question:      $seqs_rightlength\n";
print OUT "No. of these also with copy no. > $min_copyno:   $seqs_mincopy\n";
print OUT "No. of these that are also pCCMD-reads:          $seqs_used\n";
print OUT "Nucleotide composition of reads used:\n";
for my $nt qw(A C G T) {
    printf OUT ( "%10s   %5.3f\n", $nt, $nt_fraction{$nt} );
}
print OUT "\n";
print OUT "--------------------\n";
print OUT "BASE COMPOSITION FOR EACH LENGTH (find_motifs.pl)\n";
print OUT "--------------------\n";


### DERIVE MOTIFS FOR EACH DESIGNATED LENGTH AND APPEND TO OUTPUT FILE

# Hashes for motif data
my %length_to_motif;        # list of motifs found, format:  read-length => motif
my %motif_set_with_flanks;  # motif|read-length => motif_w_flank; (keeps preceding N's for anchoring)
my %motif_set_numreads;     # motif|read-length => number of reads used to determine it
my %meth_bases;             # motif|read-length => C-location, G-location (0-based)

# For each designated read length, try to derive a motif
foreach my $length_to_check (sort {$a <=> $b} keys %length_hash) {

    print OUT "\n# LENGTH = $length_to_check #\n\n";
    my $top_c_pos = $length_to_check - 17;    # 0-based location of C
    my $bot_c_pos = 16;                       # 0-based location of G

    # Reads of this length that have passed all of the above tests
    my %reads_to_check;    # set of filtered reads of this length
    foreach my $seq (keys %nonred_sequences) {
	next unless ($nonred_sequences{$seq} == $length_to_check);
       	$reads_to_check{ $seq } = 0;
    }

    # Check there are enough reads of this length to work with.
    my $number_of_reads = () = keys %reads_to_check;
    if ($number_of_reads < $MIN_READS) {
	print OUT "\tToo few unique reads of length $length_to_check.  No motif derived.\n\n";
	next;
    }

    # If so, try to derive a motif, to be stored in these variables
    my $final_motif;   # motif for $length_to_check
    my @symbols;       # string of IUPAC bases making up final motif

    # For each position along the read length, look at nt distribution of the alignment column
    printf OUT ( "%8s %12s %8s %8s %8s %8s\n\n", "Column", "Prediction", "A", "C", "G", "T" );
    for my $i (0..$length_to_check-1) {

        # For each column of the alignment, count the number of each base
	my %column_count;   # no. reads with A, C, G, T in this column; base => no. reads
	foreach my $seq (keys %reads_to_check) {
	    $column_count{ substr($seq, $i, 1) }++;
	}

        # Convert hash of absolute numbers to array of base frequencies
	my @column_freq;    # base frequencies, in order A, C, G, T
	for my $base qw(A C G T) {
	    $column_count{$base} = 0 unless ($column_count{$base});
	    push( @column_freq, $column_count{$base}/$number_of_reads );
	}

	# Using KL-divergence, compare this array to that of each IUPAC symbol
	# and choose the best match.
	my $kl_prediction = kl( \@column_freq, \%nt_fraction );

	# The invariant, methylated bases are flagged by brackets in output
	my $predicted_base;
	if ($i == $top_c_pos) {
	    $predicted_base = "[C]";  # invariant
	} elsif ($i == $bot_c_pos) {
	    $predicted_base = "[G]";  # invariant
	} else {
	    $predicted_base = " " . $kl_prediction . " ";
	}

	# Print info for this position to output file
	printf OUT ( "%8d %12s %8.3f %8.3f %8.3f %8.3f\n", $i, $predicted_base, 
		     @column_freq );

	# Add the predicted base at this position to the final motif
	push( @symbols, $predicted_base );
    }

    # Remove trailing N's from both ends of motif and string the list together
    while( $symbols[-1] eq ' N ' ) {
	pop @symbols;
    }
    my $motif_with_flank = join "", @symbols;
    while( $symbols[0] eq ' N ' ) {
	shift @symbols;
	$top_c_pos--;
	$bot_c_pos--;
    }
    $final_motif = join "", @symbols;

    # Get rid of brackets
    $final_motif =~ s/\W//g;
    $motif_with_flank =~ s/\W//g;
    print OUT "\n      Predicted motif for length $length_to_check = $final_motif\n\n";

    # Add the motif to the list of all motifs.
    # Sometimes same motif can appear at multiple lengths, so keep track of read numbers.
    $length_to_motif{$length_to_check} = $final_motif;
    $motif_set_with_flanks{"$final_motif|$length_to_check"} = $motif_with_flank;
    $motif_set_numreads{"$final_motif|$length_to_check"} = $number_of_reads;
    $meth_bases{"$final_motif|$length_to_check"} = "$top_c_pos,$bot_c_pos";
}


### ASSESS MOTIF QUALITY (SPECIFICITY) AND PRINT SUMMARY TO OUTPUT

print OUT "--------------------\n";
print OUT "SUMMARY OF MOTIFS, ITERATION $iteration (find_motifs.pl)\n";
print OUT "--------------------\n\n";

my $total_reads_start = () = keys %nonred_sequences;
print "From $total_reads_start unique reads, MOTIFS FOUND : ";

my %motif_scores;   # motif sequence => per-base specificity score
my $most_specific_motif = '';
my $most_specific_motif_with_flank = '';  # save preceding N's
my $most_specific_score = 0;
my $most_specific_read_length = 0;
my $most_specific_read_number = 0;
foreach( sort keys %length_to_motif ) {
   
    my $this_len = $_;
    my $this_motif = $length_to_motif{$this_len};

    # See what's left of the motif after removing the methylated bases and all N's
    my $for_scoring = $this_motif;   # motif after removing 1 C, 1 G and all N's
    $for_scoring =~ s/C//;
    $for_scoring =~ s/G//;
    $for_scoring =~ s/N//g;
    $for_scoring =~ s/\W//g;         # also get rid of brackets

    # Ignore motif if low specificity
    if (sp_score($for_scoring) < 10) {
	delete ($length_to_motif{$this_len});
	next;
    }

    # Per-base specificity = raw specificity score / length of the whole motif
    $motif_scores{$this_motif} = sp_score($for_scoring)/length( $this_motif );

    # Print motif to output file and to STDOUT
    print OUT "Length $this_len : $this_motif\tSCORE = $motif_scores{$this_motif}\n";
    print "$this_motif($this_len) ";

    # Keep track of most specific motif
    if ($motif_scores{$this_motif} >= $most_specific_score) {

	# This motif may already have been identified, but at a different length.  If so, choose
	#   the length associated with the highest number of reads.
	if (($most_specific_motif eq $this_motif) &&
	    ($motif_set_numreads{"$this_motif|$this_len"} < $most_specific_read_number)) {
	    next;
	}

	$most_specific_score = $motif_scores{$this_motif};
	$most_specific_motif = $this_motif;
	$most_specific_read_length = $this_len;
	$most_specific_read_number = $motif_set_numreads{"$this_motif|$this_len"};
	$most_specific_motif_with_flank = 
	    $motif_set_with_flanks{"$most_specific_motif|$most_specific_read_length"};
    }
}

# Case where no motifs are found at any length
unless (%length_to_motif) {
    print OUT "None found.\n\n\n";
    close OUT;

    # Print to STDOUT to pass back to pipeline
    print "0 ; MOST_SPECIFIC_MOTIF=0\n";
    exit;
}


### OF THOSE FOUND, PRINT DATA ON THE MOTIF MOST LIKELY TO BE "REAL"

print OUT "\nMOST SPECIFIC MOTIF = $most_specific_motif (read length ";
print OUT "$most_specific_read_length)\n\n\n";
close OUT;

# Print to STDOUT, in the format "motif(read_length,top_C_pos,bottom_C_pos)"
my $meth_base_pos = $meth_bases{"$most_specific_motif|$most_specific_read_length"};
print "|| MOST_SPECIFIC_MOTIF=$most_specific_motif($most_specific_read_length,$meth_base_pos)\n";


### PREPARE TEMPORARY FASTA FILE FOR THE NEXT ITERATION

# Output FASTA file will contain all reads from the input file EXCEPT those that contain
# the most specific motif at the proper cut distance from at least one end.

# Asymmetric sites may not be centered on the read, so determine cut distance by saving
# the flanking N's before the motif, which "anchors" it to the end of the read.  Check
# cut distance on the other end by anchoring to the start of the reverse complement.
my $motif_with_flank_regex = makeregex( $most_specific_motif_with_flank );

open( FASTA, ">$fasta_output" ) || die "Failed to create $fasta_output: $!\n";

# Test each read for the presence of the motif with >=1 correct cut
foreach my $read ( keys %seqs_for_output ) {
    my $read_rc = rev_comp( $read );  # reverse complement of read sequence

    # Check for presence of motif at exact distance from left end
    next if ($read =~ /^$motif_with_flank_regex/);

    # Check for presence of motif at exact distance from right end
    next if ($read_rc =~ /^$motif_with_flank_regex/);

    # If motif not present at either place, keep for next iteration
    print FASTA "$seqs_for_output{$read}\n$read\n";
}
close FASTA;





## Subroutine to return the minimum of two numbers
##
sub min {
    my( $num1, $num2 ) = @_;
    if ($num1 <= $num2) {
	return $num1;
    } else {
	return $num2;
    }
}



## Subroutine to convert RE sites to regex format
##
sub makeregex {
        my %code = (     # nucleotide and degenerate code
             'A' => 'A',               'T' => 'T',
             'G' => 'G',               'C' => 'C',
             'R' => '[AG]',            'Y' => '[CT]',
             'M' => '[AC]',            'K' => '[GT]',
             'S' => '[CG]',            'W' => '[AT]',
             'H' => '[ACT]',           'B' => '[CGT]',
             'V' => '[ACG]',           'D' => '[AGT]',
             'N' => '[ACGT]'
	    );

        my @site = split //, $_[0] ;   # chars of RE site from file
        my $regex;                     # RE site written as regex
        foreach my $letter (@site) {
	    $regex .= $code{$letter};  # build letter by letter
        }
        return $regex;
}


## Subroutine to perform KL-divergence vs. IUPAC symbols
##
sub kl {
    my @col = @{$_[0]};
    my %f = %{$_[1]};
    my $min_kl = 10000;
    my $min_symb = '-';

    # Frequencies of [A, C, G, T] respectively, associated with each symbol
    my %frequencies = (
	'A' => [ 1, 0, 0, 0 ],
	'C' => [ 0, 1, 0, 0 ],
	'G' => [ 0, 0, 1, 0 ],
	'T' => [ 0, 0, 0, 1 ],
	'R' => [ $f{'A'}/($f{'A'}+$f{'G'}), 0, $f{'G'}/($f{'A'}+$f{'G'}), 0 ],
	'Y' => [ 0, $f{'C'}/($f{'C'}+$f{'T'}), 0, $f{'T'}/($f{'C'}+$f{'T'}) ],
	'M' => [ $f{'A'}/($f{'A'}+$f{'C'}), $f{'C'}/($f{'A'}+$f{'C'}), 0, 0 ],
	'K' => [ 0, 0, $f{'G'}/($f{'G'}+$f{'T'}), $f{'T'}/($f{'G'}+$f{'T'}) ],
	'S' => [ 0, $f{'C'}/($f{'C'}+$f{'G'}), $f{'G'}/($f{'C'}+$f{'G'}), 0 ],
	'W' => [ $f{'A'}/($f{'A'}+$f{'T'}), 0, 0, $f{'T'}/($f{'A'}+$f{'T'}) ],
	'H' => [ $f{'A'}/($f{'A'}+$f{'C'}+$f{'T'}), $f{'C'}/($f{'A'}+$f{'C'}+$f{'T'}), 
		 0, $f{'T'}/($f{'A'}+$f{'C'}+$f{'T'}) ],
	'B' => [ 0, $f{'C'}/($f{'C'}+$f{'G'}+$f{'T'}), 
		 $f{'G'}/($f{'C'}+$f{'G'}+$f{'T'}), $f{'T'}/($f{'C'}+$f{'G'}+$f{'T'}) ],
	'V' => [ $f{'A'}/($f{'A'}+$f{'C'}+$f{'G'}), $f{'C'}/($f{'A'}+$f{'C'}+$f{'G'}), 
		 $f{'G'}/($f{'A'}+$f{'C'}+$f{'G'}), 0 ],
	'D' => [ $f{'A'}/($f{'A'}+$f{'G'}+$f{'T'}), 0, 
		 $f{'G'}/($f{'A'}+$f{'G'}+$f{'T'}), $f{'T'}/($f{'A'}+$f{'G'}+$f{'T'}) ],
	'N' => [ $f{'A'}, $f{'C'}, $f{'G'}, $f{'T'} ]
	);

    # Determine KL( column, iupac_freq ) for each IUPAC symbol
    foreach my $iupac (keys %frequencies) {
	my $kl_div;  # KL-divergence for this pair
	my @iupac_freq = @{ $frequencies{$iupac} };
	for my $b (0..3) {
	    $iupac_freq[$b] = 0.00001 if ($iupac_freq[$b]  == 0);
	    $col[$b] = 0.00001 if ($col[$b] == 0);
	        $kl_div += $col[$b] * log($col[$b]/$iupac_freq[$b]) + 
		    $iupac_freq[$b] * log($iupac_freq[$b]/$col[$b]);
	}
	$kl_div = $kl_div/2;
	if ($kl_div < $min_kl) {
	    $min_kl = $kl_div;
	    $min_symb = $iupac;
	}
    }
    return $min_symb;
}



## Subroutine to calculate specificity score of a motif, after it has already
## been stripped of all N's and the anchoring (methylated) C and G
##
sub sp_score {
    my $sequence = shift;  # Already stripped of N's and anchoring C and G
    return 0 unless ($sequence);

    my @chars = split '', $sequence;

    # Calculate specificity score
    my $specificity = 0;
    foreach my $letter (@chars) {
        $specificity += 8 if ($letter =~ /[ACGT]/);
        $specificity += 4 if ($letter =~ /[RYMKSW]/);
        $specificity += 2 if ($letter =~ /[HBVD]/);
        $specificity += 1 if ($letter =~ /N/);
    }

    return $specificity;
}


## Subroutine to return the reverse complement of a DNA sequence
##
sub rev_comp {
# Return reverse complement of a DNA sequence
    my $seq = shift;
    my $rcseq = reverse $seq;
    $rcseq = lc( $rcseq );
    $rcseq =~ s/g/C/g;
    $rcseq =~ s/c/g/g;
    $rcseq =~ s/C/c/g;
    $rcseq =~ s/a/T/g;
    $rcseq =~ s/t/a/g;
    $rcseq =~ s/T/t/g;
    $rcseq = uc( $rcseq );
    return $rcseq;
}
