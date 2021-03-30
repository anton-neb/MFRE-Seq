#!/usr/bin/perl -w

use strict;

# For all of the motifs found, determine statistics.

# Inputs are

# Given a FASTA file of reads (output from process_reads.pl) and a set of motifs, search
# each of those lengths separately for conserved motifs.  This is done by examining an ungapped
# alignment of reads column by column and comparing the nucleotide distribution at that column with
# expected background distributions of each IUPAC base symbol using KL-divergence.

# Input is FASTA file of reads (output from process_reads.pl); name of the overall output file
# from the pipeline; a string of motif sequences and associated read lengths; and a minimum 
# redundancy (output from read_stats.pl).

# The vast majority of reads generated are of the type (16,16) (exact) and (16,17) (extra base on
# one end).  As output, print statistics on the length, number, and redundancy of (16,16) and
# (16,17) reads for each motif found.  Print both to pipeline output file and to STDOUT.


### DEFINE VARIABLES

# Input FASTA file, output filename, space-separated list of motifs, and a minimum copy number.
my( $read_file, $outfile, $motif_string, $min_copyno ) = @ARGV;

# Parse motif string into array, each element of the form "motif(read_length,top_c_pos,bottom_c_pos)"
my @motif_list = split " ", $motif_string;

# Further parse motif data into hashes
my %readlength_motifs;      # motif|read-length => motif
my %readlength_top_c_positions;  # motif|read-length => 0-based position of top-strand m5C
my %readlength_bot_c_positions;  # motif|read-length => 0-based position of bottom-strand m5C (i.e., top-strand G)
my %lengths_list;   # form read_length => 0 for all read lengths and read lengths +1
foreach( @motif_list ) {
    /(.+)\((\d+),(\d+),(\d+)\)/;
    $readlength_motifs{"$1|$2"} = $1;
    $readlength_top_c_positions{"$1|$2"} = $3;
    $readlength_bot_c_positions{"$1|$2"} = $4;
    $lengths_list{$2} = 0;
    $lengths_list{$2 + 1} = 0;
}

# Parameters for deriving motifs
$min_copyno = 2 if ($min_copyno < 2);


### READ DATA FROM FILE

# Read each reference-matching entry in the FASTA file
my( $thislen, $thiscopies );         # length and number of copies of present read
my( %read_lengths, %copy_numbers );  # form read_sequence => length or copy number
open( IN, $read_file ) || die "Failed to open $read_file: $!\n";
while (<IN>) {
    chomp;

    # Get read length and number of copies from ID line
    if ( /^>/ ) {
	last unless ( /reference.*length (\d+)\|copies (\d+)$/ );
	$thislen = $1;
	$thiscopies  = $2;
	next;
    }

    # Keep only reads of lengths relevant to the motif set and of sufficient copy number
    next unless ( exists $lengths_list{$thislen} );
    next unless ( $thiscopies >= $min_copyno );

    $read_lengths{$_} = $thislen;
    $copy_numbers{$_} = $thiscopies;
}
close IN;


### GET READ STATISTICS FOR EACH MOTIF AND PRINT TO OUTPUT

open( OUT, ">>$outfile" ) || die "Failed to append to $outfile: $!\n";
print OUT "--------------------------------------\n";
print OUT "FINAL SET OF MOTIFS IDENTIFIED (final_motifs.pl)\n";
print OUT "--------------------------------------\n\n";
printf OUT ( "%20s %15s %15s %15s %15s %15s %15s\n", "MOTIF", "LEN_(16,16)", "READS_(16,16)",
	     "REDUND_(16,16)", "LEN_(16,17)", "READS_(16,17)", "REDUND_(16,17)" );

# Also print to STDOUT
printf( "%20s %15s %15s %15s %15s %15s %15s\n", "MOTIF", "LEN_(16,16)", "READS_(16,16)",
	"REDUND_(16,16)", "LEN_(16,17)", "READS_(16,17)", "REDUND_(16,17)" );

# Call (16,16) reads and (16,17) reads 00 and 01, respectively, for brevity
foreach( keys %readlength_motifs ) {
    my( $motif, $length00 ) = split /\|/, $_;
    my $length01 = $length00 + 1;
    my $num_00_reads = 0;
    my $mean_00_redundancy = 0;
    my $num_01_reads = 0;
    my $mean_01_redundancy = 0;

    # Pad the motif with 5' N's to align it with a properly cut read (where bottom m5C = 16)
    my $num_Ns = 16 - $readlength_bot_c_positions{$_};
    my $motif_aligned = ("N" x $num_Ns) . $motif;
    my $motif_aligned_regex = makeregex( $motif_aligned );

    # Test each read for the presence of the motif the proper distance from the ends
    foreach my $read (keys %read_lengths ) {

	my $thislen = $read_lengths{$read};
	my $read_rc = rev_comp( $read );  # reverse complement of read sequence

	if ($thislen == $length00) {
	    # check if a (16,16) read, only need to check one side
	    if ($read =~ /^$motif_aligned_regex/) {
		$num_00_reads++;
		$mean_00_redundancy += $copy_numbers{$read};
	    }
	} elsif ($thislen == $length01) {
	    # check if a (16,17) read, need to check both sides
	    if (($read =~ /^$motif_aligned_regex/) || ($read_rc =~ /^$motif_aligned_regex/)) {
		$num_01_reads++;
                $mean_01_redundancy += $copy_numbers{$read};
	    }
	}
    }

    # Convert absolute numbers to mean values
    if ($num_00_reads > 0) {
	$mean_00_redundancy = $mean_00_redundancy/$num_00_reads;
    }
    if ($num_01_reads > 0) {
	$mean_01_redundancy = $mean_01_redundancy/$num_01_reads;
    }

    # Print statistics on the (0,0) and (0,1) reads generated by this motif
    printf OUT ( "%20s %15d %15d %15.2f %15d %15d %15.2f\n", $motif, $length00, $num_00_reads,
            $mean_00_redundancy, $length01, $num_01_reads, $mean_01_redundancy );
    printf( "%20s %15d %15d %15.2f %15d %15d %15.2f\n", $motif, $length00, $num_00_reads,
	    $mean_00_redundancy, $length01, $num_01_reads, $mean_01_redundancy );
}
close OUT;
	    


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
