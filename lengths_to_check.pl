#!/usr/bin/perl -w

use strict;

# Given a FASTA file of reads (output from process_reads.pl) and a set of background values
# (mean number of reads, redundancy, and fraction of pCCMD reads) (output from read_stats.pl),
# determine which read lengths would be productive to mine for motifs.  These are lengths for 
# which at least one of those parameters are significantly higher than background.

# Input is FASTA file of reads (output from process_reads.pl); name of the overall output file
# from the pipeline; output values from read_stats.pl (background read number, redundancy, and
# fraction of pCCMD reads), and the iteration number from the pipeline (for annotation of the 
# output file).

# Outputs are (1) length values and recommendations appended to the overall output file from the
# pipeline, and (2) list of lengths to check printed to STDOUT to be passed back to the shell script.

# Note:  read_stats.pl is run once at the start of the pipeline to determine background values.
# Lengths_to_check.pl may be run multiple times on different sets of reads as motifs are 
# iteratively found and reads eliminated from the set.


### DEFINE VARIABLES

# Input FASTA file and (optionally) background parameters read-number, redundancy, and pCCMD-fraction
my( $read_file, $outfile, $background_readnum, $background_redund, $background_pccmd_frac, $iteration ) = @ARGV;

# Four hashes for data, format:  read length => number of copies.
# For reads matching reference, and pCCMD reads matching reference, store both
# total number of reads and number of non-redundant ("nr") reads.
my( %matching_total, %matching_nr, %pccmd_total, %pccmd_nr );

# Parameters for calculating mean background values
my $low_bound_all = 26;          # smallest read length to check for motifs
my $low_bound_background = 46;   # background calculated for 46 bases and up
my $upper_bound = 80;            # largest read length
my $num_lengths = $low_bound_background - $low_bound_all;
my $MIN_READS = 20;              # minimum reads necessary for motif-finding; see also find_motifs.pl

# Array of read lengths worth checking for motifs
my @lengths_to_check;


### READ DATA FROM FILE

# Read each reference-matching entry in the FASTA file
my( $thislen, $thiscopies );  # length and number of copies of present read
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

    # Ignore very short reads, which are more common in Ion Torrent data
    next if ($thislen < $low_bound_all);
    next if ($thislen >= $low_bound_background);

    # Add read-number values to hashes
    $matching_total{$thislen} += $thiscopies;
    $matching_nr{$thislen}++;

    # Check if pCCMD and if so, add read-number values to hashes
    my $seq = $_;
    my $top_c_pos = $thislen - 17;    # 0-based location of C in a pCCMD read
    my $bot_c_pos = 16;               # 0-based location of G in a pCCMD read
    if (substr( $seq, $top_c_pos, 1 ) eq 'C' && substr( $seq, $bot_c_pos, 1 ) eq 'G') {
        $pccmd_total{$thislen} += $thiscopies;
	$pccmd_nr{$thislen}++;
    }
}
close IN;


### CALCULATE VALUES, APPEND RESULTS TO PIPELINE OUTPUT FILE, AND DETERMINE LENGTHS TO CHECK

open( OUT, ">>$outfile" ) || die "Failed to append to $outfile: $!\n";
print OUT "--------------------------------------\n";
print OUT "ITERATION $iteration\n";
print OUT "--------------------------------------\n";
print OUT "READ ANALYSIS (lengths_to_check.pl)\n";
print OUT "--------------------\n\n";
print OUT "LENGTH\t\tMATCH\t.\t.\t\tpCCMD\t.\t.\t.\t\tCHECK\n";
print OUT "\t\tALL\tNR\tREDUND\t\tALL\t(FRAC)\tNR\t(FRAC)\n";
foreach my $len ($low_bound_all..$low_bound_background-1) {

    # Populate missing values with 0's
    unless (exists $matching_total{$len}) {
	$matching_total{$len} = 0;
	$matching_nr{$len} = 0;
    }
    unless (exists $pccmd_total{$len}) {
	$pccmd_total{$len} = 0;
	$pccmd_nr{$len} = 0;
    }

    # Don't bother printing empty rows after the first iteration
    next if ($iteration > 1 && $matching_total{$len} == 0);

    # Calculate redundancy of different types of reads
    my $redundancy = 0;          # total reads / nr reads
    my $frac_pccmd_all = 0;      # pccmd reads / total reads
    my $frac_pccmd_nr = 0;       # nr pccmd reads / nr reads
    if ($matching_nr{$len} > 0) {
	$redundancy = int($matching_total{$len}/$matching_nr{$len} + 0.5);  # round to nearest positive integer
	$frac_pccmd_all = $pccmd_total{$len}/$matching_total{$len};
	$frac_pccmd_nr  = $pccmd_nr{$len}/$matching_nr{$len};    
    }

    # Print length, numbers of reads, and redundancy for this length
    print OUT "$len\t\t$matching_total{$len}\t$matching_nr{$len}\t$redundancy\t\t$pccmd_total{$len}\t";
    printf OUT ( "%.3f\t", $frac_pccmd_all );
    print OUT "$pccmd_nr{$len}\t";
    printf OUT ( "%.3f\t\t", $frac_pccmd_nr );

    if ($len == 33) {  # this is not a valid CCMD length
	print OUT "\n";
	next;
    }

    # Determine whether this length is worth checking for motifs, i.e., if any of the parameters are 2x over
    # background value.  (This tends to be very permissive.)
    # Also require a minimum number of reads to bother checking the length.
    if ( ($matching_total{$len} >= 2*($background_readnum) || $redundancy >= 2*($background_redund)
	 || $frac_pccmd_all >= 2*($background_pccmd_frac)) && $pccmd_nr{$len} >= $MIN_READS ) {
	push( @lengths_to_check, $len );
	print OUT "* $len\n";
    } else {
	print OUT "\n";
    }
}

if (@lengths_to_check) {
    print OUT "\n";
} else {
    print OUT "\nNO LENGTHS TO CHECK.\n\n";
}
close OUT;


### PRINT RECOMMENDATIONS TO STDOUT

print "LENGTHS TO CHECK=";
if (@lengths_to_check) {
    foreach my $len (@lengths_to_check) {
	print "$len ";
    }
} else {
    print "0\n";
}
