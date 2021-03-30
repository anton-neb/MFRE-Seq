#!/usr/bin/perl -w

use strict;

# From a FASTA file of reads (output from process_reads.pl), determine the "background" values:
# mean number of reads, redundancy, and fraction of putative CCMD ("pCCMD") reads in the set.
# pCCMD reads have a C as 17th base from 3' end and G as 17th base from 5' end.

# Sole input is the FASTA file above, which must be in the format specified by process_reads.pl.
# Only the reference-matching reads from this file are used for calculation.

# Outputs are (1) file containing the background values calculated from the input file, and
# (2) background values printed to STDOUT to be passed to the next script in the pipeline. 



### DEFINE VARIABLES

# Input FASTA filename
my( $read_file, $outfile ) = @ARGV;
unless (defined $read_file && defined $outfile) {
    die "(read_stats.pl) : Need to supply both an input FASTA file and an output filename.\n";
}
# Test formatting of input file
open( IN, $read_file ) || die "Failed to open $read_file: $!\n";
my $first_line = <IN>;
close IN;
unless ($first_line =~ /\>(.*)reference(.*)copies/) {
    die "FASTA IDs not properly formatted, or else has no reference-matching reads.\n";
}


# Output variables
my( $background_readnum, $background_redund, $background_pccmd_frac );

# Four hashes for data, format:  read length => number of copies.
# For reads matching reference, and pCCMD reads matching reference, store both
# total number of reads and number of non-redundant ("nr") reads.
my( %matching_total, %matching_nr, %pccmd_total, %pccmd_nr );

# Parameters for calculating mean background values
my $low_bound_background = 46;   # background calculated for 46 bases and up
my $upper_bound = 80;            # largest read length
my $num_lengths = $upper_bound - $low_bound_background + 1;


### READ DATA FROM FASTA FILE

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

    # Keep only data for background values
    next if ($thislen < $low_bound_background);

    # Add read-number values to hashes
    $matching_total{$thislen} += $thiscopies;
    $matching_nr{$thislen}++;

    # Check if pCCMD and if so, add read-number values to hashes
    my $seq = $_;                     # DNA sequence of read
    my $top_c_pos = $thislen - 17;    # 0-based location of C in a pCCMD read
    my $bot_c_pos = 16;               # 0-based location of G in a pCCMD read
    if (substr( $seq, $top_c_pos, 1 ) eq 'C' && substr( $seq, $bot_c_pos, 1 ) eq 'G') {
        $pccmd_total{$thislen} += $thiscopies;
	$pccmd_nr{$thislen}++;
    }
}
close IN;


### CALCULATE BACKGROUND VALUES

# Get read totals
foreach my $len ($low_bound_background..$upper_bound) {

    # Populate missing values with 0's
    unless (exists $matching_total{$len}) {
	$matching_total{$len} = 0;
	$matching_nr{$len} = 0;
    }
    unless (exists $pccmd_total{$len}) {
	$pccmd_total{$len} = 0;
	$pccmd_nr{$len} = 0;
    }

    # Calculate redundancy of different types of reads
    my $redundancy = 0;          # total reads / nr reads
    my $frac_pccmd_all = 0;      # pccmd reads / total reads
    my $frac_pccmd_nr = 0;       # nr pccmd reads / nr reads ; not used further at present
    if ($matching_nr{$len} > 0) {
	$redundancy = int($matching_total{$len}/$matching_nr{$len} + 0.5);  # round to nearest positive integer
	$frac_pccmd_all = $pccmd_total{$len}/$matching_total{$len};
	$frac_pccmd_nr  = $pccmd_nr{$len}/$matching_nr{$len};    
    } else {  # No reads of this length, don't include in mean calculation
	$num_lengths--;
    }

    # Add to sum for mean
    $background_readnum += $matching_total{$len};
    $background_redund += $redundancy;
    $background_pccmd_frac += $frac_pccmd_all;
}

# Convert total numbers to mean values by dividing by number of lengths
$background_readnum = int( $background_readnum/$num_lengths );
$background_redund = int( $background_redund/$num_lengths );
$background_pccmd_frac = sprintf( "%.3f", $background_pccmd_frac/$num_lengths );


### PRINT TO OUTPUT

# Create output file
open( OUT, ">$outfile" ) || die "Failed to create $outfile: $!\n";
print OUT "----------------------------------------------------------------------------\n";
print OUT "INPUT FILE = $read_file\n";
print OUT "----------------------------------------------------------------------------\n";
print OUT "BACKGROUND PARAMETERS (read_stats.pl)\n";
print OUT "--------------------------------------\n";
print OUT "BACKGROUND LENGTH RANGE:  $low_bound_background TO $upper_bound\n\n";
print OUT "READ NUMBER:              $background_readnum\n";
print OUT "REDUNDANCY:               $background_redund\n";
print OUT "CC-READ FRACTION:         $background_pccmd_frac\n\n\n";
close OUT;

# Print also to STDOUT to pass to next pipeline module
print "BACKGROUND PARAMETERS: read_number_(all)=$background_readnum ; ";
print "redundancy=$background_redund ; pCCMD-read fraction=$background_pccmd_frac .\n";
