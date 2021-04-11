#!/usr/bin/perl -w

use strict;

# Convert a FASTA file of Illumina or Ion Torrent reads into a format suitable for
# the motif-finding pipeline.  This involves (1) collapsing duplicate reads into a
# single entry, with the number of copies in the header line, (2) mapping reads to
# a genome reference and marking all reads that are exact matches, (3) ordering the
# reads in the output file so that reference matches are listed first, and (4) removing
# any reads outside the range of 20-80 bp.

# The range parameters can be adjusted as necessary, but the motif-finding pipeline
# uses reads in the 46-80 bp range to determine background parameters and reads in the
# 20-45 bp range for motif searching.

# Inputs are (1) path to a FASTA-formatted file of sequence reads, and (2) path to a
# genome reference sequence file, also in FASTA format.  Output files are (1) a 
# "processed FASTA file" formatted as described above, and (2) a summary text file 
# describing the overall number of reads found, the matching process, and a breakdown 
# by length.  Output files will be in the same directory as the sequence read file.

# A progress report and abbreviated results are also printed to STDOUT.


### DEFINE VARIABLES

# Input FASTA file paths
my( $read_path, $ref_path ) = @ARGV;
unless (defined $read_path && defined $ref_path) {
    die "Need to supply both a sequence-read file and a reference genome file.\n";
}

# Get read filename from path
my $read_filename = $read_path;
if ($read_path =~ /.+\/(.*)/) {
    $read_filename = $1;
}
$read_filename =~ /(.+)\./;
my $base_name = $1;

# Output files
my $outfile_reads = $base_name . "_proc.fasta";
my $outfile_summary = $base_name . "_summary.txt";

# Get reference from path
my $ref_filename = $ref_path;
if ($ref_path =~ /.+\/(.*)/) {
    $ref_filename = $1;
}
my $ref_name;
if ($ref_filename =~ /(.+)\.f/) {
    $ref_name = $1;
} else {
    die "Reference should be a FASTA file (extension .fasta, .fsa, .fas, etc.)\n";
}
$ref_name = $1;

# Read length parameters
my $MIN_LENGTH = 20;
my $MAX_LENGTH = 80;

my $refseq;            # reference genome sequence

# Store non-redundant sequences from the read file, separating those in the $MIN-$MAX range.
# Some sequence output will not have any out-of-range sequences.
my( %in_range_reads_copies, %outof_range_reads_copies );     # sequence => no. copies

# Keep track of numbers of total and non-redundant reads
my( $in_range_reads_total, $in_range_reads_nr, $outof_range_reads_total, $outof_range_reads_nr ) = (0) x 4;

# Also keep track of numbers of reads of each length
my( %length_total, %length_nr );   # read-length => number of reads


### READ DATA FROM FILES

# Get reference sequence from file
print "Reading reference and sequence-read files....\n";
open( IN, "$ref_path" ) || die "Failed to open $ref_path: $!\n";
while( <IN> ) {
    next if ( /^>/ );
    chomp;
    $refseq .= $_;
}
close IN;

# Get sequence reads and length data from FASTA file
open( IN, $read_path ) || die "Failed to open $read_path: $!\n";
while( <IN> ) {
    chomp;
    next if ( /^>/ );
    my $read = $_;
    my $this_len = length( $read );
    my $read_rc = rev_comp( $read );
    $length_total{$this_len}++;

    # This process should collapse reverse-complements into one read entry
    if ($this_len >= $MIN_LENGTH && $this_len <= $MAX_LENGTH) {  # length in range
	$in_range_reads_total++;
	if (exists $in_range_reads_copies{$read}) {
	    $in_range_reads_copies{$read}++;
	} elsif (exists $in_range_reads_copies{$read_rc}) {
	    $in_range_reads_copies{$read_rc}++;
	} else {  # haven't seen this read before
	    $in_range_reads_copies{$read} = 1;
	    $length_nr{$this_len}++;
	}
    } else {  # length out of range
	$outof_range_reads_total++;
	if (exists $outof_range_reads_copies{$read}) {
            $outof_range_reads_copies{$read}++;
        } elsif (exists $outof_range_reads_copies{$read_rc}) {
            $outof_range_reads_copies{$read_rc}++;
        } else {  # haven't seen this read before
            $outof_range_reads_copies{$read} = 1;
	    $length_nr{$this_len}++;
        }
    }

}
$in_range_reads_nr = () = keys %in_range_reads_copies;
$outof_range_reads_nr = () = keys %outof_range_reads_copies;

# Print update to STDOUT
print "Range:               $MIN_LENGTH to $MAX_LENGTH\n";
print "Reads in-range:      total = $in_range_reads_total\tnr = $in_range_reads_nr\n";
print "Reads out-of-range:  total = $outof_range_reads_total\tnr = $outof_range_reads_nr\n\n";
print "Reads from selected length range below.\n\n";
print "Length\tTotal\tNR\n\n";
foreach my $lg (sort {$a <=> $b} keys %length_total) {
    print "$lg\t$length_total{$lg}\t$length_nr{$lg}\n" if ($lg > 25 && $lg < 46);
}
print "\n";


### MAP READS TO THE REFERENCE SEQUENCE

print "Mapping $in_range_reads_nr unique reads to the reference files....\n";

# Reads matching reference or not; match + nomatch = in_range for both total and nr
my( %match_reads_copies, %nomatch_reads_copies );         # read sequence => number of copies
my( $match_reads_total, $match_reads_nr, $nomatch_reads_total, $nomatch_reads_nr ) = (0) x 4;

# Split %in_range_reads_copies into %match_reads_copies and %nomatch_reads_copies
my $num_mapped = 0;
foreach( keys %in_range_reads_copies ) {
   
    # Check both strands of read against reference; don't know which orientation the hash
    # representation is in relative to the reference at this point.  In %match_reads_copies,
    # all reads will be on top strand of reference.
    my $fwd = $_;                 # original read sequence as represented in hash
    my $rc = rev_comp( $fwd );    # reverse complement of same read
    if ( $refseq =~ /$fwd/ ) {
	$match_reads_total += $in_range_reads_copies{$fwd};
	$match_reads_copies{$fwd} = $in_range_reads_copies{$fwd};
    } elsif ( $refseq =~ /$rc/ ) {
	$match_reads_total += $in_range_reads_copies{$fwd};
	$match_reads_copies{$rc} = $in_range_reads_copies{$fwd};
    } else {  # not a match to the reference
	$nomatch_reads_total += $in_range_reads_copies{$fwd};
	$nomatch_reads_copies{$fwd} = $in_range_reads_copies{$fwd};
    }
    $num_mapped++;
    print "Finished $num_mapped...\n" if ($num_mapped % 10000 == 0);
}
close IN;
$match_reads_nr  = () = keys %match_reads_copies;
$nomatch_reads_nr = () = keys %nomatch_reads_copies;

# Print update to STDOUT
print "\n";
print "Reads matching reference:  total = $match_reads_total\tnr = $match_reads_nr\n";
print "Reads not matching:        total = $nomatch_reads_total\tnr = $nomatch_reads_nr\n\n";

# Mean number of copies (= total reads / non-redundant reads)
my $in_range_meancopies = $in_range_reads_nr ? int( $in_range_reads_total/$in_range_reads_nr ) : 'n/d';
my    $match_meancopies = $match_reads_nr    ? int( $match_reads_total/$match_reads_nr )       : 'n/d';
my  $nomatch_meancopies = $nomatch_reads_nr  ? int( $nomatch_reads_total/$nomatch_reads_nr)    : 'n/d';


### PRINT OUTPUT FASTA FILE

# Only the reference-matching reads are used by the motif-finding pipeline, but all reads are printed
# to the output FASTA file.

print "Printing output files....\n";
open( OUT, ">$outfile_reads" ) || die "Failed to create $outfile_reads: $!\n";

# Print reference matches first
my( %length_match_total, %length_match_nr );   # read-length => number of reads
foreach my $seq (sort keys %match_reads_copies) {
    my $len = length( $seq );
    $length_match_total{$len} += $match_reads_copies{$seq};
    $length_match_nr{$len}++;
    print OUT ">$base_name|reference $ref_name|length $len|copies $match_reads_copies{$seq}\n$seq\n";
}

# Then print non-matches
foreach my $seq (sort keys %nomatch_reads_copies) {
    my $len = length( $seq );
    print OUT ">$base_name|no match|length $len|copies $nomatch_reads_copies{$seq}\n$seq\n";
}

# Then print out-of-range sequences
foreach my $seq (sort keys %outof_range_reads_copies) {
    my $len = length( $seq );
    print OUT ">$base_name|out of range|length $len|copies $outof_range_reads_copies{$seq}\n$seq\n";
}
close OUT;


### PRINT OUTPUT SUMMARY FILE

# Summary table for all types of reads
open( OUT, ">$outfile_summary" ) || die "Failed to create $outfile_summary: $!\n";
print OUT "# Input files #\n\n";
print OUT "\tSequence read file = $read_filename\n";
print OUT "\t    Reference file = $ref_filename\n\n\n";
print OUT "# Summary of reads in the FASTA output file $outfile_reads #\n\n";
printf( OUT "%20s %15s %20s %15s\n\n", "Read_Set", "Total_Reads",
    "Non-redundant_reads", "Mean_copies" );
printf OUT "%20s %15s %20s %15s\n", "Reference matches", $match_reads_total, $match_reads_nr,
    $match_meancopies;
printf OUT "%20s %15s %20s %15s\n", "No matches", $nomatch_reads_total, $nomatch_reads_nr,
    $nomatch_meancopies;
printf OUT "%20s %15s %20s %15s\n\n\n", "Out of range", $outof_range_reads_total, 
    $outof_range_reads_nr, "-";

# Length distribution for reference-matching reads only
print OUT "# Numbers of reads of each length in range; reference = $ref_name #\n\n";
printf( OUT "%15s %15s %15s %15s %20s\n\n", "Length", "Total", "NR", "Total_Ref_Match",
	"NR_Ref_Match");
foreach my $len (sort {$a <=> $b} keys %length_match_total) {
    printf( OUT "%15s %15s %15s %15s %20s\n", $len, $length_total{$len}, $length_nr{$len}, 
	    $length_match_total{$len}, $length_match_nr{$len} );
}
close OUT;



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
