#!/usr/bin/perl -w

use strict;

# From a FASTA file of reads (output from process_reads.pl) and a reference sequence in a 
# separate FASTA file, get frequency statistics for a given motif.  This script is
# particularly useful for degenerate motifs, which this script will parse into all 
# possible non-degenerate sequences and determine statistics separately for each.

# Inputs are the two FASTA files (reads and reference), a sequence motif, and an expected
# (16,16) read length.  (The expected readlength will enable determination of the methylated
# bases within the motif.)


### DEFINE VARIABLES

# Inputs and variable checking
my( $ref_path, $read_path, $motif, $expected_length ) = @ARGV;
unless (defined $ref_path && defined $read_path && defined $motif && defined $expected_length) {
    die "Need to supply:  path_to_reference_file, path_to_sequence_read_file, motif, read_length.\n";
}
if ($motif !~ /[ACGTRYSWKMBDHVNacgtryswkmbdhvn]/) {
    die "Impermissible characters in motif string.\n";
}
if ($expected_length =~ /\D/ || $expected_length < 20 || $expected_length > 80) {
    die "Expected read length should be a number between 20 and 80.\n";
}

# Output file name
my $read_filename = $read_path;
if ($read_path =~ /.+\/(.*)/) {
    $read_filename = $1;
}
my $base_name;
if ($read_filename =~ /(.+)\.f/) {
    $base_name = $1;
} else {
    die "Read file should be a FASTA file (extension .fasta, .fsa, .fas, etc.)\n";
}
my $outfile = $base_name . "_" . $motif . "_" . $expected_length . "_stats.txt";

# Get reference from path
my $ref_filename = $ref_path;
if ($ref_path =~ /.+\/(.*)/) {
    $ref_filename = $1;
}

# Data variables
my $refseq;            # sequence of the reference genome
my $nonred_total;      # expected number of non-redundant, reference-matching reads 20-45 bp
my %read_set;          # all reads within 1 bp of $expected_length; sequence => length


### READ DATA FROM FILES

# Get reference sequence from file
open( IN, "$ref_path" ) || die "Failed to open $ref_path©: $!\n";
while( <IN> ) {
    chomp;
    next if ( /^>/ );
    $refseq .= $_;
}
close IN;

# Read each reference-matching entry in the FASTA file around $expected_length
my $thislen;  # length of present read
open( IN, $read_path ) || die "Failed to open $read_path: $!\n";
while (<IN>) {
    chomp;

    # Get read length from ID line
    if ( /^>/ ) {
	last unless ( /reference.*length (\d+)\|copies/ );
	$thislen = $1;
	next;
    }

    # Keep only those reads within 1 bp of $expected_length
    next unless ($thislen >= $expected_length-1 && $thislen <= $expected_length+1);

    my $seq = $_;
    $read_set{ $seq } = $thislen;
}
close IN;
my $total_reads = () = keys %read_set;


### CREATE OUTPUT FILE

open( OUT, ">$outfile" ) || die "Failed to create $outfile: $!\n";
print OUT "# Data summary #\n\n";
printf( OUT "%30s    %s\n", "Read file:", $read_filename );
printf( OUT "%30s    %s\n", "Reference genome:", $ref_filename );
printf( OUT "%30s    %s\n", "Motif ((16,16) length):", "$motif ($expected_length)" );
printf( OUT "%30s %8s\n\n", "Non-redundant reads used:", $total_reads );

print OUT "# Genome data #\n";
print OUT "\n\t####  Motif $motif  (read length $expected_length)  ###########\n\n";
printf( OUT "\t%10s %15s %15s %15s %15s\n", "Motif", "Sites in",
	"Fraction", "No. Covered", "Frac. Covered" );
printf( OUT "\t%10s %15s %15s %15s %15s\n\n", "Instance", "Genome",
	"of Motif", "by Reads", "by Reads" );


### CALCULATE AND PRINT STATS FOR MOTIF AND ALL NON-DEGENERATE VERSIONS

my %read_data;   # text strings for printing to output file; submotif => text_string

# Get non-degenerate instances of the motif
my @nondegen_set = get_nondegen( $motif );
my $total_motif_sites;    # total number of sites in genome for $motif
my $total_motif_reads;    # total number of reads derived from cutting at $motif

# Analyze motif and then all non-degenerate versions
print "Getting data for motif $motif....\n";
foreach my $submotif ($motif, @nondegen_set) {
    print "\tWorking on $submotif....\n";
    my $submotif_regex = makeregex( $submotif );

    # Coordinates of all instances of $submotif in genome
    my %submotif_starts;

    # Coordinates of all possible (16,16), (16,17), and (15,16) reads for this $submotif
    my %expected_read_starts;   # expected read start => corresponding submotif start
    
    # Find all submotif and corresponding read coordinates
    while ($refseq =~ /(?=($submotif_regex))/g) {
	# Starting coordinate of submotif
	my $submotif_site = length($`);    
	$submotif_starts{ $submotif_site } = 0;

	# Starting coord of corresponding (16,16) ("perfect") read;
	# ignore cases at very ends of genome for now
	my $perfect_read_start = $submotif_site - ($expected_length - length($motif))/2;
	next if ($perfect_read_start < 0 || $perfect_read_start+$expected_length > length($refseq));

	# Account for (15,16), (16,16), and (16,17) reads
	for my $start ($perfect_read_start-1..$perfect_read_start+1) {
	    $expected_read_starts{ $start } = $submotif_site;
	}
    }
    
    # Total number of sites in genome for $submotif
    my $total_submotif_sites = () = keys %submotif_starts;
    $total_motif_sites = $total_submotif_sites if ($submotif eq $motif);
    
    # Number of sequence reads corresponding to (15,16), (16,16), and (16,17) reads for $submotif
    my $reads_with_submotif = 0;

    # Number of genome sites of $submotif for which there is >=1 read
    my $sites_covered = 0;

    # Check every read to see if it is one of the expected reads for this $submotif
    foreach my $read (keys %read_set) {

	# Work only with reads that are exact matches to reference
	next unless ($refseq =~ /$read/);

	# Work only with reads that are derived from cleavage of a $submotif instance
	my $read_start = length( $` );
	if (exists $expected_read_starts{ $read_start }) {
	    my $corresponding_submotif_start = $expected_read_starts{ $read_start };

	    # The corresponding submotif is covered by reads, so delete it if not deleted already
	    if (exists $submotif_starts{ $corresponding_submotif_start }) {
		delete ($submotif_starts{ $corresponding_submotif_start });
		$sites_covered++;
	    }
	    $reads_with_submotif++;
	}
    }
    $total_motif_reads = $reads_with_submotif if ($submotif eq $motif);

    unless ($total_motif_sites) {
	# Something went wrong, since $motif always processed before submotifs
	die "No sites for motif $motif identified in $ref_filename.\n";
    }

    # Print $submotif-centric data to output file
    if ($total_submotif_sites > 0) {
	printf( OUT "\t%10s %15s %15.3f %15s %15.3f\n", $submotif, $total_submotif_sites,
		$total_submotif_sites/$total_motif_sites, $sites_covered, $sites_covered/$total_submotif_sites );
    } else {
	printf( OUT "\t%10s %15s %15.3f %15s %15.3f\n", $submotif, $total_submotif_sites,
		$total_submotif_sites/$total_motif_sites, $sites_covered, 0 );
    }
    print OUT "\n" if ($submotif eq $motif);
    
    # Collect read-centric data to print to output file later
    $read_data{$submotif} = $reads_with_submotif . " " . ($reads_with_submotif/$total_motif_reads);
}


### PRINT READ-CENTRIC DATA TO OUTPUT FILE

print "Printing read data to file....\n";

print OUT "\n\n# Read data #\n";
print OUT "\n\t####  Motif $motif  (read length $expected_length)  ###########\n\n";
printf( OUT "\t%10s %15s %15s\n", "Motif", "No. Reads", "Frac. Reads" );
printf( OUT "\t%10s %15s %15s\n\n", "Instance", "w/ Motif", "w/ Motif" );

# Print data for motif and all non-degenerate versions
foreach my $submotif ($motif, @nondegen_set) {
    my @data = split " ", $read_data{$submotif};
    printf( OUT "\t%10s %15s %15.3f\n", $submotif, @data );
    print OUT "\n" if ($submotif eq $motif);
}

close OUT;



## Subroutine to return "non-degenerate" instances of a motif (not incl. N's)
##
sub get_nondegen {
    my $degen = $_[0];   # degenerate motif
    my @nondegen = ( 'X' );        # set of non-degenerate equivalents
    my @temp;                      # holds updated version of @nondegen

    my %code = ( # expand each symbol to degenerate possibilities                                                                        
		 A => ['A'],     T => ['T'],     C => ['C'],   G => ['G'],
		 R => ['A','G'],   Y => ['C','T'],   S => ['C','G'], 
		 W => ['A','T'],   M => ['A','C'],   K => ['G','T'],
		 H => ['A','C','T'], B => ['C','G','T'], V => ['A','C','G'], D => ['A','G','T'],
		 N => ['A','C','G','T']
	);

    # Walk along each base of $degen motif and get the base
    for my $i (0..length($degen)-1) {
	my $degen_base = substr( $degen, $i, 1 );

	# For each specific instance of this base, append to all @nondegen strings
	foreach my $j ( @{ $code{$degen_base} } ) {
	    foreach my $n (@nondegen) {
		push( @temp, $n . $j );
	    }
	}

	# Update @nondegen
	@nondegen = @temp;
	undef( @temp );
    }

    # Remove leading 'X' from each
    foreach( @nondegen ) {
	s/^.//;
    }
    return @nondegen;
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
