# MFRE-Seq
Code to accompany the publication "Genome-Wide Identification of 5-Methylcytosine Sites in Bacterial Genomes By High-Throughput Sequencing of MspJI Restriction Fragments" (PLOS One, 2021, in review).

OVERVIEW

Code comprises a preprocessing perl script (process_reads.pl), plus a pipeline with a bash driver script (mfre_pipeline.sh) and four associated perl scripts (read_stats.pl, lengths_to_check.pl, find_motifs.pl, and final_motifs.pl).  The pipeline requires a FASTA file that is in the format specified by the preprocessing script, so preprocessing must be done first.  The preprocessing step can take more than an hour for large files, while the pipeline generally runs in seconds.

Preprocessing (process_reads.pl):
  Input (1) a FASTA file of Illumina or Ion Torrent sequence reads and (2) a FASTA file of the reference genome from which the reads were derived.
  Output is a FASTA file of deduplicated reads in which all of the reads exactly matching the reference are clustered at the beginning of the file.
  This output serves as input to the MFRE-Seq motif-finding pipeline.

Pipeline (mfre_pipeline.sh), sole input is the FASTA output file from the preprocessing step:
  Determine background parameters from non-motif-containing read lengths (46-80 bp).
  While motifs continue to be found:
    From read lengths in the 20-45 bp range, identify those that might contain methylated motifs (lengths_to_check.pl).
    For each of those lengths, see if there is a conserved motif (find_motifs.pl).
    Identify the most specific of those motifs and eliminate all reads derived from cleavage of that motif (find_motifs.pl).
    Create a reduced FASTA file of reads after elimination and return to the start of the loop.
  Print a summary of all motifs found, and the number of (16,16) and (16,17) reads derived from each (final_motifs.pl).

A summary of operations and results is printed to STDOUT as the pipeline runs, and a more verbose file of output is also created.

NOTES ON INDIVIDUAL PERL SCRIPTS

lengths_to_check.pl

Any length in the 20-45 bp range where the read number, read copy-number, or fraction of pCCMD reads is 2x the background parameter determined by read_stats.pl, AND where there are at least 20 unique pCCMD reads, is recommended to check.  This tends to be very permissive, but the effort to check non-productive lengths is minimal.
After the first iteration of the pipeline, this step is mostly redundant, since those lengths not meeting all those requirements except the minimum number of unique reads will have already been eliminated.

find_motifs.pl

This script will attempt to derive a single motif for each recommended length.
Multiple motifs may be obtained, some of which can be genuine and some of which can be artifacts.
From all motifs obtained, identify the motif most likely to be real, based on (1) the specificity of the motif, and (2) when the same motif is found at different read lengths, the number of unique reads at the given length (take the higher number).
Motif specificity is calculated by an empirical scoring system in which each base of the motif is scored based on degeneracy (the less degenerate, the higher the score), and then dividing by the total number of bases in the motif.
Once a most likely motif has been found, all reads derived from that motif (i.e., cut properly on at least one side) are eliminated from the pool of reads to be considered in the next iteration of the pipeline.
Motifs are identified by examining the nt distribution at each position along the read length and comparing it to the distributions expected for each IUPAC symbol given the background nucleotide content of the reference, using KL-divergence.

final_motifs.pl

Each iteration of the pipeline except the last will have identified one motif.  For each of these motifs, determine the number of (16,16) and (16,17) reads in the original set.
Output is permissive in the sense that (1) the same motif can be associated with multiple read lengths, and (2) the same read length can be associated with multiple motifs (each identified at a different iteration).  When these cases occur, it is generally because one of the two cases is an artifact, but we allow for these results so as not to overlook potentially biologically meaningful, but rare, cases.
Two sets of results are instructive, below.

In Agrobacterium gelatinovorum, we obtained the following final motifs:

               MOTIF     LEN_(16,16)   READS_(16,16)  REDUND_(16,16)     LEN_(16,17)   READS_(16,17)  REDUND_(16,17)
              ACCGGT              32             140            6.24              33               1            4.00
              ACCGGT              30            1296            9.38              31             973            7.17
               CCWGG              29            1314            6.73              30             721            5.58
               CCWGG              31              82            5.11              32               3            6.33
               
Both motifs appear twice, apparently because MFRE cleavage left an extra base on both sides a bit more than usual, causing the program to mistake what are actually (17,17) reads for (16,16) reads.  For ACCGGT, the results at length 30 are real, while those at length 32 are spurious.  For CCWGG, the results at length 29 are real, while those at length 31 are spurious.  Both the aboslute number of apparent (16,16) reads and the extremely small number of apparent (16,17) reads are indicative of this.  (The numbers are small, because the apparent (16,17) reads are in truth (17,18) reads, which are very rare.)

In Bacillus stearothermophilus, we obtained the following final motifs:

               MOTIF     LEN_(16,16)   READS_(16,16)  REDUND_(16,16)     LEN_(16,17)   READS_(16,17)  REDUND_(16,17)
              RCCGGY              30             523            4.59              31             224            5.23
              VCCGGB              30             573            4.54              31             240            5.13

Here we have two similar motifs obtained at the same read length.  The RCCGGY motif is real and was extracted first.  Because all VCCGGB reads are also RCCGGY, the number of "additional" reads to identify VCCGGB was only 50 (573 - 523), suggesting that the additional bases at the outside positions may be a detectable off-target activity of RCCGGY methylation.  Indeed, examination of the verbose output from iteration 1, which detected RCCGGY, shows trace activity at C in the R position and G in the Y position:

                           A        C        G        T
      12           R     0.241    0.045    0.714    0.000
      13          [C]    0.000    1.000    0.000    0.000
      14           C     0.003    0.993    0.000    0.003
      15           G     0.000    0.002    0.995    0.003
      16          [G]    0.000    0.000    1.000    0.000
      17           Y     0.000    0.640    0.047    0.314

Manual inspection and careful consideration of the pipeline results is always recommended.
