#!/bin/bash

#$ -S /bin/bash
#$ -pe smp 8

set -e

echo ""
echo "Started on $(date) on $(uname -n)"

### manually entered variables
root="/home/anton/work/mfre-seq/for_github"
scripts_dir="pipeline_scripts"
temp="temp_reads.fasta"

### sole argument is the FASTA read file
READFILE=$1

### initialize for first iteration
motifs_found="INIT"
input_file="$READFILE"
output_file="$READFILE.mfre-seq.out"
iteration=1
final_motifs=""

### run initial read analysis to get background parameters
echo "--------------------"
echo "Background parameters"
echo "--------------------"
echo "Running initial sequence-read analysis..."
read_stats_output=`$root/$scripts_dir/read_stats.pl $input_file $output_file`
echo "   $read_stats_output"
echo ""

### capture background read-number, redundancy, and pCCMD-fraction parameters
rsout="read_number_\(all\)=(.+) ; redundancy=(.+) ; pCCMD-read fraction=(.+) \."
[[ "$read_stats_output" =~ $rsout ]] && rn="${BASH_REMATCH[1]}" && re="${BASH_REMATCH[2]}" && cc="${BASH_REMATCH[3]}"

### iterate until no more motifs are found
while [ $motifs_found != "0" ]; do
    # run read analysis to get background parameters
    echo "--------------------"
    echo "Iteration $iteration"
    echo "--------------------"
    echo "Determining lengths to check..."
    lengths_output=`$root/$scripts_dir/lengths_to_check.pl $input_file $output_file $rn $re $cc $iteration`
    echo "   $lengths_output"
    echo ""

    # capture lengths from output
    length_string="LENGTHS TO CHECK=(.*)$"
    [[ "$lengths_output" =~ $length_string ]] && len="${BASH_REMATCH[1]}"

    # exit loop if no lengths worth checking
    if [[ $len == '0' ]]; then
	break
    fi

    # otherwise, run motif finding
    echo "Running motif-finding..."
    find_motifs_output=`$root/$scripts_dir/find_motifs.pl $input_file $output_file $temp "$len" $re $iteration`
    echo "   $find_motifs_output"

    # capture most specific motif
    re_mot="MOST_SPECIFIC_MOTIF=(.*)$"
    [[ "$find_motifs_output" =~ $re_mot ]] && motifs_found="${BASH_REMATCH[1]}"
    if [ $motifs_found = "0" ]; then
	echo "   No motifs found; exiting."
    else
	echo "   Deleted reads derived from $motifs_found and trying again"
	final_motifs+=" $motifs_found"
    fi
    echo ""
    input_file=$temp
    iteration=$((iteration+1))

    # run a maximum of 20 iterations to avoid possible infinite loops
    if [[ "$iteration" == '21' ]]; then
	echo "Reached maximum of 20 iterations."
	break
    fi
done

final_output=`$root/$scripts_dir/final_motifs.pl $READFILE $output_file "$final_motifs" $re`
echo "--------------------"
echo "Final Motif Set"
echo "--------------------"
echo "$final_output"

# clean up; remove temp FASTA file and move output to current directory
`rm $temp`
`mv $output_file .`