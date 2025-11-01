#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file=$1

# Extract contig lengths and contig IDs
awk '/^>/ {if (seqlen){print seqlen "\t" scaffold}; scaffold=$0; seqlen=0; next} {seqlen += length($0)} END {print seqlen "\t" scaffold}' "$input_file" > contig_lengths_with_ids.txt

# Sort contig lengths in descending order
sort -nrk1,1 contig_lengths_with_ids.txt > sorted_contig_lengths_with_ids.txt

# Calculate total length of all contigs
total_length=$(awk '{sum += $1} END {print sum}' sorted_contig_lengths_with_ids.txt)

# Calculate N50 and L50
awk -v total_length=$total_length '
BEGIN {sum=0; n50=0; l50=0}
{
  sum += $1;
  l50++;
  if (sum >= total_length / 2 && n50 == 0) {
    n50 = $1;
    print "N50: " n50;
    print "L50: " l50;
    exit;
  }
}' sorted_contig_lengths_with_ids.txt

# Find the maximum length
max_length=$(sort -nr contig_lengths_with_ids.txt | head -n 1)

# Output the longest contig length
echo "The length of the longest contig is: $max_length"

# Output of the genome size
echo "The Genome size: $total_length"
