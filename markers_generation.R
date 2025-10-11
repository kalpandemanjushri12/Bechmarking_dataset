# Install and load required packages
install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)

# Set working directory
setwd("/home/Desktop")
getwd()

# Define the path to your FASTA file
fasta_file <- "chrN.fa"

# Read the FASTA file
chrN_seq <- readDNAStringSet(fasta_file)

# Extract the sequence
seq <- chrN_seq[[1]]

# Get the sequence length
seq_length <- length(seq)
cat("Length of the chromosome N sequence:", seq_length, "\n")

# Define marker size and step size
marker_size <- 2000
step_size <- 100000

# Create a list to store markers
markers <- list()

# Extract markers every 100,000 base pairs, each 2000 bases long
for (start_pos in seq(1, seq_length - marker_size + 1, by = step_size)) {
  marker <- subseq(seq, start = start_pos, end = start_pos + marker_size - 1)
  markers[[length(markers) + 1]] <- list(
    start = start_pos,
    end = start_pos + marker_size - 1,
    sequence = marker
  )
}

cat("Number of markers extracted:", length(markers), "\n")

# Save markers in FASTA format
output_fasta_file <- "chrN_100kb_2kblen_markers.fa"
fasta_lines <- character()

for (i in seq_along(markers)) {
  header <- paste0(">chrN:", markers[[i]]$start, "-", markers[[i]]$end, " marker", i)
  sequence <- as.character(markers[[i]]$sequence)
  fasta_lines <- c(fasta_lines, header, sequence)
}

writeLines(fasta_lines, output_fasta_file)
cat("Markers saved to:", output_fasta_file, "\n")

