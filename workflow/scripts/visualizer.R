library(ggplot2)
library(GenomicFeatures)
library(GenomicAlignments)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)

# Function to calculate contig lengths
calc_contig_lengths <- function(contigs_file) {
  contigs <- readDNAStringSet(contigs_file)  # Read DNA sequences from the file
  lengths <- width(contigs)  # Calculate the lengths of each contig
  return(lengths)  # Return the lengths
}

# Function to calculate GC content for contigs
calc_gc_content <- function(contigs_file) {
  contigs <- readDNAStringSet(contigs_file)  # Read DNA sequences from the file
  gc_content <- letterFrequency(contigs, letters = c("G", "C"), as.prob = TRUE)  # Calculate GC content as a proportion
  gc_content <- rowSums(gc_content)  # Sum the GC content for each contig
  return(gc_content)  # Return the GC content
}

# Function to plot contig lengths
plot_contig_lengths <- function(lengths) {
  p <- ggplot(data.frame(Length = lengths), aes(x = Length)) +  # Create a data frame and set up the plot
    geom_histogram(binwidth = 1000, fill = "cornflowerblue", color = "black") +  # Add histogram with specified bin width and colors
    theme_minimal() +  # Use a minimal theme
    xlab("Contig Length") + ylab("Frequency") +  # Label the axes
    ggtitle("Distribution of Contig Lengths")  # Add a title
  return(p)  # Return the plot
}

# Function to plot GC content
plot_gc_content <- function(gc_content) {
  p <- ggplot(data.frame(GC = gc_content), aes(x = GC)) +  # Create a data frame and set up the plot
    geom_histogram(binwidth = 0.01, fill = "forestgreen", color = "black") +  # Add histogram with specified bin width and colors
    theme_minimal() +  # Use a minimal theme
    xlab("GC Content") + ylab("Frequency") +  # Label the axes
    ggtitle("Distribution of GC Content across Contigs")  # Add a title
  return(p)  # Return the plot
}

# Function to parse the provided Prodigal output format
parse_prodigal_output <- function(filepath) {
    content <- readLines(filepath)  # Read lines from the Prodigal output file
    seqnames <- character()  # Initialize vector for sequence names
    starts <- integer()  # Initialize vector for start positions
    ends <- integer()  # Initialize vector for end positions
    strands <- character()  # Initialize vector for strands
    ids <- character()  # Initialize vector for IDs
    seq_counter <- 1  # Initialize a counter for sequences

    for (i in seq_along(content)) {
        line <- content[i]  # Get the current line

        if (grepl("^\\s+CDS", line)) {  # Check if the line contains CDS information
            location <- gsub(".*CDS\\s+", "", line)  # Extract the location information
            location <- gsub("[<>]", "", location)  # Remove partial indicators
            strand <- ifelse(grepl("complement", location), "-", "+")  # Determine the strand
            location <- gsub("complement\\(|\\)", "", location)  # Clean up the location string
            start_end <- unlist(strsplit(location, "\\.."))  # Split the location into start and end
            note_line <- content[i + 1]  # Get the next line which contains the ID
            id <- gsub('.*ID="([^"]+).*', '\\1', note_line)  # Extract the ID
            seqname <- paste0("Sequence", seq_counter)  # Create a sequence name
            seq_counter <- seq_counter + 1  # Increment the sequence counter

            seqnames <- c(seqnames, seqname)  # Add to sequence names vector
            starts <- c(starts, as.integer(start_end[1]))  # Add to starts vector
            ends <- c(ends, as.integer(start_end[2]))  # Add to ends vector
            strands <- c(strands, strand)  # Add to strands vector
            ids <- c(ids, id)  # Add to IDs vector
        }
    }

    gr <- GRanges(
        seqnames = Rle(seqnames),
        ranges = IRanges(start = starts, end = ends),
        strand = strands,
        ID = ids
    )  # Create a GRanges object with the extracted information

    return(gr)  # Return the GRanges object
}

# Function to plot gene lengths from a GRanges object
plot_genes_from_gr <- function(gr) {
    gene_lengths <- width(gr)  # Get the lengths of the genes

    p <- ggplot(data.frame(Length=gene_lengths), aes(x=Length)) +  # Create a data frame and set up the plot
        geom_histogram(binwidth=50, fill="blue", color="black") +  # Add histogram with specified bin width and colors
        theme_minimal() +  # Use a minimal theme
        xlab("Gene Length") + ylab("Frequency") +  # Label the axes
        ggtitle("Distribution of Gene Lengths from Prodigal Output")  # Add a title
    return(p)  # Return the plot
}

# Function to plot protein lengths from a file
plot_proteins <- function(proteins_faa_file) {
    proteins <- readAAStringSet(proteins_faa_file)  # Read amino acid sequences from the file
    protein_lengths <- width(proteins)  # Calculate the lengths of each protein

    p <- ggplot(data.frame(Length=protein_lengths), aes(x=Length)) +  # Create a data frame and set up the plot
        geom_histogram(binwidth=10, fill="skyblue", color="black") +  # Add histogram with specified bin width and colors
        theme_minimal() +  # Use a minimal theme
        xlab("Protein Length") + ylab("Frequency") +  # Label the axes
        ggtitle("Distribution of Protein Lengths from Prodigal Output")  # Add a title
    return(p)  # Return the plot
}

# Function to read and bin coverage data
read_and_bin_coverage <- function(coverage_file, bin_size = 1000) {
  coverage_data <- fread(coverage_file, col.names = c("Ref", "Position", "Coverage"))  # Read coverage data from the file

  # Creating bins for positions
  coverage_data[, Bin := ceiling(Position / bin_size)]  # Create bins based on position

  # Summarizing coverage within bins
  binned_coverage <- coverage_data[, .(AverageCoverage = mean(Coverage)), by = .(Bin)]  # Calculate average coverage for each bin

  # Calculating the position for plotting as the center of each bin
  binned_coverage[, Position := (Bin - 0.5) * bin_size]  # Calculate the bin centers

  return(binned_coverage)  # Return the binned coverage data
}

# Function to plot the binned coverage data
plot_binned_coverage <- function(binned_coverage) {
  p <- ggplot(binned_coverage, aes(x = Position, y = AverageCoverage)) +  # Set up the plot with binned coverage data
    geom_line() +  # Add a line plot
    theme_minimal() +  # Use a minimal theme
    xlab("Position") + ylab("Average Coverage") +  # Label the axes
    ggtitle("Binned Genomic Coverage")  # Add a title

  return(p)  # Return the plot
}

# Main execution
genes_gff_file <- snakemake@input$genes_gff  # Get the input GFF file from Snakemake
proteins_faa_file <- snakemake@input$proteins_faa  # Get the input FAA file from Snakemake
contigs_file <- snakemake@input$contigs  # Get the input contigs file from Snakemake
output_pdf <- snakemake@output$visualization  # Get the output PDF file path from Snakemake

# Calculate metrics
contig_lengths <- calc_contig_lengths(contigs_file)  # Calculate the lengths of the contigs
gc_content <- calc_gc_content(contigs_file)  # Calculate the GC content of the contigs

# Open a PDF device to save all plots
pdf(output_pdf, width = 8, height = 10)  # Open the PDF device

# Existing plots
gr <- parse_prodigal_output(genes_gff_file)  # Parse the Prodigal output
p1 <- plot_proteins(proteins_faa_file)  # Plot the protein lengths
p2 <- plot_genes_from_gr(gr)  # Plot the gene lengths
print(p1)  # Print the protein lengths plot
print(p2)  # Print the gene lengths plot

# Include new plots for contig analysis
p_length <- plot_contig_lengths(contig_lengths)  # Plot the contig lengths
p_gc <- plot_gc_content(gc_content)  # Plot the GC content
print(p_length)  # Print the contig lengths plot
print(p_gc)  # Print the GC content plot

# Close the PDF device
dev.off()