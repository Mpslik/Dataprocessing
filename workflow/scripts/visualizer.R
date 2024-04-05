library(ggplot2)
library(GenomicFeatures)
library(GenomicAlignments)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)

# Function to calculate contig lengths
calc_contig_lengths <- function(contigs_file) {
  contigs <- readDNAStringSet(contigs_file)
  lengths <- width(contigs)
  return(lengths)
}

# Function to calculate GC content for contigs
calc_gc_content <- function(contigs_file) {
  contigs <- readDNAStringSet(contigs_file)
  gc_content <- letterFrequency(contigs, letters = c("G", "C"), as.prob = TRUE)
  gc_content <- rowSums(gc_content)
  return(gc_content)
}

# Function to plot contig lengths
plot_contig_lengths <- function(lengths) {
  p <- ggplot(data.frame(Length = lengths), aes(x = Length)) +
    geom_histogram(binwidth = 1000, fill = "cornflowerblue", color = "black") +
    theme_minimal() +
    xlab("Contig Length") + ylab("Frequency") +
    ggtitle("Distribution of Contig Lengths")
  return(p)
}

# Function to plot GC content
plot_gc_content <- function(gc_content) {
  p <- ggplot(data.frame(GC = gc_content), aes(x = GC)) +
    geom_histogram(binwidth = 0.01, fill = "forestgreen", color = "black") +
    theme_minimal() +
    xlab("GC Content") + ylab("Frequency") +
    ggtitle("Distribution of GC Content across Contigs")
  return(p)
}
# Function to parse the provided Prodigal output format
parse_prodigal_output <- function(filepath) {
    content <- readLines(filepath)
    seqnames <- character()
    starts <- integer()
    ends <- integer()
    strands <- character()
    ids <- character()
    seq_counter <- 1

    for (i in seq_along(content)) {
        line <- content[i]

        if (grepl("^\\s+CDS", line)) {
            location <- gsub(".*CDS\\s+", "", line)
            location <- gsub("[<>]", "", location)  # Remove partial indicators
            strand <- ifelse(grepl("complement", location), "-", "+")
            location <- gsub("complement\\(|\\)", "", location)
            start_end <- unlist(strsplit(location, "\\.."))
            note_line <- content[i + 1]
            id <- gsub('.*ID="([^"]+).*', '\\1', note_line)
            seqname <- paste0("Sequence", seq_counter)
            seq_counter <- seq_counter + 1

            seqnames <- c(seqnames, seqname)
            starts <- c(starts, as.integer(start_end[1]))
            ends <- c(ends, as.integer(start_end[2]))
            strands <- c(strands, strand)
            ids <- c(ids, id)
        }
    }

    gr <- GRanges(
        seqnames = Rle(seqnames),
        ranges = IRanges(start = starts, end = ends),
        strand = strands,
        ID = ids
    )

    return(gr)
}

plot_genes_from_gr <- function(gr) {
    gene_lengths <- width(gr)

    p <- ggplot(data.frame(Length=gene_lengths), aes(x=Length)) +
        geom_histogram(binwidth=50, fill="blue", color="black") +
        theme_minimal() +
        xlab("Gene Length") + ylab("Frequency") +
        ggtitle("Distribution of Gene Lengths from Prodigal Output")
    return(p)
}

# Adjusted plot_proteins function to return plot object instead of plotting directly
plot_proteins <- function(proteins_faa_file) {
    proteins <- readAAStringSet(proteins_faa_file)
    protein_lengths <- width(proteins)

    p <- ggplot(data.frame(Length=protein_lengths), aes(x=Length)) +
        geom_histogram(binwidth=10, fill="skyblue", color="black") +
        theme_minimal() +
        xlab("Protein Length") + ylab("Frequency") +
        ggtitle("Distribution of Protein Lengths from Prodigal Output")
    return(p)
}

# Adjusted plot_coverage function to calculate and return plot object
read_and_bin_coverage <- function(coverage_file, bin_size = 1000) {
  coverage_data <- fread(coverage_file, col.names = c("Ref", "Position", "Coverage"))

  # Creating bins for positions
  coverage_data[, Bin := ceiling(Position / bin_size)]

  # Summarizing coverage within bins
  binned_coverage <- coverage_data[, .(AverageCoverage = mean(Coverage)), by = .(Bin)]

  # Calculating the position for plotting as the center of each bin
  binned_coverage[, Position := (Bin - 0.5) * bin_size]

  return(binned_coverage)
}

# Function to plot the binned coverage data
plot_binned_coverage <- function(binned_coverage) {
  p <- ggplot(binned_coverage, aes(x = Position, y = AverageCoverage)) +
    geom_line() +
    theme_minimal() +
    xlab("Position") + ylab("Average Coverage") +
    ggtitle("Binned Genomic Coverage")

  return(p)
}


# Main execution
genes_gff_file <- snakemake@input$genes_gff
proteins_faa_file <- snakemake@input$proteins_faa
contigs_file <- snakemake@input$contigs  # Ensure this matches the actual input
output_pdf <- snakemake@output$visualization

# Calculate metrics
contig_lengths <- calc_contig_lengths(contigs_file)
gc_content <- calc_gc_content(contigs_file)

# Open a PDF device to save all plots
pdf(output_pdf, width = 8, height = 10)

# Existing plots
gr <- parse_prodigal_output(genes_gff_file)
p1 <- plot_proteins(proteins_faa_file)
p2 <- plot_genes_from_gr(gr)
print(p1)
print(p2)

# Include new plots for contig analysis
p_length <- plot_contig_lengths(contig_lengths)
p_gc <- plot_gc_content(gc_content)
print(p_length)
print(p_gc)

# Close the PDF device
dev.off()