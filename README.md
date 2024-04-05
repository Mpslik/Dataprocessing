# Genomic Data Analysis Project

## Overview
This project focuses on the analysis of genomic data to identify genes, proteins, and other genomic features from metagenomic samples. Utilizing a suite of bioinformatics tools and R for data visualization, the pipeline processes raw sequencing reads, performs assembly, gene prediction, functional annotation, and visualizes key metrics such as contig lengths, GC content, and coverage.

## Getting Started

### Prerequisites
- Conda environment
- Python 3.9
- R 4.1.0
- Snakemake
- Bioinformatics tools: BWA, Samtools, MegaHit, Trimmomatic, Prodigal
- NCBI datasets CLI

### Installation
Clone the repository:
```bash
git clone https://github.com/Mpslik/Dataprocessing
cd Dataprocessing/
```
Set up the Conda environment:

conda env create -f environment.yml
conda activate ngs-analysis-core

markdown
Copy code

Download and prepare the reference genome and raw reads:

- Refer to the `download_genome_data` and `download_reads` rules in the Snakefile for instructions.

### Usage
Configure your analysis:

- Edit the `config/config.yaml` file to specify your samples, paths to raw reads, and other configuration options.

Run the pipeline:
```bash
snakemake --cores [number_of_cores] -s workflow/Snakefile
```
markdown
Copy code

- Replace `[number_of_cores]` with the number of CPU cores you wish to use.


## Pipeline Steps
1. **Read Preprocessing:** Quality control and trimming of raw reads.
2. **Genome Assembly:** Assembling reads into contigs.
3. **Alignment:** Mapping reads back to the reference genome.
4. **Gene Prediction:** Identifying gene structures in the assembled genome.
5. **Functional Annotation:** Annotating predicted genes with functional information.
6. **Visualization:** Generating visual summaries of key metrics.
![DAG Visualization](dag.png "Directed Acyclic Graph")

### Visualization
The visualization step produces a PDF report containing:
- Distribution of contig lengths
- GC content across contigs
- Coverage depth distribution
- Additional plots as configured

See the `visualize_results` rule in the Snakefile for more details.

## License
This project is licensed under the GPU3 - see the LICENSE.md file for details.

## Acknowledgments
- Tools and libraries: Snakemake, BWA, Samtools, and others.
- Data sources: NCBI datasets.
Now, the markdown is correctly formatted without breaking the code block structure. Remember to replace [repository name] and [number_of_cores] with your specific details.








