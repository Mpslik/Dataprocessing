# Genomic Data Analysis Project

## Overview
### this project
This project focuses on the analysis of genomic data to identify genes, proteins, and other genomic features from metagenomic samples. Utilizing a suite of bioinformatics tools and R for data visualization, the pipeline processes raw sequencing reads, performs assembly, gene prediction, functional annotation, and visualizes key metrics such as contig lengths, GC content, and coverage.

### Original pipeline:
The article [2] describes the development of a pipeline designed to accelerate the analysis of next-generation sequencing (NGS) data. This pipeline allows for the reliable and clear presentation of results for identifying both new and known viruses from associated or environmental samples, with an emphasis on virus discovery. The data [3] used for this purpose was derived from diarrheic American mink (Neovison vison).

![Original Diagram]( misc/Original_pipeline_diagram.jpeg )


### differences: 
the original article includes more steps than the current pipeline here. This project runs up to the assembly of the contigs and then performs an analyses step with prodigal. 
the change from MGA to prodigal occurred due to an error running the MGA on this combination of packages and version numbers and would configure in this project pipeline. 
The rest of the steps that branch off on the original diagram from the contigs file where also omitted due to time constraints of this school project.
and for the project a visualisation step is included that is not present in the original pipeline 


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

See the `visualize_results` rule in the Snakefile for more details.

## License
This project is licensed under the GPL-3.0 license - see the LICENSE.md file for details.

## Acknowledgments / Resources 
- [1] Tools and libraries: Snakemake, BWA, Samtools, and others.

- Data sources:
  - [2] **Article**: [Link to Article](https://academic.oup.com/ve/article/6/2/veaa091/6017186?login=false#373674411)
  - [3] **Data**: [Link to Data Repository](https://bitbucket.org/plyusnin/lazypipe/src/master)





