# CDAC_Assignment
# Gene Phylogeny Pipeline

The `gene_phylo_pipeline.py` script is designed for comprehensive gene analysis and phylogenetic tree construction. It processes genomic data to identify, extract, and analyze paralogous gene sequences across multiple genomes. This script automates the tasks of BLAST searches, sequence extraction, CDS extraction, and sequence alignment.

## Features

- Genome Parsing: Reads and processes genome contig `.fasta` files.
- BLAST Search: Identifies target genes in genomic sequences with configurable parameters.
- Result Filtering: Filters BLAST results based on identity and query coverage.
- Sequence Extraction: Extracts paralogous gene sequences from genome contigs.
- CDS Extraction: Uses Prodigal to extract coding sequences from contigs.
- Sequence Alignment: Aligns gene sequences using ClustalW.

## Prerequisites

Ensure the following tools are installed and accessible in your environment:

- BLAST (`blastn`)
- Prodigal (`prodigal`)
- ClustalW (`clustalw`)
- Python 3.8+ with the following packages:
  - `biopython`
  - `pandas`

## Installation

1. Install Required Python Packages:
   `pip install biopython pandas`

2. Install BLAST:
   Follow instructions from the [BLAST website](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to download and install.

3. Install Prodigal:
   Download from the [Prodigal website](https://github.com/hyattpd/Prodigal) and follow the installation instructions.

4. Install ClustalW:
   Download from the [ClustalW website](http://www.clustal.org/clustal2/) and follow the installation instructions.

## Usage

1. Prepare Input Data:
   - Place genome contig `.fasta` files in the `genomes` directory.
   - Place gene FASTA files in the `Genes` directory.

2. Run the Script:
   `python gene_phylo_pipeline.py`

3. Output:
   - BLAST Results: Saved in the `results` directory.
   - Filtered BLAST Results: `filtered_combined_blast_results.csv` in the `results` directory.
   - Extracted Sequences: Saved in the `phylogeny` directory.
   - Aligned Sequences: Aligned outputs in the `phylogeny` directory.

## Script Details

- Parsing Genome Contigs: Reads `.fasta` files from the `genomes` directory and extracts sequences.
- BLAST Search: Executes `blastn` to identify target genes in the genomes.
- Filtering Results: Retains BLAST results with ≥ 95% identity and full query coverage.
- Extracting Sequences: Extracts paralogous sequences based on filtered BLAST results.
- CDS Extraction: Runs Prodigal to extract coding sequences from genome contigs.
- Aligning Sequences: Aligns gene sequences using ClustalW.
