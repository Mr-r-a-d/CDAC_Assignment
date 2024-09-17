import os
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import subprocess
from Bio.Align.Applications import ClustalwCommandline

# Define base directory
base_dir = Path("/home/enigma101/CDAC/Test")

# Define paths
results_dir = base_dir / "results"
phylogeny_dir = results_dir / "phylogeny"
genome_dir = base_dir / "genomes"
genes_dir = base_dir / "Genes"

# Create results and phylogeny directories if they don't exist
results_dir.mkdir(parents=True, exist_ok=True)
phylogeny_dir.mkdir(parents=True, exist_ok=True)

# Step 1: Parse Genome Contigs
genome_files = [file.name for file in genome_dir.iterdir() if file.suffix == ".fasta"]
contigs = []
for genome_file in genome_files:
    genome_file_path = genome_dir / genome_file
    try:
        contigs.extend(list(SeqIO.parse(genome_file_path, "fasta")))
    except FileNotFoundError:
        print(f"Error: File not found - {genome_file_path}")
    except Exception as e:
        print(f"Error: Failed to parse file - {genome_file_path}: {str(e)}")

# Step 2: Verify Gene Names in FASTA Files
gene_fasta_files = list(genes_dir.glob("*.fasta"))
gene_names = [Path(file).stem for file in gene_fasta_files]
print(f"Gene names in FASTA files: {gene_names}")

# Define genes_of_interest
genes_of_interest = set(gene_names)  # Initialize the variable

# Step 3: Perform BLAST Search for Genes of Interest
query_files = gene_fasta_files
subject_files = [genome_dir / file for file in genome_files]
output_files = [results_dir / f"blast_results_{Path(subject_file).stem}_{Path(query_file).stem}.txt"
                for subject_file in subject_files
                for query_file in query_files]

for subject_file in subject_files:
    for query_file in query_files:
        output_file = results_dir / f"blast_results_{Path(subject_file).stem}_{Path(query_file).stem}.txt"
        blast_command = f"blastn -query {query_file} -subject {subject_file} -out {output_file} -outfmt 6 -perc_identity 95"
        result = subprocess.run(blast_command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error running BLAST on {subject_file} with query {query_file}: {result.stderr}")
        else:
            print(f"BLAST completed for {subject_file} with query {query_file}")

# Step 4: Parse and Filter BLAST Results
blast_results = []
for output_file in results_dir.glob("blast_results_*.txt"):
    if output_file.stat().st_size > 0:
        blast_result = pd.read_csv(output_file, sep='\t', header=None)
        blast_result.columns = ['query', 'subject', 'identity', 'alignment_length', 'mismatches', 'gap_opens', 
                                 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
        blast_results.append(blast_result)
    else:
        print(f"Warning: No BLAST results found in file {output_file}")

# Combine all filtered results into a single DataFrame
combined_results = pd.concat(blast_results, ignore_index=True)

# Filter by identity >= 95% and query coverage >= 100%
filtered_results = combined_results[(
    combined_results['identity'] >= 95) & 
    (combined_results['alignment_length'] == (combined_results['q_end'] - combined_results['q_start'] + 1)) & 
    (combined_results['query'].isin(genes_of_interest))]

# Save the combined and filtered BLAST results
filtered_results.to_csv(results_dir / "filtered_combined_blast_results.csv", index=False)

def extract_sequences_from_genome(genome_file, blast_results_df, gene_names, output_dir):
    genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    
    for gene_name in gene_names:
        gene_blast_results = blast_results_df[blast_results_df['query'] == gene_name]
        
        for _, row in gene_blast_results.iterrows():
            subject_id = row['subject']
            
            # First, try exact match of subject ID
            if subject_id in genome_sequences:
                contig_seq = genome_sequences[subject_id].seq
            else:
                # Clean the subject ID by stripping any suffix after "_"
                cleaned_subject_id = subject_id.split('_')[0]
                
                if cleaned_subject_id in genome_sequences:
                    contig_seq = genome_sequences[cleaned_subject_id].seq
                else:
                    print(f"Warning: Neither full nor cleaned Subject ID {subject_id} found in genome sequences for {genome_file}")
                    continue

            # Extract sequence
            start = min(row['s_start'], row['s_end']) - 1  # Ensure valid start position
            end = max(row['s_start'], row['s_end'])        # Ensure valid end position

            if start < end and end <= len(contig_seq):
                homologous_seq = contig_seq[start:end]

                with open(output_dir / f"{gene_name}_{subject_id}.fasta", "w") as handle:
                    handle.write(f">{gene_name}_{subject_id}\n")
                    handle.write(str(homologous_seq) + "\n")

                print(f"Extracted sequence for {gene_name} from {subject_id} in {genome_file}")
            else:
                print(f"Warning: Invalid range {start+1}-{end} for {subject_id} in {genome_file}. Skipping...")

# Process each genome and extract sequences based on BLAST results
for genome_file in genome_files:
    genome_file_path = genome_dir / genome_file
    extract_sequences_from_genome(genome_file_path, filtered_results, genes_of_interest, phylogeny_dir)

print("Homologous regions extraction completed.")

# Step 6: Extract CDS with Prodigal
cds_files = [phylogeny_dir / f"cds_output_{Path(file).stem}.faa" for file in genome_files]
for contig_file, cds_file in zip(subject_files, cds_files):
    prodigal_command = f"prodigal -i {contig_file} -a {cds_file}"
    result = subprocess.run(prodigal_command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running Prodigal on {contig_file}: {result.stderr}")
    else:
        print(f"Prodigal CDS extraction completed for {contig_file}")

# Step 7: Align sequences
# Define the list of genes and their corresponding fasta files
genes = ['gapA', 'infB', 'mdh', 'pgi', 'phoE', 'rpoB', 'tonB']

# Loop over each gene and align sequences
for gene in genes:
    # Get all files for this gene
    gene_files = list(phylogeny_dir.glob(f"{gene}_*.fasta"))
    
    # Combine the sequences for this gene into a single fasta file
    combined_fasta = phylogeny_dir / f"{gene}_combined.fasta"
    
    with open(combined_fasta, "w") as output_handle:
        for file in gene_files:
            for record in SeqIO.parse(file, "fasta"):
                SeqIO.write(record, output_handle, "fasta")
    
    # Align the sequences using ClustalW
    aligned_output = phylogeny_dir / f"aligned_{gene}.aln"
    clustalw_cline = ClustalwCommandline(infile=str(combined_fasta), outfile=str(aligned_output))
    stdout, stderr = clustalw_cline()
    
    if stderr:
        print(f"Error in alignment for {gene}: {stderr}")
    else:
        print(f"Alignment completed for {gene}")

print("Alignment completed.")
