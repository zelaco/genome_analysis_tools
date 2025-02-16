# Genome Analysis Tools

This repository contains a collection of bioinformatics tools for genome analysis. I will keep adding more scripts based on the workflows I used to analyse genomes during my PhD.

## Tools

1. **Genome Statistics**: Calculate various genome metrics from FASTA files.
2. **Retrieve Antibiotic Resistance Genes(ARG) sequences**: Reteieve the gene sequences of the output result of ABRicate analysis to a FASTA file.

## Installation

Clone the repository:
```bash
git clone https://github.com/zelaco/genome_analysis_tools.git
cd genome_analysis_tools
```

## Install the dependencies:

```bash
pip install -r requirements.txt
```

## Usage

**Genome Statistics**
```bash
python genome_statistics/genome_statistics.py --input-dir path/to/fasta-files --output-file path/to/output.csv
```

**Retrieve Antibiotic Resistance Genes(ARG) sequences**
```bash
python retrieve_arg_seqs/seq_retrieve.py --excel path/to/excel-file --fasta-dir path/to/genome-containing-directory --nucleotide-output path/to/nucleotide_sequences.fasta --protein-output path/to/protein_sequences.fasta
```
