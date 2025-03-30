# Genome Analysis Tools

This repository contains a collection of bioinformatics tools and scripts for genome analysis after assembly/annotation. The goal is to provide a combination of workflows and analysis that were useful to me in the past, as well as some things I want to improve at.

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
python genome_statistics/genome_statistics.py --input-dir path/to/fasta-files \
                                              --output-file path/to/output.csv
```

**Retrieve Antibiotic Resistance Genes (ARG) sequences based on the output file from ABRicate**
```bash
python retrieve_arg_seqs/seq_retrieve.py --excel path/to/excel-file \
                                         --fasta-dir path/to/genome-containing-directory \
                                         --nucleotide-output path/to/nucleotide_sequences.fasta \
                                         --protein-output path/to/protein_sequences.fasta
```
