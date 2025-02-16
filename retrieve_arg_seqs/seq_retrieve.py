import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import os
import logging
import argparse

# Configure logging
logging.basicConfig(filename='sequence_extraction.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def extract_sequences(excel_file, fasta_dir, nucleotide_output, protein_output):
    # Load the Excel file
    try:
        data = pd.read_excel(excel_file)
    except Exception as e:
        logging.error(f"Error reading Excel file: {e}")
        return

    # Open output FASTA files
    try:
        with open(nucleotide_output, "w") as nuc_out, open(protein_output, "w") as prot_out:
            # Process each line in the Excel file
            for index, row in data.iterrows():
                # Extract details from the row
                file_name = row.get("#FILE")
                sequence_id = row.get("SEQUENCE")
                start = int(row.get("START"))
                end = int(row.get("END"))
                strand = row.get("STRAND")
                gene = row.get("GENE")

                # Check if necessary columns are present
                if pd.isna([file_name, sequence_id, start, end, strand, gene]).any():
                    logging.warning(f"Missing data in row {index}. Skipping.")
                    continue

                fasta_file = os.path.join(fasta_dir, f"{file_name}.fasta")  # Dynamically select the FASTA file

                # Check if the FASTA file exists
                if not os.path.exists(fasta_file):
                    logging.warning(f"FASTA file {fasta_file} not found. Skipping.")
                    continue

                # Parse the FASTA file and search for the sequence
                fasta_sequences = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}

                if sequence_id in fasta_sequences:
                    genome_seq = fasta_sequences[sequence_id].seq
                    extracted_seq = genome_seq[start-1:end]  # Adjust for 1-based indexing

                    # Reverse complement if the strand is "-"
                    if strand == "-":
                        extracted_seq = extracted_seq.reverse_complement()

                    # Construct the FASTA header
                    fasta_header = f">{file_name}_{gene}"

                    # Write nucleotide sequence to the nucleotide FASTA file
                    nuc_out.write(f"{fasta_header}\n{extracted_seq}\n")

                    # Translate to protein and write to the protein FASTA file
                    protein_seq = extracted_seq.translate()
                    prot_out.write(f"{fasta_header}\n{protein_seq}\n")
                else:
                    logging.warning(f"Sequence ID {sequence_id} not found in {fasta_file}.")
    except Exception as e:
        logging.error(f"Error writing output files: {e}")

def main():
    parser = argparse.ArgumentParser(description="Extract nucleotide and protein sequences from FASTA files based on an Excel file.")
    parser.add_argument("--excel", required=True, help="Path to the input Excel file.")
    parser.add_argument("--fasta-dir", required=True, help="Directory containing FASTA files.")
    parser.add_argument("--nucleotide-output", required=True, help="Output file for nucleotide sequences.")
    parser.add_argument("--protein-output", required=True, help="Output file for protein sequences.")

    args = parser.parse_args()

    extract_sequences(args.excel, args.fasta_dir, args.nucleotide_output, args.protein_output)

if __name__ == "__main__":
    main()
