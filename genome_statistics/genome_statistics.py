from Bio import SeqIO
import os
import csv
import numpy as np
import argparse
from typing import List, Tuple
import logging
from concurrent.futures import ProcessPoolExecutor

def calculate_genome_metrics(fasta_file: str) -> Tuple[int, float, int, int, int, int, float]:
    """
    Calculate various genome metrics from a FASTA file.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        Tuple[int, float, int, int, int, int, float]: Genome metrics.
    """
    contig_lengths = []
    gc_count = 0
    total_length = 0

    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = record.seq
            length = len(seq)
            contig_lengths.append(length)
            total_length += length
            gc_count += seq.count('G') + seq.count('C')

        if total_length > 0:
            gc_content = (gc_count / total_length) * 100
        else:
            gc_content = 0

        contig_lengths.sort(reverse=True)
        cumulative_lengths = np.cumsum(contig_lengths)
        half_genome_size = total_length / 2
        n50 = next(length for length, cum_len in zip(contig_lengths, cumulative_lengths) if cum_len >= half_genome_size)

        num_contigs = len(contig_lengths)
        longest_contig = contig_lengths[0] if contig_lengths else 0
        shortest_contig = contig_lengths[-1] if contig_lengths else 0
        average_contig_length = total_length / num_contigs if num_contigs > 0 else 0

        return total_length, gc_content, num_contigs, n50, longest_contig, shortest_contig, average_contig_length

    except Exception as e:
        logging.error(f"Error processing file {fasta_file}: {e}")
        return 0, 0, 0, 0, 0, 0, 0

def process_fasta_files_in_directory(directory: str) -> List[Tuple[str, int, float, int, int, int, int, float]]:
    """
    Process all FASTA files in a directory to calculate genome metrics.

    Args:
        directory (str): Path to the directory containing FASTA files.

    Returns:
        List[Tuple[str, int, float, int, int, int, int, float]]: List of genome metrics for each file.
    """
    results = []
    fasta_files = [f for f in os.listdir(directory) if f.endswith((".fasta", ".fa"))]

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(calculate_genome_metrics, os.path.join(directory, filename)) for filename in fasta_files]
        for future, filename in zip(futures, fasta_files):
            metrics = future.result()
            isolate_name = os.path.splitext(filename)[0]
            results.append((isolate_name, *metrics))

    return results

def save_results_to_csv(results: List[Tuple[str, int, float, int, int, int, int, float]], output_file: str) -> None:
    """
    Save the results to a CSV file.

    Args:
        results (List[Tuple[str, int, float, int, int, int, int, float]]): List of genome metrics.
        output_file (str): Path to the output CSV file.
    """
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Isolate Name', 'Genome Size', 'GC Content', 'Number of Contigs', 'N50', 'Longest Contig', 'Shortest Contig', 'Average Contig Length']
        writer = csv.writer(csvfile)
        writer.writerow(fieldnames)
        for result in results:
            writer.writerow(result)

def main():
    parser = argparse.ArgumentParser(description="Calculate genome metrics from FASTA files.")
    parser.add_argument('--input-dir', required=True, help="Directory containing FASTA files")
    parser.add_argument('--output-file', required=True, help="Output CSV file path")
    parser.add_argument('--log-file', help="Log file path")
    args = parser.parse_args()

    logging.basicConfig(filename=args.log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    results = process_fasta_files_in_directory(args.input_dir)
    save_results_to_csv(results, args.output_file)
    logging.info(f"Results saved to {args.output_file}")

if __name__ == "__main__":
    main()
