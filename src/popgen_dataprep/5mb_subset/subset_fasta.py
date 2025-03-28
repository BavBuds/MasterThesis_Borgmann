import os
import logging
import time
from Bio import SeqIO
from pathlib import Path

# Define paths
OUTPUT_DIR = "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/dataprep/5mb_subset"

# Define input FASTA files
input_files = {
    "HT": "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Hexaplex_ONT_FLYE/assembly.fasta",
    "BB": "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/assembly__captus-asm/01_assembly/assembly.fasta",
    "HT2": "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Hexaplex_Captus/02_assemblies/SRR28865916__captus-asm/01_assembly/assembly.fasta"
}

# Set up logging
def setup_logging():
    """Configure logging to output to both console and a file."""
    log_dir = os.path.join(OUTPUT_DIR, "logs")
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"subset_fasta_{timestamp}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return log_file

def get_fasta_stats(fasta_file):
    """Get statistics for a FASTA file."""
    if not os.path.exists(fasta_file):
        return None
    
    contigs = list(SeqIO.parse(fasta_file, "fasta"))
    total_length = sum(len(contig.seq) for contig in contigs)
    contig_lengths = [len(contig.seq) for contig in contigs]
    
    file_size = os.path.getsize(fasta_file)
    human_file_size = f"{file_size / (1024*1024):.2f} MB"
    
    stats = {
        "file_path": fasta_file,
        "file_size_bytes": file_size,
        "file_size": human_file_size,
        "contig_count": len(contigs),
        "total_bases": total_length,
        "longest_contig": max(contig_lengths) if contig_lengths else 0,
        "shortest_contig": min(contig_lengths) if contig_lengths else 0,
        "avg_contig_length": int(sum(contig_lengths) / len(contig_lengths)) if contig_lengths else 0,
        "n50": calculate_n50(contig_lengths) if contig_lengths else 0
    }
    return stats

def calculate_n50(contig_lengths):
    """Calculate N50 value for a list of contig lengths."""
    sorted_lengths = sorted(contig_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    running_sum = 0
    for length in sorted_lengths:
        running_sum += length
        if running_sum >= total_length / 2:
            return length
    return 0

def subset_fasta(input_fasta, output_fasta, max_size=5000000):
    """Subset the largest contigs until reaching at least 5MB (even if slightly over)."""
    # Check if the output file already exists
    if os.path.exists(output_fasta):
        logging.info(f"✅ Skipping {output_fasta}, file already exists.")
        return get_fasta_stats(output_fasta)
    
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    
    # Get input file stats
    input_stats = get_fasta_stats(input_fasta)
    logging.info(f"Processing {input_fasta}")
    logging.info(f"  - Input contains {input_stats['contig_count']} contigs with {input_stats['total_bases']:,} bases")
    
    contigs = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Sort contigs by length (longest first)
    contigs.sort(key=lambda x: len(x.seq), reverse=True)
    
    total_length = 0
    selected_contigs = []
    
    # Select contigs until we reach the size limit
    for contig in contigs:
        if total_length + len(contig.seq) > max_size:
            # If we're still under 5MB, add one more contig to ensure we go over
            if total_length < max_size:
                selected_contigs.append(contig)
                total_length += len(contig.seq)
                logging.debug(f"Added contig {contig.id} to reach threshold ({len(contig.seq):,} bp)")
            break
        selected_contigs.append(contig)
        total_length += len(contig.seq)
        logging.debug(f"Added contig {contig.id} ({len(contig.seq):,} bp)")
    
    # Write to new FASTA file
    SeqIO.write(selected_contigs, output_fasta, "fasta")
    
    # Get output file stats
    output_stats = get_fasta_stats(output_fasta)
    
    logging.info(f"✅ Created {output_fasta}")
    logging.info(f"  - Selected {len(selected_contigs)} contigs ({total_length:,} bases)")
    
    return output_stats

def print_stats_table(stats_collection):
    """Print a formatted table of statistics for all processed files."""
    if not stats_collection:
        return
    
    logging.info("\n" + "="*80)
    logging.info("SUMMARY OF PROCESSED FILES")
    logging.info("="*80)
    
    headers = ["Sample", "File Size", "Contigs", "Total Bases", "Longest", "N50"]
    row_format = "{:<6} {:<12} {:<8} {:<15} {:<12} {:<12}"
    
    logging.info(row_format.format(*headers))
    logging.info("-"*80)
    
    for sample, stats in stats_collection.items():
        if stats:
            row = [
                sample,
                stats["file_size"],
                f"{stats['contig_count']}",
                f"{stats['total_bases']:,}",
                f"{stats['longest_contig']:,}",
                f"{stats['n50']:,}"
            ]
            logging.info(row_format.format(*row))
    
    logging.info("="*80)

def main():
    """Main function to process all input files."""
    log_file = setup_logging()
    logging.info(f"Starting FASTA subset process. Log file: {log_file}")
    
    # Process only missing files
    output_stats = {}
    
    for sample, input_fasta in input_files.items():
        output_fasta = os.path.join(OUTPUT_DIR, f"{sample}_assembly_5mb_subset.fasta")
        logging.info(f"\nProcessing sample: {sample}")
        
        try:
            stats = subset_fasta(input_fasta, output_fasta)
            output_stats[sample] = stats
        except Exception as e:
            logging.error(f"Error processing {sample}: {str(e)}")
            output_stats[sample] = None
    
    # Print summary table
    print_stats_table(output_stats)
    
    logging.info(f"\nProcess completed. Log file saved to: {log_file}")

if __name__ == "__main__":
    main()