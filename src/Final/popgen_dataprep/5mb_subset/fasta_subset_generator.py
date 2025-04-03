#!/usr/bin/env python3
"""
FASTA Subset Generator

This script extracts subsets of large genome assemblies by selecting the largest contigs
until reaching a specified size threshold. It creates smaller, more manageable FASTA
files while preserving the most significant contigs from the original assemblies.

The tool calculates and reports key assembly statistics such as:
- N50 values
- Contig counts and lengths
- Total sequence size
- File sizes

Author: Max Borgmann
Date: March 2025
"""

import os
import sys
import logging
import time
import argparse
import statistics
from Bio import SeqIO
from pathlib import Path


class FastaSubsetGenerator:
    """
    A tool for creating size-limited subsets of genome assemblies.
    
    This class implements methods to select the largest contigs from FASTA files
    up to a specified size threshold, generating smaller files for testing or
    focused analysis while preserving important genomic regions.
    
    Attributes:
        args: Command-line arguments
        logger: Configured logger instance
        input_files: Dictionary mapping sample IDs to input file paths
    """
    
    def __init__(self):
        """Initialize the FASTA subset generator with command-line arguments."""
        self.args = self._parse_arguments()
        self.logger = self._setup_logging()
        
        # Create main output directory if it doesn't exist
        os.makedirs(self.args.output_dir, exist_ok=True)
        
        # Create a folder for sample FASTAs inside the output directory
        self.samples_output_dir = os.path.join(self.args.output_dir, "Input_fastas")
        os.makedirs(self.samples_output_dir, exist_ok=True)
        
        # Initialize input files from command line or use defaults
        self.input_files = self._initialize_input_files()
        
        self.logger.info("FASTA Subset Generator initialized")
        self.logger.info(f"Output directory: {self.args.output_dir}")
        self.logger.info(f"Sample FASTA files will be saved in: {self.samples_output_dir}")
        self.logger.info(f"Size threshold: {self.args.size_threshold:,} bp")
        
    def _parse_arguments(self):
        """Parse command-line arguments for the FASTA subset generator."""
        parser = argparse.ArgumentParser(
            description="FASTA Subset Generator - Create smaller FASTA files from large assemblies",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        
        parser.add_argument("-o", "--output_dir", required=True,
                            help="Output directory for subset FASTA files")
        parser.add_argument("-s", "--size_threshold", type=int, default=5000000,
                            help="Size threshold in base pairs (bp)")
        parser.add_argument("-f", "--force", action="store_true",
                            help="Force overwrite existing output files")
        parser.add_argument("-c", "--config",
                            help="Optional JSON config file with input files")
        parser.add_argument("-i", "--inputs", nargs="+",
                            help="Input files in format 'SAMPLE_ID:PATH_TO_FASTA'")
        parser.add_argument("-q", "--quiet", action="store_true",
                            help="Reduce console output verbosity")
        
        return parser.parse_args()
        
    def _setup_logging(self):
        """Configure logging to output to both console and a file."""
        log_dir = os.path.join(self.args.output_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)
        
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        log_file = os.path.join(log_dir, f"subset_fasta_{timestamp}.log")
        
        # Set up logger
        logger = logging.getLogger('fasta_subset')
        logger.setLevel(logging.DEBUG)
        
        # Clear any existing handlers
        if logger.hasHandlers():
            logger.handlers.clear()
        
        # File handler always logs everything
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
        
        # Console handler level depends on quiet flag
        console_handler = logging.StreamHandler()
        console_level = logging.WARNING if self.args.quiet else logging.INFO
        console_handler.setLevel(console_level)
        console_formatter = logging.Formatter("[%(levelname)s] %(message)s")
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        
        logger.info(f"Logging initialized. Log file: {log_file}")
        return logger
        
    def _initialize_input_files(self):
        """
        Initialize input files from arguments, config file, or use defaults.
        
        Returns:
            dict: Dictionary mapping sample IDs to input file paths
        """
        # Default input files if nothing specified
        default_files = {
            "HT": "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Hexaplex_ONT_FLYE/assembly.fasta",
            "BB": "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/assembly__captus-asm/01_assembly/assembly.fasta",
            "HT2": "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Hexaplex_Captus/02_assemblies/SRR28865916__captus-asm/01_assembly/assembly.fasta"
        }
        
        # If custom inputs are provided, parse them
        input_files = {}
        if self.args.inputs:
            for input_spec in self.args.inputs:
                if ":" in input_spec:
                    sample_id, file_path = input_spec.split(":", 1)
                    input_files[sample_id] = file_path
                else:
                    self.logger.warning(f"Ignoring malformed input specification: {input_spec}")
                    self.logger.warning("Format should be 'SAMPLE_ID:PATH_TO_FASTA'")
            
            # If no valid inputs found, fallback to defaults
            if not input_files:
                self.logger.warning("No valid input files specified, using defaults")
                input_files = default_files
        else:
            self.logger.info("Using default input files")
            input_files = default_files
            
        # TODO: Add JSON config file parsing if self.args.config is provided
        
        # Verify that input files exist
        for sample_id, file_path in input_files.items():
            if not os.path.exists(file_path):
                self.logger.warning(f"Input file not found for {sample_id}: {file_path}")
            else:
                self.logger.debug(f"Found input file for {sample_id}: {file_path}")
        
        self.logger.info(f"Processing {len(input_files)} input files")
        return input_files
    
    def get_fasta_stats(self, fasta_file):
        """
        Get comprehensive statistics for a FASTA file.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            dict: Dictionary with FASTA statistics or None if file not found
        """
        if not os.path.exists(fasta_file):
            self.logger.warning(f"Cannot compute stats, file not found: {fasta_file}")
            return None
        
        self.logger.debug(f"Computing statistics for: {fasta_file}")
        
        try:
            contigs = list(SeqIO.parse(fasta_file, "fasta"))
            total_length = sum(len(contig.seq) for contig in contigs)
            contig_lengths = [len(contig.seq) for contig in contigs]
            
            file_size = os.path.getsize(fasta_file)
            human_file_size = f"{file_size / (1024*1024):.2f} MB"
            
            mean_gc, stdev_gc = self.calculate_gc_content(fasta_file)
            
            stats = {
                "file_path": fasta_file,
                "file_size_bytes": file_size,
                "file_size": human_file_size,
                "contig_count": len(contigs),
                "total_bases": total_length,
                "longest_contig": max(contig_lengths) if contig_lengths else 0,
                "shortest_contig": min(contig_lengths) if contig_lengths else 0,
                "avg_contig_length": int(sum(contig_lengths) / len(contig_lengths)) if contig_lengths else 0,
                "n50": self.calculate_n50(contig_lengths) if contig_lengths else 0,
                "gc_mean": mean_gc,
                "gc_stdev": stdev_gc
            }
            
            self.logger.debug(f"Statistics computed: {stats['contig_count']} contigs, {stats['total_bases']:,} bases")
            return stats
            
        except Exception as e:
            self.logger.error(f"Error computing statistics for {fasta_file}: {str(e)}")
            return None

    def calculate_gc_content(self, fasta_file):
        """Calculate mean and standard deviation of GC content across contigs."""
        gc_contents = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq).upper()
            gc_count = seq.count("G") + seq.count("C")
            if len(seq) > 0:
                gc_contents.append((gc_count / len(seq)) * 100)
        
        if gc_contents:
            mean_gc = statistics.mean(gc_contents)
            stdev_gc = statistics.stdev(gc_contents) if len(gc_contents) > 1 else 0.0
        else:
            mean_gc = 0.0
            stdev_gc = 0.0
        
        return round(mean_gc, 2), round(stdev_gc, 2)

    def calculate_n50(self, contig_lengths):
        """
        Calculate N50 value for a list of contig lengths.
        
        Args:
            contig_lengths: List of contig lengths
            
        Returns:
            int: N50 value
        """
        sorted_lengths = sorted(contig_lengths, reverse=True)
        total_length = sum(sorted_lengths)
        running_sum = 0
        
        for length in sorted_lengths:
            running_sum += length
            if running_sum >= total_length / 2:
                return length
                
        return 0
    def subset_fasta(self, input_fasta, output_fasta, max_size=None):
        """
        Subset the largest contigs until reaching the size threshold.
        
        Args:
            input_fasta: Path to input FASTA file
            output_fasta: Path to output FASTA file
            max_size: Maximum size in base pairs (default: from args)
            
        Returns:
            dict: Dictionary with statistics for output FASTA or None if failed
        """
        if max_size is None:
            max_size = self.args.size_threshold
            
        if os.path.exists(output_fasta) and not self.args.force:
            self.logger.info(f"✅ Skipping {output_fasta}, file already exists. Use --force to overwrite.")
            return self.get_fasta_stats(output_fasta)
        
        os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
        
        input_stats = self.get_fasta_stats(input_fasta)
        if not input_stats:
            self.logger.error(f"Could not read input file: {input_fasta}")
            return None
            
        self.logger.info(f"Processing {input_fasta}")
        self.logger.info(f"  - Input contains {input_stats['contig_count']} contigs with {input_stats['total_bases']:,} bases")
        
        try:
            contigs = list(SeqIO.parse(input_fasta, "fasta"))
            contigs.sort(key=lambda x: len(x.seq), reverse=True)
            
            total_length = 0
            selected_contigs = []
            
            for contig in contigs:
                if total_length + len(contig.seq) > max_size:
                    if total_length < max_size:
                        selected_contigs.append(contig)
                        total_length += len(contig.seq)
                        self.logger.debug(f"Added contig {contig.id} to reach threshold ({len(contig.seq):,} bp)")
                    break
                
                selected_contigs.append(contig)
                total_length += len(contig.seq)
                self.logger.debug(f"Added contig {contig.id} ({len(contig.seq):,} bp)")
            
            SeqIO.write(selected_contigs, output_fasta, "fasta")
            output_stats = self.get_fasta_stats(output_fasta)
            
            self.logger.info(f"✅ Created {output_fasta}")
            self.logger.info(f"  - Selected {len(selected_contigs)} contigs ({total_length:,} bases)")
            
            return output_stats
            
        except Exception as e:
            self.logger.error(f"Error subsetting FASTA file {input_fasta}: {str(e)}")
            return None


    def print_stats_table(self, original_stats, subset_stats):
        """
        Print a formatted comparison table of original and subset FASTA statistics.

        Args:
            original_stats: Dict of original full assembly stats
            subset_stats: Dict of subset FASTA stats
        """
        self.logger.info("\n" + "="*140)
        self.logger.info("COMPARISON OF ORIGINAL VS SUBSET ASSEMBLIES")
        self.logger.info("="*140)

        headers = [
            "Sample", "Contigs (orig)", "Contigs (sub)",
            "N50 (orig)", "N50 (sub)", "GC% (orig)", "GC% (sub)",
            "Bases (orig)", "Bases (sub)"
        ]
        row_fmt = "{:<8} {:<15} {:<14} {:<12} {:<10} {:<12} {:<12} {:<14} {:<14}"
        self.logger.info(row_fmt.format(*headers))
        self.logger.info("-"*140)

        for sample in original_stats:
            orig = original_stats.get(sample)
            sub = subset_stats.get(sample)

            if not orig or not sub:
                continue

            row = [
                sample,
                f"{orig['contig_count']}",
                f"{sub['contig_count']}",
                f"{orig['n50']:,}",
                f"{sub['n50']:,}",
                f"{orig['gc_mean']}±{orig['gc_stdev']}",
                f"{sub['gc_mean']}±{sub['gc_stdev']}",
                f"{orig['total_bases']:,}",
                f"{sub['total_bases']:,}"
            ]
            self.logger.info(row_fmt.format(*row))

        self.logger.info("="*140)


    def run(self):
        """
        Run the FASTA subset generator for all input files and compare to original assemblies.
        
        Returns:
            int: 0 if successful, 1 if any errors occurred
        """
        start_time = time.time()
        self.logger.info("Starting FASTA subset generation process")
        
        original_stats = {}
        subset_stats = {}
        error_count = 0
        
        for sample, input_fasta in self.input_files.items():
            self.logger.info(f"\nProcessing sample: {sample}")
            sample_output_dir = os.path.join(self.samples_output_dir, sample)
            os.makedirs(sample_output_dir, exist_ok=True)
            output_fasta = os.path.join(sample_output_dir, f"{sample}_assembly_5mb_subset.fasta")

            try:
                # Stats for original
                original_stats[sample] = self.get_fasta_stats(input_fasta)

                # Generate subset + stats
                subset_stats[sample] = self.subset_fasta(input_fasta, output_fasta)

                if not subset_stats[sample]:
                    error_count += 1

            except Exception as e:
                self.logger.error(f"Unhandled error processing {sample}: {str(e)}")
                import traceback
                self.logger.debug(traceback.format_exc())
                original_stats[sample] = None
                subset_stats[sample] = None
                error_count += 1

        self.print_stats_table(original_stats, subset_stats)
        duration = time.time() - start_time
        self.logger.info(f"\nProcess completed in {duration:.2f} seconds")
        self.logger.info(f"Processed {len(self.input_files)} files with {error_count} errors")

        return 0 if error_count == 0 else 1



if __name__ == "__main__":
    try:
        generator = FastaSubsetGenerator()
        exit_code = generator.run()
        sys.exit(exit_code)
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)