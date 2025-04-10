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
- GC content statistics

Advanced filtering capabilities include:
- Size-ranked selection (largest contigs first)
- GC content filtering (retaining contigs within ±2.5 SD of the assembly mean)
- Coverage filtering (excluding contigs with outlier coverage)

Author: Max Borgmann
Date: March 2025
"""

import os
import sys
import logging
import time
import argparse
import statistics
import json
from Bio import SeqIO


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
        coverage_files: Dictionary mapping sample IDs to coverage file paths
        is_long_read: Dictionary mapping sample IDs to boolean indicating long read data
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
        self.logger.info(f"GC content difference threshold: {self.args.gc_threshold}%")
        
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
        parser.add_argument("-g", "--gc_threshold", type=float, default=2.0,
                            help="Maximum GC%% difference threshold")
        parser.add_argument("-f", "--force", action="store_true",
                            help="Force overwrite existing output files")
        parser.add_argument("-c", "--config",
                            help="Optional JSON config file with input files")
        parser.add_argument("-i", "--inputs", nargs="+",
                            help="Input files in format 'SAMPLE_ID:PATH_TO_FASTA'")
        parser.add_argument("--coverage_files", nargs="+",
                            help="Coverage files in format 'SAMPLE_ID:PATH_TO_COVERAGE'")
        parser.add_argument("--long_read", nargs="+",
                            help="Specify which samples use long read technology (higher coverage threshold)")
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
        
        # Check if config file is provided
        if self.args.config and os.path.exists(self.args.config):
            try:
                with open(self.args.config, 'r') as f:
                    config = json.load(f)
                
                # Extract input files from config
                input_files = config.get('input_files', {})
                
                # Extract coverage files and long read flags if present
                self.coverage_files = config.get('coverage_files', {})
                self.is_long_read = config.get('is_long_read', {})
                
                self.logger.info(f"Loaded configuration from {self.args.config}")
                
            except Exception as e:
                self.logger.error(f"Error reading config file: {str(e)}")
                input_files = {}
        else:
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
            
        # Verify that input files exist
        for sample_id, file_path in input_files.items():
            if not os.path.exists(file_path):
                self.logger.warning(f"Input file not found for {sample_id}: {file_path}")
            else:
                self.logger.debug(f"Found input file for {sample_id}: {file_path}")
        
        # Initialize coverage files if provided via command line
        self.coverage_files = getattr(self, 'coverage_files', {})
        if self.args.coverage_files:
            for cov_spec in self.args.coverage_files:
                if ":" in cov_spec:
                    sample_id, file_path = cov_spec.split(":", 1)
                    self.coverage_files[sample_id] = file_path
                    if not os.path.exists(file_path):
                        self.logger.warning(f"Coverage file not found for {sample_id}: {file_path}")
        
        # Initialize long read flags if provided via command line
        self.is_long_read = getattr(self, 'is_long_read', {})
        if self.args.long_read:
            for sample_id in self.args.long_read:
                self.is_long_read[sample_id] = True
        
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
    
    def calculate_gc_for_contig(self, seq):
        """Calculate GC content for a single contig sequence."""
        seq = str(seq).upper()
        gc_count = seq.count("G") + seq.count("C")
        if len(seq) > 0:
            return (gc_count / len(seq)) * 100
        return 0.0

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
    
    def compare_gc_content(self, original_stats, subset_stats, threshold=2.0):
        """
        Compare GC content between original and subset assemblies.
        
        Args:
            original_stats: Statistics of the original assembly
            subset_stats: Statistics of the subset assembly
            threshold: Maximum acceptable GC% difference (default: 2.0%)
            
        Returns:
            bool: True if GC difference is within threshold, False otherwise
        """
        gc_diff = abs(original_stats['gc_mean'] - subset_stats['gc_mean'])
        
        if gc_diff > threshold:
            self.logger.warning(f"GC content differs by {gc_diff:.2f}% (original: {original_stats['gc_mean']}%, "
                               f"subset: {subset_stats['gc_mean']}%), exceeding threshold of {threshold}%")
            return False
        else:
            self.logger.info(f"GC content difference is {gc_diff:.2f}% (within threshold of {threshold}%)")
            return True
    
    def filter_by_gc_content(self, contigs, assembly_mean_gc, assembly_stdev_gc, sd_threshold=2.5):
        """
        Filter contigs by GC content, retaining only those within ±SD threshold
        of the assembly mean.
        
        Args:
            contigs: List of SeqIO contig records
            assembly_mean_gc: Mean GC% of the assembly
            assembly_stdev_gc: Standard deviation of GC% in the assembly
            sd_threshold: Number of standard deviations to use as threshold
            
        Returns:
            list: Filtered list of contigs
        """
        filtered_contigs = []
        lower_bound = assembly_mean_gc - (sd_threshold * assembly_stdev_gc)
        upper_bound = assembly_mean_gc + (sd_threshold * assembly_stdev_gc)
        
        self.logger.info(f"Filtering contigs by GC content: {lower_bound:.2f}% to {upper_bound:.2f}%")
        
        for contig in contigs:
            gc_percent = self.calculate_gc_for_contig(contig.seq)
            
            if lower_bound <= gc_percent <= upper_bound:
                filtered_contigs.append(contig)
            else:
                self.logger.debug(f"Excluded contig {contig.id} with GC% of {gc_percent:.2f}%")
        
        self.logger.info(f"GC content filter retained {len(filtered_contigs)} of {len(contigs)} contigs")
        return filtered_contigs
    
    def filter_by_coverage(self, contigs, is_long_read=False, coverage_file=None):
        """
        Filter contigs by coverage, excluding those with outlier coverage.
        
        Args:
            contigs: List of SeqIO contig records
            is_long_read: Whether the data is from long read (LR) technology
            coverage_file: Optional path to a file with coverage data (format: contig_id\tcoverage)
            
        Returns:
            list: Filtered list of contigs
        """
        # Default coverage thresholds
        coverage_threshold = 500 if is_long_read else 250
        
        # Coverage data dictionary
        coverage_data = {}
        
        # If coverage file is provided, read it
        if coverage_file and os.path.exists(coverage_file):
            try:
                with open(coverage_file, 'r') as f:
                    for line in f:
                        if line.strip():
                            parts = line.strip().split()
                            if len(parts) >= 2:
                                contig_id = parts[0]
                                coverage = float(parts[1])
                                coverage_data[contig_id] = coverage
                self.logger.info(f"Loaded coverage data for {len(coverage_data)} contigs from {coverage_file}")
            except Exception as e:
                self.logger.error(f"Error reading coverage file {coverage_file}: {str(e)}")
        else:
            # Try to extract coverage from contig headers
            # This is format-dependent and may not work for all assemblies
            self.logger.warning("No coverage file provided, attempting to extract from contig headers")
            
            for contig in contigs:
                # Example pattern: looking for "cov" or "coverage" in the header
                header = contig.description
                coverage = None
                
                # Look for common coverage annotations in headers
                if "cov=" in header:
                    try:
                        cov_part = header.split("cov=")[1].split()[0]
                        coverage = float(cov_part)
                    except ValueError:
                        pass
                elif "coverage=" in header:
                    try:
                        cov_part = header.split("coverage=")[1].split()[0]
                        coverage = float(cov_part)
                    except ValueError:
                        pass
                
                if coverage is not None:
                    coverage_data[contig.id] = coverage
            
            if not coverage_data:
                self.logger.warning("Could not extract coverage information from contig headers")
                self.logger.warning("Coverage filtering will not be applied")
                return contigs
        
        # Filter contigs by coverage
        filtered_contigs = []
        for contig in contigs:
            coverage = coverage_data.get(contig.id)
            
            # If we don't have coverage data for this contig, keep it
            if coverage is None:
                filtered_contigs.append(contig)
                continue
            
            if coverage <= coverage_threshold:
                filtered_contigs.append(contig)
            else:
                self.logger.debug(f"Excluded contig {contig.id} with coverage {coverage:.2f}×")
        
        self.logger.info(f"Coverage filter retained {len(filtered_contigs)} of {len(contigs)} contigs")
        return filtered_contigs

    def subset_fasta(self, input_fasta, output_fasta, max_size=None, gc_threshold=2.0, 
                      apply_gc_filter=False, apply_cov_filter=False, coverage_file=None, 
                      is_long_read=False, iteration_history=None):
        """
        Subset the largest contigs until reaching the size threshold with optional GC and coverage filtering.
        
        Args:
            input_fasta: Path to input FASTA file
            output_fasta: Path to output FASTA file
            max_size: Maximum size in base pairs (default: from args)
            gc_threshold: Maximum acceptable GC% difference between input and output
            apply_gc_filter: Whether to apply GC content filtering
            apply_cov_filter: Whether to apply coverage filtering
            coverage_file: Optional path to a file with coverage data
            is_long_read: Whether the data is from long read (LR) technology
            iteration_history: List to track filtering iterations
            
        Returns:
            tuple: (output_stats, iteration_history) with statistics and filtering history
        """
        if max_size is None:
            max_size = self.args.size_threshold
            
        if iteration_history is None:
            iteration_history = []
            
        # Track applied filters in this iteration
        current_filters = []
        if apply_gc_filter:
            current_filters.append("GC Content")
        if apply_cov_filter:
            current_filters.append("Coverage")
        if not current_filters:
            current_filters.append("Size only")
            
        if os.path.exists(output_fasta) and not self.args.force:
            self.logger.info(f"✅ Skipping {output_fasta}, file already exists. Use --force to overwrite.")
            output_stats = self.get_fasta_stats(output_fasta)
            # Add empty iteration history for existing files
            if not iteration_history:
                iteration_history.append({
                    "filters": ["Size only (existing file)"],
                    "stats": output_stats
                })
            return output_stats, iteration_history
        
        os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
        
        input_stats = self.get_fasta_stats(input_fasta)
        if not input_stats:
            self.logger.error(f"Could not read input file: {input_fasta}")
            return None, iteration_history
            
        self.logger.info(f"Processing {input_fasta}")
        self.logger.info(f"  - Input contains {input_stats['contig_count']} contigs with {input_stats['total_bases']:,} bases")
        self.logger.info(f"  - Using filters: {', '.join(current_filters)}")
        
        try:
            # Load all contigs from the input file
            contigs = list(SeqIO.parse(input_fasta, "fasta"))
            
            # Apply GC filtering if requested
            if apply_gc_filter:
                contigs = self.filter_by_gc_content(
                    contigs, 
                    input_stats['gc_mean'], 
                    input_stats['gc_stdev']
                )
                
            # Apply coverage filtering if requested
            if apply_cov_filter:
                contigs = self.filter_by_coverage(
                    contigs,
                    is_long_read=is_long_read,
                    coverage_file=coverage_file
                )
            
            # Sort contigs by length for size-based selection
            contigs.sort(key=lambda x: len(x.seq), reverse=True)
            
            # Select largest contigs up to size threshold
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
            
            # Write the output FASTA
            SeqIO.write(selected_contigs, output_fasta, "fasta")
            output_stats = self.get_fasta_stats(output_fasta)
            
            self.logger.info(f"✅ Created {output_fasta}")
            self.logger.info(f"  - Selected {len(selected_contigs)} contigs ({total_length:,} bases)")
            
            # Add this iteration to history
            iteration_history.append({
                "filters": current_filters,
                "stats": output_stats
            })
            
            # Check if GC content differs significantly from the original
            if not apply_gc_filter and not apply_cov_filter:
                if not self.compare_gc_content(input_stats, output_stats, threshold=gc_threshold):
                    self.logger.warning("Significant GC content difference detected. Refiltering with additional criteria...")
                    
                    # Remove the current output file
                    os.remove(output_fasta)
                    
                    # Recursively call with additional filters
                    return self.subset_fasta(
                        input_fasta=input_fasta,
                        output_fasta=output_fasta,
                        max_size=max_size,
                        gc_threshold=gc_threshold,
                        apply_gc_filter=True,
                        apply_cov_filter=apply_cov_filter if coverage_file else False,
                        coverage_file=coverage_file,
                        is_long_read=is_long_read,
                        iteration_history=iteration_history
                    )
            
            return output_stats, iteration_history
            
        except Exception as e:
            self.logger.error(f"Error subsetting FASTA file {input_fasta}: {str(e)}")
            import traceback
            self.logger.debug(traceback.format_exc())
            return None, iteration_history

    def write_gc_report(self, original_stats, subset_stats, filtering_history):
        """
        Write a detailed GC content report to a text file.
        
        Args:
            original_stats: Dictionary of original assembly statistics
            subset_stats: Dictionary of subset assembly statistics
            filtering_history: Dictionary tracking filtering iterations by sample
        """
        report_path = os.path.join(self.args.output_dir, "gc_content_report.txt")
        
        with open(report_path, "w") as f:
            f.write("=" * 80 + "\n")
            f.write("GC CONTENT REPORT - FASTA SUBSET GENERATOR\n")
            f.write(f"Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 80 + "\n\n")
            
            for sample in original_stats:
                orig = original_stats.get(sample)
                sub = subset_stats.get(sample)
                history = filtering_history.get(sample, [])
                
                if not orig or not sub:
                    continue
                    
                f.write(f"SAMPLE: {sample}\n")
                f.write("-" * 40 + "\n")
                
                # Original stats
                f.write("Original Assembly:\n")
                f.write(f"  - File: {orig['file_path']}\n")
                f.write(f"  - Contigs: {orig['contig_count']}\n")
                f.write(f"  - Total Size: {orig['total_bases']:,} bp\n")
                f.write(f"  - GC Content: {orig['gc_mean']}% ± {orig['gc_stdev']}%\n\n")
                
                # Final subset stats
                f.write("Final Subset Assembly:\n")
                f.write(f"  - File: {sub['file_path']}\n")
                f.write(f"  - Contigs: {sub['contig_count']}\n")
                f.write(f"  - Total Size: {sub['total_bases']:,} bp\n")
                f.write(f"  - GC Content: {sub['gc_mean']}% ± {sub['gc_stdev']}%\n")
                f.write(f"  - GC Difference: {abs(orig['gc_mean'] - sub['gc_mean']):.2f}%\n\n")
                
                # Filtering iterations if any
                if history:
                    f.write("Filtering History:\n")
                    for i, iteration in enumerate(history, 1):
                        filters = iteration.get("filters", [])
                        stats = iteration.get("stats", {})
                        
                        if not stats:
                            continue
                            
                        filter_str = ", ".join(filters) if filters else "Size only"
                        f.write(f"  Iteration {i} ({filter_str}):\n")
                        f.write(f"    - Contigs: {stats.get('contig_count', 'N/A')}\n")
                        f.write(f"    - Total Size: {stats.get('total_bases', 'N/A'):,} bp\n")
                        f.write(f"    - GC Content: {stats.get('gc_mean', 'N/A')}% ± {stats.get('gc_stdev', 'N/A')}%\n")
                        f.write(f"    - GC Difference: {abs(orig['gc_mean'] - stats.get('gc_mean', 0)):.2f}%\n\n")
                else:
                    f.write("No additional filtering iterations were required.\n\n")
                
                f.write("\n" + "=" * 80 + "\n\n")
        
        self.logger.info(f"✅ GC content report written to: {report_path}")


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
            "Bases (orig)", "Bases (sub)", "GC Diff", "Filters Applied"
        ]
        row_fmt = "{:<8} {:<15} {:<14} {:<12} {:<10} {:<12} {:<12} {:<14} {:<14} {:<8} {:<15}"
        self.logger.info(row_fmt.format(*headers))
        self.logger.info("-"*140)

        for sample in original_stats:
            orig = original_stats.get(sample)
            sub = subset_stats.get(sample)

            if not orig or not sub:
                continue
                
            gc_diff = abs(orig['gc_mean'] - sub['gc_mean'])
            filters = []
            if gc_diff > self.args.gc_threshold:
                filters.append("GC")
            if self.coverage_files.get(sample):
                filters.append("COV")
            
            filters_str = ",".join(filters) if filters else "None"

            row = [
                sample,
                f"{orig['contig_count']}",
                f"{sub['contig_count']}",
                f"{orig['n50']:,}",
                f"{sub['n50']:,}",
                f"{orig['gc_mean']}±{orig['gc_stdev']}",
                f"{sub['gc_mean']}±{sub['gc_stdev']}",
                f"{orig['total_bases']:,}",
                f"{sub['total_bases']:,}",
                f"{gc_diff:.2f}%",
                filters_str
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
        filtering_history = {}
        error_count = 0
        
        for sample, input_fasta in self.input_files.items():
            self.logger.info(f"\nProcessing sample: {sample}")
            sample_output_dir = os.path.join(self.samples_output_dir, sample)
            os.makedirs(sample_output_dir, exist_ok=True)
            output_fasta = os.path.join(sample_output_dir, f"{sample}_assembly_5mb_subset.fasta")

            try:
                # Get coverage file if available
                coverage_file = self.coverage_files.get(sample)
                is_long_read = self.is_long_read.get(sample, False)
                
                # Stats for original
                original_stats[sample] = self.get_fasta_stats(input_fasta)

                # Generate subset + stats with history tracking
                subset_stats_result, history = self.subset_fasta(
                    input_fasta=input_fasta, 
                    output_fasta=output_fasta,
                    gc_threshold=self.args.gc_threshold,
                    coverage_file=coverage_file,
                    is_long_read=is_long_read,
                    iteration_history=[]
                )
                
                subset_stats[sample] = subset_stats_result
                filtering_history[sample] = history

                if not subset_stats[sample]:
                    error_count += 1

            except Exception as e:
                self.logger.error(f"Unhandled error processing {sample}: {str(e)}")
                import traceback
                self.logger.debug(traceback.format_exc())
                original_stats[sample] = None
                subset_stats[sample] = None
                error_count += 1

        # Print comparison table
        self.print_stats_table(original_stats, subset_stats)
        
        # Generate GC report
        self.write_gc_report(original_stats, subset_stats, filtering_history)
        
        duration = time.time() - start_time
        self.logger.info(f"\nProcess completed in {duration:.2f} seconds")
        self.logger.info(f"Processed {len(self.input_files)} files with {error_count} errors")
        self.logger.info(f"GC content report available at: {os.path.join(self.args.output_dir, 'gc_content_report.txt')}")

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