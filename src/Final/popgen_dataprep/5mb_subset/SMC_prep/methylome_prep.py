#!/usr/bin/env python3

import sys
import os
import gzip
import argparse
import subprocess
import tempfile
import io
import pysam
from collections import defaultdict

def extract_methylation_modkit(bam_file, sample_id, output_dir, min_prob=0.8, min_cov=10):
    """
    Extract methylation data from ONT BAM file using modkit
    Returns a bedGraph file path with methylation data
    """
    # Check if modkit is installed
    try:
        subprocess.run(["modkit", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        sys.stderr.write("modkit not found. Please install it with: pip install ont-modkit\n")
        return None
    
    # Create output files
    bedgraph_file = os.path.join(output_dir, f"{sample_id}_methyl.bedgraph")
    
    # Run modkit to extract methylation information
    cmd = f"modkit pileup {bam_file} {bedgraph_file} --mod m --cpg --only-tabs --min-probability {min_prob} --min-coverage {min_cov}"
    sys.stderr.write(f"Running: {cmd}\n")
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        return bedgraph_file
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running modkit: {e}\n")
        return None

def extract_methylation_pysam(bam_file, sample_id, output_dir, min_prob=0.8, min_cov=10):
    """
    Extract methylation data from ONT BAM file using pysam
    Fallback method when modkit is not available
    """
    sys.stderr.write("Using pysam to extract methylation data (slower but doesn't require modkit)\n")
    
    bedgraph_file = os.path.join(output_dir, f"{sample_id}_methyl.bedgraph")
    
    try:
        import pysam
    except ImportError:
        sys.stderr.write("pysam not found. Please install it with: pip install pysam\n")
        return None
    
    # Dictionary to store methylation counts at each position
    methylation_counts = defaultdict(lambda: {'meth': 0, 'total': 0})
    
    # Parse the BAM file and extract methylation information
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                # Skip if no methylation info
                if not read.has_tag('MM') or not read.has_tag('ML'):
                    continue
                
                # Extract methylation tags
                mod_types = read.get_tag('MM').split(';')
                mod_probs = read.get_tag('ML').split(';')
                
                # Skip if no modifications
                if not mod_types or len(mod_types) < 1:
                    continue
                
                # Process each modification
                for i, mod_string in enumerate(mod_types):
                    if not mod_string or ',' not in mod_string:
                        continue
                    
                    mod_type, *positions = mod_string.split(',')
                    if mod_type != 'm':  # Only process 5mC modifications
                        continue
                    
                    # Get probabilities for this modification type
                    if i < len(mod_probs):
                        probs = [float(p)/255 for p in mod_probs[i].split(',')]
                    else:
                        continue
                    
                    # Process each modified position
                    for j, pos_delta in enumerate(positions):
                        if j >= len(probs):
                            break
                        
                        # Calculate genomic position
                        if int(pos_delta) >= 0:
                            genomic_pos = read.reference_start + int(pos_delta)
                            chrom = read.reference_name
                            
                            # Only include if probability exceeds threshold
                            if probs[j] >= min_prob:
                                methylation_counts[(chrom, genomic_pos)]['meth'] += 1
                            methylation_counts[(chrom, genomic_pos)]['total'] += 1
    except Exception as e:
        sys.stderr.write(f"Error processing BAM file with pysam: {e}\n")
        return None
    
    # Write the bedGraph file
    try:
        with open(bedgraph_file, 'w') as f:
            f.write("track type=bedGraph name=Methylation description=Methylation\n")
            
            for (chrom, pos), counts in sorted(methylation_counts.items()):
                if counts['total'] >= min_cov:
                    meth_percent = (counts['meth'] / counts['total']) * 100
                    f.write(f"{chrom}\t{pos}\t{pos+1}\t{meth_percent:.2f}\n")
        
        return bedgraph_file
    except Exception as e:
        sys.stderr.write(f"Error writing bedGraph file: {e}\n")
        return None

def read_bedgraph(filename, min_cov=10):
    """Read methylation data from a bedGraph file"""
    methylation_data = {}
    
    # Open file (handle both plain text and gzipped files)
    if filename.endswith('.gz'):
        f = io.TextIOWrapper(gzip.open(filename, 'r'))
    else:
        f = open(filename, 'r')
    
    # Parse the file
    for line in f:
        if line.startswith('#') or line.startswith('track'):
            continue
            
        fields = line.strip().split('\t')
        if len(fields) < 4:
            continue
            
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        try:
            meth_percent = float(fields[3])
        except ValueError:
            continue
        
        # Use the start position (for CpG sites this is typically a specific position)
        pos = start
        
        # Determine methylation state
        if meth_percent > 80:
            state = 'M'  # Methylated
        elif meth_percent < 20:
            state = 'D'  # Unmethylated
        else:
            state = 'C'  # Intermediate
        
        # Store the methylation state
        methylation_data[(chrom, pos)] = state
    
    f.close()
    return methylation_data

class MaskIterator:
    def __init__(self, filename, negative=False):
        if filename.endswith(".gz"):
            self.file = io.TextIOWrapper(gzip.open(filename, "r"))
        else:
            self.file = open(filename, "r")
        self.eof = False
        self.lastPos = 1
        self.negative = negative
        self.readLine()
  
    def readLine(self):
        try:
            line = next(self.file)
            fields = line.strip().split()
            if len(fields) == 2:
                self.start = int(fields[0])
                self.end = int(fields[1])
            else:
                self.start = int(fields[1]) + 1
                self.end = int(fields[2])
        except StopIteration:
            self.eof = True
  
    def getVal(self, pos):
        if pos < self.lastPos:
            return False  # Skip positions that go backwards
            
        self.lastPos = pos
        while not self.eof and pos > self.end:
            self.readLine()
        if self.eof:
            return None
        if pos >= self.start and pos <= self.end:
            return True if not self.negative else False
        else:
            return False if not self.negative else True

class MergedMask:
    def __init__(self, mask_iterators):
        self.maskIterators = mask_iterators
  
    def getVal(self, pos):
        return all((m.getVal(pos) for m in self.maskIterators))

def convert_ont_methylation(bam_files, sample_ids, masks=None, negative_masks=None, min_cov=10, min_prob=0.8, chr_name=None):
    """Convert ONT methylation BAM files to eSMC2 input format"""
    
    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Extract methylation data from each BAM file to bedGraph format
        bedgraph_files = []
        valid_sample_ids = []
        
        for i, bam_file in enumerate(bam_files):
            sample_id = sample_ids[i]
            sys.stderr.write(f"Processing ONT methylation BAM for sample {sample_id}: {bam_file}\n")
            
            # Try using modkit first, fall back to pysam if needed
            bedgraph_file = extract_methylation_modkit(bam_file, sample_id, temp_dir, min_prob, min_cov)
            if not bedgraph_file:
                bedgraph_file = extract_methylation_pysam(bam_file, sample_id, temp_dir, min_prob, min_cov)
            
            if bedgraph_file and os.path.exists(bedgraph_file):
                bedgraph_files.append(bedgraph_file)
                valid_sample_ids.append(sample_id)
            else:
                sys.stderr.write(f"Failed to extract methylation data for sample {sample_id}\n")
        
        if not bedgraph_files:
            sys.stderr.write("No methylation data could be extracted from any BAM file\n")
            return
        
        # Create mask iterators
        mask_iterators = []
        if masks:
            for f in masks:
                sys.stderr.write(f"Adding mask: {f}\n")
                mask_iterators.append(MaskIterator(f))
        if negative_masks:
            for nm in negative_masks:
                sys.stderr.write(f"Adding negative mask: {nm}\n")
                mask_iterators.append(MaskIterator(nm, True))
        
        merged_mask = MergedMask(mask_iterators) if mask_iterators else None
        
        # Read methylation data from each bedGraph file
        all_methylation_data = {}
        for i, bedgraph_file in enumerate(bedgraph_files):
            sample_id = valid_sample_ids[i]
            sys.stderr.write(f"Reading methylation bedGraph for sample {sample_id}: {bedgraph_file}\n")
            
            methylation_data = read_bedgraph(bedgraph_file, min_cov)
            
            # Combine methylation data from this sample with others
            for pos, state in methylation_data.items():
                if pos not in all_methylation_data:
                    all_methylation_data[pos] = ['C'] * len(bedgraph_files)  # Default to 'C'
                all_methylation_data[pos][i] = state
        
        # Sort positions by chromosome and position
        sorted_positions = sorted(all_methylation_data.keys())
        if not sorted_positions:
            sys.stderr.write("No methylation positions found in any BAM file\n")
            return
        
        # Write the output
        prev_pos = 0
        prev_chrom = None
        
        for chrom, pos in sorted_positions:
            # Apply mask if provided
            if merged_mask and not merged_mask.getVal(pos):
                continue
                
            # Calculate distance
            if chrom != prev_chrom:
                distance = pos  # First position in chromosome
            else:
                distance = pos - prev_pos
            
            # Override chromosome name if requested
            output_chrom = chr_name if chr_name is not None else chrom
            
            # Create the methylation state string for all samples
            state_string = ''.join(all_methylation_data[(chrom, pos)])
            
            # Output the methylation data line
            print(f"{output_chrom}\t{pos}\t{distance}\t{state_string}")
            
            prev_pos = pos
            prev_chrom = chrom

def main():
    parser = argparse.ArgumentParser(description="Convert ONT methylation BAM files to eSMC2 input format")
    parser.add_argument("--bam", required=True, action="append", help="ONT methylation BAM file (can be given multiple times)")
    parser.add_argument("--sample", required=True, action="append", help="Sample ID for each BAM file (must match order of --bam arguments)")
    parser.add_argument("--mask", action="append", help="Apply mask in BED format (regions to include)")
    parser.add_argument("--negative_mask", action="append", help="Apply negative mask in BED format (regions to exclude)")
    parser.add_argument("--min_coverage", type=int, default=10, help="Minimum read coverage to include a methylation site")
    parser.add_argument("--min_probability", type=float, default=0.8, help="Minimum probability threshold for methylation calls (0-1)")
    parser.add_argument("--chr", help="Override chromosome names in output")
    
    args = parser.parse_args()
    
    # Validate input arguments
    if len(args.bam) != len(args.sample):
        sys.stderr.write("Error: Number of BAM files must match number of sample IDs\n")
        return 1
    
    convert_ont_methylation(
        args.bam,
        args.sample,
        args.mask,
        args.negative_mask,
        args.min_coverage,
        args.min_probability,
        args.chr
    )
    
    return 0

if __name__ == "__main__":
    sys.exit(main())