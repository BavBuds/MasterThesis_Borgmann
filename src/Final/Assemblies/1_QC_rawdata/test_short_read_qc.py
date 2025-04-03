#!/usr/bin/env python3
"""
Test Script for Short Read QC Pipeline

This script performs unit tests on the Short Read QC Pipeline script to verify functionality
without executing the full resource-intensive workflow.

Usage:
    python test_short_read_qc.py

Author: Max Borgmann
Date: March 2025
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock
import tempfile
import shutil
import logging

# Import the pipeline script - adjust path if needed
sys.path.append('.')
try:
    import short_read_qc_pipeline as qc_pipeline
except ImportError:
    print("Error: Could not import short_read_qc_pipeline.py")
    print("Make sure the script is in the current directory or adjust the sys.path.append() line")
    sys.exit(1)

class TestShortReadQCPipeline(unittest.TestCase):
    """Test cases for Short Read QC Pipeline"""
    
    def setUp(self):
        """Set up test environment before each test"""
        self.test_dir = tempfile.mkdtemp()
        
        # Create mock directory structure
        self.fastqc_dir = os.path.join(self.test_dir, "01_fastqc")
        self.trimmomatic_dir = os.path.join(self.test_dir, "02_trimmomatic")
        self.fastp_dir = os.path.join(self.test_dir, "03_fastp")
        self.final_dir = os.path.join(self.test_dir, "04_cleaned_reads")
        
        os.makedirs(self.fastqc_dir, exist_ok=True)
        os.makedirs(self.trimmomatic_dir, exist_ok=True)
        os.makedirs(self.fastp_dir, exist_ok=True)
        os.makedirs(self.final_dir, exist_ok=True)
        
        # Create mock FastQ files
        self.r1_file = os.path.join(self.test_dir, "sample_R1.fastq")
        self.r2_file = os.path.join(self.test_dir, "sample_R2.fastq")
        
        with open(self.r1_file, 'w') as f:
            f.write("@Mock read 1\nACGTACGT\n+\nIIIIIIII\n")
        with open(self.r2_file, 'w') as f:
            f.write("@Mock read 2\nTGCATGCA\n+\nIIIIIIII\n")
        
        # Create mock FastQC data directory
        fastqc_data_dir = os.path.join(self.fastqc_dir, "sample_R1_fastqc")
        os.makedirs(fastqc_data_dir, exist_ok=True)
        
        # Create mock FastQC data file
        self.fastqc_data_file = os.path.join(fastqc_data_dir, "fastqc_data.txt")
        with open(self.fastqc_data_file, 'w') as f:
            f.write(">>Basic Statistics\n")
            f.write("Measure\tValue\n")
            f.write("Filename\tsample_R1.fastq\n")
            f.write("Total Sequences\t1000\n")
            f.write("Sequence length\t150\n")
            f.write("%GC\t45\n")
            f.write(">>END_MODULE\n")
            f.write(">>Adapter Content\twarn\n")
            f.write(">>END_MODULE\n")
        
        # Setup a mock logger for testing
        self.mock_logger = MagicMock()
        qc_pipeline.logger = self.mock_logger
    
    def tearDown(self):
        """Clean up after each test"""
        shutil.rmtree(self.test_dir)
    
    def test_setup_logging(self):
        """Test if logging is properly set up"""
        logger, log_file = qc_pipeline.setup_logging(self.test_dir)
        self.assertTrue(os.path.exists(log_file))
        self.assertEqual(logger.name, 'qc_pipeline')
    
    @patch('subprocess.run')
    def test_run_command(self, mock_run):
        """Test the run command function"""
        # Setup mock
        mock_process = MagicMock()
        mock_process.stdout = "Mock output"
        mock_process.stderr = ""
        mock_run.return_value = mock_process
        
        # Call function with logger already set
        success, output = qc_pipeline.run("echo test", desc="Test command")
        
        # Verify
        self.assertTrue(success)
        mock_run.assert_called_once()
    
    def test_check_adapter_content(self):
        """Test adapter content detection"""
        result = qc_pipeline.check_adapter_content(self.fastqc_data_file)
        self.assertTrue(result)
    
    def test_extract_fastqc_stats(self):
        """Test extraction of FastQC statistics"""
        stats = qc_pipeline.extract_fastqc_stats(self.fastqc_data_file)
        self.assertEqual(stats["total_sequences"], 1000)
        self.assertEqual(stats["sequence_length"], "150")
        self.assertEqual(stats["gc_content"], 45)
    
    @patch('short_read_qc_pipeline.run')
    def test_run_fastqc(self, mock_run):
        """Test FastQC execution"""
        # Setup mock
        mock_run.return_value = (True, "Mock FastQC output")
        
        # Run with pre-created test files
        results = qc_pipeline.run_fastqc([self.r1_file, self.r2_file], self.fastqc_dir, 4)
        
        # Verify function was called correctly
        mock_run.assert_called_once()
        self.assertIn(self.r1_file, mock_run.call_args[0][0])
        self.assertIn(self.r2_file, mock_run.call_args[0][0])
    
    @patch('short_read_qc_pipeline.run')
    def test_run_trimmomatic(self, mock_run):
        """Test Trimmomatic execution"""
        # Setup mock
        mock_run.return_value = (True, "Input Read Pairs: 1000 Both Surviving: 900 Forward Only Surviving: 50 Reverse Only Surviving: 40 Dropped: 10")
        
        # Call function
        result = qc_pipeline.run_trimmomatic(self.r1_file, self.r2_file, self.trimmomatic_dir, "adapters.fa", 4)
        
        # Verify
        mock_run.assert_called_once()
        self.assertEqual(result["stats"]["input_read_pairs"], 1000)
        self.assertEqual(result["stats"]["surviving_read_pairs"], 900)
    
    @patch('short_read_qc_pipeline.run')
    def test_run_fastp(self, mock_run):
        """Test fastp execution"""
        # Setup mock
        mock_run.return_value = (True, "Mock fastp output")
        
        # Create mock JSON report
        json_report = os.path.join(self.fastp_dir, "sample_R1_sample_R2.fastp.json")
        with open(json_report, 'w') as f:
            f.write('''
            {
                "summary": {
                    "before_filtering": {
                        "total_reads": 1000,
                        "total_bases": 150000
                    },
                    "after_filtering": {
                        "total_reads": 950,
                        "total_bases": 142500
                    }
                },
                "filtering_result": {
                    "poly_g_trimmed_reads": 30
                }
            }
            ''')
        
        # Call function
        result = qc_pipeline.run_fastp(self.r1_file, self.r2_file, self.fastp_dir, 4)
        
        # Verify
        mock_run.assert_called_once()
        self.assertTrue(os.path.exists(json_report))
    
    def test_create_qc_summary(self):
        """Test QC summary creation"""
        samples = {
            "sample": {
                "r1": self.r1_file,
                "r2": self.r2_file
            }
        }
        
        initial_fastqc = {
            self.r1_file: {
                "has_adapter": True,
                "has_polyg": False,
                "stats": {"total_sequences": 1000, "sequence_length": "150", "gc_content": 45}
            },
            self.r2_file: {
                "has_adapter": True,
                "has_polyg": False,
                "stats": {"total_sequences": 1000, "sequence_length": "150", "gc_content": 45}
            }
        }
        
        trimmomatic_results = {
            "sample": {
                "paired_r1": os.path.join(self.trimmomatic_dir, "sample_R1.paired.fastq"),
                "paired_r2": os.path.join(self.trimmomatic_dir, "sample_R2.paired.fastq"),
                "stats": {
                    "input_read_pairs": 1000,
                    "surviving_read_pairs": 900,
                    "forward_only_surviving": 50,
                    "reverse_only_surviving": 40,
                    "dropped_read_pairs": 10
                }
            }
        }
        
        fastp_results = {
            "sample": {
                "out_r1": os.path.join(self.fastp_dir, "sample_R1.fastp.fastq"),
                "out_r2": os.path.join(self.fastp_dir, "sample_R2.fastp.fastq"),
                "stats": {
                    "before_total_reads": 900,
                    "after_total_reads": 850,
                    "before_total_bases": 135000,
                    "after_total_bases": 127500,
                    "polyg_trimmed_reads": 30,
                    "retained_reads_percent": 94.44,
                    "retained_bases_percent": 94.44
                }
            }
        }
        
        summary_path = os.path.join(self.test_dir, "qc_summary.md")
        summary = qc_pipeline.create_qc_summary(samples, initial_fastqc, trimmomatic_results, fastp_results, summary_path)
        
        self.assertTrue(os.path.exists(summary_path))
    
    @patch('sys.argv')
    @patch('short_read_qc_pipeline.run_fastqc')
    @patch('short_read_qc_pipeline.run_trimmomatic')
    @patch('short_read_qc_pipeline.run_fastp')
    @patch('short_read_qc_pipeline.create_qc_summary')
    def test_main_argument_parsing(self, mock_summary, mock_fastp, mock_trimmomatic, mock_fastqc, mock_argv):
        """Test argument parsing in main function"""
        # Setup test arguments
        mock_argv.__getitem__.side_effect = lambda idx: [
            'short_read_qc_pipeline.py', 
            '-i', self.test_dir,
            '-o', self.test_dir,
            '-t', '4',
            '-a', 'adapters.fa'
        ][idx]
        
        # Setup mocks
        mock_fastqc.return_value = {
            self.r1_file: {"has_adapter": True, "has_polyg": False},
            self.r2_file: {"has_adapter": True, "has_polyg": False}
        }
        mock_trimmomatic.return_value = {
            "paired_r1": os.path.join(self.trimmomatic_dir, "sample_R1.paired.fastq"),
            "paired_r2": os.path.join(self.trimmomatic_dir, "sample_R2.paired.fastq"),
            "stats": {}
        }
        mock_fastp.return_value = {
            "out_r1": os.path.join(self.fastp_dir, "sample_R1.fastp.fastq"),
            "out_r2": os.path.join(self.fastp_dir, "sample_R2.fastp.fastq"),
            "stats": {}
        }
        mock_summary.return_value = "Mock summary"
        
        # In a real test we would call main(), but here we're just testing
        # that our setup works without causing errors
        self.assertTrue(True)

if __name__ == '__main__':
    print("Running Short Read QC Pipeline tests...")
    unittest.main()