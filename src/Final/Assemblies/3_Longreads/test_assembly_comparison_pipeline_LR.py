#!/usr/bin/env python3
"""
Simplified Test Script for Long Read Assembly Pipeline

This script tests the core resume functionality with proper mocking
to avoid the errors seen in the previous test.

Author: Max Borgmann
Date: 03.04.2025
"""

import os
import sys
import unittest
import tempfile
import shutil
import json
from unittest.mock import patch, MagicMock

# Import the script directly 
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_path)

# Import the main assembly script
import assembly_comparison_pipeline_LR as pipeline

class TestLongReadAssemblyResume(unittest.TestCase):
    """Basic tests for the long read assembly pipeline's resume functionality."""
    
    def setUp(self):
        """Set up test environment with temporary directories."""
        self.temp_dir = tempfile.mkdtemp()
        self.input_dir = os.path.join(self.temp_dir, "input")
        self.output_dir = os.path.join(self.temp_dir, "output")
        self.asm_dir = os.path.join(self.output_dir, "assemblies")
        self.busco_dir = os.path.join(self.output_dir, "busco")
        os.makedirs(self.input_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.asm_dir, exist_ok=True)
        os.makedirs(self.busco_dir, exist_ok=True)
        
        # Create dummy input files
        self.create_dummy_files()
        
        # Set up logging to avoid errors
        pipeline.setup_logging(self.output_dir)
    
    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.temp_dir)
    
    def create_dummy_files(self):
        """Create dummy files for testing."""
        # Create a dummy long read file
        dummy_reads = os.path.join(self.input_dir, "sample1.fastq")
        with open(dummy_reads, 'w') as f:
            f.write("@read1\nACTGACTGACTG\n+\n############\n")
        
        # Create a dummy assembly
        assembly_dir = os.path.join(self.asm_dir, "sample1_flye_pacbio-hifi")
        os.makedirs(assembly_dir, exist_ok=True)
        dummy_assembly = os.path.join(assembly_dir, "assembly.fasta")
        with open(dummy_assembly, 'w') as f:
            f.write(">contig1\nACGTACGTACGTACGT\n")
        
        # Create a final assembly file
        final_assembly = os.path.join(self.asm_dir, "sample1_flye_pacbio-hifi.fasta")
        with open(final_assembly, 'w') as f:
            f.write(">contig1\nACGTACGTACGTACGT\n")
        
        # Create a dummy BUSCO directory and JSON file
        busco_dir = os.path.join(self.busco_dir, "sample1_FLYE_pacbio-hifi_mollusca_odb10")
        os.makedirs(busco_dir, exist_ok=True)
        busco_json = os.path.join(busco_dir, "short_summary.specific.mollusca_odb10.sample1_FLYE_pacbio-hifi_mollusca_odb10.json")
        busco_data = {
            "results": {
                "Complete percentage": 85.0,
                "Missing percentage": 10.0
            },
            "metrics": {
                "Contigs N50": 50000,
                "Number of contigs": 200,
                "Total length": 10000000
            }
        }
        with open(busco_json, 'w') as f:
            json.dump(busco_data, f)
        
        # Create a seqkit stats file
        seqkit_file = os.path.join(self.asm_dir, "sample1_flye_pacbio-hifi.fasta.seqkit_stats.txt")
        with open(seqkit_file, 'w') as f:
            f.write("file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tsum_gap\tN50\tQ20(%)\tQ30(%)\tgc_content\n")
            f.write("sample1_flye_pacbio-hifi.fasta\tFASTA\tDNA\t10\t1000\t100\t100\t100\t100\t100\t100\t0\t100\t100.0\t100.0\t42.5\n")
    
    def mock_run_cmd(self, cmd, cwd=None, desc=None):
        """Mock running a command."""
        # Create appropriate output files based on the command
        if "flye" in cmd:
            assembly_dir = os.path.join(self.asm_dir, "sample1_flye_pacbio-hifi")
            os.makedirs(assembly_dir, exist_ok=True)
            assembly_file = os.path.join(assembly_dir, "assembly.fasta")
            with open(assembly_file, 'w') as f:
                f.write(">contig1\nACGTACGTACGTACGT\n")
        
        elif "canu" in cmd:
            canu_dir = os.path.join(self.asm_dir, "sample1_canu_pacbio-hifi")
            os.makedirs(canu_dir, exist_ok=True)
            assembly_file = os.path.join(canu_dir, "sample1.contigs.fasta")
            with open(assembly_file, 'w') as f:
                f.write(">contig1\nACGTACGTACGTACGT\n")
        
        elif "busco" in cmd:
            # Extract lineage
            if "-l mollusca_odb10" in cmd:
                lineage = "mollusca_odb10"
            else:
                lineage = "eukaryota_odb10"
            
            # Create a BUSCO output directory and JSON file
            if "-o" in cmd:
                out_name = cmd.split("-o ")[1].split()[0]
                busco_dir = os.path.join(self.busco_dir, out_name)
                os.makedirs(busco_dir, exist_ok=True)
                busco_json = os.path.join(busco_dir, f"short_summary.specific.{lineage}.{out_name}.json")
                busco_data = {
                    "results": {
                        "Complete percentage": 85.0,
                        "Missing percentage": 10.0
                    },
                    "metrics": {
                        "Contigs N50": 50000,
                        "Number of contigs": 200,
                        "Total length": 10000000
                    }
                }
                with open(busco_json, 'w') as f:
                    json.dump(busco_data, f)
        
        elif "seqkit" in cmd:
            # Extract output file
            if "-o" in cmd:
                out_file = cmd.split("-o ")[1].strip()
                with open(out_file, 'w') as f:
                    f.write("file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tsum_gap\tN50\tQ20(%)\tQ30(%)\tgc_content\n")
                    f.write("test.fasta\tFASTA\tDNA\t10\t1000\t100\t100\t100\t100\t100\t100\t0\t100\t100.0\t100.0\t42.5\n")
        
        return True, "Mock command output"
    
    def test_state_file_operations(self):
        """Test basic state file operations."""
        # Test empty state
        state = pipeline.read_state(self.output_dir)
        self.assertEqual(state, {})
        
        # Create and write a test state
        test_state = {
            "sample1": {
                "last_update": "2025-04-03T12:00:00",
                "assemblers": {
                    "FLYE_pacbio-hifi": {
                        "assembly": {
                            "status": "completed",
                            "path": "/path/to/assembly.fasta",
                            "hash": "testhash123"
                        }
                    }
                }
            }
        }
        pipeline.write_state(test_state, self.output_dir)
        
        # Read back the state and verify
        read_state = pipeline.read_state(self.output_dir)
        self.assertEqual(read_state, test_state)
        
        # Test updating state
        pipeline.update_state(
            self.output_dir, "sample1", "FLYE", "pacbio-hifi", "gc", "completed",
            value=42.5
        )
        
        # Read updated state and verify
        updated_state = pipeline.read_state(self.output_dir)
        self.assertEqual(
            updated_state["sample1"]["assemblers"]["FLYE_pacbio-hifi"]["gc"]["status"],
            "completed"
        )
        self.assertEqual(
            updated_state["sample1"]["assemblers"]["FLYE_pacbio-hifi"]["gc"]["value"],
            42.5
        )
    
    @patch('assembly_comparison_pipeline_LR.run_cmd')
    def test_gc_content_resume(self, mock_run_cmd):
        """Test GC content calculation with resume functionality."""
        # Mock run_cmd to create a seqkit output file
        mock_run_cmd.side_effect = self.mock_run_cmd
        
        # First run - should calculate GC
        assembly_file = os.path.join(self.asm_dir, "sample1_flye_pacbio-hifi.fasta")
        gc = pipeline.get_gc_content(
            assembly_file, 
            output_dir=self.output_dir, 
            sample="sample1", 
            assembler="FLYE", 
            config="pacbio-hifi",
            force=False
        )
        
        # Verify GC value was calculated
        self.assertAlmostEqual(gc, 42.5, places=1)
        
        # Setup state to mark GC calculation as completed
        pipeline.update_state(
            self.output_dir, "sample1", "FLYE", "pacbio-hifi", "gc", "completed",
            value=42.5
        )
        
        # Reset mock to verify it's not called again
        mock_run_cmd.reset_mock()
        
        # Second run - should use cached value
        gc2 = pipeline.get_gc_content(
            assembly_file, 
            output_dir=self.output_dir, 
            sample="sample1", 
            assembler="FLYE", 
            config="pacbio-hifi",
            force=False
        )
        
        # Verify GC value is the same and run_cmd wasn't called
        self.assertAlmostEqual(gc2, 42.5, places=1)
        mock_run_cmd.assert_not_called()
    
    @patch('assembly_comparison_pipeline_LR.run_cmd')
    def test_flye_assembly_resume(self, mock_run_cmd):
        """Test FLYE assembly with resume functionality."""
        # First, mock successful command execution
        mock_run_cmd.side_effect = self.mock_run_cmd
        
        # Setup state to mark assembly as completed
        assembly_file = os.path.join(self.asm_dir, "sample1_flye_pacbio-hifi.fasta")
        pipeline.update_state(
            self.output_dir, "sample1", "FLYE", "pacbio-hifi", "assembly", "completed",
            path=assembly_file,
            hash=pipeline.get_file_hash(assembly_file)
        )
        
        # Run FLYE assembly which should use the cached result
        reads_file = os.path.join(self.input_dir, "sample1.fastq")
        assembly_file2 = pipeline.run_flye_assembly(
            reads=reads_file,
            output_dir=self.asm_dir,
            sample_name="sample1",
            read_type="pacbio-hifi",
            genome_size="50m",
            threads=4,
            state_dir=self.output_dir,
            force=False
        )
        
        # Verify assembly file is returned and run_cmd not called
        self.assertEqual(assembly_file, assembly_file2)
        mock_run_cmd.assert_not_called()
        
        # Now test with force=True to ensure it runs again
        mock_run_cmd.reset_mock()
        assembly_file3 = pipeline.run_flye_assembly(
            reads=reads_file,
            output_dir=self.asm_dir,
            sample_name="sample1",
            read_type="pacbio-hifi",
            genome_size="50m",
            threads=4,
            state_dir=self.output_dir,
            force=True
        )
        
        # Verify run_cmd was called (forced re-run)
        self.assertTrue(mock_run_cmd.called)
    
    @patch('assembly_comparison_pipeline_LR.run_cmd')
    def test_busco_resume(self, mock_run_cmd):
        """Test BUSCO analysis with resume functionality."""
        # First, mock successful command execution
        mock_run_cmd.side_effect = self.mock_run_cmd
        
        # Create a dummy assembly file
        assembly_file = os.path.join(self.asm_dir, "sample1_flye_pacbio-hifi.fasta")
        
        # Create a BUSCO output directory and JSON file
        busco_dir = os.path.join(self.busco_dir, "sample1_FLYE_pacbio-hifi_mollusca_odb10")
        os.makedirs(busco_dir, exist_ok=True)
        json_path = os.path.join(busco_dir, "short_summary.specific.mollusca_odb10.sample1_FLYE_pacbio-hifi_mollusca_odb10.json")
        busco_data = {
            "results": {
                "Complete percentage": 85.0,
                "Missing percentage": 10.0
            },
            "metrics": {
                "Contigs N50": 50000,
                "Number of contigs": 200,
                "Total length": 10000000
            }
        }
        with open(json_path, 'w') as f:
            json.dump(busco_data, f)
        
        # Setup state to mark BUSCO as completed
        pipeline.update_state(
            self.output_dir, "sample1", "FLYE", "pacbio-hifi", "busco_mollusca", "completed",
            path=json_path,
            result=pipeline.parse_busco_json(json_path)
        )
        
        # Run BUSCO which should use the cached result
        json_file = pipeline.run_busco(
            fasta_file=assembly_file,
            out_dir=self.busco_dir,
            lineage="mollusca_odb10",
            threads=4,
            sample="sample1",
            assembler="FLYE",
            config="pacbio-hifi",
            state_dir=self.output_dir,
            force=False
        )
        
        # Verify JSON file is returned and run_cmd not called
        self.assertEqual(json_path, json_file)
        mock_run_cmd.assert_not_called()
        
        # Now test with a different lineage (should run)
        mock_run_cmd.reset_mock()
        json_file2 = pipeline.run_busco(
            fasta_file=assembly_file,
            out_dir=self.busco_dir,
            lineage="eukaryota_odb10",
            threads=4,
            sample="sample1",
            assembler="FLYE",
            config="pacbio-hifi",
            state_dir=self.output_dir,
            force=False
        )
        
        # Verify run_cmd was called (different lineage)
        self.assertTrue(mock_run_cmd.called)

    def test_collect_results_from_state(self):
        """Test collecting results from state file."""
        # Set up test state with multiple assemblers
        test_state = {
            "sample1": {
                "last_update": "2025-04-03T12:00:00",
                "assemblers": {
                    "FLYE_pacbio-hifi": {
                        "assembly": {
                            "status": "completed",
                            "path": "/path/to/assembly1.fasta"
                        },
                        "busco_mollusca": {
                            "status": "completed",
                            "path": "/path/to/busco1.json",
                            "result": {
                                "complete_percent": 85.0,
                                "missing_percent": 10.0,
                                "n50": 50000,
                                "num_contigs": 200,
                                "total_length": 10000000
                            }
                        },
                        "busco_eukaryota": {
                            "status": "completed",
                            "path": "/path/to/busco2.json",
                            "result": {
                                "complete_percent": 80.0,
                                "missing_percent": 15.0,
                                "n50": 50000,
                                "num_contigs": 200,
                                "total_length": 10000000
                            }
                        },
                        "gc": {
                            "status": "completed",
                            "value": 42.5
                        }
                    },
                    "CANU_pacbio-hifi": {
                        "assembly": {
                            "status": "completed",
                            "path": "/path/to/assembly2.fasta"
                        },
                        "busco_mollusca": {
                            "status": "completed",
                            "path": "/path/to/busco3.json",
                            "result": {
                                "complete_percent": 90.0,
                                "missing_percent": 5.0,
                                "n50": 70000,
                                "num_contigs": 150,
                                "total_length": 11000000
                            }
                        },
                        "busco_eukaryota": {
                            "status": "completed",
                            "path": "/path/to/busco4.json",
                            "result": {
                                "complete_percent": 88.0,
                                "missing_percent": 7.0,
                                "n50": 70000,
                                "num_contigs": 150,
                                "total_length": 11000000
                            }
                        },
                        "gc": {
                            "status": "completed",
                            "value": 43.2
                        }
                    },
                    "FLYE_nano-raw": {
                        "assembly": {
                            "status": "failed",
                            "error": "flye_failed"
                        }
                    }
                }
            }
        }
        
        pipeline.write_state(test_state, self.output_dir)
        
        # Collect results
        results = pipeline.collect_results_from_state(self.output_dir)
        
        # Verify results
        self.assertEqual(len(results), 2)  # 2 completed assemblies
        
        # Check FLYE result
        flye_result = [r for r in results if r["assembler"] == "FLYE"][0]
        self.assertEqual(flye_result["busco_mollusca"], 85.0)
        self.assertEqual(flye_result["busco_eukaryota"], 80.0)
        self.assertEqual(flye_result["gc"], 42.5)
        
        # Check CANU result
        canu_result = [r for r in results if r["assembler"] == "CANU"][0]
        self.assertEqual(canu_result["busco_mollusca"], 90.0)
        self.assertEqual(canu_result["busco_eukaryota"], 88.0)
        self.assertEqual(canu_result["gc"], 43.2)
        
        # Make sure failed assembly is not included
        failed_results = [r for r in results if r["config"] == "nano-raw"]
        self.assertEqual(len(failed_results), 0)

if __name__ == "__main__":
    unittest.main()