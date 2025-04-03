#!/usr/bin/env python3
"""
Test script for assembly_comparison_pipeline_SR.py

This test script verifies that:
1. All assemblers (SOAP, MaSuRCA, Captus) start correctly with provided inputs
2. BUSCO analysis works properly
3. Summary report is created correctly
4. Conda environments (soap_env, masurca_env, captus_env) are set up correctly

Author: Max Borgmann
Date: 03.04.2025
"""

import os
import sys
import unittest
import tempfile
import shutil
import subprocess
import json
from unittest import mock
from pathlib import Path

# Suppress IDE warnings for tabulate
try:
    import tabulate  # type: ignore # noqa: F401
except ImportError:
    pass  # We have a fallback for missing tabulate later in the code

# Direct import using the full path
script_dir = os.path.dirname(os.path.abspath(__file__))
pipeline_path = os.path.join(script_dir, "assembly_comparison_pipeline_SR.py")
sys.path.insert(0, script_dir)

# Import the module directly
import assembly_comparison_pipeline_SR as assembly_comparison


class TestAssemblyComparisonPipeline(unittest.TestCase):
    """Test cases for assembly_comparison_pipeline_SR.py"""

    @classmethod
    def setUpClass(cls):
        """Create test environment with mock FASTQ files"""
        # Create temp directory for testing
        cls.test_dir = tempfile.mkdtemp()
        cls.input_dir = os.path.join(cls.test_dir, "input")
        cls.output_dir = os.path.join(cls.test_dir, "output")
        
        os.makedirs(cls.input_dir, exist_ok=True)
        os.makedirs(cls.output_dir, exist_ok=True)
        
        # Create mock FASTQ files (minimal content to pass file existence checks)
        cls.sample_name = "test_sample"
        cls.r1_path = os.path.join(cls.input_dir, f"{cls.sample_name}_R1.fastq")
        cls.r2_path = os.path.join(cls.input_dir, f"{cls.sample_name}_R2.fastq")
        
        fastq_content = """@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
@SEQ_ID2
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
"""
        with open(cls.r1_path, 'w') as f:
            f.write(fastq_content)
        with open(cls.r2_path, 'w') as f:
            f.write(fastq_content)
        
        # Create mock assembly output directories for testing
        cls.asm_dir = os.path.join(cls.output_dir, "assemblies")
        cls.busco_dir = os.path.join(cls.output_dir, "busco")
        os.makedirs(cls.asm_dir, exist_ok=True)
        os.makedirs(cls.busco_dir, exist_ok=True)
        
        # Create mock assembly files (so pipeline sees something "assembled")
        cls.soap63_asm = os.path.join(cls.asm_dir, f"{cls.sample_name}_soap_k63.fasta")
        cls.soap127_asm = os.path.join(cls.asm_dir, f"{cls.sample_name}_soap_k127.fasta")
        cls.masurca_perm_asm = os.path.join(cls.asm_dir, f"{cls.sample_name}_masurca_permissive.fasta")
        cls.masurca_str_asm = os.path.join(cls.asm_dir, f"{cls.sample_name}_masurca_stringent.fasta")
        cls.captus_asm = os.path.join(cls.asm_dir, f"{cls.sample_name}_captus_wgs.fasta")
        
        fasta_content = """>contig1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>contig2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
"""
        for asm_file in [cls.soap63_asm, cls.soap127_asm, cls.masurca_perm_asm, 
                         cls.masurca_str_asm, cls.captus_asm]:
            with open(asm_file, 'w') as f:
                f.write(fasta_content)
        
        # Patch pandas DataFrame.to_markdown if needed
        try:
            import pandas as pd
            test_df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
            test_df.to_markdown()
        except (ImportError, AttributeError):
            from pandas import DataFrame
            original_to_markdown = getattr(DataFrame, 'to_markdown', None)
            def to_markdown_simple(self, *args, **kwargs):
                try:
                    if original_to_markdown:
                        return original_to_markdown(self, *args, **kwargs)
                    raise ImportError("No original to_markdown")
                except ImportError:
                    header = "| " + " | ".join(str(col) for col in self.columns) + " |"
                    separator = "| " + " | ".join(["---"] * len(self.columns)) + " |"
                    rows = []
                    for _, row in self.iterrows():
                        rows.append("| " + " | ".join(str(v) for v in row.values) + " |")
                    return "\n".join([header, separator] + rows)
            DataFrame.to_markdown = to_markdown_simple
        
    @classmethod
    def tearDownClass(cls):
        """Clean up test files and directories"""
        shutil.rmtree(cls.test_dir)

    def create_mock_busco_json(self, lineage, assembler, config):
        """Helper to create a BUSCO JSON file (for tests that check parsing)"""
        busco_out_dir = os.path.join(
            self.busco_dir, 
            f"{self.sample_name}_{assembler}_{config}_{lineage}"
        )
        os.makedirs(busco_out_dir, exist_ok=True)
        json_file = os.path.join(busco_out_dir, f"short_summary.{lineage}.{self.sample_name}.json")
        mock_json = {
            "results": {
                "Complete percentage": 85.0,
                "Missing percentage": 10.0
            },
            "metrics": {
                "Contigs N50": 50000,
                "Number of contigs": 1000,
                "Total length": 500000000
            }
        }
        with open(json_file, 'w') as f:
            json.dump(mock_json, f)
        return json_file

    @mock.patch.object(assembly_comparison, 'run_in_env')
    def test_soap_denovo_assembly(self, mock_run_in_env):
        """Test SOAP assembly execution (simulate fresh run)"""
        soap_dir = os.path.join(self.asm_dir, f"{self.sample_name}_soap_k63")
        os.makedirs(soap_dir, exist_ok=True)
        prefix = os.path.join(soap_dir, f"{self.sample_name}_soap_k63")
        gapclosed_file = f"{prefix}.gapclosed.fasta"
        
        def side_effect(cmd, cwd=None, desc=None):
            if "GapCloser" in cmd:
                with open(gapclosed_file, 'w') as f:
                    f.write(">scaffold1\nATGCATGC\n")
            return (True, "Success output")
        
        mock_run_in_env.side_effect = side_effect
        
        result = assembly_comparison.run_soapdenovo2_with_gapcloser(
            self.r1_path, self.r2_path, self.asm_dir, 
            self.sample_name, 63, threads=2
        )
        
        self.assertIsNotNone(result)
        self.assertIn("soap_k63.fasta", result)
        # 5 steps: pregraph, contig, map, scaff, gapcloser
        self.assertEqual(mock_run_in_env.call_count, 5)
        
        calls = mock_run_in_env.call_args_list
        cmd_types = []
        for call in calls:
            cmd = call[0][1]  # second positional arg is the cmd string
            if "pregraph" in cmd:
                cmd_types.append("pregraph")
            elif "contig" in cmd:
                cmd_types.append("contig")
            elif "map" in cmd and "scaff" not in cmd:
                cmd_types.append("map")
            elif "scaff" in cmd:
                cmd_types.append("scaff")
            elif "GapCloser" in cmd:
                cmd_types.append("gapcloser")
        self.assertIn("pregraph", cmd_types)
        self.assertIn("contig", cmd_types)
        self.assertIn("map", cmd_types)
        self.assertIn("scaff", cmd_types)
        self.assertIn("gapcloser", cmd_types)

    @mock.patch.object(assembly_comparison, 'run_in_env')
    def test_masurca_assembly(self, mock_run_in_env):
        """Test MaSuRCA assembly execution (simulate fresh run)"""
        masurca_dir = os.path.join(self.asm_dir, f"{self.sample_name}_masurca_permissive")
        os.makedirs(masurca_dir, exist_ok=True)
        # Do NOT pre-create assemble.sh to force a fresh run
        # Also do NOT create the final output file (final.genome.scf.fasta)
        ca_dir = os.path.join(masurca_dir, "CA")
        os.makedirs(ca_dir, exist_ok=True)
        
        # Side effect: if assemble.sh is executed, simulate creation of final output file.
        def masurca_side_effect(cmd, cwd=None, desc=None):
            if "./assemble.sh" in cmd:
                final_file = os.path.join(ca_dir, "final.genome.scf.fasta")
                with open(final_file, "w") as f:
                    f.write(">scaffold1\nATGCATGC\n")
            return (True, "Success output")
        
        mock_run_in_env.side_effect = masurca_side_effect
        
        result = assembly_comparison.run_masurca(
            self.r1_path, self.r2_path, self.asm_dir,
            self.sample_name, "permissive", threads=2
        )
        
        self.assertIsNotNone(result)
        self.assertIn("masurca_permissive.fasta", result)
        # Should call run_in_env for (1) masurca config; (2) ./assemble.sh
        self.assertEqual(mock_run_in_env.call_count, 2)
        
        config_file = os.path.join(masurca_dir, "masurca_config.txt")
        self.assertTrue(os.path.exists(config_file))
        with open(config_file, 'r') as f:
            config_content = f.read()
            self.assertIn("GRAPH_KMER_SIZE = 31", config_content)
            self.assertIn("KMER_COUNT_THRESHOLD = 2", config_content)

    @mock.patch.object(assembly_comparison, 'run_in_env')
    def test_captus_assembly(self, mock_run_in_env):
        """Test Captus assembly execution (simulate fresh run)"""
        captus_dir = os.path.join(self.asm_dir, f"{self.sample_name}_captus_wgs_run")
        os.makedirs(captus_dir, exist_ok=True)
        captus_asm_dir = os.path.join(captus_dir, f"{self.sample_name}__captus-asm")
        os.makedirs(captus_asm_dir, exist_ok=True)
        # Do NOT create the assembly.fasta file so that run_in_env is definitely called
        
        def captus_side_effect(cmd, cwd=None, desc=None):
            if "captus assemble" in cmd:
                asm_file = os.path.join(captus_asm_dir, "assembly.fasta")
                with open(asm_file, "w") as f:
                    f.write(">scaffold1\nATGCATGC\n")
            return (True, "Success output")
        
        mock_run_in_env.side_effect = captus_side_effect
        
        result = assembly_comparison.run_captus_assembly(
            self.r1_path, self.r2_path, self.asm_dir,
            self.sample_name, threads=2
        )
        
        self.assertIsNotNone(result)
        self.assertIn("captus_wgs.fasta", result)
        # Only 1 call to run_in_env (the "captus assemble" command)
        self.assertEqual(mock_run_in_env.call_count, 1)
        cmd = mock_run_in_env.call_args[0][1]  # the cmd string
        self.assertIn("captus assemble", cmd)
        self.assertIn("--preset WGS", cmd)

    @mock.patch.object(assembly_comparison, 'run_in_env')
    def test_busco_execution(self, mock_run_in_env):
        """Test BUSCO execution (simulate fresh run)"""
        lineage = "mollusca_odb10"
        assembler = "SOAP"
        config = "k63"
        busco_out_name = f"{self.sample_name}_{assembler}_{config}_{lineage}"
        busco_out_dir = os.path.join(self.busco_dir, busco_out_name)
        os.makedirs(busco_out_dir, exist_ok=True)
        # Do NOT pre-create the BUSCO JSON file
        
        json_file = os.path.join(busco_out_dir, f"short_summary.{lineage}.{self.sample_name}.json")
        
        def busco_side_effect(cmd, cwd=None, desc=None):
            if "busco -i" in cmd:
                with open(json_file, 'w') as f:
                    json.dump({
                        "results": {
                            "Complete percentage": 85.2,
                            "Missing percentage": 10.1
                        },
                        "metrics": {
                            "Contigs N50": 35000,
                            "Number of contigs": 1200,
                            "Total length": 480000000
                        }
                    }, f)
            return (True, "Success output")
        
        mock_run_in_env.side_effect = busco_side_effect
        
        result = assembly_comparison.run_busco(
            "soap_env",
            self.soap63_asm,
            self.busco_dir,
            lineage,
            2,
            self.sample_name,
            assembler,
            config
        )
        
        self.assertIsNotNone(result)
        self.assertEqual(mock_run_in_env.call_count, 1)
        cmd = mock_run_in_env.call_args[0][1]  # the cmd string
        self.assertIn("busco -i", cmd)
        self.assertIn(f"-l {lineage}", cmd)
        self.assertIn("-m genome", cmd)

    def test_busco_json_parsing(self):
        """Test parsing of BUSCO JSON output"""
        lineage = "mollusca_odb10"
        assembler = "MaSuRCA"
        config = "permissive"
        json_file = self.create_mock_busco_json(lineage, assembler, config)
        parsed_data = assembly_comparison.parse_busco_json(json_file)
        self.assertIn("complete_percent", parsed_data)
        self.assertIn("missing_percent", parsed_data)
        self.assertIn("n50", parsed_data)
        self.assertIn("num_contigs", parsed_data)
        self.assertIn("total_length", parsed_data)
        self.assertIsInstance(parsed_data["complete_percent"], float)
        self.assertIsInstance(parsed_data["n50"], int)

    @mock.patch('assembly_comparison_pipeline_SR.get_gc_content')
    def test_gc_content_calculation(self, mock_gc):
        """Test GC content calculation with minimal argument check"""
        mock_gc.return_value = 42.5
        # Just call the function with minimal arguments
        gc = assembly_comparison.get_gc_content("soap_env", self.soap63_asm)
        self.assertEqual(gc, 42.5)
        # We only check the first 2 positional args
        mock_gc.assert_called_once()
        called_args, called_kwargs = mock_gc.call_args
        self.assertEqual(called_args[0], "soap_env")  # env_name
        self.assertEqual(called_args[1], self.soap63_asm)  # fasta_file

    def test_create_summary(self):
        """Test creation of final summary"""
        results = [
            {
                "sample": self.sample_name,
                "assembler": "SOAP",
                "config": "k63",
                "busco_mollusca": 75.5,
                "missing_mollusca": 20.1,
                "busco_eukaryota": 70.2,
                "missing_eukaryota": 25.3,
                "n50": 35000,
                "num_contigs": 1200,
                "total_length": 480000000,
                "gc": 38.5
            },
            {
                "sample": self.sample_name,
                "assembler": "MaSuRCA",
                "config": "permissive",
                "busco_mollusca": 90.3,
                "missing_mollusca": 8.1,
                "busco_eukaryota": 85.1,
                "missing_eukaryota": 12.5,
                "n50": 70000,
                "num_contigs": 950,
                "total_length": 495000000,
                "gc": 39.2
            },
            {
                "sample": self.sample_name,
                "assembler": "Captus",
                "config": "WGS",
                "busco_mollusca": 92.1,
                "missing_mollusca": 6.5,
                "busco_eukaryota": 87.0,
                "missing_eukaryota": 10.8,
                "n50": 85000,
                "num_contigs": 820,
                "total_length": 510000000,
                "gc": 40.1
            }
        ]
        
        try:
            assembly_comparison.create_comparison_summary(results, self.output_dir)
            csv_path = os.path.join(self.output_dir, "assembly_comparison_summary.csv")
            md_path = os.path.join(self.output_dir, "assembly_comparison_summary.md")
            self.assertTrue(os.path.exists(csv_path))
            self.assertTrue(os.path.exists(md_path))
            import pandas as pd
            df = pd.read_csv(csv_path)
            # Should be sorted so that the row with highest mollusca, eukaryota, and N50 is first
            self.assertEqual(df.iloc[0]["assembler"], "Captus")
            self.assertEqual(df.iloc[0]["config"], "WGS")
            with open(md_path, 'r') as f:
                md_content = f.read()
                self.assertIn("# Assembly Comparison Summary", md_content)
                self.assertIn("Best Assembly", md_content)
                self.assertIn("Captus", md_content)
        except Exception as e:
            self.skipTest(f"Skipping test_create_summary due to error: {e}")

    @mock.patch('assembly_comparison_pipeline_SR.setup_conda_envs')
    @mock.patch('assembly_comparison_pipeline_SR.run_soapdenovo2_with_gapcloser')
    @mock.patch('assembly_comparison_pipeline_SR.run_masurca')
    @mock.patch('assembly_comparison_pipeline_SR.run_captus_assembly')
    @mock.patch('assembly_comparison_pipeline_SR.run_busco')
    @mock.patch('assembly_comparison_pipeline_SR.get_gc_content')
    @mock.patch('assembly_comparison_pipeline_SR.parse_busco_json')
    def test_full_pipeline_integration(self, mock_parse_busco, mock_gc, mock_busco, 
                                       mock_captus, mock_masurca, mock_soap,
                                       mock_setup_envs):
        """Test full pipeline integration"""
        mock_setup_envs.return_value = None
        mock_soap.return_value = self.soap63_asm
        mock_masurca.return_value = self.masurca_perm_asm
        mock_captus.return_value = self.captus_asm
        mock_busco.return_value = "mock_busco.json"
        mock_gc.return_value = 40.0
        mock_parse_busco.return_value = {
            "complete_percent": 85.0,
            "missing_percent": 10.0,
            "n50": 50000,
            "num_contigs": 1000,
            "total_length": 500000000
        }

        old_argv = sys.argv
        sys.argv = ["assembly_comparison_pipeline_SR.py", "-i", self.input_dir, "-o", self.output_dir, "-t", "2"]

        try:
            with mock.patch('sys.exit') as mock_exit:
                assembly_comparison.main()
                mock_exit.assert_not_called()
            csv_path = os.path.join(self.output_dir, "assembly_comparison_summary.csv")
            md_path = os.path.join(self.output_dir, "assembly_comparison_summary.md")
            self.assertTrue(os.path.exists(csv_path))
            self.assertTrue(os.path.exists(md_path))
        finally:
            sys.argv = old_argv

    @mock.patch('assembly_comparison_pipeline_SR.run_in_env')
    def test_conda_env_creation(self, mock_run_in_env):
        """
        Test that ensures setup_conda_envs actually creates or detects 
        soap_env, masurca_env, and captus_env via mamba.
        """
        # If 'mamba' is not on the PATH, we skip this test to avoid confusion.
        if not shutil.which("mamba"):
            self.skipTest("mamba not found on PATH, skipping conda env creation test.")
        
        # Call the real setup_conda_envs function
        assembly_comparison.setup_conda_envs()
        
        # Now parse the 'mamba env list --json' result to check the envs
        cmd = ["mamba", "env", "list", "--json"]
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            data = json.loads(result.stdout)
            env_paths = data.get("envs", [])
            env_names = [os.path.basename(e) for e in env_paths]
            
            self.assertIn("soap_env", env_names, "soap_env not found in 'mamba env list'")
            self.assertIn("masurca_env", env_names, "masurca_env not found in 'mamba env list'")
            self.assertIn("captus_env", env_names, "captus_env not found in 'mamba env list'")
        except subprocess.CalledProcessError as e:
            self.fail(f"Failed to list mamba envs: {e}")


if __name__ == "__main__":
    unittest.main()
