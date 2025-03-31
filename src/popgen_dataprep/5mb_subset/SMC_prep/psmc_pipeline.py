#!/usr/bin/env python3
"""
PSMC Preparation Pipeline

This script prepares genomic data for Pairwise Sequentially Markovian Coalescent (PSMC)
analysis to infer historical effective population sizes. It processes filtered VCF files,
generates a diploid consensus FASTA using IUPAC codes (via bcftools consensus),
converts it to PSMC input format using seqtk and fq2psmcfa, and runs the PSMC analysis
with bootstrapping.

Author: Max Borgmann
Date: March 2025
"""

import os
import sys
import argparse
import subprocess
import concurrent.futures
import logging
import time
import shutil
import glob
from pathlib import Path

class PSMCPipeline:
    """
    Pipeline for preparing and running PSMC analysis on genomic data.
    
    The workflow includes:
    1. Converting filtered VCF to a diploid consensus FASTA with IUPAC codes,
       then converting it to PSMCFA using seqtk and fq2psmcfa.
    2. Running PSMC analysis.
    3. Performing bootstrapping.
    4. Generating plots (PDF/EPS only).
    
    Attributes:
        args: Command-line arguments
        logger: Logger for pipeline progress
    """
    
    def __init__(self):
        """Initialize the PSMC pipeline with command-line arguments."""
        self.args = self._parse_arguments()
        self.logger, self.log_path = self._setup_logging()
        
        # Create output directories
        os.makedirs(self.args.output_dir, exist_ok=True)
        self.logger.info("PSMC Preparation Pipeline initialized")
        self.logger.info(f"Log file: {self.log_path}")
        
        # Check for required tools
        self._check_dependencies()
    
    def _parse_arguments(self):
        """Parse command-line arguments."""
        parser = argparse.ArgumentParser(
            description="PSMC Preparation Pipeline for demographic history inference",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        
        parser.add_argument("-i", "--input_dir", required=True,
                            help="Directory containing filtered VCF files")
        parser.add_argument("-r", "--ref_dir", required=True,
                            help="Directory containing reference FASTA files")
        parser.add_argument("-o", "--output_dir", required=True,
                            help="Output directory for PSMC files")
        parser.add_argument("-s", "--samples", nargs="+", required=True,
                            help="Sample names to process")
        parser.add_argument("-t", "--threads", type=int, default=30,
                            help="Number of threads to use")
        parser.add_argument("-b", "--bootstraps", type=int, default=100,
                            help="Number of bootstrap replicates")
        parser.add_argument("-p", "--pattern", default="4+25*2+4+6",
                            help="PSMC pattern parameter")
        parser.add_argument("--min_depth", type=int, default=10,
                            help="Minimum depth for filtering")
        parser.add_argument("--max_depth", type=int, default=100,
                            help="Maximum depth for filtering")
        parser.add_argument("--min_qual", type=int, default=0,
                            help="Minimum quality score for filtering")
        parser.add_argument("--generation_time", type=float, default=1.0,
                            help="Generation time in years for the species")
        parser.add_argument("--mutation_rate", type=float, default=0.9e-8,
                            help="Mutation rate per base per generation")
        parser.add_argument("--psmc_n", type=int, default=25,
                            help="PSMC -N parameter (number of iterations)")
        parser.add_argument("--psmc_t", type=int, default=15,
                            help="PSMC -t parameter (maximum 2N0 coalescent time)")
        parser.add_argument("--psmc_r", type=int, default=5,
                            help="PSMC -r parameter (initial theta/rho ratio)")
        parser.add_argument("--plot_options", default="-p",
                            help="Extra options for psmc_plot.pl (e.g. -p for log-scaled x-axis)")
        parser.add_argument("--y_scale", type=int, default=50000,
                            help="Default Y-axis scale for the plot (pY parameter)")
        parser.add_argument("--mode", choices=["full", "prepare", "analysis"], default="full",
                            help="Run mode: 'full' for complete pipeline, 'prepare' for input only, 'analysis' for PSMC only")
        parser.add_argument("--skip_plotting", action="store_true",
                            help="Skip plotting and consider pipeline done if combined PSMC is created")
        parser.add_argument("--force", action="store_true",
                            help="Force overwrite of existing output files")
        
        return parser.parse_args()
    
    def _setup_logging(self):
        """Configure logging to file and console."""
        os.makedirs(self.args.output_dir, exist_ok=True)
        log_file_path = os.path.join(self.args.output_dir, "psmc_pipeline.log")
        
        logger = logging.getLogger('psmc_pipeline')
        logger.setLevel(logging.DEBUG)
        
        if logger.hasHandlers():
            logger.handlers.clear()
        
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.DEBUG)
        
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_formatter = logging.Formatter('%(levelname)s: %(message)s')
        
        file_handler.setFormatter(file_formatter)
        console_handler.setFormatter(console_formatter)
        
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
        
        return logger, log_file_path
    
    def _check_dependencies(self):
        """Check if required tools are installed."""
        required_tools = [
            "bcftools", "samtools", "fq2psmcfa", "psmc", 
            "splitfa", "psmc_plot.pl", "seqtk"
        ]
        
        self.logger.info("Checking for required tools...")
        missing_tools = []
        
        for tool in required_tools:
            try:
                subprocess.run(f"which {tool}", shell=True, check=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.logger.debug(f"Tool found: {tool}")
            except subprocess.CalledProcessError:
                missing_tools.append(tool)
        
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            sys.exit(1)
        
        try:
            subprocess.run("which vcfutils.pl", shell=True, check=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.logger.debug("vcfutils.pl found")
        except subprocess.CalledProcessError:
            self.logger.error("vcfutils.pl not found in PATH; required for generating consensus from VCF")
            sys.exit(1)
        
        self.logger.info("All required tools are available")
    
    def run_command(self, cmd, workdir=None, check=True, env=None):
        """Run a shell command with logging."""
        self.logger.info(f"Running: {cmd}")
        start_time = time.time()
        
        if env is None:
            env = os.environ.copy()
        
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                check=check,
                cwd=workdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                env=env
            )
            
            if result.stdout.strip():
                self.logger.debug(f"Command stdout: {result.stdout.strip()}")
            if result.stderr.strip() and "uninitialized value" not in result.stderr:
                self.logger.debug(f"Command stderr: {result.stderr.strip()}")
            
            duration = time.time() - start_time
            self.logger.debug(f"Command completed in {duration:.2f}s")
            return True
        
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed with exit code {e.returncode}")
            if e.stdout:
                self.logger.error(f"stdout: {e.stdout.strip()}")
            if e.stderr:
                self.logger.error(f"stderr: {e.stderr.strip()}")
            return False

    def vcf_to_psmcfa(self, sample):
        """
        1) bcftools consensus --> IUPAC-coded FASTA

        2) fq2psmcfa -
        """
        sample_dir = os.path.join(self.args.output_dir, sample)
        os.makedirs(sample_dir, exist_ok=True)

        input_vcf = os.path.join(self.args.input_dir, sample, f"{sample}_filtered.vcf.gz")
        ref_fa    = os.path.join(self.args.ref_dir, "Input_fastas", sample, f"{sample}_assembly_5mb_subset.fasta")
        iupac_fa  = os.path.join(sample_dir, f"{sample}_iupac.fa")
        psmcfa_file = os.path.join(sample_dir, f"{sample}.psmcfa")

        if not os.path.exists(input_vcf):
            self.logger.error(f"VCF not found: {input_vcf}")
            return None
        if not os.path.exists(ref_fa):
            self.logger.error(f"Reference FASTA not found: {ref_fa}")
            return None

        # Create IUPAC-coded FASTA
        self.logger.info(f"Generating IUPAC-coded consensus for {sample}")
        cmd_consensus = f"bcftools consensus --iupac-codes -f {ref_fa} {input_vcf} > {iupac_fa}"
        if not self.run_command(cmd_consensus):
            self.logger.error("Consensus FASTA generation failed")
            return None
        
        # Check how many heterozygous bases we got
        self.logger.info("Counting ambiguous (IUPAC) bases:")
        self.run_command(f"grep '[RYMKSW]' {iupac_fa} | wc -l")

        # Convert to PSMCFA
        self.logger.info(f"Generating PSMCFA for {sample}")
        cmd_psmcfa = f"cat {iupac_fa} | fq2psmcfa - > {psmcfa_file}"
        if not self.run_command(cmd_psmcfa):
            self.logger.error("PSMCFA generation failed")
            return None
        
        return psmcfa_file
    
    def run_psmc(self, sample, psmcfa_file):
        """Run PSMC on the generated .psmcfa."""
        sample_dir = os.path.join(self.args.output_dir, sample)
        psmc_file  = os.path.join(sample_dir, f"{sample}.psmc")
        
        if os.path.exists(psmc_file) and not self.args.force:
            self.logger.info(f"PSMC file exists, skipping: {psmc_file}")
            return psmc_file
        
        self.logger.info(f"Running PSMC for {sample}")
        psmc_cmd = (
            f"psmc -N{self.args.psmc_n} -t{self.args.psmc_t} -r{self.args.psmc_r} "
            f"-p \"{self.args.pattern}\" -o {psmc_file} {psmcfa_file}"
        )
        if not self.run_command(psmc_cmd):
            self.logger.error("PSMC run failed")
            return None
        
        if not os.path.exists(psmc_file) or os.path.getsize(psmc_file) == 0:
            self.logger.error(f"PSMC output empty: {psmc_file}")
            return None
        
        # Quick look at the file
        self.run_command(f"head -n 10 {psmc_file}")
        return psmc_file

    def split_psmcfa(self, sample, psmcfa_file):
        """Split .psmcfa into smaller chunks for bootstrap sampling."""
        sample_dir = os.path.join(self.args.output_dir, sample)
        split_file = os.path.join(sample_dir, f"{sample}.split.psmcfa")

        if os.path.exists(split_file) and not self.args.force:
            self.logger.info(f"Split file exists, skipping: {split_file}")
            return split_file

        self.logger.info(f"Splitting PSMCFA for bootstraps: {sample}")
        split_cmd = f"splitfa {psmcfa_file} > {split_file}"
        if not self.run_command(split_cmd):
            self.logger.error("splitfa failed")
            return None
        
        if not os.path.exists(split_file) or os.path.getsize(split_file) == 0:
            self.logger.error(f"Split file empty: {split_file}")
            return None
        
        return split_file
    
    def run_bootstrap(self, sample, split_file, idx):
        """Run a single bootstrap replicate."""
        sample_dir = os.path.join(self.args.output_dir, sample)
        boot_dir   = os.path.join(sample_dir, "bootstraps")
        os.makedirs(boot_dir, exist_ok=True)

        boot_file  = os.path.join(boot_dir, f"bootstrap.{idx}.psmc")
        if os.path.exists(boot_file) and not self.args.force:
            self.logger.debug(f"Bootstrap {idx} exists, skipping.")
            return boot_file
        
        cmd_boot = (
            f"psmc -N{self.args.psmc_n} -t{self.args.psmc_t} -r{self.args.psmc_r} "
            f"-b -p \"{self.args.pattern}\" -o {boot_file} {split_file}"
        )
        self.logger.debug(f"Running bootstrap {idx} for {sample}")
        success = self.run_command(cmd_boot, check=False)
        if not success:
            self.logger.warning(f"Bootstrap {idx} failed for {sample}")
            return None
        return boot_file

    def run_bootstraps(self, sample, split_file):
        """Run multiple bootstraps in parallel."""
        self.logger.info(f"Running {self.args.bootstraps} bootstraps for {sample}")
        boot_files = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.args.threads) as exe:
            futures = {
                exe.submit(self.run_bootstrap, sample, split_file, i): i
                for i in range(1, self.args.bootstraps + 1)
            }
            for fut in concurrent.futures.as_completed(futures):
                idx = futures[fut]
                try:
                    res = fut.result()
                    if res:
                        boot_files.append(res)
                    else:
                        self.logger.warning(f"Bootstrap {idx} returned None for {sample}")
                except Exception as e:
                    self.logger.error(f"Bootstrap {idx} raised exception: {e}")
        
        self.logger.info(f"Completed {len(boot_files)}/{self.args.bootstraps} bootstraps")
        return boot_files
    
    def combine_results(self, sample, psmc_file, boot_files):
        """Concatenate original + bootstrap .psmc into a single file for plotting."""
        sample_dir   = os.path.join(self.args.output_dir, sample)
        combined_psmc = os.path.join(sample_dir, f"{sample}_combined.psmc")
        
        if os.path.exists(combined_psmc) and not self.args.force:
            self.logger.info(f"Combined file exists, skipping: {combined_psmc}")
            return combined_psmc
        
        self.logger.info(f"Combining PSMC + bootstraps for {sample}")
        cmd_cat = f"cat {psmc_file} {' '.join(boot_files)} > {combined_psmc}"
        if not self.run_command(cmd_cat):
            self.logger.error("Combining .psmc files failed")
            return None
        
        if not os.path.exists(combined_psmc) or os.path.getsize(combined_psmc) == 0:
            self.logger.error(f"Combined file empty: {combined_psmc}")
            return None
        
        return combined_psmc

    def _auto_scale_psmc(self, psmc_file):
        """
        Parse TR/RS lines from .psmc to find maximum Ne,
        and pick a Y scale based on the data. Also clamp it if it’s huge.
        """
        bin_size = 100
        mu       = self.args.mutation_rate

        theta0 = None
        max_lambda = 0.0

        try:
            with open(psmc_file, "r") as f:
                for line in f:
                    line = line.strip()
                    # TR line e.g. "TR 0.002034 0.002034"
                    if line.startswith("TR"):
                        parts = line.split()
                        val = float(parts[1])
                        if val > 0:
                            theta0 = val
                    # RS line e.g. "RS 0 0.000000 1.000000 ..."
                    if line.startswith("RS"):
                        parts = line.split()
                        if len(parts) >= 4:
                            lam = float(parts[3])
                            if lam > max_lambda:
                                max_lambda = lam
        except Exception as e:
            self.logger.warning(f"auto_scale parse error: {e}")
            return self.args.y_scale

        if not theta0 or theta0 <= 0:
            self.logger.warning("No valid theta0 found; using user y_scale.")
            return self.args.y_scale
        
        # Calculate N0 = theta0 / (4 * mu) / bin_size
        N0    = (theta0 / (4.0 * mu)) / bin_size
        max_N = N0 * max_lambda

        # PSMC’s Y scale is in "1e4" units. For example, if max_N=1e6, that’s Y=100.
        # We add ~20% buffer.
        y_auto = max_N / 1e4 * 1.2
        y_auto_int = int(y_auto + 1)

        # Clamp the value if it’s too big.
        if y_auto_int > 500:
            y_auto_int = 150
        
        self.logger.info(f"Auto-chosen Y scale: {y_auto_int} (N0={N0:.2f}, max_lambda={max_lambda:.2f})")
        return y_auto_int

    # ----------------------------------------------------------------
    # UPDATED PLOT FUNCTION
    # ----------------------------------------------------------------
    def plot_results(self, sample, combined_file):
        """Generate PSMC plots (PDF/EPS only)."""
        if self.args.skip_plotting:
            self.logger.info(f"Plotting skipped for {sample}")
            return True, []
        
        sample_dir  = os.path.join(self.args.output_dir, sample)
        plot_prefix = os.path.join(sample_dir, sample)

        # If plots already exist (PDF or EPS) and we're not forcing a re-run, skip
        existing_plots = sorted(
            glob.glob(f"{plot_prefix}*.pdf") + glob.glob(f"{plot_prefix}*.eps")
        )
        if existing_plots and not self.args.force:
            self.logger.info(f"Plot files exist, skipping: {', '.join(existing_plots)}")
            return True, existing_plots

        # Auto-scale Y
        y_auto = self._auto_scale_psmc(combined_file)
        final_y = y_auto

        self.logger.info(f"Plotting with Y scale = {final_y}")
        
        # Invoke psmc_plot.pl
        plot_cmd = (
            f"psmc_plot.pl "
            f"-u {self.args.mutation_rate} "
            f"-g {self.args.generation_time} "
            f"-Y{final_y} "
            f"{self.args.plot_options} "
            f"{plot_prefix} {combined_file}"
        )
        
        success = self.run_command(plot_cmd)
        if not success:
            self.logger.error("psmc_plot command failed")
        
        # Collect any generated PDF/EPS files that start with the same prefix
        generated_plots = sorted(
            glob.glob(f"{plot_prefix}*.pdf") + glob.glob(f"{plot_prefix}*.eps")
        )
        if not generated_plots:
            self.logger.error("No PDF/EPS plots were created!")
            return False, []
        
        for gp in generated_plots:
            self.logger.info(f"Plot generated: {gp}")
        
        return True, generated_plots
    # ----------------------------------------------------------------

    def process_sample(self, sample):
        """Complete workflow for one sample."""
        self.logger.info(f"Processing sample: {sample}")
        results = {
            "sample": sample,
            "psmcfa_file": None,
            "psmc_file": None,
            "combined_file": None,
            "plot_files": [],
            "status": "Failed",
            "error": None
        }
        
        try:
            # If analysis mode, skip VCF->PSMCFA
            if self.args.mode == "analysis":
                psmcfa_file = os.path.join(self.args.output_dir, sample, f"{sample}.psmcfa")
                if not os.path.exists(psmcfa_file):
                    msg = f"PSMCFA not found for analysis: {psmcfa_file}"
                    self.logger.error(msg)
                    results["error"] = msg
                    return results
                self.logger.info(f"Starting from existing psmcfa: {psmcfa_file}")
                results["psmcfa_file"] = psmcfa_file
            else:
                # Generate PSMCFA
                psmcfa_file = self.vcf_to_psmcfa(sample)
                if not psmcfa_file:
                    msg = f"Failed to create psmcfa for {sample}"
                    self.logger.error(msg)
                    results["error"] = msg
                    return results
                results["psmcfa_file"] = psmcfa_file
            
            # If prepare only, stop
            if self.args.mode == "prepare":
                self.logger.info("Prepare-only mode: done after PSMCFA.")
                results["status"] = "Completed - PSMCFA only"
                return results

            # 2) Run PSMC
            psmc_file = self.run_psmc(sample, psmcfa_file)
            if not psmc_file:
                msg = f"PSMC run failed for {sample}"
                self.logger.error(msg)
                results["error"] = msg
                return results
            results["psmc_file"] = psmc_file

            # 3) Bootstrapping
            split_file = self.split_psmcfa(sample, psmcfa_file)
            if not split_file:
                msg = f"split_psmcfa failed for {sample}"
                self.logger.error(msg)
                results["error"] = msg
                return results

            boot_files = self.run_bootstraps(sample, split_file)
            if not boot_files:
                msg = f"No bootstrap results for {sample}"
                self.logger.error(msg)
                results["error"] = msg
                return results

            # 4) Combine
            combined_file = self.combine_results(sample, psmc_file, boot_files)
            if not combined_file:
                msg = f"Combining results failed for {sample}"
                self.logger.error(msg)
                results["error"] = msg
                return results
            results["combined_file"] = combined_file
            
            # 5) Plot
            if self.args.skip_plotting:
                self.logger.info("Skipping plotting per user request.")
                results["status"] = "Completed - No plots"
                return results

            plot_ok, plot_files = self.plot_results(sample, combined_file)
            results["plot_files"] = plot_files
            if plot_ok and plot_files:
                self.logger.info(f"{sample} completed successfully.")
                results["status"] = "Completed"
            else:
                self.logger.warning(f"{sample} plot missing or failed.")
                results["status"] = "Completed - Plot failures"
                results["error"] = "Plot generation failed but combined file was created"
            
            return results
        
        except Exception as e:
            msg = f"Unexpected error for {sample}: {e}"
            self.logger.error(msg)
            import traceback
            self.logger.debug(traceback.format_exc())
            results["error"] = msg
            return results

    def run(self):
        """Run pipeline for all samples."""
        self.logger.info("Starting PSMC pipeline")
        self.logger.info(f"Mode: {self.args.mode}")
        self.logger.info(f"Samples: {', '.join(self.args.samples)}")
        self.logger.info(f"Mutation rate: {self.args.mutation_rate}, Generation time: {self.args.generation_time}")

        results = []
        for sample in self.args.samples:
            t0 = time.time()
            out = self.process_sample(sample)
            out["duration"] = f"{time.time()-t0:.2f} seconds"
            results.append(out)
            self.logger.info(f"Sample {sample} => {out['status']} in {out['duration']}")

        # Summary
        completed = sum(1 for r in results if r["status"].startswith("Completed"))
        self.logger.info("\n==== Pipeline Summary ====")
        for r in results:
            smp = r["sample"]
            stat = r["status"]
            dur = r.get("duration","?")
            nplots = len(r["plot_files"])
            self.logger.info(f"{smp}: {stat}, {nplots} plot(s), {dur}")
            if r["error"]:
                self.logger.info(f"  Error: {r['error']}")
        
        self.logger.info(f"Processed {len(results)} sample(s). {completed} completed successfully.")
        return completed == len(results)

if __name__ == "__main__":
    pipeline = PSMCPipeline()
    ok = pipeline.run()
    sys.exit(0 if ok else 1)
