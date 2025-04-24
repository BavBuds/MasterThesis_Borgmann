#!/usr/bin/env python3
"""
Simplified Coverage‑&‑Variant QC pipeline
────────────────────────────────────────
• Maps reads **against the full assembly** for every sample
• Uses the full BAM directly with the subset FASTA for variant calling
• Generates metrics both before and after 3x median coverage filtering
• Includes all samples: BB, HT, HT2, and HT2_ON_HTref (cross-mapping)

Author: Max Borgmann (Modified)
"""

from __future__ import annotations

import argparse
import logging
import subprocess
from pathlib import Path
from typing import Dict, List

import numpy as np


# ────────────────────────────────
# CLI & logging helpers
# ────────────────────────────────
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Depth / SNP QC on subset contigs, mapping to full assemblies"
    )
    p.add_argument(
        "-o",
        "--base",
        default=(
            "/data/proj2/home/students/m.borgmann/Master_thesis/"
            "data/processed/Popgen_analysis/dataprep/5mb_subset"
        ),
        help="Base output directory (same one used by FASTA‑subset script)",
    )
    p.add_argument("-t", "--threads", type=int, default=20, help="CPU threads")
    return p.parse_args()


def setup_logging(outdir: Path) -> logging.Logger:
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / "pipeline.log"

    lg = logging.getLogger("qc")
    lg.setLevel(logging.DEBUG)
    lg.handlers.clear()

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    lg.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    lg.addHandler(ch)

    lg.info("Logging → %s", log_file)
    return lg


def run(cmd: str, lg: logging.Logger, cwd: Path | None = None) -> bool:
    """Run a shell command and report its result."""
    lg.info("↪ %s", cmd)
    try:
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return True
    except subprocess.CalledProcessError as err:  # noqa: PERF203
        lg.error("command failed (%s)", err.returncode)
        lg.debug(err.stdout)
        lg.debug(err.stderr)
        return False


# ────────────────────────────────
# sample configuration
# ────────────────────────────────
BASE_FASTA_DIR = Path(
    "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies"
)

SAMPLES: Dict[str, Dict] = {
    # long‑read ONT assembly
    "HT": {
        "assembly_full": BASE_FASTA_DIR
        / "Hexaplex_ONT_FLYE/Polished_assembly/Final/Hexaplex_assembly_ONT_FLYE_polished.fasta",
        "assembly_subset": "Input_fastas/HT/HT_assembly_5mb_subset.fasta",
        "reads": (
            "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/Hexaplex_trunculus_nanopore/20241217_DNA_Tellier_HT_31_PB/20241217_1155_2A_PAY72199_a0bfda36/fastq_pass/merged_HT_reads.fastq"
        ),
        "mapper": "minimap2",
        "datatype": "ONT",
    },
    # short‑read assemblies
    "HT2": {
        "assembly_full": BASE_FASTA_DIR
        / "Hexaplex_Captus/02_assemblies/SRR28865916__captus-asm/01_assembly/assembly.fasta",
        "assembly_subset": "Input_fastas/HT2/HT2_assembly_5mb_subset.fasta",
        "reads_1": (
            "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/"
            "Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/"
            "SRR28865916/SRR28865916_R1.fastq"
        ),
        "reads_2": (
            "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/"
            "Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/"
            "SRR28865916/SRR28865916_R2.fastq"
        ),
        "mapper": "bwa‑mem2",
        "datatype": "ILLUMINA",
    },
    "BB": {
        "assembly_full": BASE_FASTA_DIR
        / "Bolinus_Captus/SRR28863561__captus-asm/01_assembly/assembly.fasta",
        "assembly_subset": "Input_fastas/BB/BB_assembly_5mb_subset.fasta",
        "reads_1": (
            "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/"
            "Reference_data/WGS_projects/PRJNA1106534_Bolinus_brandaris/"
            "SRR28863561/SRR28863561_R1.fastq"
        ),
        "reads_2": (
            "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/"
            "Reference_data/WGS_projects/PRJNA1106534_Bolinus_brandaris/"
            "SRR28863561/SRR28863561_R2.fastq"
        ),
        "mapper": "bwa‑mem2",
        "datatype": "ILLUMINA",
    },
    # NEW cross‑mapping: short‑read HT2 to long‑read HT reference
    "HT2_ON_HTref": {
        "assembly_full": BASE_FASTA_DIR
        / "Hexaplex_ONT_FLYE/Polished_assembly/Final/Hexaplex_assembly_ONT_FLYE_polished.fasta",
        "assembly_subset": "Input_fastas/HT/HT_assembly_5mb_subset.fasta",
        "reads_1": (
            "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/"
            "Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/"
            "SRR28865916/SRR28865916_R1.fastq"
        ),
        "reads_2": (
            "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/"
            "Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/"
            "SRR28865916/SRR28865916_R2.fastq"
        ),
        "mapper": "bwa‑mem2",
        "datatype": "ILLUMINA",
    },
}


# ────────────────────────────────
# helper utilities
# ────────────────────────────────
def map_reads(
    sample: str,
    meta: Dict,
    sam_path: Path,
    threads: int,
    lg: logging.Logger,
) -> bool:
    """
    Map reads and write *SAM* to sam_path (Illumina) or pipe straight to BAM (ONT).

    For ONT we add **-L** so oversize CIGARs go to the CG tag, keeping the SAM
    compatible with samtools' long‑CIGAR mode later.
    """
    if meta["mapper"] == "minimap2":  # ONT
        cmd = (
            f"minimap2 -t {threads} -ax map-ont -L {meta['assembly_full']} "
            f"{meta['reads']} > {sam_path}"
        )
        return run(cmd, lg)

    if meta["mapper"] == "bwa‑mem2":  # Illumina
        if not run(f"bwa-mem2 index {meta['assembly_full']}", lg):
            return False
        cmd = (
            f"bwa-mem2 mem -t {threads} {meta['assembly_full']} "
            f"{meta['reads_1']} {meta['reads_2']} > {sam_path}"
        )
        return run(cmd, lg)

    lg.error("Unknown mapper for %s", sample)
    return False


def extract_depth_for_subset(
    full_bam: Path, subset_fasta: Path, depth_txt: Path, threads: int, lg: logging.Logger
) -> bool:
    """Extract depth data by analyzing entire BAM and then filtering for subset contigs."""
    # 1. Check if files exist
    if not full_bam.exists():
        lg.error(f"BAM file does not exist: {full_bam}")
        return False
    if not subset_fasta.exists():
        lg.error(f"Subset FASTA file does not exist: {subset_fasta}")
        return False
    
    # 2. Create FASTA index if needed and get contig names
    if not run(f"samtools faidx {subset_fasta}", lg):
        lg.error(f"Failed to index FASTA file: {subset_fasta}")
        return False
        
    # Read contig names from FASTA index
    fasta_contigs = []
    with open(f"{subset_fasta}.fai", "r") as fai:
        for line in fai:
            contig = line.split()[0]
            fasta_contigs.append(contig)
    
    if not fasta_contigs:
        lg.error("No contigs found in FASTA index")
        return False
    
    lg.info(f"Found {len(fasta_contigs)} contigs in FASTA subset")
    
    # Get the size of the BAM file to estimate processing time
    bam_size = full_bam.stat().st_size / (1024*1024*1024)  # Size in GB
    lg.info(f"Processing BAM file of {bam_size:.2f} GB - this may take a while...")
    
    # Run samtools depth on entire BAM but limit output to our contigs using grep
    # This is more efficient than running depth on whole BAM and filtering afterwards
    contig_pattern = "|".join(fasta_contigs)
    cmd = f"samtools depth -@ {threads} {full_bam} | grep -E '^({contig_pattern})' > {depth_txt}"
    
    lg.info("Extracting depth for subset contigs from full BAM...")
    try:
        # Use subprocess.call directly for this piped command
        result = subprocess.call(cmd, shell=True)
        if result == 0:
            lg.info("Depth extraction successful")
            
            # Check if output file has content
            if depth_txt.exists() and depth_txt.stat().st_size > 0:
                # Get line count to verify
                line_count = int(subprocess.check_output(f"wc -l < {depth_txt}", shell=True).strip())
                lg.info(f"Extracted {line_count} depth data points")
                return True
            else:
                lg.warning("Depth file exists but is empty")
        else:
            lg.error(f"Depth extraction command failed with code {result}")
    except Exception as e:
        lg.error(f"Error during depth extraction: {e}")
    
    # If we reach here, the command failed or produced empty output
    lg.warning("Full BAM depth extraction failed, trying one contig at a time...")
    
    # Try extracting depth one contig at a time
    with open(depth_txt, "w") as outfile:
        for contig in fasta_contigs:
            lg.info(f"Extracting depth for contig: {contig}")
            temp_file = depth_txt.parent / f"temp_{contig}.depth"
            try:
                # Use -a flag to output all positions, including zero-coverage
                single_cmd = f"samtools depth -a {full_bam} -r {contig} > {temp_file}"
                subprocess.call(single_cmd, shell=True)
                
                if temp_file.exists() and temp_file.stat().st_size > 0:
                    with open(temp_file, "r") as infile:
                        outfile.write(infile.read())
                    lg.info(f"Successfully extracted depth for {contig}")
                else:
                    lg.warning(f"No depth data for {contig}")
                
                # Clean up temp file
                if temp_file.exists():
                    temp_file.unlink()
            except Exception as e:
                lg.error(f"Error extracting depth for {contig}: {e}")
    
    # Check if we got any data
    if depth_txt.exists() and depth_txt.stat().st_size > 0:
        lg.info("Successfully extracted depth data")
        return True
    
    lg.error("All depth extraction methods failed")
    return False


def calculate_percent_above_10x(depths: np.ndarray) -> float:
    """Calculate percentage of bases with coverage ≥ a threshold."""
    if depths.size == 0:
        return 0.0
    return (depths >= 10).sum() / depths.size * 100.0


# ─────────────────────────────────────────────────────────────────────────────
# Write QC report
# ─────────────────────────────────────────────────────────────────────────────
def write_qc_report(results: list[dict], qc_dir: Path, lg: logging.Logger) -> None:
    """Create a text table mirroring the old pipeline's QC_report.txt."""
    if not results:
        return
    
    # Generate two reports: before and after filtering
    for report_type in ["before_filtering", "after_filtering"]:
        hdr = (
            "Sample".ljust(12)
            + "Mean Cov".rjust(10)
            + "Median Cov".rjust(12)
            + "SD Cov".rjust(12)
            + "Dyn Cap".rjust(10)
            + "% ≥10x".rjust(10)
            + "Het Sites".rjust(12)
            + "Het (%)".rjust(10)
            + "Missing".rjust(10)
            + "Miss (%)".rjust(10)
            + "Status".rjust(25)
            + "\n"
            + "-" * 130
        )
        lines: List[str] = [hdr]
        
        for r in results:
            # Skip samples that don't have both before and after stats
            if f"{report_type}_stats" not in r:
                continue
                
            stats = r[f"{report_type}_stats"]
            line = (
                r["sample"].ljust(12)
                + f"{stats['mean_cov']:.2f}".rjust(10)
                + f"{stats['median_cov']:.2f}".rjust(12)
                + f"{stats['sd_cov']:.2f}".rjust(12)
                + f"{stats['depth_cap']}".rjust(10)
                + f"{stats['pct_10x']:.2f}".rjust(10)
                + f"{stats['het']}".rjust(12)
                + f"{stats['het_pct']:.2f}".rjust(10)
                + "0".rjust(10)  # Missing
                + "0.00".rjust(10)  # Miss %
                + "Completed successfully".rjust(25)
            )
            lines.append(line)
        
        report = "\n".join(lines)
        report_path = qc_dir / f"QC_report_{report_type}.txt"
        report_path.write_text(report, encoding="utf-8")
        lg.info("%s QC report written → %s", report_type, report_path)
    
    # Also create a LaTeX-formatted table for potential inclusion in the paper
    latex_report_path = qc_dir / "QC_report_latex.txt"
    latex_lines = [
        "\\begin{table}[H]",
        "\\centering",
        "\\caption{Coverage and variant statistics across all samples before and after implementing 3$\\times$ median coverage filtering.}",
        "\\label{tab:CoverageStatistics_comparison}",
        "\\footnotesize",
        "\\renewcommand{\\arraystretch}{1.2}",
        "\\setlength{\\tabcolsep}{4pt}",
        "\\begin{tabular}{lccccccc}",
        "\\toprule",
        "\\textbf{Sample} & \\textbf{Filter} & \\textbf{Mean Cov} & \\textbf{Median Cov} & \\textbf{\\% $\\geq$10x} & \\textbf{Het Sites} & \\textbf{Het (\\%)} & \\textbf{Depth Cap} \\\\",
        "\\midrule"
    ]
    
    for r in results:
        if "before_filtering_stats" not in r or "after_filtering_stats" not in r:
            continue
            
        before = r["before_filtering_stats"]
        after = r["after_filtering_stats"]
        
        latex_lines.extend([
            f"{r['sample']} & Before & {before['mean_cov']:.2f} & {before['median_cov']:.2f} & {before['pct_10x']:.2f} & {before['het']} & {before['het_pct']:.2f} & N/A \\\\",
            f"{r['sample']} & After & {after['mean_cov']:.2f} & {after['median_cov']:.2f} & {after['pct_10x']:.2f} & {after['het']} & {after['het_pct']:.2f} & {after['depth_cap']} \\\\"
        ])
    
    latex_lines.extend([
        "\\bottomrule",
        "\\end{tabular}",
        "\\end{table}"
    ])
    
    latex_report = "\n".join(latex_lines)
    latex_report_path.write_text(latex_report, encoding="utf-8")
    lg.info("LaTeX-formatted QC report written → %s", latex_report_path)


# ─────────────────────────────────────────────────────────────────────────────
# main entry‑point
# ─────────────────────────────────────────────────────────────────────────────
def main() -> None:  # noqa: C901 (large function kept single‑file)
    args = parse_args()
    base = Path(args.base)
    qc_dir = base / "QC"
    tmp_dir = qc_dir / "Intermediate_data"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    lg = setup_logging(qc_dir)
    threads = args.threads

    all_results: List[dict] = []

    # --------------------------------------------------------------------- #
    for sample, meta_raw in SAMPLES.items():
        lg.info("\n🟢  SAMPLE %s", sample)
        meta = meta_raw.copy()
        meta["assembly_subset"] = base / meta["assembly_subset"]

        work = tmp_dir / sample
        work.mkdir(exist_ok=True)

        sam = work / f"{sample}.sam"
        raw_bam = work / f"{sample}_raw.bam"
        filtered_bam = work / f"{sample}_filtered.bam"
        full_bam = work / f"{sample}_full.sorted.bam"
        depth_txt = work / f"{sample}_depth.txt"
        
        # New files for before/after comparison
        raw_vcf_nofilter = work / f"{sample}_raw_nofilter.vcf.gz"
        filt_vcf_nofilter = work / f"{sample}_filtered_nofilter.vcf.gz"
        raw_vcf_withfilter = work / f"{sample}_raw_withfilter.vcf.gz"
        filt_vcf_withfilter = work / f"{sample}_filtered_withfilter.vcf.gz"

        # 1 ─ mapping
        if not raw_bam.exists():
            if not map_reads(sample, meta, sam, threads, lg):
                lg.error("Mapping failed – skipping sample")
                continue
            if not run(
                f"samtools view -@ {threads} -O BAM {sam} -o {raw_bam}",
                lg,
            ):
                continue
            sam.unlink(missing_ok=True)

        # 2 ─ duplicate/secondary filtering
        if not filtered_bam.exists():
            if meta["datatype"] == "ONT":
                if not run(
                    f"samtools view -@ {threads} -F 0x900 -b {raw_bam} -o {filtered_bam}",
                    lg,
                ):
                    continue
            else:
                tmp = work / "tmp"
                steps = [
                    f"samtools sort -@ {threads} -n {raw_bam} -o {tmp}.n",
                    f"samtools fixmate -@ {threads} -m {tmp}.n {tmp}.fx",
                    f"samtools sort -@ {threads} {tmp}.fx -o {tmp}.cs",
                    f"samtools markdup -@ {threads} {tmp}.cs {tmp}.md",
                    (
                        "samtools view -@ {t} -b -F 0x400 -F 0x900 {i} -o {o}".format(
                            t=threads, i=f"{tmp}.md", o=filtered_bam
                        )
                    ),
                ]
                if not all(run(cmd, lg) for cmd in steps):
                    continue
                for f in work.glob("tmp*"):
                    f.unlink(missing_ok=True)

        # 3 ─ full‑assembly BAM
        if not full_bam.exists():
            if not run(f"samtools sort -@ {threads} {filtered_bam} -o {full_bam}", lg):
                continue
            run(f"samtools index -@ {threads} {full_bam}", lg)

        # 4 ─ Extract depth for subset contigs only
        if not depth_txt.exists():
            if not extract_depth_for_subset(full_bam, meta["assembly_subset"], depth_txt, threads, lg):
                lg.error("Failed to extract depth for %s", sample)
                continue
        
        # Check if depth file actually has content
        try:
            depths = np.loadtxt(depth_txt, usecols=2, dtype=int)
            if depths.size == 0:
                lg.error("No depth data in file for %s - check BAM file and FASTA subset", sample)
                continue
        except Exception as e:
            lg.error("Error loading depth data for %s: %s", sample, e)
            continue
            
        # Calculate metrics before filtering
        median_cov = float(np.median(depths))
        mean_cov = float(depths.mean())
        sd_cov = float(np.std(depths))
        pct_10x = calculate_percent_above_10x(depths)
        depth_cap = int(3 * median_cov)
        
        lg.info("[BEFORE] median=%.1f×  mean=%.1f×  sd=%.1f×  ≥10x=%.2f%%", 
                median_cov, mean_cov, sd_cov, pct_10x)
        
        # Store the before-filtering metrics
        before_stats = {
            "median_cov": median_cov,
            "mean_cov": mean_cov,
            "sd_cov": sd_cov,
            "pct_10x": pct_10x,
            "depth_cap": 0,  # No cap applied yet
        }
        
        # 5a ─ variant calling WITHOUT coverage cap (before filtering)
        if not raw_vcf_nofilter.exists():
            # Direct use of full BAM with subset FASTA - no subsetting needed
            mpile = (
                f"bcftools mpileup -Ou -f {meta['assembly_subset']} "
                f"--threads {threads} -a AD,DP,SP {full_bam}"
            )
            call = f"bcftools call --ploidy 2 -m -A -Oz -o {raw_vcf_nofilter}"
            if not run(f"{mpile} | {call}", lg):
                continue
            run(f"bcftools index --threads {threads} {raw_vcf_nofilter}", lg)

        if not filt_vcf_nofilter.exists():
            if not run(
                "bcftools filter -e 'QUAL<30 || INFO/DP<10' "
                f"{raw_vcf_nofilter} -Oz -o {filt_vcf_nofilter}",
                lg,
            ):
                continue
            run(f"bcftools index --threads {threads} {filt_vcf_nofilter}", lg)

        # Calculate variant stats before filtering
        het_before = int(
            subprocess.check_output(
                f"bcftools view -g het {filt_vcf_nofilter} | wc -l", shell=True
            )
        )
        sites_before = int(
            subprocess.check_output(
                f"bcftools view -H {filt_vcf_nofilter} | wc -l", shell=True
            )
        )
        het_pct_before = 0.0 if sites_before == 0 else het_before / sites_before * 100.0
        lg.info("[BEFORE] het=%d  sites=%d  het%%=%.2f", het_before, sites_before, het_pct_before)
        
        # Add variant stats to before-filtering metrics
        before_stats.update({
            "het": het_before,
            "het_pct": het_pct_before,
            "sites": sites_before,
        })
        
        # 5b ─ variant calling WITH coverage cap (after filtering)
        if not raw_vcf_withfilter.exists():
            # Direct use of full BAM with subset FASTA, but add depth cap
            mpile = (
                f"bcftools mpileup -Ou -f {meta['assembly_subset']} "
                f"-d {depth_cap} -a AD,DP,SP --threads {threads} {full_bam}"
            )
            call = f"bcftools call --ploidy 2 -m -A -Oz -o {raw_vcf_withfilter}"
            if not run(f"{mpile} | {call}", lg):
                continue
            run(f"bcftools index --threads {threads} {raw_vcf_withfilter}", lg)

        if not filt_vcf_withfilter.exists():
            if not run(
                "bcftools filter -e 'QUAL<30 || INFO/DP<10' "
                f"{raw_vcf_withfilter} -Oz -o {filt_vcf_withfilter}",
                lg,
            ):
                continue
            run(f"bcftools index --threads {threads} {filt_vcf_withfilter}", lg)

        # Calculate variant stats after filtering
        het_after = int(
            subprocess.check_output(
                f"bcftools view -g het {filt_vcf_withfilter} | wc -l", shell=True
            )
        )
        sites_after = int(
            subprocess.check_output(
                f"bcftools view -H {filt_vcf_withfilter} | wc -l", shell=True
            )
        )
        het_pct_after = 0.0 if sites_after == 0 else het_after / sites_after * 100.0
        lg.info("[AFTER] het=%d  sites=%d  het%%=%.2f", het_after, sites_after, het_pct_after)
        
        # Calculate and store metrics with filtering
        # For filtered metrics, we need to create a capped depth array
        capped_depths = np.minimum(depths, depth_cap)
        mean_cov_after = float(capped_depths.mean())
        # Median is unchanged by capping at 3x median
        pct_10x_after = calculate_percent_above_10x(capped_depths)
        
        lg.info("[AFTER] median=%.1f×  mean=%.1f×  ≥10x=%.2f%%", 
                median_cov, mean_cov_after, pct_10x_after)
        
        # Store the after-filtering metrics
        after_stats = {
            "median_cov": median_cov,  # Median stays the same
            "mean_cov": mean_cov_after,
            "sd_cov": sd_cov,          # SD is from original calculation
            "pct_10x": pct_10x_after,
            "depth_cap": depth_cap,
            "het": het_after,
            "het_pct": het_pct_after,
            "sites": sites_after,
        }
        
        # Store results for report
        all_results.append({
            "sample": sample,
            "before_filtering_stats": before_stats,
            "after_filtering_stats": after_stats,
        })

    # ----------------------------------------------------------------- #
    # Write the QC reports (before and after filtering)
    write_qc_report(all_results, qc_dir, lg)


if __name__ == "__main__":
    main()