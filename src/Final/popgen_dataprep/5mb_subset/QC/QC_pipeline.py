#!/usr/bin/env python3
"""
Simplified Coverage-&-Variant QC pipeline
─────────────────────────────────────────
• Maps reads against the full assembly for every sample
• Calls variants only on a 5 Mb subset (via --regions-file)
• Computes coverage / variant metrics before & after a 3×-median depth cap
• Memory-aware: chunked I/O, adaptive thread counts, soft RSS limits
"""

from __future__ import annotations

import argparse
import gc
import logging
import random
import resource
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

# ────────────────────────────────
# psutil (optional, graceful fallback)
# ────────────────────────────────
try:
    import psutil
except ImportError:                     # minimal stub
    class _VM: total = 100 << 30; available = 80 << 30; percent = 20
    class _Proc:
        def __init__(self, pid): ...
        def memory_info(self):
            class _MI: rss = 0
            return _MI()
    class _Ps:
        Process = staticmethod(lambda pid: _Proc(pid))
        virtual_memory = staticmethod(lambda: _VM)
    psutil = _Ps()                      # type: ignore


# ────────────────────────────────
# CLI & logging helpers
# ────────────────────────────────
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Depth/SNP QC on subset contigs")
    p.add_argument(
        "-o", "--base",
        default="/data/proj2/home/students/m.borgmann/Master_thesis/"
                "data/processed/Popgen_analysis/dataprep/5mb_subset",
        help="Base output directory (same used by FASTA-subset step)",
    )
    p.add_argument("-t", "--threads", type=int, default=20, help="CPU threads")
    p.add_argument(
        "-m", "--max-memory", type=int, default=80,
        help="Maximum memory (GB) for external tools",
    )
    p.add_argument("--sequential", action="store_true",
                   help="Process samples sequentially to save memory")
    p.add_argument("--chunk-size", type=int, default=1_000_000,
                   help="Lines per chunk when parsing depth")
    return p.parse_args()


def setup_logging(outdir: Path) -> logging.Logger:
    outdir.mkdir(parents=True, exist_ok=True)
    lg = logging.getLogger("qc")
    lg.setLevel(logging.DEBUG)
    lg.handlers.clear()

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(outdir / "pipeline.log"); fh.setFormatter(fmt); fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler();                        ch.setFormatter(fmt); ch.setLevel(logging.INFO)
    lg.addHandler(fh); lg.addHandler(ch)
    lg.info("Logging → %s", outdir / "pipeline.log")
    return lg


# ────────────────────────────────
# generic helpers
# ────────────────────────────────
def set_soft_limit(gb: int) -> None:
    try:
        resource.setrlimit(resource.RLIMIT_AS, (gb << 30, resource.RLIM_INFINITY))
    except (ValueError, resource.error):
        pass


def reset_soft_limit() -> None:
    try:
        resource.setrlimit(resource.RLIMIT_AS,
                           (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except (ValueError, resource.error):
        pass


def run(cmd: str, lg: logging.Logger,
        cwd: Path | None = None,
        max_memory_gb: Optional[int] = None) -> bool:
    """Run a shell command with optional soft memory limit."""
    lg.info("↪ %s", cmd)
    limited = False
    if max_memory_gb and any(tok in cmd for tok in ("bcftools", "minimap2", "bwa-mem2")):
        set_soft_limit(int(max_memory_gb * 0.9)); limited = True
    try:
        subprocess.run(cmd, shell=True, cwd=cwd, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return True
    except subprocess.CalledProcessError as e:
        lg.error("command failed (%d)", e.returncode)
        lg.debug(e.stdout); lg.debug(e.stderr)
        return False
    finally:
        if limited:
            reset_soft_limit()
            lg.debug("Reset memory limit")


# ────────────────────────────────
#  Reference / region sanity checks
# ────────────────────────────────
def check_region_compatibility(ref_fa: Path, regions: Path, lg: logging.Logger) -> bool:
    """Return False if any contig in regions is missing from ref_fa.fai"""
    fai_ctg = {l.split('\t', 1)[0] for l in open(f"{ref_fa}.fai")}
    for line in open(regions):
        ctg = line.rstrip().split('\t', 1)[0]
        if ctg not in fai_ctg:
            lg.error("Contig %s present in regions-file but missing from %s", ctg, ref_fa)
            return False
    return True


def check_bam_reference_compatibility(bam: Path, ref_fa: Path, lg: logging.Logger) -> bool:
    fai_ctg = {l.split('\t', 1)[0] for l in open(f"{ref_fa}.fai")}
    bam_hdr = subprocess.check_output(f"samtools view -H {bam}", shell=True, text=True)
    for l in bam_hdr.splitlines():
        if l.startswith("@SQ"):
            ctg = l.split("\t")[1][3:]            # SN:<name>
            if ctg not in fai_ctg:
                lg.error("BAM refers to contig %s not present in %s", ctg, ref_fa)
                return False
    return True


# ────────────────────────────────
#  get_opt_threads – adaptive thread helper
# ────────────────────────────────
def get_opt_threads(tool: str, max_thr: int, max_mem_gb: int) -> int:
    if tool == "bwa-mem2": return min(max_thr, max(1, max_mem_gb // 4))
    if tool == "minimap2": return min(max_thr, max(1, max_mem_gb // 2))
    if tool == "bcftools": return min(max_thr, max(1, max_mem_gb // 8))
    return max_thr


# ────────────────────────────────
#  Variant-calling wrapper
# ────────────────────────────────
def _ref_from_mpileup(cmd: str) -> Path:
    parts = shlex.split(cmd)
    try:
        return Path(parts[parts.index("-f") + 1])
    except (ValueError, IndexError):
        raise RuntimeError("Could not parse reference FASTA from mpileup cmd")


def run_variant_pipe(mpile_cmd: str, call_cmd: str, lg: logging.Logger,
                     regions_file: Optional[Path] = None,
                     max_memory_gb: Optional[int] = None) -> bool:
    ref_fa = _ref_from_mpileup(mpile_cmd)

    if regions_file:
        if not check_region_compatibility(ref_fa, regions_file, lg):
            lg.error("Region / reference mismatch – aborting variant call")
            return False
        mpile_cmd = mpile_cmd.replace("mpileup ",
                                      f"mpileup --regions-file {regions_file} ", 1)

    full_cmd = f"{mpile_cmd} | {call_cmd}"
    lg.info("↪ %s", full_cmd)

    if max_memory_gb:
        set_soft_limit(int(max_memory_gb * 0.9))
    try:
        proc = subprocess.run(full_cmd, shell=True, text=True,
                              capture_output=True)
        if proc.returncode != 0:
            lg.error("variant pipeline exited with %d", proc.returncode)
            lg.error("stderr:\n%s", proc.stderr.strip() or "<EMPTY>")
            return False
        return True
    finally:
        if max_memory_gb:
            reset_soft_limit()
            lg.debug("Reset memory limit")


# ────────────────────────────────
#  Reliable heterozygosity counter
# ────────────────────────────────
def count_hets_and_sites(vcf_gz: Path, lg: logging.Logger) -> tuple[int, int]:
    cmd_het   = f"bcftools view -H -g het  -v snps {vcf_gz} | wc -l"
    cmd_sites = f"bcftools view -H -v snps {vcf_gz} | wc -l"
    het   = int(subprocess.check_output(cmd_het,   shell=True))
    sites = int(subprocess.check_output(cmd_sites, shell=True))
    lg.debug("HET=%d  SITES=%d  %.2f%%", het, sites, het*100/sites if sites else 0)
    return het, sites


# ────────────────────────────────
# FASTA utilities
# ────────────────────────────────
def get_fasta_length(fasta_path: Path, lg: logging.Logger) -> int:
    """Get the total length of a FASTA file by summing contig lengths from .fai file."""
    fai_path = f"{fasta_path}.fai"
    if not Path(fai_path).exists():
        lg.info("Creating FASTA index for %s", fasta_path)
        subprocess.run(f"samtools faidx {fasta_path}", shell=True, check=True)
    
    total_length = 0
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                total_length += int(parts[1])
    return total_length


# ────────────────────────────────
# mapping helper
# ────────────────────────────────
def map_reads(sample: str, meta: Dict, sam_path: Path,
              threads: int, max_gb: int, lg: logging.Logger) -> bool:
    if sam_path.exists() and sam_path.stat().st_size > 0:
        lg.info("SAM exists – skip mapping"); return True
    if sam_path.exists(): sam_path.unlink()

    if meta["mapper"] == "minimap2":
        thr = get_opt_threads("minimap2", threads, max_gb)
        mem_per_thr_mb = int((max_gb * 1024) / (thr * 1.5))
        cmd = (f"minimap2 -t {thr} -K {mem_per_thr_mb}M -ax map-ont -L "
               f"{meta['assembly_full']} {meta['reads']} > {sam_path}")
        return run(cmd, lg, max_memory_gb=max_gb)

    if meta["mapper"] == "bwa-mem2":
        thr = get_opt_threads("bwa-mem2", threads, max_gb)
        if not Path(f"{meta['assembly_full']}.bwt.2bit.64").exists():
            lg.info("Building bwa-mem2 index…")
            if not run(f"bwa-mem2 index {meta['assembly_full']}", lg, max_memory_gb=max_gb):
                return False
        cmd = (f"bwa-mem2 mem -M -t {thr} -K 10000000 {meta['assembly_full']} "
               f"{meta['reads_1']} {meta['reads_2']} > {sam_path}")
        return run(cmd, lg, max_memory_gb=max_gb)

    lg.error("Unknown mapper for %s", sample); return False


def calculate_tstv(vcf_gz: Path, lg: logging.Logger) -> float:
    """
    Calculate the transition/transversion ratio (Ts/Tv) for SNPs in the VCF.
    """
    cmd = ["bcftools", "view", "-H", "-v", "snps", str(vcf_gz)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    ts = 0
    tv = 0
    for line in p.stdout:
        fields = line.split('\t')
        ref, alt = fields[3], fields[4]
        # transitions
        if (ref, alt) in (("A","G"),("G","A"),("C","T"),("T","C")):
            ts += 1
        else:
            tv += 1
    p.wait()
    if tv == 0:
        return float("nan")
    return ts / tv


# ────────────────────────────────
# depth utilities
# ────────────────────────────────
def process_depths_in_chunks(depth_file: Path, lg: logging.Logger,
                             chunk_size: int = 1_000_000,
                             max_samples: int = 1_000_000) -> dict:
    if not depth_file.exists():
        lg.error("Depth file missing: %s", depth_file)
        return {"total_positions": 0}

    tot, s, sq, ge10, sample = 0, 0, 0, 0, []
    with open(depth_file) as fh:
        while (lines := [fh.readline().strip() for _ in range(chunk_size)]):
            depths = [int(l.split()[2]) for l in lines if l]
            if not depths: break
            n = len(depths)
            tot += n; s += sum(depths); sq += sum(d*d for d in depths)
            ge10 += sum(d >= 10 for d in depths)
            if len(sample) < max_samples:
                need = min(max_samples - len(sample), n)
                sample.extend(random.sample(depths, need) if need < n else depths)

    if tot == 0: return {"total_positions": 0}
    mean = s / tot
    med  = float(np.median(sample))
    var  = sq / tot - mean*mean
    sd   = np.sqrt(var) if var > 0 else 0.0
    return {"mean": mean, "median": med, "sd": sd,
            "pct_10x": ge10 / tot * 100, "total_positions": tot}


def extract_depth_for_subset(full_bam: Path, subset_fa: Path, depth_txt: Path,
                             threads: int, max_gb: int, lg: logging.Logger) -> bool:
    if depth_txt.exists() and depth_txt.stat().st_size > 0:
        lg.info("Depth already extracted"); return True

    run(f"samtools faidx {subset_fa}", lg, max_memory_gb=max_gb)
    contigs = [l.split()[0] for l in open(f"{subset_fa}.fai")]
    run(f"samtools index -@ {threads} {full_bam}", lg, max_memory_gb=max_gb)

    tmp_dir = depth_txt.parent / "tmp_depth"; tmp_dir.mkdir(exist_ok=True)
    for i, c in enumerate(contigs, 1):
        tmp = tmp_dir / f"{c}.depth"
        if tmp.exists() and tmp.stat().st_size > 0:
            continue
        lg.info("Depth %s (%d/%d)", c, i, len(contigs))
        if not run(f"samtools depth -@ {threads} -a {full_bam} -r {c} > {tmp}",
                   lg, max_memory_gb=max_gb):
            return False

    with open(depth_txt, "w") as out:
        for c in contigs:
            tmp = tmp_dir / f"{c}.depth"
            with open(tmp) as fh:
                out.writelines(fh)
            tmp.unlink(missing_ok=True)
    tmp_dir.rmdir()
    lg.info("Depth written (%.1f MB)", depth_txt.stat().st_size / 1e6)
    return True


def calculate_filtered_depth_stats(depth_file: Path, cap: int, lg: logging.Logger,
                                   chunk_size: int = 1_000_000) -> dict:
    tot, s, ge10 = 0, 0, 0
    with open(depth_file) as fh:
        while (lines := [fh.readline().strip() for _ in range(chunk_size)]):
            depths = [min(int(l.split()[2]), cap) for l in lines if l]
            if not depths: break
            tot += len(depths); s += sum(depths); ge10 += sum(d >= 10 for d in depths)
    if not tot: return {"total_positions": 0}
    return {"mean": s / tot, "pct_10x": ge10 / tot * 100, "total_positions": tot}



# ────────────────────────────────
# QC report writer
# ────────────────────────────────
def write_qc_report(res: List[dict], qc_dir: Path, lg: logging.Logger) -> None:
    if not res: return

    for phase, tag in [("before_filtering_stats","before"),
                       ("after_filtering_stats", "after")]:
        hdr = (
            "Sample".ljust(12) +
            "Mean".rjust(10) + "Median".rjust(10) + "SD".rjust(10) +
            "Cap".rjust(8)  + "%≥10x".rjust(8)   +
            "Het".rjust(10) + "Het%".rjust(8)    +
            "Ts/Tv".rjust(8)
        )
        lines = [hdr, "-"*len(hdr)]
        for r in res:
            s = r[phase]
            lines.append(
                r["sample"].ljust(12) +
                f"{s['mean_cov']:.2f}".rjust(10)   +
                f"{s['median_cov']:.2f}".rjust(10) +
                f"{s['sd_cov']:.2f}".rjust(10)     +
                f"{s['depth_cap']}".rjust(8)       +
                f"{s['pct_10x']:.2f}".rjust(8)     +
                f"{s['het']}".rjust(10)            +
                f"{s['het_pct']:.2f}".rjust(8)     +  # Use higher precision for small percentages
                f"{s['tstv']:.2f}".rjust(8)
            )
        (qc_dir / f"QC_{tag}.txt").write_text("\n".join(lines), encoding="utf8")
        lg.info("Wrote QC_%s.txt", tag)




# ────────────────────────────────
# sample configuration
# ────────────────────────────────
BASE_FASTA = Path("/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies")

SAMPLES: Dict[str, Dict] = {
    "HT": {
        "assembly_full": BASE_FASTA / "Hexaplex_ONT_FLYE/Polished_assembly/Final/Hexaplex_assembly_ONT_FLYE_polished.fasta",
        "assembly_subset": "Input_fastas/HT/HT_assembly_5mb_subset.fasta",
        "reads": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/Hexaplex_trunculus_nanopore/20241217_DNA_Tellier_HT_31_PB/20241217_1155_2A_PAY72199_a0bfda36/fastq_pass/merged_HT_reads.fastq",
        "mapper": "minimap2", 
        "datatype": "ONT",
    },
    "HT2": {
        "assembly_full": BASE_FASTA / "Hexaplex_Captus/02_assemblies/SRR28865916__captus-asm/01_assembly/assembly.fasta",
        "assembly_subset": "Input_fastas/HT2/HT2_assembly_5mb_subset.fasta",
        "reads_1": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/SRR28865916/SRR28865916_R1.fastq",
        "reads_2": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/SRR28865916/SRR28865916_R2.fastq",
        "mapper": "bwa-mem2", 
        "datatype": "ILLUMINA",
    },
    "BB": {
        "assembly_full": BASE_FASTA / "Bolinus_Captus/SRR28863561__captus-asm/01_assembly/assembly.fasta",
        "assembly_subset": "Input_fastas/BB/BB_assembly_5mb_subset.fasta",
        "reads_1": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106534_Bolinus_brandaris/SRR28863561/SRR28863561_R1.fastq",
        "reads_2": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106534_Bolinus_brandaris/SRR28863561/SRR28863561_R2.fastq",
        "mapper": "bwa-mem2", 
        "datatype": "ILLUMINA",
    },
    "HT2_ON_HTref": {
        "assembly_full": BASE_FASTA / "Hexaplex_ONT_FLYE/Polished_assembly/Final/Hexaplex_assembly_ONT_FLYE_polished.fasta",
        "assembly_subset": "Input_fastas/HT/HT_assembly_5mb_subset.fasta",
        "reads_1": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/SRR28865916/SRR28865916_R1.fastq",
        "reads_2": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/SRR28865916/SRR28865916_R2.fastq",
        "mapper": "bwa-mem2", 
        "datatype": "ILLUMINA",
    },
}


# ────────────────────────────────
# main
# ────────────────────────────────
def main() -> None:
    args = parse_args()
    base = Path(args.base)
    qc_dir = base / "QC"
    tmp_dir = qc_dir / "Intermediate_data"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    lg = setup_logging(qc_dir)

    threads, max_gb, chunk = args.threads, args.max_memory, args.chunk_size
    lg.info("System memory: %.1f GB total, %.1f GB available",
            psutil.virtual_memory().total / 1e9,
            psutil.virtual_memory().available / 1e9)

    results: List[dict] = []

    for sample, meta0 in SAMPLES.items():
        lg.info("\n🟢  SAMPLE %s", sample)
        meta = meta0.copy()
        meta["assembly_subset"] = base / meta["assembly_subset"]

        # Get subset FASTA length for proper heterozygosity calculation
        subset_length = get_fasta_length(meta["assembly_subset"], lg)
        lg.info("Subset FASTA total length: %d bp", subset_length)

        work = tmp_dir / sample
        work.mkdir(exist_ok=True)
        sam          = work / f"{sample}.sam"
        raw_bam      = work / f"{sample}_raw.bam"
        filtered_bam = work / f"{sample}_filtered.bam"
        full_bam     = work / f"{sample}_full.sorted.bam"
        depth_txt    = work / f"{sample}_depth.txt"

        raw_vcf_uncapped   = work / f"{sample}_raw_uncapped.vcf.gz"
        filt_vcf_uncapped  = work / f"{sample}_filtered_uncapped.vcf.gz"
        raw_vcf_capped     = work / f"{sample}_raw_capped.vcf.gz"
        filt_vcf_capped    = work / f"{sample}_filtered_capped.vcf.gz"

        subset_list = work / f"{sample}_subset_contigs.txt"
        if not subset_list.exists():
            run(f"samtools faidx {meta['assembly_subset']}", lg, max_memory_gb=max_gb)
            subprocess.run(f"cut -f1 {meta['assembly_subset']}.fai > {subset_list}",
                           shell=True, check=True)

        # 1 ─ mapping
        if not raw_bam.exists():
            if not map_reads(sample, meta, sam, threads, max_gb, lg):
                lg.error("mapping failed")
                continue
            run(f"samtools view -@ {threads} -O BAM {sam} -o {raw_bam}", lg)
            sam.unlink(missing_ok=True)

        # 2 ─ duplicate / secondary filtering
        if not filtered_bam.exists():
            if meta["datatype"] == "ONT":
                run(f"samtools view -@ {threads} -F 0x900 -b {raw_bam} -o {filtered_bam}", lg)
            else:
                tmp = work / "tmp"
                steps = [
                    f"samtools sort -@ {threads} -n {raw_bam} -o {tmp}.n",
                    f"samtools fixmate -@ {threads} -m {tmp}.n {tmp}.fx",
                    f"samtools sort -@ {threads} {tmp}.fx -o {tmp}.cs",
                    f"samtools markdup -@ {threads} {tmp}.cs {tmp}.md",
                    f"samtools view -@ {threads} -b -F 0x400 -F 0x900 {tmp}.md -o {filtered_bam}",
                ]
                for cmd in steps:
                    if not run(cmd, lg):
                        break
                for f in work.glob("tmp*"):
                    f.unlink(missing_ok=True)

        # 3 ─ coordinate-sort full BAM
        if not full_bam.exists():
            run(f"samtools sort -@ {threads} {filtered_bam} -o {full_bam}", lg)
            run(f"samtools index -@ {threads} {full_bam}", lg)

        # 4 ─ depth extraction
        if not depth_txt.exists():
            if not extract_depth_for_subset(full_bam, meta["assembly_subset"],
                                            depth_txt, threads, max_gb, lg):
                continue

        # 5 ─ depth stats & cap
        stats = process_depths_in_chunks(depth_txt, lg, chunk_size=chunk)
        if stats["total_positions"] == 0:
            lg.error("no depth")
            continue
        cap = int(3 * stats["median"])
        lg.info("Depth cap = %d×", cap)

        bcft_thr = get_opt_threads("bcftools", threads, max_gb)

        # 6a ─ variants BEFORE cap
        if not raw_vcf_uncapped.exists():
            mpile = (f"bcftools mpileup -Ou -f {meta['assembly_full']} "
                     f"-d 1000 -a AD,DP,SP --threads {bcft_thr} {full_bam}")
            call  = f"bcftools call -m -A --ploidy 2 -Oz -o {raw_vcf_uncapped}"
            if not run_variant_pipe(mpile, call, lg,
                                    regions_file=subset_list,
                                    max_memory_gb=max_gb):
                continue
            run(f"bcftools index --threads {bcft_thr} {raw_vcf_uncapped}", lg)

        if not filt_vcf_uncapped.exists():
            run(f"bcftools view -m2 -M2 -v snps {raw_vcf_uncapped} | "
                f"bcftools filter -e 'QUAL<30 || INFO/DP<10' -Oz -o {filt_vcf_uncapped}", lg)
            run(f"bcftools index --threads {bcft_thr} {filt_vcf_uncapped}", lg)

        # heterozygosity & Ts/Tv before cap
        het_before, sites_before = count_hets_and_sites(filt_vcf_uncapped, lg)
        # Calculate heterozygosity as percentage of subset FASTA size
        het_pct_before = (het_before / subset_length * 100) if subset_length else 0
        tstv_before   = calculate_tstv(filt_vcf_uncapped, lg)

        # 6b ─ variants AFTER cap
        if not raw_vcf_capped.exists():
            mpile = (f"bcftools mpileup -Ou -f {meta['assembly_full']} "
                     f"-d {cap} -a AD,DP,SP --threads {bcft_thr} {full_bam}")
            call  = f"bcftools call -m -A --ploidy 2 -Oz -o {raw_vcf_capped}"
            if not run_variant_pipe(mpile, call, lg,
                                    regions_file=subset_list,
                                    max_memory_gb=max_gb):
                continue
            run(f"bcftools index --threads {bcft_thr} {raw_vcf_capped}", lg)

        if not filt_vcf_capped.exists():
            run(f"bcftools view -m2 -M2 -v snps {raw_vcf_capped} | "
                f"bcftools filter -e 'QUAL<30 || INFO/DP<10' -Oz -o {filt_vcf_capped}", lg)
            run(f"bcftools index --threads {bcft_thr} {filt_vcf_capped}", lg)

        # heterozygosity & Ts/Tv after cap
        het_after, sites_after = count_hets_and_sites(filt_vcf_capped, lg)
        # Calculate heterozygosity as percentage of subset FASTA size
        het_pct_after = (het_after / subset_length * 100) if subset_length else 0
        tstv_after   = calculate_tstv(filt_vcf_capped, lg)

        # 7 ─ filtered depth metrics
        filt_stats  = calculate_filtered_depth_stats(depth_txt, cap, lg, chunk)
        mean_after  = filt_stats["mean"]
        pct10_after = filt_stats["pct_10x"]

        # collect all results
        results.append({
            "sample": sample,
            "before_filtering_stats": {
                "median_cov":  stats["median"],
                "mean_cov":    stats["mean"],
                "sd_cov":      stats["sd"],
                "pct_10x":     stats["pct_10x"],
                "depth_cap":   0,
                "het":         het_before,
                "het_pct":     het_pct_before,
                "tstv":        tstv_before,
            },
            "after_filtering_stats": {
                "median_cov":  stats["median"],
                "mean_cov":    mean_after,
                "sd_cov":      stats["sd"],
                "pct_10x":     pct10_after,
                "depth_cap":   cap,
                "het":         het_after,
                "het_pct":     het_pct_after,
                "tstv":        tstv_after,
            },
        })

        if args.sequential:
            lg.info("GC…")
            gc.collect()

    write_qc_report(results, qc_dir, lg)
    lg.info("Pipeline finished ✔")



if __name__ == "__main__":
    main()