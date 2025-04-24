#!/usr/bin/env python3
"""
bed_to_smc.py  –  Convert ONT bedMethyl pile-ups to the SMCm per-site table
────────────────────────────────────────────────────────────────────────────
* Reads a (possibly gzipped) bedMethyl file from Dorado/nf-core nanopore
* Retrieves the surrounding trinucleotide from the reference FASTA
* Applies a one-sided binomial test using the user-supplied false-positive
  methyl-call rate (p0) and Benjamini–Hochberg FDR (α = 0.01)
* Outputs the 10-column table expected by SMCm / the Tellier-lab papers:
      seqnames  start  strand  context  counts.methylated  counts.total
      posteriorMax  status  rc.meth.lvl  context.trinucleotide
────────────────────────────────────────────────────────────────────────────
author :  Max Borgmann & contributors  •  MIT licence  •  2025-04-23
"""

import argparse
import csv
import gzip

import numpy as np  # noqa: F401
from Bio import SeqIO
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests

# ───────────────────────── helper functions ──────────────────────────────
def get_trinuc(seq_dict, chrom, pos0, strand):
    """Return the 3-mer (5'→3') centred on cytosine at pos0 (0-based)."""
    seq = seq_dict[chrom].seq
    if strand == '+':
        tri = seq[pos0 : pos0 + 3]
    else:
        if pos0 < 2:
            return None
        tri = seq[pos0 - 2 : pos0 + 1].reverse_complement()
    tri = str(tri).upper()
    return None if 'N' in tri else tri


def jeffreys_shrink(k, n):
    """Bayesian shrinkage of mC fraction with Jeffreys(½,½) prior."""
    return (k + 0.5) / (n + 1.0)


# ─────────────────────────────── main ────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description="Convert bedMethyl to SMCm table")
    p.add_argument("-i", "--bed", required=True, help="input .bed[.gz]")
    p.add_argument("-f", "--fasta", required=True, help="reference FASTA")
    p.add_argument("-o", "--out", required=True, help="output .txt")
    p.add_argument("--min_cov", type=int, default=5,
                   help="discard sites with coverage <N (default: 1)")
    p.add_argument("--alpha", type=float, default=0.01,
                   help="FDR threshold (default: 0.01)")
    p.add_argument("--p0", type=float, default=0.005,
                   help="false-positive methyl-call rate for binomial H0 "
                        "(≈ non-conversion rate; default: 0.005)")
    p.add_argument("--contigs", nargs="*", default=None,
                   help="optional whitelist of contigs to keep")
    args = p.parse_args()

    # Load reference genome
    seqs = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    # Pass 1 – parse bed file and collect candidate rows
    rows = []
    pvals = []
    with gzip.open(args.bed, "rt") if args.bed.endswith(".gz") else open(args.bed) as bed:
        for ln in bed:
            if ln.startswith("#") or not ln.strip():
                continue
            c = ln.rstrip("\n").split("\t")
            chrom = c[0]
            if args.contigs and chrom not in args.contigs:
                continue
            if chrom not in seqs:
                continue

            start0 = int(c[1])
            strand = c[5]
            total = int(c[9])
            if total < args.min_cov:
                continue
            methylated = int(c[12])
            context = c[3].split(",")[1]
            tri = get_trinuc(seqs, chrom, start0, strand)
            if tri is None:
                continue

            # one-sided (greater) binomial P-value
            pval = binom.sf(methylated - 1, total, args.p0)
            pvals.append(pval)
            rows.append([chrom, start0 + 1, strand, context,
                         methylated, total, tri])

    if not rows:
        raise SystemExit("No sites passed the filters – nothing written.")

    # FDR correction
    rej, qvals, _, _ = multipletests(pvals, alpha=args.alpha, method="fdr_bh")

    # Pass 2 – write table
    with open(args.out, "w", newline="") as fo:
        wr = csv.writer(fo, delimiter="\t")
        wr.writerow(["seqnames", "start", "strand", "context",
                     "counts.methylated", "counts.total",
                     "posteriorMax", "status", "rc.meth.lvl",
                     "context.trinucleotide"])

        for rec, q, is_meth in zip(rows, qvals, rej):
            chrom, pos1, strand, ctx, meth, tot, tri = rec
            post_max = 1.0 - q        # upper posterior ≈ 1-q
            status = "M" if is_meth else "U"
            rc_lvl = jeffreys_shrink(meth, tot)
            wr.writerow([chrom, pos1, strand, ctx,
                         meth, tot,
                         f"{post_max:.4f}", status,
                         f"{rc_lvl:.4f}", tri])


if __name__ == "__main__":
    main()
