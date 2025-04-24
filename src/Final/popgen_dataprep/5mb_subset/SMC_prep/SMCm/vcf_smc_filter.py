#!/usr/bin/env python3
# vcf_smc_filter.py – filter VCF for SMCm / eSMC2
#   • keeps only biallelic SNPs that show an ALT allele in the sample
#   • writes bgzip‑compressed output + .tbi index
#   • no quality / depth filters unless you ask for them

import argparse, gzip, io
import pysam
from pathlib import Path

def load_contigs(fai):
    with open(fai) as fh:
        return {ln.split('\t',1)[0] for ln in fh}

def wanted(f, contigs):
    chrom, ref, alt, gt = f[0], f[3], f[4], f[9].split(':')[0]
    if contigs and chrom not in contigs:
        return False
    if alt == '.' or ',' in alt or len(ref)!=1 or len(alt)!=1:
        return False            # monomorphic, multiallelic, or not SNP
    if gt.replace('|','/').startswith('0/0'):
        return False            # homo‑ref
    return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-i','--vcf_in',  required=True)
    ap.add_argument('-o','--vcf_out', required=True)
    ap.add_argument('--sample', default='SAMPLE')
    ap.add_argument('--fasta',  help='FASTA (needs .fai) to restrict contigs')
    args = ap.parse_args()

    contigs = load_contigs(args.fasta + '.fai') if args.fasta else None
    out_path = Path(args.vcf_out).with_suffix('.gz')

    with pysam.BGZFile(out_path, 'w') as bgzf:              # binary mode
        out = io.TextIOWrapper(bgzf, encoding='utf-8')      # text wrapper

        opener = gzip.open if args.vcf_in.endswith('.gz') else open
        with opener(args.vcf_in, 'rt') as inp:
            for ln in inp:
                if ln.startswith('##'):
                    out.write(ln)
                elif ln.startswith('#CHROM'):
                    header = ln.rstrip('\n').split('\t')
                    header[9] = args.sample
                    out.write('\t'.join(header) + '\n')
                else:
                    f = ln.rstrip('\n').split('\t')
                    if wanted(f, contigs):
                        out.write(ln)

        out.flush()  # make sure wrapper pushes data to BGZFile

    pysam.tabix_index(str(out_path), preset='vcf', force=True)

if __name__ == '__main__':
    main()
