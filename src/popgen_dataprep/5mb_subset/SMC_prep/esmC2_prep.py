#!/usr/bin/env python3
import os
import subprocess

# Directories (adjust these paths as needed)
REFERENCE_DIR = "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/dataprep/5mb_subset/reference"
VCF_DIR = "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/dataprep/5mb_subset/VCFs"
ESMC2_OUT_DIR = "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/dataprep/5mb_subset/eSMC2_preparation"

# File names
reference_fasta = os.path.join(REFERENCE_DIR, "reference.fasta")
# Assume HaplotypeCaller produced these per-sample GVCFs for a given chromosome (e.g., Chr2)
gvcf_HT = os.path.join(VCF_DIR, "HT.Chr2.g.vcf.gz")
gvcf_HT2 = os.path.join(VCF_DIR, "HT2.Chr2.g.vcf.gz")
gvcf_BB = os.path.join(VCF_DIR, "BB.Chr2.g.vcf.gz")

# Define cohorts:
# Group1: HT and HT2 combined
# Group2: BB (single sample)
cohorts = {
    "HT_combined": [("HT", gvcf_HT), ("HT2", gvcf_HT2)],
    "BB": [("BB", gvcf_BB)]
}

# Ensure output base directory exists
os.makedirs(ESMC2_OUT_DIR, exist_ok=True)

def run(cmd, workdir=None):
    print(f"🔧 Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True, cwd=workdir)

# Loop over cohorts
for cohort_name, sample_list in cohorts.items():
    print(f"\n🚀 Processing cohort: {cohort_name}")
    cohort_dir = os.path.join(ESMC2_OUT_DIR, cohort_name)
    os.makedirs(cohort_dir, exist_ok=True)

    # 1. Combine GVCFs using GATK CombineGVCFs
    combined_gvcf = os.path.join(cohort_dir, f"{cohort_name}.Chr2.g.vcf.gz")
    combine_cmd = f"gatk CombineGVCFs -R {reference_fasta}"
    for sample, gvcf in sample_list:
        combine_cmd += f" --variant {gvcf}"
    combine_cmd += f" -O {combined_gvcf}"
    run(combine_cmd)

    # 2. Genotype GVCFs with GenotypeGVCFs
    allsites_vcf = os.path.join(cohort_dir, f"{cohort_name}.allsites.geno.vcf.gz")
    genotype_cmd = f'gatk --java-options "-Xmx4g" GenotypeGVCFs -R {reference_fasta} -V {combined_gvcf} --include-non-variant-sites -O {allsites_vcf}'
    run(genotype_cmd)

    # 3. Create SNP-only VCF (final.filtered.vcf.gz)
    final_filtered_vcf = os.path.join(cohort_dir, f"{cohort_name}.final.filtered.vcf.gz")
    snp_view_cmd = (
        f'bcftools view {allsites_vcf} '
        f'--genotype ^miss '
        f'--apply-filters .,PASS '
        f'--include \'TYPE="snp"\' '
        f'-Oz -o {final_filtered_vcf}'
    )
    run(snp_view_cmd)

    # 4. Create callable sites VCF (mask VCF)
    allsites_filtered_vcf = os.path.join(cohort_dir, f"{cohort_name}.final.allsites.geno.filtered.vcf.gz")
    mask_vcf_cmd = (
        f'bcftools view {allsites_vcf} '
        f'--genotype ^miss '
        f'--apply-filters .,PASS '
        f'--include \'TYPE="snp" && INFO/DP > 30 || TYPE="ref" && INFO/DP > 30\' '
        f'-Oz -o {allsites_filtered_vcf}'
    )
    run(mask_vcf_cmd)

    # 5. Convert callable sites VCF to BED and merge intervals to create the mask file
    mask_bed = os.path.join(cohort_dir, f"{cohort_name}.final.mask.bed")
    # Convert to BED using vcf2bed (vcfutils.pl is assumed available in your PATH)
    vcf2bed_cmd = f"zcat {allsites_filtered_vcf} | vcf2bed > {mask_bed}"
    run(vcf2bed_cmd)
    # Merge adjacent intervals with bedtools
    merged_mask_bed = os.path.join(cohort_dir, f"{cohort_name}.final.mask.merged.bed")
    merge_cmd = f"bedtools merge -i {mask_bed} > {merged_mask_bed}"
    run(merge_cmd)
    # Gzip the final mask BED
    run(f"gzip -f {merged_mask_bed}")
    final_mask = merged_mask_bed + ".gz"

    # 6. Split the final.filtered VCF into individual sample VCFs.
    # Get the list of sample names from the VCF:
    sample_list_cmd = f"bcftools query -l {final_filtered_vcf}"
    sample_names = subprocess.check_output(sample_list_cmd, shell=True).decode().strip().split()
    per_sample_vcfs = []
    for s in sample_names:
        out_vcf = os.path.join(cohort_dir, f"{cohort_name}.{s}.vcf.gz")
        extract_cmd = f"bcftools view --samples {s} {final_filtered_vcf} --output-type z --output-file {out_vcf}"
        run(extract_cmd)
        per_sample_vcfs.append(out_vcf)

    # 7. Download (if needed) and run generate_multihetsep.py to generate the Multihetsep file.
    gen_mhs_script = os.path.join(cohort_dir, "generate_multihetsep.py")
    if not os.path.exists(gen_mhs_script):
        run(f"wget -O {gen_mhs_script} https://raw.githubusercontent.com/stschiff/msmc-tools/master/generate_multihetsep.py")
        run(f"chmod u+x {gen_mhs_script}")
    multihetsep_out = os.path.join(cohort_dir, f"{cohort_name}.Chr02.mhs")
    # Build the command using the mask and the per-sample VCFs
    mhs_cmd = f"python3 {gen_mhs_script} --mask={final_mask} " + " ".join(per_sample_vcfs) + f" > {multihetsep_out}"
    run(mhs_cmd)

    print(f"✅ Finished eSMC2 preparation for cohort: {cohort_name}")

print("\n✅ eSMC2 preparation completed for all cohorts.")
