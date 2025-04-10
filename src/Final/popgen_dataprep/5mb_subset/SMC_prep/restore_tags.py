import pysam
import sys

unmapped_bam_path = sys.argv[1]  # Dorado BAM (unmapped, methylation tags)
aligned_bam_path = sys.argv[2]   # BAM after minimap2 (aligned)
output_bam_path = sys.argv[3]    # Output BAM with methylation tags restored

# Load unmapped reads and store tags
unmapped_bam = pysam.AlignmentFile(unmapped_bam_path, "rb")
tags_dict = {}
for read in unmapped_bam.fetch(until_eof=True):
    tags_dict[read.query_name] = read.get_tags()
unmapped_bam.close()

# Open aligned BAM and prepare output
aligned_bam = pysam.AlignmentFile(aligned_bam_path, "rb")
out_bam = pysam.AlignmentFile(output_bam_path, "wb", template=aligned_bam)

# Transfer methylation tags back
for read in aligned_bam.fetch():
    tags = tags_dict.get(read.query_name, [])
    methylation_tags = [(tag, val) for tag, val in tags if tag in ["MM", "ML"]]
    read.set_tags(read.get_tags() + methylation_tags)
    out_bam.write(read)

aligned_bam.close()
out_bam.close()
