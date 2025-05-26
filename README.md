# smartseq3-mutation-pipeline
A streamlined workflow to detect and analyze mutations in specific genes (e.g., Braf, Trp53) from Smart-seq3 plate-level BAM files using GATK for variant calling and R for downstream analysis and visualization. Designed for mouse data aligned to GRCm38.

**#Input Requirement**

Smart-seq3 aligned BAM file (.bam) and its index (.bai)

Reference genome: Ensembl GRCm38 (no chr prefix)

Gene annotation: optional for visualization

# Output

braf_tp53.variants.vcf: contains variant calls

Visualizations: made in R using ggplot2 or GenomicRanges

**#Tips**

Ensure reference and BAM use same contig naming (no chr prefix mismatch)

You can batch this for multiple BAMs

Filter .vcf by GT to classify WT vs mutant
