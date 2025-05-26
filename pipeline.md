# 🧬 Smart-seq3 Mutation Analysis Pipeline (BAM -> VCF)

This pipeline guides you through analyzing **mutational status of specific genes** (e.g., *Braf*, *Trp53*) from **Smart-seq3 BAM files** using **GATK (in bash)** and **R** for interpretation.

---

## 📁 1. Input Requirements
- Smart-seq3 aligned BAM file (`.bam`) and its index (`.bai`)
- Reference genome: Ensembl **GRCm38 (no chr prefix)**
- Gene annotation: optional for visualization

---

## ⚙️ 2. Environment Setup (WSL/Linux)
```bash
# Install required tools
sudo apt update && sudo apt install -y samtools unzip openjdk-17-jre wget

# Install GATK
cd ~
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
sudo ln -s ~/gatk-4.5.0.0/gatk /usr/local/bin/gatk
```

---

## 🧬 3. Prepare Reference Files
```bash
# Create reference directory
mkdir -p ~/gatk_ref && cd ~/gatk_ref

# Download Ensembl GRCm38 primary assembly FASTA
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

# Unzip and index
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
samtools faidx Mus_musculus.GRCm38.dna.primary_assembly.fa

gatk CreateSequenceDictionary \
  -R Mus_musculus.GRCm38.dna.primary_assembly.fa \
  -O Mus_musculus.GRCm38.dna.primary_assembly.dict
```

---

## 🔍 4. Filter BAM for Gene Regions of Interest
```bash
# Example: Braf (chr6:37497000-37507000), Trp53 (chr11:69566000-69586000)
samtools view -b input.bam \
  6:37497000-37507000 11:69566000-69586000 \
  -o braf_tp53.bam

samtools index braf_tp53.bam
```

---

## 🧼 5. Preprocess BAM for GATK
```bash
# Split spliced reads
gatk SplitNCigarReads \
  -R Mus_musculus.GRCm38.dna.primary_assembly.fa \
  -I braf_tp53.bam \
  -O braf_tp53.split.bam

# Add read groups if missing
gatk AddOrReplaceReadGroups \
  -I braf_tp53.split.bam \
  -O braf_tp53.rg.bam \
  -RGID plate1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM sample1

samtools index braf_tp53.rg.bam
```

---

## 🔬 6. Call Variants with GATK
```bash
gatk HaplotypeCaller \
  -R Mus_musculus.GRCm38.dna.primary_assembly.fa \
  -I braf_tp53.rg.bam \
  -O braf_tp53.variants.vcf \
  -L 6:37497000-37507000 \
  -L 11:69566000-69586000 \
  --dont-use-soft-clipped-bases \
  --standard-min-confidence-threshold-for-calling 20.0
```

---

## 📊 7. Analyze Variants in R
```r
library(VariantAnnotation)

vcf <- readVcf("C:/your/path/braf_tp53.variants.vcf", genome = "GRCm38")
vr <- rowRanges(vcf)

# Check Braf V637E (mouse equivalent of human BRAF V600E)
braf_v637e <- vr[seqnames(vr) == "6" & start(vr) == 37498645]

# Check Trp53 region variants
trp53 <- vr[seqnames(vr) == "11" & start(vr) > 69566000 & start(vr) < 69586000]

# Genotypes
geno(vcf)$GT[rownames(vcf) %in% names(trp53), ]
```

---

## 📁 Output
- `braf_tp53.variants.vcf`: contains variant calls
- Visualizations: made in R using `ggplot2` or `GenomicRanges`

---

## ✅ Tips
- Ensure reference and BAM use same contig naming (no `chr` prefix mismatch)
- You can batch this for multiple BAMs
- Filter `.vcf` by `GT` to classify WT vs mutant

---

## 📌 To Do
- Automate in Snakemake or shell script
- Extend to full-exome or CNV analysis with `GATK4` and `AllelicCNV`

---
© 2025 – Ready to upload to GitHub. Add README + sample data + optional RMarkdown notebook for visualization.
