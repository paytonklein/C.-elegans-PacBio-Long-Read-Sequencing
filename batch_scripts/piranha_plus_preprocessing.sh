#!/bin/bash
#SBATCH --job-name=piranha
#SBATCH --output=piranha.out
#SBATCH --error=piranha.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G

module load conda/latest
conda activate piranha_env

cd /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments/BAMs

# -----------------------------
# Step 1: Filter to standard chromosomes
# -----------------------------
for f in IP_hnRNP-H1 Input_hnRNP-HI IP_IgG Input_IgG
do
    samtools view -h ${f}.trimmed_Aligned.sortedByCoord.out.bam \
    | grep -E '^@|^chr[0-9XYM]' \
    | samtools view -b > ${f}.filtered.bam
done

# -----------------------------
# Step 2: Convert BAM → BED
# -----------------------------
for f in IP_hnRNP-H1 Input_hnRNP-HI IP_IgG Input_IgG
do
    bedtools bamtobed -i ${f}.filtered.bam > ${f}.bed
done

# -----------------------------
# Step 3: Run Piranha (hnRNP-H1)
# -----------------------------
Piranha -b 20 -s \
-c Input_hnRNP-HI.bed \
IP_hnRNP-H1.bed \
> hnRNP-H1_peaks.bed

# -----------------------------
# Step 4: Run Piranha (IgG)
# -----------------------------
Piranha -b 20 -s \
-c Input_IgG.bed \
IP_IgG.bed \
> IgG_peaks.bed