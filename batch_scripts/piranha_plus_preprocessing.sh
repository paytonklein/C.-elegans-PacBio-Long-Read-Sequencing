#!/bin/bash
#SBATCH --job-name=piranha
#SBATCH --output=piranha.out
#SBATCH --error=piranha.err
#SBATCH --time=02:00:00
#SBATCH --mem=40G

module load conda/latest
conda activate binding_sites
echo "Active environment name: $CONDA_DEFAULT_ENV"

cd /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments/BAMs
pwd

# -----------------------------
# Step 1: Filter to standard chromosomes
# -----------------------------
echo "Filtering BAMs"
for f in IP_hnRNP-H1 Input_hnRNP-HI IP_IgG Input_IgG
do
    samtools view -h ${f}.trimmed_Aligned.sortedByCoord.out.bam \
    | awk '$1 ~ /^@/ || $3 ~ /^chr([0-9]+|X|Y|M)$/' \
    | samtools view -b > ${f}.filtered.bam
done

# -----------------------------
# Step 2: BAM → CLEAN BED6 (CRITICAL FIX)
# -----------------------------
echo "Converting to clean BED6"
for f in IP_hnRNP-H1 Input_hnRNP-HI IP_IgG Input_IgG
do
    bedtools bamtobed -i ${f}.filtered.bam \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,".",$5,$6}' \
    | sort -k1,1 -k2,2n > ${f}.clean.bed
done

# -----------------------------
# Switch environment for Piranha
# -----------------------------
conda deactivate
conda activate piranha_env
echo "Active environment name: $CONDA_DEFAULT_ENV"

# -----------------------------
# Step 3: Run Piranha
# -----------------------------
echo "Running Piranha (hnRNP-H1)"
Piranha -b 20 -s -a 0 \
-c Input_hnRNP-HI.clean.bed \
IP_hnRNP-H1.clean.bed \
> hnRNP-H1_peaks.bed

echo "Running Piranha (IgG)"
Piranha -b 20 -s -a 0 \
-c Input_IgG.clean.bed \
IP_IgG.clean.bed \
> IgG_peaks.bed

echo "Done"