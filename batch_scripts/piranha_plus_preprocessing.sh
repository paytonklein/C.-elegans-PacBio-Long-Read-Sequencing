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
# Step 1: Filter BAMs
# -----------------------------
echo "Filtering BAMs"
for f in IP_hnRNP-H1 Input_hnRNP-HI IP_IgG Input_IgG
do
    samtools view -h ${f}.trimmed_Aligned.sortedByCoord.out.bam \
    | awk '$1 ~ /^@/ || $3 ~ /^chr([0-9]+|X|Y|M)$/' \
    | samtools view -b > ${f}.filtered.bam

    # (optional but recommended fix for idxstats warning)
    samtools index ${f}.filtered.bam
done

# -----------------------------
# Step 2: Create genome file
# -----------------------------
echo "Creating genome file"
samtools idxstats IP_hnRNP-H1.filtered.bam \
| awk '$1!="*" {print $1"\t"$2}' > hg38.genome

# -----------------------------
# Step 3: Make windows (FIXED)
# -----------------------------
echo "Creating 100bp windows (fixed for memory stability)"
bedtools makewindows -g hg38.genome -w 100 > windows.bed

# -----------------------------
# Step 4: Compute coverage per bin
# -----------------------------
echo "Generating binned coverage"
for f in IP_hnRNP-H1 Input_hnRNP-HI IP_IgG Input_IgG
do
    bedtools coverage -a windows.bed -b ${f}.filtered.bam \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,".",$5,"+"}' \
    > ${f}.binned.bed
done

# -----------------------------
# Switch environment
# -----------------------------
conda deactivate
conda activate piranha_env
echo "Active environment name: $CONDA_DEFAULT_ENV"

# -----------------------------
# Step 5: Run Piranha
# -----------------------------
echo "Running Piranha (hnRNP-H1)"
Piranha -s -a 0 \
-c Input_hnRNP-HI.binned.bed \
IP_hnRNP-H1.binned.bed \
> hnRNP-H1_peaks.bed

echo "Running Piranha (IgG)"
Piranha -s -a 0 \
-c Input_IgG.binned.bed \
IP_IgG.binned.bed \
> IgG_peaks.bed

echo "Done"