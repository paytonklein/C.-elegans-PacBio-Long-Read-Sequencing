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
echo "starting filter to standard chromosomes"
for f in IP_hnRNP-H1 Input_hnRNP-HI IP_IgG Input_IgG
do
    samtools view -h ${f}.trimmed_Aligned.sortedByCoord.out.bam \
    | awk '$1 ~ /^@/ || $3 ~ /^chr([0-9]+|X|Y|M)$/' \
    | samtools view -b > ${f}.filtered.bam
done
echo "done chromosome filtering w samtools"

conda deactivate
conda activate piranha_env
echo "Active environment name: $CONDA_DEFAULT_ENV"

echo "start converting from bam -> bed"

# -----------------------------
# Step 2: Convert BAM → BED
# -----------------------------
for f in IP_hnRNP-H1 Input_hnRNP-HI IP_IgG Input_IgG
do
    bedtools bamtobed -i ${f}.filtered.bam > ${f}.bed
done

echo "done converting to bed files"

echo "sanity checks: (hopefully return nothing)"
cut -f1 IP_hnRNP-H1.bed | sort | uniq > ip.chroms
cut -f1 Input_hnRNP-HI.bed | sort | uniq > input.chroms
diff ip.chroms input.chroms

echo "sorting bed files" 
# -----------------------------
# Step 2.5: Sort BED files (CRITICAL)
# -----------------------------
for f in IP_hnRNP-H1 Input_hnRNP-HI IP_IgG Input_IgG
do
    sort -k1,1 -k2,2n ${f}.bed > ${f}.sorted.bed
done

echo "run Piranha for hnRNP-H1"
# -----------------------------
# Step 3: Run Piranha (hnRNP-H1)
# -----------------------------
Piranha -b 20 -s -a 0 \
-c Input_hnRNP-HI.sorted.bed \
IP_hnRNP-H1.sorted.bed \
> hnRNP-H1_peaks.bed

echo "run Piranha for IgG"
# -----------------------------
# Step 4: Run Piranha (IgG)
# -----------------------------
Piranha -b 20 -s -a 0 \
-c Input_IgG.sorted.bed \
IP_IgG.sorted.bed \
> IgG_peaks.bed

echo "done!"