#!/bin/bash
#SBATCH --job-name=piranha
#SBATCH --output=piranha.out
#SBATCH --error=piranha.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G

module load conda/latest
conda activate piranha_env

cd /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments/BAMs

# Convert BAM → BED
bedtools bamtobed -i IP_hnRNP-H1.trimmed_Aligned.sortedByCoord.out.bam > IP_hnRNP-H1.bed
bedtools bamtobed -i Input_hnRNP-HI.trimmed_Aligned.sortedByCoord.out.bam > Input_hnRNP-HI.bed

bedtools bamtobed -i IP_IgG.trimmed_Aligned.sortedByCoord.out.bam > IP_IgG.bed
bedtools bamtobed -i Input_IgG.trimmed_Aligned.sortedByCoord.out.bam > Input_IgG.bed

# Run Piranha
Piranha -b 20 -s -c Input_hnRNP-HI.bed IP_hnRNP-H1.bed > hnRNP-H1_peaks.bed
Piranha -b 20 -s -c Input_IgG.bed IP_IgG.bed > IgG_peaks.bed