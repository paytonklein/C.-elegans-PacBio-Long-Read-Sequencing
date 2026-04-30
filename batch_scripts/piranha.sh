#!/bin/bash
#SBATCH --job-name=piranha
#SBATCH --output=piranha.out
#SBATCH --error=piranha.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

module load conda/latest
conda activate piranha_env

cd /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments/BAMs

Piranha -b 20 -s \
-c Input_hnRNP-HI.trimmed_Aligned.sortedByCoord.out.bam \
IP_hnRNP-H1.trimmed_Aligned.sortedByCoord.out.bam \
> hnRNP-H1_peaks.bed

Piranha -b 20 -s \
-c Input_IgG.trimmed_Aligned.sortedByCoord.out.bam \
IP_IgG.trimmed_Aligned.sortedByCoord.out.bam \
> IgG_peaks.bed