#!/bin/bash
#SBATCH --job-name=piranha
#SBATCH --output=piranha.out
#SBATCH --error=piranha.err
#SBATCH --time=02:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4

module load conda/latest
conda activate piranha_env

cd /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments/BAMs

# -----------------------------
# hnRNP-H1 IP vs Input
# -----------------------------
Piranha \
  -b 20 \
  -s \
  -a \
  Input_hnRNP-HI.trimmed_Aligned.sortedByCoord.out.bam \
  IP_hnRNP-H1.trimmed_Aligned.sortedByCoord.out.bam \
  > hnRNP-H1_peaks.bed
# -----------------------------
# IgG control (optional comparison)
# -----------------------------
Piranha \
  -b 20 \
  -s \
  -a \
  Input_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  IP_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  > IgG_peaks.bed