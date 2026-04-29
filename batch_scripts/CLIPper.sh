#!/bin/bash
#SBATCH --job-name=clipper
#SBATCH --output=clipper.out
#SBATCH --error=clipper.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G

module load conda/latest
conda activate clipper_env

clipper \
  -b IP_hnRNP-H1.trimmed_Aligned.sortedByCoord.out.bam \
  -s hg38 \
  --min-counts 5 \
  --threshold-method binomial \
  -o IP_hnRNP-H1_clipper_peaks.bed

clipper \
  -b IP_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  -s hg38 \
  --min-counts 5 \
  --threshold-method binomial \
  -o IP_IgG_clipper_peaks.bed