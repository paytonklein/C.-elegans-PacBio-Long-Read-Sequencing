#!/bin/bash
#SBATCH --job-name=macs2
#SBATCH --output=macs2.out
#SBATCH --error=macs2.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G

module load conda/latest
conda activate binding_sites   # or wherever MACS2 is installed

cd /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments/BAMs

# hnRNP-H1 peak calling
macs2 callpeak \
  -t IP_hnRNP-H1.trimmed_Aligned.sortedByCoord.out.bam \
  -c Input_hnRNP-HI.trimmed_Aligned.sortedByCoord.out.bam \
  -f BAM \
  -g hs \
  -n hnRNP-H1 \
  --outdir macs2_hnRNP-H1 \
  -q 0.05

# IgG control (optional)
macs2 callpeak \
  -t IP_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  -c Input_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  -f BAM \
  -g hs \
  -n IgG \
  --outdir macs2_IgG \
  -q 0.05