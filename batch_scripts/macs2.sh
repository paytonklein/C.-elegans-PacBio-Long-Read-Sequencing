#!/bin/bash
#SBATCH --job-name=macs2
#SBATCH --output=macs2.out
#SBATCH --error=macs2.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G

module load conda/latest

conda activate macs2_clean_env
echo "Active environment: $CONDA_DEFAULT_ENV"

cd /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments/BAMs

# -----------------------------
# Run MACS2 peak calling
# -----------------------------
echo "Running MACS2 (hnRNP-H1)"
macs2 callpeak \
  -t IP_hnRNP-H1.trimmed_Aligned.sortedByCoord.out.bam \
  -c Input_hnRNP-HI.trimmed_Aligned.sortedByCoord.out.bam \
  -f BAM \
  -g hs \
  -n hnRNP-H1 \
  --outdir macs2_hnRNP-H1

echo "Running MACS2 (IgG control)"
macs2 callpeak \
  -t IP_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  -c Input_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  -f BAM \
  -g hs \
  -n IgG \
  --outdir macs2_IgG

echo "Done"