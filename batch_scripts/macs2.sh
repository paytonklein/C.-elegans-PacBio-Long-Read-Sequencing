#!/bin/bash
#SBATCH --job-name=macs3
#SBATCH --output=macs3.out
#SBATCH --error=macs3.err
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

module load conda/latest

conda activate macs3_env
echo "Active environment: $CONDA_DEFAULT_ENV"

cd /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments/BAMs

mkdir -p macs3_hnRNP-H1 macs3_IgG

# -----------------------------
# Run MACS3 peak calling
# -----------------------------
echo "Running MACS3 (hnRNP-H1)"

macs3 callpeak \
  -t IP_hnRNP-H1.trimmed_Aligned.sortedByCoord.out.bam \
  -c Input_hnRNP-HI.trimmed_Aligned.sortedByCoord.out.bam \
  -f BAM \
  -g hs \
  -n hnRNP-H1 \
  --outdir macs3_hnRNP-H1 \
  --nomodel \
  --extsize 75 \
  -q 0.05

echo "Running MACS3 (IgG control)"

macs3 callpeak \
  -t IP_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  -c Input_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  -f BAM \
  -g hs \
  -n IgG \
  --outdir macs3_IgG \
  --nomodel \
  --extsize 75 \
  -q 0.05

echo "Done"