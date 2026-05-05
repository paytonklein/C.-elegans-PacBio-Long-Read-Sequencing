#!/bin/bash
#SBATCH --partition=uri-cpu,cpu
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --job-name=bam2bw
#SBATCH --output=bam2bw_%J.out

# Go to working directory
cd /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments/BAMs

# Activate conda environment
module load conda/latest
conda activate binding_sites
echo "Active environment: $CONDA_DEFAULT_ENV"

# Output directory
OUTPUT_DIR=bigwigs
mkdir -p $OUTPUT_DIR

# -----------------------------
# Run bamCoverage for each sample
# -----------------------------

bamCoverage -b Input_hnRNP-HI.trimmed_Aligned.sortedByCoord.out.bam \
  -o $OUTPUT_DIR/Input_hnRNP-HI.bw \
  --normalizeUsing CPM \
  --binSize 10

bamCoverage -b IP_hnRNP-H1.trimmed_Aligned.sortedByCoord.out.bam \
  -o $OUTPUT_DIR/IP_hnRNP-H1.bw \
  --normalizeUsing CPM \
  --binSize 10

bamCoverage -b Input_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  -o $OUTPUT_DIR/Input_IgG.bw \
  --normalizeUsing CPM \
  --binSize 10

bamCoverage -b IP_IgG.trimmed_Aligned.sortedByCoord.out.bam \
  -o $OUTPUT_DIR/IP_IgG.bw \
  --normalizeUsing CPM \
  --binSize 10

echo "Done generating bigWigs"