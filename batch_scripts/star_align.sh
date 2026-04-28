#!/bin/bash
#SBATCH --job-name=STAR_align
#SBATCH --output=star_%A_%a.out
#SBATCH --error=star_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-4

set -euo pipefail

# -----------------------------
# Paths
# -----------------------------
SAMPLE_LIST=/home/payton_klein_uri_edu/HNRNPH1/samples.txt
FASTQ_DIR=/home/payton_klein_uri_edu/HNRNPH1
INDEX=/scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_index
OUTDIR=/scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_alignments

mkdir -p $OUTDIR

# Get sample
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
BASE=$(basename "$SAMPLE" .fastq)

echo "Processing sample: $SAMPLE"

# Run STAR
STAR \
--runThreadN $SLURM_CPUS_PER_TASK \
--genomeDir $INDEX \
--readFilesIn $FASTQ_DIR/$SAMPLE \
--outFileNamePrefix $OUTDIR/${BASE}_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--limitGenomeGenerateRAM 30000000000