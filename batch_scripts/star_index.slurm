#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --output=star_index.out
#SBATCH --error=star_index.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G

module load STAR
# or: conda activate star_env

STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /path/to/STAR_index \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/annotation.gtf \
  --sjdbOverhang 100