#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --output=star_index.out
#SBATCH --error=star_index.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G

module load conda/latest
conda activate binding_sites

STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/STAR_index \
  --genomeFastaFiles /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile /scratch4/workspace/payton_klein_uri_edu-CMB320/hg38/gencode.v38.annotation.gtf \
  --sjdbOverhang 74