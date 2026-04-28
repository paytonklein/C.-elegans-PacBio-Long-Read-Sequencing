#!/bin/bash
#SBATCH --job-name=STAR_align
#SBATCH --output=star_%A_%a.out
#SBATCH --error=star_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-4