#!/bin/bash
#SBATCH --job-name=pbsv_call
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam/slurm_logs/call_%j.out
#SBATCH --error=/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam/slurm_logs/call_%j.err

# activate the conda environment
module load conda/latest
conda activate pacbiosv

# some more paths
which pbsv
which samtools
echo "PATH=$PATH"

# file/directory paths
BASE="/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam"
REF="/work/pi_nhowlett_uri_edu/jessie/reference/WBce1235_genomic.fa"
SVSIG_DIR="${BASE}/pbsv_svsig"
OUTVCF="${BASE}/pbsv_variants/all_samples.pbsv.vcf"

mkdir -p "${BASE}/pbsv_variants"

echo "Collecting svsig files..."

# automatically grab all svsig files
SVSIG_FILES=("$SVSIG_DIR"/*.svsig.gz)
echo "Running pbsv call..."

# debugging statements to see what is being passed to each variable

# checking the reference file
echo "REF argument (raw):"
printf '%q\n' "$REF"

if [ ! -f "$REF" ]; then
    echo "ERROR: Reference FASTA not found!"
    exit 1
fi

ls -lh "$REF"

head -10 "$REF" | cat -A

# checking the svsig_files array built from the directory
echo "SVSIG_FILES array:"
for f in "${SVSIG_FILES[@]}"; do
    printf '%q\n' "$f"
done

echo "Number of SVSIG files: ${#SVSIG_FILES[@]}"


set -x # debugging caller to output the exact command after pass throughs
# call variants with pbsv
pbsv call --ccs -j ${SLURM_CPUS_PER_TASK} "$REF" "${SVSIG_FILES[@]}" "$OUTVCF"

echo "Finished pbsv call"