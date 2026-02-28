#!/bin/bash
#SBATCH --job-name=pbsv_discover
#SBATCH --array=1-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam/slurm_logs/discover_%A_%a.out
#SBATCH --error=/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam/slurm_logs/discover_%A_%a.err

# activate the conda env
module load conda/latest
conda activate pacbiosv

# file/directory paths
BASE="/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam"
SAMPLE_LIST="${BASE}/samples.txt"
ALIGN_DIR="${BASE}/pbsv_aligned"
DISCOVER_DIR="${BASE}/pbsv_svsig"
mkdir -p $DISCOVER_DIR

TANDEM_REPEATS="/path/to/tandem_repeats.bed"    # provides repeat annotation - recommended

# get the sample prefix for each sample
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST)
BAM="${ALIGN_DIR}/${SAMPLE}.aligned.bam"
SVSIG="${DISCOVER_DIR}/${SAMPLE}.svsig.gz"

# run pbsv discover to get all of the signatures for the SVs
pbsv discover $BAM $SVSIG --tandem-repeats $TANDEM_REPEATS
tabix -c '#' -s 3 -b 4 -e 4 $   # optional but recommended if we need to make a pbsv random call 