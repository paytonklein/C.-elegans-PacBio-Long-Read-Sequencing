#!/bin/bash
#SBATCH --job-name=pbmm2_align
#SBATCH --array=1-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/align_%A_%a.out
#SBATCH --error=logs/align_%A_%a.err

# activate the conda env
module load conda/latest
conda activate pbsv

# set file/directory locations
BASE="/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam"
REF="/work/pi_nhowlett_uri_edu/Celegans_PacBio_CNV_2025/GCF_000002985.6_WBcel235_genomic.fna"
SAMPLE_LIST="${BASE}/samples.txt"

OUTDIR="${BASE}/pbsv_aligned"
mkdir -p $OUTDIR

# get sample name from txt file
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST)

SAMPLE_DIR="${BASE}/${SAMPLE}"
FASTQ="${SAMPLE_DIR}/${SAMPLE}.fastq.gz"
OUTPUT_BAM="${OUTDIR}/${SAMPLE}.aligned.bam"

# make sure fastq is found in each <sample_prefix> directory
if [[ ! -f "$FASTQ" ]]; then
    echo "ERROR: FASTQ not found:"
    echo "$FASTQ"
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "FASTQ: $FASTQ"

# run alignment with pbmm2
pbmm2 align $REF \
    $FASTQ \
    $OUTPUT_BAM \
    --preset CCS \
    --sort \
    --sample $SAMPLE \
    --rg "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
    -j ${SLURM_CPUS_PER_TASK}

echo "Finished $SAMPLE"