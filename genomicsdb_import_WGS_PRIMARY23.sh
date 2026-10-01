#!/usr/bin/env bash
#SBATCH --job-name=GDB_WGS_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --array=1-5%2
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=logs/UVC/primary23/genomicsdb_P23_%A_%a.out
#SBATCH --error=logs/UVC/primary23/genomicsdb_P23_%A_%a.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
REF="${ROOT}/reference/TAIR10/TAIR10_chr_all.fa"
MAP="${ROOT}/metadata/wgs_gvcf_sample_map_primary23.tsv"

CHR="${SLURM_ARRAY_TASK_ID}"
WORKSPACE="${ROOT}/results/wgs_variants/genomicsdb_primary23/chr${CHR}"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load gatk/4.6.1.0

export TMPDIR="${SLURM_TMPDIR}"

echo "========================================"
echo "PRIMARY23 GENOMICSDB IMPORT"
echo "CHROMOSOME: ${CHR}"
echo "WORKSPACE:  ${WORKSPACE}"
echo "HOST:       $(hostname)"
echo "START:      $(date)"
echo "TMPDIR:     ${TMPDIR}"
echo "========================================"

echo "===== INPUT VALIDATION ====="

for f in \
    "$REF" \
    "${REF}.fai" \
    "${ROOT}/reference/TAIR10/TAIR10_chr_all.dict" \
    "$MAP"
do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing or empty input: $f"
        exit 1
    fi

    printf "PASS\t%s\t%s bytes\n" \
        "$f" "$(stat -Lc '%s' "$f")"
done

record_count=$(wc -l < "$MAP")
unique_samples=$(cut -f1 "$MAP" | sort -u | wc -l)
unique_paths=$(cut -f2 "$MAP" | sort -u | wc -l)

echo "MAP RECORDS:       ${record_count}"
echo "UNIQUE SAMPLES:    ${unique_samples}"
echo "UNIQUE GVCF PATHS: ${unique_paths}"

if [[ "$record_count" -ne 23 ||
      "$unique_samples" -ne 23 ||
      "$unique_paths" -ne 23 ]]
then
    echo "ERROR: PRIMARY23 sample map validation failed"
    exit 1
fi

while IFS=$'\t' read -r sample gvcf; do
    if [[ ! -s "$gvcf" || ! -s "${gvcf}.tbi" ]]; then
        echo "ERROR: missing gVCF or index for ${sample}"
        exit 1
    fi
done < "$MAP"

echo "ALL 23 GVCFS AND INDEXES: PASS"

if [[ -e "$WORKSPACE" ]]; then
    echo "ERROR: workspace already exists: $WORKSPACE"
    echo "No existing workspace was changed."
    exit 1
fi

echo "===== GENOMICSDB IMPORT ====="

gatk \
    --java-options "-Xmx18g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
    GenomicsDBImport \
    -R "$REF" \
    --sample-name-map "$MAP" \
    --genomicsdb-workspace-path "$WORKSPACE" \
    -L "$CHR" \
    --reader-threads 4 \
    --tmp-dir "$SLURM_TMPDIR"

echo "===== OUTPUT VALIDATION ====="

if [[ ! -d "$WORKSPACE" ]]; then
    echo "ERROR: workspace was not created: $WORKSPACE"
    exit 1
fi

file_count=$(find "$WORKSPACE" -type f | wc -l)
workspace_size=$(du -sh "$WORKSPACE" | awk '{print $1}')

if [[ "$file_count" -eq 0 ]]; then
    echo "ERROR: workspace contains no files"
    exit 1
fi

echo "WORKSPACE FILES: ${file_count}"
echo "WORKSPACE SIZE:  ${workspace_size}"
echo "GENOMICSDB IMPORT VALIDATED: CHR${CHR}"

echo "PRIMARY18 DATABASE RETAINED"
echo "ALL GVCF INPUTS RETAINED"
echo "PRIMARY23 GENOMICSDB COMPLETE: CHR${CHR}"
echo "FINISH: $(date)"
