#!/usr/bin/env bash
#SBATCH --job-name=Dist_WGS_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/primary23/genetic_distance_P23_%j.out
#SBATCH --error=logs/UVC/primary23/genetic_distance_P23_%j.err

set -euxo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
PGEN="${ROOT}/results/wgs_variants/PRIMARY23_PCA/PRIMARY23_PCA_input"
GROUP_FILE="${ROOT}/metadata/wgs_primary23_groups.tsv"
OUTDIR="${ROOT}/results/wgs_variants/PRIMARY23_relatedness"
EXPORT="${OUTDIR}/PRIMARY23_5119_SNP_genotypes"
PREFIX="PRIMARY23_5119_SNP"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load plink/2.0.0-a.6.32
module load python/3.11

echo "========================================"
echo "PRIMARY23 GENETIC-DISTANCE ANALYSIS"
echo "START: $(date)"
echo "HOST:  $(hostname)"
echo "========================================"

echo "===== INPUT VALIDATION ====="

for file in \
    "${PGEN}.pgen" \
    "${PGEN}.pvar" \
    "${PGEN}.psam" \
    "$GROUP_FILE" \
    "${ROOT}/scripts/calculate_PRIMARY23_genetic_distances.py"
do
    test -s "$file"
    echo "PASS  $file"
done

sample_count=$(
    awk 'NR>1 {n++} END {print n+0}' "${PGEN}.psam"
)

variant_count=$(
    awk '$1 !~ /^#/ {n++} END {print n+0}' "${PGEN}.pvar"
)

echo "SAMPLES:  $sample_count"
echo "VARIANTS: $variant_count"

if [[ "$sample_count" -ne 23 ]]; then
    echo "ERROR: expected 23 samples"
    exit 1
fi

if [[ "$variant_count" -ne 5119 ]]; then
    echo "ERROR: expected 5119 variants"
    exit 1
fi

echo "===== EXPORT ADDITIVE GENOTYPES ====="

plink2 \
  --pfile "$PGEN" \
  --allow-extra-chr \
  --export A \
  --threads 2 \
  --out "$EXPORT"

test -s "${EXPORT}.raw"

echo "===== CALCULATE PAIRWISE DISTANCES ====="

python3 scripts/calculate_PRIMARY23_genetic_distances.py \
  "${EXPORT}.raw" \
  "$GROUP_FILE" \
  "$OUTDIR" \
  "$PREFIX"

echo "===== OUTPUT VALIDATION ====="

PAIRWISE="${OUTDIR}/${PREFIX}_pairwise_distances.tsv"
MATRIX="${OUTDIR}/${PREFIX}_IBS_distance_matrix.tsv"
NEAREST="${OUTDIR}/${PREFIX}_nearest_neighbors.tsv"
UVC="${OUTDIR}/${PREFIX}_UVC_ranked_neighbors.tsv"

for file in "$PAIRWISE" "$MATRIX" "$NEAREST" "$UVC"; do
    test -s "$file"
    echo "PASS  $file"
done

pair_count=$(
    awk 'NR>1 {n++} END {print n+0}' "$PAIRWISE"
)

nearest_count=$(
    awk 'NR>1 {n++} END {print n+0}' "$NEAREST"
)

uvc_count=$(
    awk 'NR>1 {n++} END {print n+0}' "$UVC"
)

echo "PAIRWISE RECORDS: $pair_count"
echo "NEAREST RECORDS:  $nearest_count"
echo "UVC RANK RECORDS: $uvc_count"

if [[ "$pair_count" -ne 253 ]]; then
    echo "ERROR: expected 253 pairwise comparisons"
    exit 1
fi

if [[ "$nearest_count" -ne 23 ]]; then
    echo "ERROR: expected 23 nearest-neighbor records"
    exit 1
fi

if [[ "$uvc_count" -ne 110 ]]; then
    echo "ERROR: expected 110 UVC-neighbor records"
    exit 1
fi

echo
echo "===== UVC TOP FIVE GENETIC NEIGHBORS ====="

awk -F '\t' '
    NR==1 || $2<=5
' "$UVC" | column -t

echo
echo "===== ALL-SAMPLE NEAREST NEIGHBORS ====="

column -t "$NEAREST"

echo
echo "PRIMARY23 GENETIC-DISTANCE ANALYSIS: COMPLETE"
echo "FINISH: $(date)"
