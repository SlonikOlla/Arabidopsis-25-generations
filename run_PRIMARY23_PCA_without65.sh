#!/usr/bin/env bash
#SBATCH --job-name=PCA_no65_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/primary23/PCA_no65_P23_%j.out
#SBATCH --error=logs/UVC/primary23/PCA_no65_P23_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
PCADIR="${ROOT}/results/wgs_variants/PRIMARY23_PCA"
INPUT="${PCADIR}/PRIMARY23_PCA_input"
SNPS="${PCADIR}/PRIMARY23_PCA_5kb_SNPs_v2.txt"
FREQ="${PCADIR}/PRIMARY23_all23_5kb_frequencies.afreq"
REMOVE="${PCADIR}/remove_SRR17867665_FID_IID.txt"
OUTPUT="${PCADIR}/PRIMARY23_without_SRR17867665_5kb_v4"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load plink/2.0.0-a.6.32

printf '#FID\tIID\nSRR17867665\tSRR17867665\n' > "$REMOVE"

echo "===== REMOVAL FILE ====="
cat "$REMOVE"

echo "===== SENSITIVITY PCA ====="

plink2 \
  --pfile "$INPUT" \
  --allow-extra-chr \
  --extract "$SNPS" \
  --remove "$REMOVE" \
  --read-freq "$FREQ" \
  --pca 10 allele-wts \
  --threads 2 \
  --out "$OUTPUT"

test -s "${OUTPUT}.eigenvec"
test -s "${OUTPUT}.eigenval"

sample_count=$(
  awk 'NR>1 {n++} END {print n+0}' \
    "${OUTPUT}.eigenvec"
)

echo "SENSITIVITY PCA RECORDS: $sample_count"

if [[ "$sample_count" -ne 22 ]]; then
    echo "ERROR: expected 22 samples"
    exit 1
fi

if grep -q 'SRR17867665' "${OUTPUT}.eigenvec"; then
    echo "ERROR: SRR17867665 remains in sensitivity PCA"
    exit 1
fi

echo "SRR17867665 ABSENT: PASS"

echo "===== EIGENVALUES ====="
nl -ba "${OUTPUT}.eigenval"

echo "PRIMARY23 SENSITIVITY PCA: COMPLETE"
