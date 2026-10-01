#!/usr/bin/env bash
#SBATCH --job-name=PCA5kb_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/primary23/PCA5kb_P23_%j.out
#SBATCH --error=logs/UVC/primary23/PCA5kb_P23_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
PCADIR="${ROOT}/results/wgs_variants/PRIMARY23_PCA"
INPUT="${PCADIR}/PRIMARY23_PCA_input"
SNPS="${PCADIR}/PRIMARY23_PCA_5kb_SNPs_v2.txt"
ALL="${PCADIR}/PRIMARY23_all23_5kb_v2"
NO65="${PCADIR}/PRIMARY23_without_SRR17867665_5kb_v2"
REMOVE="${PCADIR}/remove_SRR17867665_5kb_v2.txt"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load plink/2.0.0-a.6.32

echo "===== INPUT CHECK ====="
test -s "${INPUT}.pgen"
test -s "${INPUT}.pvar"
test -s "${INPUT}.psam"
test -s "$SNPS"

echo "THINNED SNPS: $(wc -l < "$SNPS")"

echo "===== PCA: ALL 23 SAMPLES ====="

plink2 \
  --pfile "$INPUT" \
  --allow-extra-chr \
  --extract "$SNPS" \
  --pca 10 allele-wts \
  --threads 2 \
  --out "$ALL"

test -s "${ALL}.eigenvec"
test -s "${ALL}.eigenval"

echo "ALL-SAMPLE PCA RECORDS: $(awk 'NR>1 {n++} END {print n+0}' "${ALL}.eigenvec")"

echo "===== PCA: WITHOUT SRR17867665 ====="

printf '#IID\nSRR17867665\n' > "$REMOVE"

plink2 \
  --pfile "$INPUT" \
  --allow-extra-chr \
  --extract "$SNPS" \
  --remove "$REMOVE" \
  --pca 10 allele-wts \
  --threads 2 \
  --out "$NO65"

test -s "${NO65}.eigenvec"
test -s "${NO65}.eigenval"

echo "SENSITIVITY PCA RECORDS: $(awk 'NR>1 {n++} END {print n+0}' "${NO65}.eigenvec")"

echo "===== EIGENVALUES: ALL 23 ====="
nl -ba "${ALL}.eigenval"

echo "===== EIGENVALUES: WITHOUT SRR17867665 ====="
nl -ba "${NO65}.eigenval"

echo "PRIMARY23 PCA: COMPLETE"
