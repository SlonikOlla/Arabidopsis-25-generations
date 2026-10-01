#!/usr/bin/env bash
#SBATCH --job-name=PCA5kb_P23_v3
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/primary23/PCA5kb_P23_v3_%j.out
#SBATCH --error=logs/UVC/primary23/PCA5kb_P23_v3_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
PCADIR="${ROOT}/results/wgs_variants/PRIMARY23_PCA"
INPUT="${PCADIR}/PRIMARY23_PCA_input"
SNPS="${PCADIR}/PRIMARY23_PCA_5kb_SNPs_v2.txt"

FREQ="${PCADIR}/PRIMARY23_all23_5kb_frequencies"
ALL="${PCADIR}/PRIMARY23_all23_5kb_v3"
NO65="${PCADIR}/PRIMARY23_without_SRR17867665_5kb_v3"
REMOVE="${PCADIR}/remove_SRR17867665_v3.txt"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load plink/2.0.0-a.6.32

echo "===== INPUT CHECK ====="

for file in \
  "${INPUT}.pgen" \
  "${INPUT}.pvar" \
  "${INPUT}.psam" \
  "$SNPS"
do
    test -s "$file"
    echo "PASS  $file"
done

echo "THINNED SNPS: $(wc -l < "$SNPS")"

echo "===== CALCULATE ALL-23 ALLELE FREQUENCIES ====="

plink2 \
  --pfile "$INPUT" \
  --allow-extra-chr \
  --extract "$SNPS" \
  --freq \
  --threads 2 \
  --out "$FREQ"

test -s "${FREQ}.afreq"

echo "FREQUENCY RECORDS: $(awk 'NR>1 {n++} END {print n+0}' "${FREQ}.afreq")"

echo "===== PCA: ALL 23 SAMPLES ====="

plink2 \
  --pfile "$INPUT" \
  --allow-extra-chr \
  --extract "$SNPS" \
  --read-freq "${FREQ}.afreq" \
  --pca 10 allele-wts \
  --threads 2 \
  --out "$ALL"

test -s "${ALL}.eigenvec"
test -s "${ALL}.eigenval"

all_count=$(awk 'NR>1 {n++} END {print n+0}' "${ALL}.eigenvec")
echo "ALL-SAMPLE PCA RECORDS: $all_count"
test "$all_count" -eq 23

echo "===== PCA: WITHOUT SRR17867665 ====="

printf '#IID\nSRR17867665\n' > "$REMOVE"

plink2 \
  --pfile "$INPUT" \
  --allow-extra-chr \
  --extract "$SNPS" \
  --remove "$REMOVE" \
  --read-freq "${FREQ}.afreq" \
  --pca 10 allele-wts \
  --threads 2 \
  --out "$NO65"

test -s "${NO65}.eigenvec"
test -s "${NO65}.eigenval"

no65_count=$(awk 'NR>1 {n++} END {print n+0}' "${NO65}.eigenvec")
echo "SENSITIVITY PCA RECORDS: $no65_count"
test "$no65_count" -eq 22

echo "===== EIGENVALUES: ALL 23 ====="
nl -ba "${ALL}.eigenval"

echo "===== EIGENVALUES: WITHOUT SRR17867665 ====="
nl -ba "${NO65}.eigenval"

echo "===== FINAL VALIDATION ====="
echo "THINNED SNPS:            $(wc -l < "$SNPS")"
echo "PRIMARY PCA SAMPLES:     $all_count"
echo "SENSITIVITY PCA SAMPLES: $no65_count"
echo "SAME ALL-23 ALLELE FREQUENCIES USED FOR BOTH ANALYSES"
echo "PRIMARY23 PCA V3: COMPLETE"
