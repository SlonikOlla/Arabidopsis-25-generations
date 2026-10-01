#!/usr/bin/env bash
#SBATCH --job-name=PCAfull_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/UVC/primary23/PCA_full_%j.out
#SBATCH --error=logs/UVC/primary23/PCA_full_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
PCADIR="${ROOT}/results/wgs_variants/PRIMARY23_PCA"

PLINK_PREFIX="${PCADIR}/PRIMARY23_PCA_input"
THIN_LIST="${PCADIR}/PRIMARY23_PCA_5kb_SNPs_v2.txt"
FREQ="${PCADIR}/PRIMARY23_all23_5kb_frequencies.afreq"
REMOVE65="${PCADIR}/remove_SRR17867665_FID_IID.txt"

PCA_ALL="${PCADIR}/PRIMARY23_all23_5kb_full"
PCA_NO65="${PCADIR}/PRIMARY23_without_SRR17867665_5kb_full"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load plink/2.0.0-a.6.32

echo "========================================"
echo "PRIMARY23 FULL PCA"
echo "START: $(date)"
echo "HOST:  $(hostname)"
echo "========================================"

echo "===== SOFTWARE ====="
plink2 --version

echo "===== INPUT VALIDATION ====="

for f in \
    "${PLINK_PREFIX}.pgen" \
    "${PLINK_PREFIX}.pvar" \
    "${PLINK_PREFIX}.psam" \
    "$THIN_LIST" \
    "$FREQ" \
    "$REMOVE65"
do
    test -s "$f"
    printf "PASS\t%s\t%s bytes\n" "$f" "$(stat -c '%s' "$f")"
done

thin_count=$(wc -l < "$THIN_LIST")
freq_count=$(awk 'NR>1 {n++} END {print n+0}' "$FREQ")

echo "THINNED SNPS:     $thin_count"
echo "FREQUENCY RECORDS: $freq_count"

test "$thin_count" -eq 608
test "$freq_count" -eq 608

echo "===== FULL PCA: ALL 23 SAMPLES ====="

plink2 \
    --pfile "$PLINK_PREFIX" \
    --allow-extra-chr \
    --extract "$THIN_LIST" \
    --read-freq "$FREQ" \
    --pca 22 allele-wts \
    --threads 2 \
    --out "$PCA_ALL"

test -s "${PCA_ALL}.eigenvec"
test -s "${PCA_ALL}.eigenval"
test -s "${PCA_ALL}.eigenvec.allele"

all_samples=$(
    awk 'NR>1 {n++} END {print n+0}' \
        "${PCA_ALL}.eigenvec"
)

all_pcs=$(wc -l < "${PCA_ALL}.eigenval")

echo "ALL-SAMPLE PCA RECORDS: $all_samples"
echo "ALL-SAMPLE PCS:         $all_pcs"

test "$all_samples" -eq 23
test "$all_pcs" -eq 22

echo "===== FULL PCA: WITHOUT SRR17867665 ====="

plink2 \
    --pfile "$PLINK_PREFIX" \
    --allow-extra-chr \
    --extract "$THIN_LIST" \
    --remove "$REMOVE65" \
    --read-freq "$FREQ" \
    --pca 21 allele-wts \
    --threads 2 \
    --out "$PCA_NO65"

test -s "${PCA_NO65}.eigenvec"
test -s "${PCA_NO65}.eigenval"
test -s "${PCA_NO65}.eigenvec.allele"

no65_samples=$(
    awk 'NR>1 {n++} END {print n+0}' \
        "${PCA_NO65}.eigenvec"
)

no65_pcs=$(wc -l < "${PCA_NO65}.eigenval")

echo "SENSITIVITY PCA RECORDS: $no65_samples"
echo "SENSITIVITY PCS:         $no65_pcs"

test "$no65_samples" -eq 22
test "$no65_pcs" -eq 21

if awk 'NR>1 {print $1; print $2}' \
    "${PCA_NO65}.eigenvec" |
    grep -Fxq 'SRR17867665'
then
    echo "ERROR: SRR17867665 remains in sensitivity PCA"
    exit 1
else
    echo "SRR17867665 ABSENT: PASS"
fi

echo "===== VARIANCE ACCOUNTING ====="

awk '
{
    eigenvalue[NR]=$1
    total+=$1
}
END {
    printf "component\teigenvalue\tvariance_percent\n"
    for (i=1; i<=NR; i++) {
        printf "PC%d\t%.10f\t%.6f\n",
               i, eigenvalue[i],
               100*eigenvalue[i]/total
    }
}
' "${PCA_ALL}.eigenval" \
> "${PCA_ALL}.variance.tsv"

awk '
{
    eigenvalue[NR]=$1
    total+=$1
}
END {
    printf "component\teigenvalue\tvariance_percent\n"
    for (i=1; i<=NR; i++) {
        printf "PC%d\t%.10f\t%.6f\n",
               i, eigenvalue[i],
               100*eigenvalue[i]/total
    }
}
' "${PCA_NO65}.eigenval" \
> "${PCA_NO65}.variance.tsv"

test -s "${PCA_ALL}.variance.tsv"
test -s "${PCA_NO65}.variance.tsv"

echo "===== ALL-23 VARIANCE: FIRST 10 PCS ====="
head -11 "${PCA_ALL}.variance.tsv"

echo "===== SENSITIVITY VARIANCE: FIRST 10 PCS ====="
head -11 "${PCA_NO65}.variance.tsv"

echo "PRIMARY23 FULL PCA: COMPLETE"
echo "FINISH: $(date)"
