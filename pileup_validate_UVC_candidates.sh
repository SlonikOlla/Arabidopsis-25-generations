#!/usr/bin/env bash
#SBATCH --job-name=UVC_pileup
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/mutation_validation/UVC_pileup_%j.out
#SBATCH --error=logs/UVC/mutation_validation/UVC_pileup_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
REF="${ROOT}/reference/TAIR10/TAIR10_chr_all.fa"
BED="${ROOT}/results/wgs_variants/PRIMARY23_mutation_analysis/PRIMARY23_UVC_sample_specific_candidates.bed"
OUTDIR="${ROOT}/results/wgs_variants/PRIMARY23_mutation_analysis"
OUT="${OUTDIR}/PRIMARY23_UVC_candidate_independent_pileup.vcf.gz"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load bcftools/1.22

bams=(
    results/wgs_alignment/UVC/F25UV_7_WGS.markdup.bam
    results/wgs_alignment/UVC/F25UV_8_WGS.markdup.bam
    results/wgs_alignment/UVC/F25UV_10_WGS.markdup.bam
    results/wgs_alignment/UVC/F25UV_11_WGS.markdup.bam
    results/wgs_alignment/UVC/F25UV_13_WGS.markdup.bam
    results/wgs_sv/bam/SRR17867641.bam
    results/wgs_sv/bam/SRR17867643.bam
    results/wgs_sv/bam/SRR17867645.bam
)

echo "===== INPUT VALIDATION ====="

test -s "$REF"
test -s "${REF}.fai"
test -s "$BED"

for bam in "${bams[@]}"; do
    test -s "$bam"
    test -s "${bam}.bai"
    echo "PASS  $bam"
done

echo "===== INDEPENDENT MPILEUP AND CALLING ====="

bcftools mpileup \
    --threads 4 \
    -f "$REF" \
    -R "$BED" \
    -q 30 \
    -Q 20 \
    -d 10000 \
    -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR \
    -Ou \
    "${bams[@]}" |
bcftools call \
    --threads 4 \
    -m \
    -A \
    -Oz \
    -o "$OUT"

bcftools index --threads 4 -t "$OUT"

test -s "$OUT"
test -s "${OUT}.tbi"
gzip -t "$OUT"

echo "===== OUTPUT VALIDATION ====="
echo "SAMPLES: $(bcftools query -l "$OUT" | wc -l)"
echo "RECORDS: $(bcftools view -H "$OUT" | wc -l)"

echo "===== SAMPLE ORDER ====="
bcftools query -l "$OUT"

echo "UVC CANDIDATE INDEPENDENT PILEUP: COMPLETE"
echo "OUTPUT: $OUT"
