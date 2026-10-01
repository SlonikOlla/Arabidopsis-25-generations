#!/usr/bin/env bash
#SBATCH --job-name=FilterMut_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/UVC/mutation_analysis/filter_mutation_%j.out
#SBATCH --error=logs/UVC/mutation_analysis/filter_mutation_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
REF="${ROOT}/reference/TAIR10/TAIR10_chr_all.fa"
RAW="${ROOT}/results/wgs_variants/joint_primary23/PRIMARY23_chr1-5.vcf.gz"
OUTDIR="${ROOT}/results/wgs_variants/PRIMARY23_mutation_analysis"

FLAGGED="${OUTDIR}/PRIMARY23_SNP_INDEL_hardfilter_flagged.vcf.gz"
PASS="${OUTDIR}/PRIMARY23_SNP_INDEL_PASS_biallelic.vcf.gz"
MASKED="${OUTDIR}/PRIMARY23_SNP_INDEL_PASS_biallelic_GTmasked.vcf.gz"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load gatk/4.6.1.0
module load bcftools/1.22

echo "===== INPUT VALIDATION ====="
test -s "$REF"
test -s "${REF}.fai"
test -s "$RAW"
test -s "${RAW}.tbi"
echo "INPUTS: PASS"

echo "===== SNP AND INDEL HARD FILTERING ====="

gatk --java-options \
  "-Xmx12g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
  VariantFiltration \
  -R "$REF" \
  -V "$RAW" \
  -O "$FLAGGED" \
  --filter-name "QUAL30" \
  --filter-expression "QUAL < 30.0" \
  --filter-name "SNP_QD2" \
  --filter-expression "vc.isSNP() && QD < 2.0" \
  --filter-name "SNP_SOR3" \
  --filter-expression "vc.isSNP() && SOR > 3.0" \
  --filter-name "SNP_FS60" \
  --filter-expression "vc.isSNP() && FS > 60.0" \
  --filter-name "SNP_MQ40" \
  --filter-expression "vc.isSNP() && MQ < 40.0" \
  --filter-name "SNP_MQRankSum" \
  --filter-expression "vc.isSNP() && MQRankSum < -12.5" \
  --filter-name "SNP_ReadPosRankSum" \
  --filter-expression "vc.isSNP() && ReadPosRankSum < -8.0" \
  --filter-name "INDEL_QD2" \
  --filter-expression "vc.isIndel() && QD < 2.0" \
  --filter-name "INDEL_FS200" \
  --filter-expression "vc.isIndel() && FS > 200.0" \
  --filter-name "INDEL_SOR10" \
  --filter-expression "vc.isIndel() && SOR > 10.0" \
  --filter-name "INDEL_ReadPosRankSum" \
  --filter-expression "vc.isIndel() && ReadPosRankSum < -20.0"

test -s "$FLAGGED"
gzip -t "$FLAGGED"

echo "===== RETAINING PASS BIALLELIC SNPs AND INDELs ====="

bcftools view \
  --threads 4 \
  -f PASS \
  -m2 \
  -M2 \
  -v snps,indels \
  -Oz \
  -o "$PASS" \
  "$FLAGGED"

bcftools index --threads 4 -t "$PASS"

echo "===== MASKING LOW-QUALITY GENOTYPES ====="
echo "Criteria: DP < 8 or GQ < 20"

bcftools filter \
  --threads 4 \
  -S . \
  -e 'FMT/DP<8 || FMT/GQ<20' \
  -Oz \
  -o "$MASKED" \
  "$PASS"

bcftools index --threads 4 -t "$MASKED"

echo "===== OUTPUT VALIDATION ====="

for file in "$FLAGGED" "$PASS" "$MASKED"; do
    test -s "$file"
    test -s "${file}.tbi"
    gzip -t "$file"
    echo "PASS  $file"
done

printf "PASS SNPs:\t"
bcftools view -H -v snps "$PASS" | wc -l

printf "PASS INDELs:\t"
bcftools view -H -v indels "$PASS" | wc -l

printf "MASKED RECORDS:\t"
bcftools view -H "$MASKED" | wc -l

echo "PRIMARY23 MUTATION FILTERING: COMPLETE"
