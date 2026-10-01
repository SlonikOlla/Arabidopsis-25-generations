#!/usr/bin/env bash
#SBATCH --job-name=FilterPCA_P23
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/UVC/primary23/filter_PCA_P23_%j.out
#SBATCH --error=logs/UVC/primary23/filter_PCA_P23_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
REF="${ROOT}/reference/TAIR10/TAIR10_chr_all.fa"
RAW="${ROOT}/results/wgs_variants/joint_primary23/PRIMARY23_chr1-5.vcf.gz"
FILTERDIR="${ROOT}/results/wgs_variants/PRIMARY23_filtered"
PCADIR="${ROOT}/results/wgs_variants/PRIMARY23_PCA"

SNP_RAW="${FILTERDIR}/PRIMARY23_SNP_raw.vcf.gz"
SNP_FLAGGED="${FILTERDIR}/PRIMARY23_SNP_hardfilter_flagged.vcf.gz"
SNP_PASS="${FILTERDIR}/PRIMARY23_SNP_PASS_biallelic.vcf.gz"
SNP_MASKED="${FILTERDIR}/PRIMARY23_SNP_PASS_biallelic_GTmasked.vcf.gz"
SNP_TAGGED="${FILTERDIR}/PRIMARY23_SNP_PASS_biallelic_GTmasked_tags.vcf.gz"
PCA_VCF="${FILTERDIR}/PRIMARY23_PCA_SNPs.vcf.gz"

PLINK_PREFIX="${PCADIR}/PRIMARY23_PCA_input"
LD_PREFIX="${PCADIR}/PRIMARY23_PCA_LD"
PCA_ALL="${PCADIR}/PRIMARY23_all23"
PCA_NO65="${PCADIR}/PRIMARY23_without_SRR17867665"
REMOVE65="${PCADIR}/remove_SRR17867665.txt"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load gatk/4.6.1.0
module load bcftools/1.22
module load plink/2.0.0-a.6.32

export JAVA_TOOL_OPTIONS="-Xmx12g"

echo "========================================"
echo "PRIMARY23 SNP FILTERING AND PCA"
echo "START: $(date)"
echo "HOST:  $(hostname)"
echo "========================================"

echo "===== SOFTWARE ====="
gatk --version
bcftools --version | head -1
plink2 --version | head -1

echo "===== INPUT VALIDATION ====="

test -s "$REF"
test -s "${REF}.fai"
test -s "$RAW"
test -s "${RAW}.tbi"

gzip -t "$RAW"

raw_samples=$(bcftools query -l "$RAW" | wc -l)
raw_variants=$(bcftools view -H "$RAW" | wc -l)

echo "RAW SAMPLES:  $raw_samples"
echo "RAW VARIANTS: $raw_variants"

if [[ "$raw_samples" -ne 23 ]]; then
    echo "ERROR: expected 23 samples"
    exit 1
fi

if [[ "$raw_variants" -ne 51173 ]]; then
    echo "ERROR: expected 51173 raw variants"
    exit 1
fi

echo "INPUT VALIDATION: PASS"

echo "===== EXTRACTING SNP RECORDS ====="

gatk --java-options \
  "-Xmx12g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
  SelectVariants \
  -R "$REF" \
  -V "$RAW" \
  --select-type-to-include SNP \
  -O "$SNP_RAW"

test -s "$SNP_RAW"
gzip -t "$SNP_RAW"

echo "RAW SNP RECORDS: $(bcftools view -H "$SNP_RAW" | wc -l)"

echo "===== GATK SNP HARD FILTERS ====="

gatk --java-options \
  "-Xmx12g -Djava.io.tmpdir=${SLURM_TMPDIR}" \
  VariantFiltration \
  -R "$REF" \
  -V "$SNP_RAW" \
  -O "$SNP_FLAGGED" \
  --filter-name "QD2" \
  --filter-expression "QD < 2.0" \
  --filter-name "QUAL30" \
  --filter-expression "QUAL < 30.0" \
  --filter-name "SOR3" \
  --filter-expression "SOR > 3.0" \
  --filter-name "FS60" \
  --filter-expression "FS > 60.0" \
  --filter-name "MQ40" \
  --filter-expression "MQ < 40.0" \
  --filter-name "MQRankSum-12.5" \
  --filter-expression "MQRankSum < -12.5" \
  --filter-name "ReadPosRankSum-8" \
  --filter-expression "ReadPosRankSum < -8.0"

test -s "$SNP_FLAGGED"
gzip -t "$SNP_FLAGGED"

echo "===== SITE-FILTER INVENTORY ====="

bcftools query -f '%FILTER\n' "$SNP_FLAGGED" |
sort |
uniq -c |
sort -nr \
  > "${FILTERDIR}/PRIMARY23_SNP_hardfilter_inventory.tsv"

cat "${FILTERDIR}/PRIMARY23_SNP_hardfilter_inventory.tsv"

echo "===== RETAINING PASS BIALLELIC SNPS ====="

bcftools view \
  --threads 4 \
  -f PASS \
  -m2 \
  -M2 \
  -v snps \
  -Oz \
  -o "$SNP_PASS" \
  "$SNP_FLAGGED"

bcftools index --threads 4 -t "$SNP_PASS"

test -s "$SNP_PASS"
test -s "${SNP_PASS}.tbi"
gzip -t "$SNP_PASS"

pass_count=$(bcftools view -H "$SNP_PASS" | wc -l)

echo "PASS BIALLELIC SNP RECORDS: $pass_count"

if [[ "$pass_count" -eq 0 ]]; then
    echo "ERROR: no SNPs passed site filtering"
    exit 1
fi

echo "===== MASKING LOW-QUALITY GENOTYPES ====="
echo "Criteria: DP < 8 or GQ < 20"

bcftools filter \
  --threads 4 \
  -S . \
  -e 'FMT/DP<8 || FMT/GQ<20' \
  -Oz \
  -o "$SNP_MASKED" \
  "$SNP_PASS"

bcftools index --threads 4 -t "$SNP_MASKED"

test -s "$SNP_MASKED"
test -s "${SNP_MASKED}.tbi"
gzip -t "$SNP_MASKED"

echo "===== RECALCULATING COHORT TAGS ====="

bcftools +fill-tags \
  "$SNP_MASKED" \
  -Oz \
  -o "$SNP_TAGGED" \
  -- \
  -t AC,AN,AF,MAF,NS,F_MISSING

bcftools index --threads 4 -t "$SNP_TAGGED"

test -s "$SNP_TAGGED"
test -s "${SNP_TAGGED}.tbi"
gzip -t "$SNP_TAGGED"

echo "===== BUILDING PCA SNP SET ====="
echo "Criteria: missingness <= 20%; MAF >= 0.05"

bcftools view \
  --threads 4 \
  -i 'INFO/F_MISSING<=0.20 && INFO/MAF>=0.05' \
  -Oz \
  -o "$PCA_VCF" \
  "$SNP_TAGGED"

bcftools index --threads 4 -t "$PCA_VCF"

test -s "$PCA_VCF"
test -s "${PCA_VCF}.tbi"
gzip -t "$PCA_VCF"

pca_vcf_samples=$(bcftools query -l "$PCA_VCF" | wc -l)
pca_vcf_variants=$(bcftools view -H "$PCA_VCF" | wc -l)

echo "PCA VCF SAMPLES:  $pca_vcf_samples"
echo "PCA VCF VARIANTS: $pca_vcf_variants"

if [[ "$pca_vcf_samples" -ne 23 ]]; then
    echo "ERROR: PCA VCF does not contain 23 samples"
    exit 1
fi

if [[ "$pca_vcf_variants" -lt 100 ]]; then
    echo "ERROR: too few SNPs for PCA"
    exit 1
fi

echo "===== PLINK CONVERSION ====="

plink2 \
  --vcf "$PCA_VCF" \
  --double-id \
  --allow-extra-chr \
  --set-missing-var-ids '@:#:$r:$a' \
  --make-pgen \
  --out "$PLINK_PREFIX"

test -s "${PLINK_PREFIX}.pgen"
test -s "${PLINK_PREFIX}.pvar"
test -s "${PLINK_PREFIX}.psam"

echo "===== LD PRUNING ====="

plink2 \
  --pfile "$PLINK_PREFIX" \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.2 \
  --out "$LD_PREFIX"

test -s "${LD_PREFIX}.prune.in"

pruned_count=$(wc -l < "${LD_PREFIX}.prune.in")

echo "LD-PRUNED SNP COUNT: $pruned_count"

if [[ "$pruned_count" -lt 100 ]]; then
    echo "ERROR: too few LD-pruned SNPs for PCA"
    exit 1
fi

echo "===== PCA: ALL 23 SAMPLES ====="

plink2 \
  --pfile "$PLINK_PREFIX" \
  --allow-extra-chr \
  --extract "${LD_PREFIX}.prune.in" \
  --pca 10 allele-wts \
  --out "$PCA_ALL"

test -s "${PCA_ALL}.eigenvec"
test -s "${PCA_ALL}.eigenval"
test -s "${PCA_ALL}.eigenvec.allele"

all_pca_samples=$(
  awk 'NR>1 {n++} END {print n+0}' \
    "${PCA_ALL}.eigenvec"
)

echo "ALL-SAMPLE PCA RECORDS: $all_pca_samples"

if [[ "$all_pca_samples" -ne 23 ]]; then
    echo "ERROR: all-sample PCA does not contain 23 samples"
    exit 1
fi

echo "===== PCA SENSITIVITY: EXCLUDE SRR17867665 ====="

printf "SRR17867665\tSRR17867665\n" > "$REMOVE65"

plink2 \
  --pfile "$PLINK_PREFIX" \
  --allow-extra-chr \
  --extract "${LD_PREFIX}.prune.in" \
  --remove "$REMOVE65" \
  --pca 10 allele-wts \
  --out "$PCA_NO65"

test -s "${PCA_NO65}.eigenvec"
test -s "${PCA_NO65}.eigenval"
test -s "${PCA_NO65}.eigenvec.allele"

no65_pca_samples=$(
  awk 'NR>1 {n++} END {print n+0}' \
    "${PCA_NO65}.eigenvec"
)

echo "NO-SRR17867665 PCA RECORDS: $no65_pca_samples"

if [[ "$no65_pca_samples" -ne 22 ]]; then
    echo "ERROR: sensitivity PCA does not contain 22 samples"
    exit 1
fi

echo "===== FINAL VALIDATION ====="

echo "RAW JOINT VARIANTS:       $raw_variants"
echo "PASS BIALLELIC SNPS:      $pass_count"
echo "PCA-ELIGIBLE SNPS:        $pca_vcf_variants"
echo "LD-PRUNED PCA SNPS:       $pruned_count"
echo "ALL-SAMPLE PCA SAMPLES:   $all_pca_samples"
echo "SENSITIVITY PCA SAMPLES:  $no65_pca_samples"

echo "RAW JOINT VCF RETAINED"
echo "SRR17867665 RETAINED IN PRIMARY PCA"
echo "PRIMARY23 FILTERING AND PCA: COMPLETE"
echo "FINISH: $(date)"
