#!/usr/bin/env bash
#SBATCH --job-name=AnnotUVC24
#SBATCH --account=def-vasilisa_cpu
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/UVC/annotation/annotate_UVC24_%j.out
#SBATCH --error=logs/UVC/annotation/annotate_UVC24_%j.err

set -euo pipefail

ROOT="/scratch/vasilisa/arabidopsis_25gen"
INDIR="${ROOT}/results/wgs_variants/PRIMARY23_mutation_analysis"
INPUT="${INDIR}/PRIMARY23_UVC_independently_confirmed_SNPs.tsv"
PREFIX="${INDIR}/PRIMARY23_UVC_confirmed_24"
REF="${ROOT}/reference/TAIR10/TAIR10_chr_all.fa"

BUILD="${ROOT}/tools/snpeff_araport11"
CONFIG="${BUILD}/snpEff.config"

cd "$ROOT"

module --force purge
module load StdEnv/2023
module load bcftools/1.22
module load snpeff/5.2
SNPEFF_JAR="${EBROOTSNPEFF}/snpEff.jar"

unset JAVA_TOOL_OPTIONS

test -s "$INPUT"
test -s "$REF"
test -s "${BUILD}/data/Araport11/snpEffectPredictor.bin"

echo "===== BUILDING 24-SNP VCF ====="

{
    echo '##fileformat=VCFv4.2'
    echo '##reference=TAIR10'
    echo '##INFO=<ID=LINEAGE,Number=1,Type=String,Description="UVC genetic lineage">'
    echo '##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Focal UVC sample">'

    awk 'BEGIN{FS="\t"}
         {print "##contig=<ID="$1",length="$2">"}' \
      "${REF}.fai"

    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'

    awk '
    BEGIN {FS=OFS="\t"}
    NR>1 {
        printf "%s\t%s\tUVC%03d\t%s\t%s\t.\tPASS\tLINEAGE=%s;SAMPLE=%s\n",
               $1,$2,NR-1,$3,$4,$7,$8
    }
    ' "$INPUT"
} > "${PREFIX}.unsorted.vcf"

bcftools sort \
  -Oz \
  -o "${PREFIX}.sorted.vcf.gz" \
  "${PREFIX}.unsorted.vcf"

echo "===== REFERENCE-ALLELE VALIDATION ====="

bcftools norm \
  -f "$REF" \
  -c e \
  -Oz \
  -o "${PREFIX}.vcf.gz" \
  "${PREFIX}.sorted.vcf.gz"

bcftools index -t "${PREFIX}.vcf.gz"

records=$(bcftools view -H "${PREFIX}.vcf.gz" | wc -l)

echo "VCF RECORDS: $records"

if [[ "$records" -ne 24 ]]; then
    echo "ERROR: expected 24 confirmed SNPs"
    exit 1
fi

echo "REFERENCE ALLELES: PASS"

echo "===== SNPEFF ARAPORT11 ANNOTATION ====="

java -Xmx6g -jar "$SNPEFF_JAR" ann \
  -c "$CONFIG" \
  -dataDir "${BUILD}/data" \
  -ud 2000 \
  -stats "${PREFIX}_snpEff_summary.html" \
  -csvStats "${PREFIX}_snpEff_summary.csv" \
  Araport11 \
  "${PREFIX}.vcf.gz" \
  > "${PREFIX}_Araport11_annotated.vcf"

bgzip -f "${PREFIX}_Araport11_annotated.vcf"
tabix -f -p vcf "${PREFIX}_Araport11_annotated.vcf.gz"

echo "===== ANNOTATION VALIDATION ====="

annotated_records=$(
  bcftools view -H "${PREFIX}_Araport11_annotated.vcf.gz" |
  wc -l
)

records_with_ann=$(
  bcftools query -f '%INFO/ANN\n' \
    "${PREFIX}_Araport11_annotated.vcf.gz" |
  awk '$0!="." && $0!="" {n++} END {print n+0}'
)

echo "ANNOTATED VCF RECORDS: $annotated_records"
echo "RECORDS WITH ANN:      $records_with_ann"

if [[ "$annotated_records" -ne 24 || "$records_with_ann" -ne 24 ]]; then
    echo "ERROR: incomplete SnpEff annotation"
    exit 1
fi

bcftools query \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/LINEAGE\t%INFO/SAMPLE\t%INFO/ANN\n' \
  "${PREFIX}_Araport11_annotated.vcf.gz" \
  > "${PREFIX}_Araport11_ANN.tsv"

echo "===== OUTPUT FILES ====="

ls -lh \
  "${PREFIX}.vcf.gz" \
  "${PREFIX}.vcf.gz.tbi" \
  "${PREFIX}_Araport11_annotated.vcf.gz" \
  "${PREFIX}_Araport11_annotated.vcf.gz.tbi" \
  "${PREFIX}_Araport11_ANN.tsv" \
  "${PREFIX}_snpEff_summary.html" \
  "${PREFIX}_snpEff_summary.csv"

echo "UVC CONFIRMED SNP ANNOTATION: COMPLETE"
