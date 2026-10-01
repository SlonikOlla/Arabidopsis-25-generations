
# Workflow map

## WGS

1. Trim and align UVC WGS reads.
2. Run HaplotypeCaller and joint genotyping.
3. Apply site filters and genotype masking.
4. Audit recurrent heterozygous calls.
5. Build homozygous-only, thinned marker set for PCA and IBS distances.
6. Call matched-background UVC candidates.
7. Validate candidates independently and audit specificity.
8. Summarize spectra, context, annotations, and candidate priorities.

## WGBS

1. Trim, align, deduplicate, and extract methylation calls.
2. Build 100-bp and 500-bp methylation summaries.
3. Apply sample QC and define the primary cohort.
4. Construct context-specific PCA/correlation matrices.
5. Classify background-stratified strict methylation states.
6. Fit background-adjusted treatment models.
7. Test interaction sensitivity and feature enrichment.

## Cross-stress comparison

1. Harmonize eligible CG loci across heat, cold, and UVC analyses.
2. Intersect switch sets and audit direction concordance.
3. Test overlap enrichment using eligible-locus universes.
4. Annotate recurrent loci and genes.
