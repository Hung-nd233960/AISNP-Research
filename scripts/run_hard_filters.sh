#!/bin/bash
# ==============================================================================
# HARD FILTERS SCRIPT
# ==============================================================================
# Run strict quality control filters:
#   1. SNP-only (remove indels, CNVs)
#   2. Biallelic (max 2 alleles)
#   3. Sample subsetting (optional)
#   4. MAF filtering (remove rare variants)
#
# These filters should be applied consistently to ALL datasets.
# ==============================================================================

set -e  # Exit on error

# Configuration
THREADS=8
INPUT_VCF="1000genomes/main_vcf/ALL_merged.vcf.gz"
SAMPLE_LIST="1000genomes/EAS_subpopulation_samples_list.csv"
OUTPUT_DIR="1000genomes/output"
MIN_AF=0.0016  # 1/612 alleles for 306 samples

# Create output directory
mkdir -p $OUTPUT_DIR

echo "=========================================="
echo "RUNNING HARD FILTERS"
echo "=========================================="
echo "  Input VCF: $INPUT_VCF"
echo "  Sample list: $SAMPLE_LIST"
echo "  Min AF: $MIN_AF"
echo "  Threads: $THREADS"
echo "=========================================="

# ------------------------------------------------------------------------------
# HARD FILTER 1: SNP-only + Biallelic + Sample Subsetting
# ------------------------------------------------------------------------------
echo ""
echo ">>> STEP 1: SNP-only and Biallelic Filter"
echo "    Removing indels, CNVs, and multi-allelic sites"
echo ""

plink2 --vcf $INPUT_VCF \
       --snps-only \
       --max-alleles 2 \
       --keep $SAMPLE_LIST \
       --make-pgen \
       --out $OUTPUT_DIR/EAS_AND_SNP_filtered_data \
       --threads $THREADS

echo ""
echo ">>> STEP 1 Complete"
echo "    Output: $OUTPUT_DIR/EAS_AND_SNP_filtered_data.*"

# ------------------------------------------------------------------------------
# Calculate frequency statistics before MAF filter
# ------------------------------------------------------------------------------
echo ""
echo ">>> Calculating pre-MAF frequency statistics..."
echo ""

plink2 --pfile $OUTPUT_DIR/EAS_AND_SNP_filtered_data \
       --freq \
       --out $OUTPUT_DIR/EAS_AND_SNP_filtered_data_info \
       --threads $THREADS

# ------------------------------------------------------------------------------
# HARD FILTER 2: Minor Allele Frequency (MAF)
# ------------------------------------------------------------------------------
echo ""
echo ">>> STEP 2: MAF Filter"
echo "    Minimum allele frequency: $MIN_AF"
echo ""

plink2 --pfile $OUTPUT_DIR/EAS_AND_SNP_filtered_data \
       --min-af $MIN_AF \
       --make-pgen \
       --out $OUTPUT_DIR/EAS_AND_SNP_filtered_data_MAF_filtered \
       --threads $THREADS

echo ""
echo ">>> STEP 2 Complete"
echo "    Output: $OUTPUT_DIR/EAS_AND_SNP_filtered_data_MAF_filtered.*"

# ------------------------------------------------------------------------------
# Calculate frequency statistics after MAF filter
# ------------------------------------------------------------------------------
echo ""
echo ">>> Calculating post-MAF frequency statistics..."
echo ""

plink2 --pfile $OUTPUT_DIR/EAS_AND_SNP_filtered_data_MAF_filtered \
       --freq \
       --out $OUTPUT_DIR/EAS_AND_SNP_filtered_data_MAF_filtered_info \
       --threads $THREADS

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo ""
echo "=========================================="
echo "HARD FILTERS COMPLETE"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - $OUTPUT_DIR/EAS_AND_SNP_filtered_data.* (after SNP/biallelic filter)"
echo "  - $OUTPUT_DIR/EAS_AND_SNP_filtered_data_MAF_filtered.* (after MAF filter)"
echo ""
echo "Frequency statistics:"
echo "  - $OUTPUT_DIR/EAS_AND_SNP_filtered_data_info.afreq"
echo "  - $OUTPUT_DIR/EAS_AND_SNP_filtered_data_MAF_filtered_info.afreq"
echo ""

# Count variants at each stage
echo "Variant counts:"
echo -n "  After SNP/biallelic: "
grep -v "^#" $OUTPUT_DIR/EAS_AND_SNP_filtered_data.pvar | wc -l

echo -n "  After MAF filter:    "
grep -v "^#" $OUTPUT_DIR/EAS_AND_SNP_filtered_data_MAF_filtered.pvar | wc -l

echo ""
echo "Next step: Run situational filters (HWE, LD pruning)"
echo "=========================================="
