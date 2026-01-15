#!/bin/bash
# ==============================================================================
# SITUATIONAL FILTERS SCRIPT
# ==============================================================================
# Run context-dependent filters:
#   1. Hardy-Weinberg Equilibrium (HWE)
#   2. Variant ID standardization
#   3. Linkage Disequilibrium (LD) pruning
#   4. FST calculation
#   5. PCA
#
# These filters may be adjusted based on study design and research questions.
# ==============================================================================

set -e  # Exit on error

# Configuration
THREADS=8
INPUT_PFILE="1000genomes/output/EAS_AND_SNP_filtered_data_MAF_filtered"
OUTPUT_DIR="1000genomes/output"

# Situational parameters (adjust as needed)
HWE_THRESHOLD="1e-6"      # Relaxed for population genetics
LD_WINDOW="1000kb"        # Aggressive for PCA
LD_STEP="1"
LD_R2="0.1"

echo "=========================================="
echo "RUNNING SITUATIONAL FILTERS"
echo "=========================================="
echo "  Input: $INPUT_PFILE"
echo "  HWE threshold: $HWE_THRESHOLD (SITUATIONAL)"
echo "  LD window: $LD_WINDOW (SITUATIONAL)"
echo "  LD R²: $LD_R2 (SITUATIONAL)"
echo "=========================================="
echo ""
echo "NOTE: These are SITUATIONAL filters."
echo "      Adjust parameters based on your study design."
echo ""

# ------------------------------------------------------------------------------
# SITUATIONAL FILTER 1: Hardy-Weinberg Equilibrium
# ------------------------------------------------------------------------------
echo ""
echo ">>> STEP 1: Calculate HWE Statistics"
echo ""

plink2 --pfile $INPUT_PFILE \
       --hardy \
       --out $OUTPUT_DIR/EAS_AND_SNP_filtered_data_MAF_filtered_hardy \
       --threads $THREADS

echo ""
echo ">>> STEP 2: Apply HWE Filter"
echo "    Threshold: P < $HWE_THRESHOLD"
echo ""

plink2 --pfile $INPUT_PFILE \
       --hwe $HWE_THRESHOLD \
       --make-pgen \
       --out $OUTPUT_DIR/EAS_SNP_MAF_HWE_filtered \
       --threads $THREADS

# ------------------------------------------------------------------------------
# SITUATIONAL: Set Variant IDs
# ------------------------------------------------------------------------------
echo ""
echo ">>> STEP 3: Set Variant IDs"
echo "    Format: chrom:pos[b37]ref,alt"
echo ""

plink2 --pfile $OUTPUT_DIR/EAS_SNP_MAF_HWE_filtered \
       --set-all-var-ids @:#[b37]\$r,\$a \
       --make-pgen \
       --out $OUTPUT_DIR/EAS_SNP_MAF_HWE_filtered_unique_ids \
       --threads $THREADS

# ------------------------------------------------------------------------------
# SITUATIONAL FILTER 2: LD Pruning (Multiple Configurations)
# ------------------------------------------------------------------------------
echo ""
echo ">>> STEP 4: Calculate LD Pruning Lists"
echo "    Testing multiple LD configurations..."
echo ""

# Aggressive pruning for PCA
echo "    Config 1: 1000kb window, step=1, r²=0.1 (aggressive)"
plink2 --pfile $OUTPUT_DIR/EAS_SNP_MAF_HWE_filtered_unique_ids \
       --indep-pairwise 1000kb 1 0.1 \
       --out $OUTPUT_DIR/EAS_SNP_MAF_HWE_LD_PRUNED_1000kb_1_0.1 \
       --threads $THREADS

# Moderate pruning
echo "    Config 2: 500kb window, step=1, r²=0.2 (moderate)"
plink2 --pfile $OUTPUT_DIR/EAS_SNP_MAF_HWE_filtered_unique_ids \
       --indep-pairwise 500kb 1 0.2 \
       --out $OUTPUT_DIR/EAS_SNP_MAF_HWE_LD_PRUNED_500kb_1_0.2 \
       --threads $THREADS

echo ""
echo ">>> STEP 5: Apply LD Pruning"
echo "    Using aggressive pruning (1000kb, r²=0.1)"
echo ""

plink2 --pfile $OUTPUT_DIR/EAS_SNP_MAF_HWE_filtered_unique_ids \
       --extract $OUTPUT_DIR/EAS_SNP_MAF_HWE_LD_PRUNED_1000kb_1_0.1.prune.in \
       --make-pgen \
       --out $OUTPUT_DIR/EAS_FINAL_DATA_FOR_FST \
       --threads $THREADS

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo ""
echo "=========================================="
echo "SITUATIONAL FILTERS COMPLETE"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - HWE stats: $OUTPUT_DIR/EAS_AND_SNP_filtered_data_MAF_filtered_hardy.hardy"
echo "  - After HWE: $OUTPUT_DIR/EAS_SNP_MAF_HWE_filtered.*"
echo "  - Unique IDs: $OUTPUT_DIR/EAS_SNP_MAF_HWE_filtered_unique_ids.*"
echo "  - LD pruned: $OUTPUT_DIR/EAS_FINAL_DATA_FOR_FST.*"
echo ""
echo "LD prune lists:"
echo "  - $OUTPUT_DIR/EAS_SNP_MAF_HWE_LD_PRUNED_1000kb_1_0.1.prune.in (aggressive)"
echo "  - $OUTPUT_DIR/EAS_SNP_MAF_HWE_LD_PRUNED_500kb_1_0.2.prune.in (moderate)"
echo ""

# Count variants at each stage
echo "Variant counts:"
echo -n "  After HWE:       "
grep -v "^#" $OUTPUT_DIR/EAS_SNP_MAF_HWE_filtered.pvar | wc -l

echo -n "  After LD prune:  "
grep -v "^#" $OUTPUT_DIR/EAS_FINAL_DATA_FOR_FST.pvar | wc -l

echo ""
echo "Next step: Add population labels and calculate FST"
echo "=========================================="
