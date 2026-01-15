#!/bin/bash
# =============================================================================
# VCF Preprocessing Pipeline
# 
# Normalizes, indexes, and merges 1000 Genomes VCF files into PLINK2 format.
# =============================================================================

set -euo pipefail

# Configuration
VCF_DIR="${1:-1000genomes}"
NORM_DIR="${VCF_DIR}/normed"
PLINK_DIR="${VCF_DIR}/plink"
MERGE_LIST="${PLINK_DIR}/merge_list.txt"
OUTPUT_PREFIX="${2:-allchr}"

# Create directories
mkdir -p "$NORM_DIR"
mkdir -p "$PLINK_DIR"

echo "=== VCF Preprocessing Pipeline ==="
echo "VCF Directory: $VCF_DIR"
echo "Output Prefix: $OUTPUT_PREFIX"
echo ""

# -----------------------------------------------------------------------------
# Step 1: Normalize VCFs (split multiallelics, left-align indels)
# -----------------------------------------------------------------------------
normalize_vcfs() {
    echo "Step 1: Normalizing VCFs..."
    
    for vcf in "$VCF_DIR"/ALL.chr*.vcf.gz; do
        [ -f "$vcf" ] || continue
        
        chr=$(basename "$vcf" .vcf.gz)
        norm_vcf="$NORM_DIR/${chr}.norm.vcf.gz"
        
        if [ -f "$norm_vcf" ]; then
            echo "  Skipping $chr (already normalized)"
            continue
        fi
        
        echo "  Normalizing: $chr"
        bcftools norm -m -both "$vcf" -Oz -o "$norm_vcf"
        tabix -p vcf "$norm_vcf"
    done
    
    echo "  Done normalizing."
}

# -----------------------------------------------------------------------------
# Step 2: Convert VCFs to PLINK2 format
# -----------------------------------------------------------------------------
convert_to_plink() {
    echo "Step 2: Converting to PLINK2 format..."
    
    for norm_vcf in "$NORM_DIR"/ALL.chr*.norm.vcf.gz; do
        [ -f "$norm_vcf" ] || continue
        
        chr=$(basename "$norm_vcf" .norm.vcf.gz)
        out_prefix="$PLINK_DIR/$chr"
        
        if [ -f "${out_prefix}.pgen" ]; then
            echo "  Skipping $chr (already converted)"
            continue
        fi
        
        echo "  Converting: $chr"
        plink2 --vcf "$norm_vcf" \
               --make-pgen \
               --out "$out_prefix"
    done
    
    echo "  Done converting."
}

# -----------------------------------------------------------------------------
# Step 3: Merge all chromosomes
# -----------------------------------------------------------------------------
merge_chromosomes() {
    echo "Step 3: Merging chromosomes..."
    
    # Find first chromosome file
    first_pgen=$(ls "$PLINK_DIR"/ALL.chr1*.pgen 2>/dev/null | head -1)
    
    if [ -z "$first_pgen" ]; then
        echo "  Error: No chromosome 1 pgen file found"
        exit 1
    fi
    
    first_prefix="${first_pgen%.pgen}"
    
    # Create merge list (all chromosomes except the first)
    rm -f "$MERGE_LIST"
    
    for pgen in "$PLINK_DIR"/ALL.chr*.pgen; do
        prefix="${pgen%.pgen}"
        if [ "$prefix" != "$first_prefix" ]; then
            echo "$prefix" >> "$MERGE_LIST"
        fi
    done
    
    # Count files to merge
    n_files=$(wc -l < "$MERGE_LIST" 2>/dev/null || echo "0")
    echo "  Merging $((n_files + 1)) chromosome files..."
    
    # Merge
    plink2 --pfile "$first_prefix" \
           --pmerge-list "$MERGE_LIST" \
           --set-all-var-ids '@:#:\$r:\$a' \
           --new-id-max-allele-len 662 \
           --make-pgen \
           --out "$PLINK_DIR/$OUTPUT_PREFIX"
    
    echo "  Done merging."
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main() {
    # Check dependencies
    command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found"; exit 1; }
    command -v plink2 >/dev/null 2>&1 || { echo "Error: plink2 not found"; exit 1; }
    command -v tabix >/dev/null 2>&1 || { echo "Error: tabix not found"; exit 1; }
    
    # Run pipeline
    normalize_vcfs
    convert_to_plink
    merge_chromosomes
    
    echo ""
    echo "=== Pipeline Complete ==="
    echo "Output: $PLINK_DIR/$OUTPUT_PREFIX.{pgen,pvar,psam}"
}

# Run if executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
