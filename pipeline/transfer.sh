#!/bin/bash
# =============================================================================
# transfer_cyp2d6.sh  — run this ON THE HPC before copying to laptop
# Copies only the two final CYP2D6 analysis VCFs to a transfer folder,
# cleans up duplicate junk files, and confirms what to copy.
# =============================================================================

TRANSFER_DIR=/etc/ace-data/home/rnagaba/reset/transfer_cyp2d6
mkdir -p "$TRANSFER_DIR"

AFR=/etc/ace-data/home/rnagaba/reset/variants/AFR/codeine/CYP2D6
EUR=/etc/ace-data/home/rnagaba/reset/variants/EUR/codeine/CYP2D6

echo "=== Cleaning up duplicate files ==="
# Remove junk files created by previous runs
for DIR in "$AFR" "$EUR"; do
    rm -fv "$DIR"/*.vep.analysis.vcf.gz
    rm -fv "$DIR"/*.vep.analysis.vcf.gz.csi
    rm -fv "$DIR"/*.vep.vep.analysis.vcf.gz
    rm -fv "$DIR"/*.vep.vep.analysis.vcf.gz.csi
done

echo ""
echo "=== Copying final analysis files to transfer folder ==="
cp -v "$AFR/AFR_CYP2D6.analysis.vcf.gz"     "$TRANSFER_DIR/"
cp -v "$AFR/AFR_CYP2D6.analysis.vcf.gz.csi" "$TRANSFER_DIR/"
cp -v "$EUR/EUR_CYP2D6.analysis.vcf.gz"     "$TRANSFER_DIR/"
cp -v "$EUR/EUR_CYP2D6.analysis.vcf.gz.csi" "$TRANSFER_DIR/"

echo ""
echo "=== Files ready to transfer to laptop ==="
ls -lh "$TRANSFER_DIR"

echo ""
echo "=== Copy command to run on your LAPTOP (not HPC) ==="
echo "scp rnagaba@kla-ac-hpc-02:$TRANSFER_DIR/* 'D:/RITAH/SCHOOL/Beast Mode/analysis/'"
