#!/bin/bash
# Download official HPO genes_to_phenotype.txt from HPO website
# This script downloads the full HPO dataset for production use
# The included genes_to_phenotype.txt is a curated subset for development/testing

set -e

HPO_DIR="$(dirname "$0")/hpo"
mkdir -p "$HPO_DIR"

echo "Downloading HPO genes_to_phenotype.txt..."
wget -q -O "$HPO_DIR/genes_to_phenotype_full.txt" \
    "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"

if [ -f "$HPO_DIR/genes_to_phenotype_full.txt" ]; then
    LINES=$(wc -l < "$HPO_DIR/genes_to_phenotype_full.txt")
    echo "Downloaded successfully: $LINES lines"
    echo "File: $HPO_DIR/genes_to_phenotype_full.txt"
    echo ""
    echo "To use the full HPO dataset, update .env:"
    echo "  CARRIER_HPO_GENES_TO_PHENOTYPE=$HPO_DIR/genes_to_phenotype_full.txt"
else
    echo "ERROR: Download failed"
    exit 1
fi
