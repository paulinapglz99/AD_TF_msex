#!/bin/bash
# ==============================================================================
# One-Command Pipeline: From Raw Data to Results
# ==============================================================================
# Usage: ./run_full_pipeline.sh /path/to/mayo/data cores
# Example: ./run_full_pipeline.sh /STORAGE/genut/alzheimer/data_by_counts/Mayo_counts/counts_by_tissue 8
# ==============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ==============================================================================
# ARGUMENT VALIDATION
# ==============================================================================

usage() {
    echo "Usage: $0 <source_data_dir> <cores>"
    echo ""
    echo "Arguments:"
    echo "  source_data_dir  Path to your counts_by_tissue directory"
    echo "  cores            Number of CPU cores to use (default: 4)"
    echo ""
    echo "Example:"
    echo "  $0 /STORAGE/genut/alzheimer/data_by_counts/Mayo_counts/counts_by_tissue 8"
}

if [ $# -lt 1 ]; then
    echo -e "${RED}Error: Missing arguments${NC}"
    usage
    exit 1
fi

SOURCE_DIR="$1"
CORES="${2:-4}"

# Validate source directory
if [ ! -d "$SOURCE_DIR" ]; then
    echo -e "${RED}Error: Source directory not found: $SOURCE_DIR${NC}"
    exit 1
fi

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}SEX-STRATIFIED RNA-SEQ PIPELINE${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo -e "${GREEN}Configuration:${NC}"
echo "  Source data: $SOURCE_DIR"
echo "  Cores: $CORES"
echo "  Current directory: $(pwd)"
echo ""

# ==============================================================================
# STEP 1: DETECT AVAILABLE TISSUES
# ==============================================================================

echo -e "${BLUE}[1/5] Detecting tissues in source directory...${NC}"

TISSUES=()
TISSUE_COUNT=0

for tissue_dir in "$SOURCE_DIR"/*; do
    if [ -d "$tissue_dir" ]; then
        tissue=$(basename "$tissue_dir")
        
        # Check if it has counts and metadata
        counts_file=$(find "$tissue_dir" -type f \( -name "*counts*.rds" -o -name "*counts*.csv" \) 2>/dev/null | head -1)
        metadata_file=$(find "$tissue_dir" -type f \( -name "*metadata*" \) 2>/dev/null | head -1)
        
        if [ ! -z "$counts_file" ] && [ ! -z "$metadata_file" ]; then
            TISSUES+=("$tissue")
            TISSUE_COUNT=$((TISSUE_COUNT + 1))
            echo -e "  ${GREEN}✓${NC} $tissue"
        fi
    fi
done

if [ $TISSUE_COUNT -eq 0 ]; then
    echo -e "${RED}Error: No tissues with both counts and metadata found${NC}"
    exit 1
fi

echo -e "${GREEN}Found $TISSUE_COUNT tissues${NC}"
echo ""

# ==============================================================================
# STEP 2: CREATE DIRECTORY STRUCTURE
# ==============================================================================

echo -e "${BLUE}[2/5] Setting up directory structure...${NC}"

mkdir -p data/counts data/metadata results logs scripts

echo -e "  ${GREEN}✓${NC} Created: data/counts"
echo -e "  ${GREEN}✓${NC} Created: data/metadata"
echo -e "  ${GREEN}✓${NC} Created: results"
echo ""

# ==============================================================================
# STEP 3: ADAPT DATA FORMAT (R SCRIPT)
# ==============================================================================

echo -e "${BLUE}[3/5] Adapting your data to pipeline format...${NC}"

# Create temporary R script for data adaptation
cat > /tmp/adapt_data_auto.R << 'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Get arguments from command line
args <- commandArgs(trailingOnly = TRUE)
source_dir <- args[1]
tissue_name <- args[2]

tissue_dir <- file.path(source_dir, tissue_name)

# Standardize tissue name for output
brain_area <- tolower(tissue_name)
brain_area <- gsub("_+", "_", brain_area)  # clean underscores

# Find files
counts_files <- list.files(tissue_dir, pattern = "count", 
                          full.names = TRUE, ignore.case = TRUE)
counts_files <- counts_files[!grepl("metadata", counts_files, ignore.case = TRUE)]
counts_file <- counts_files[1]

metadata_files <- list.files(tissue_dir, pattern = "metadata", 
                            full.names = TRUE, ignore.case = TRUE)
metadata_file <- metadata_files[1]

# Load counts
if (grepl("\\.rds$", counts_file, ignore.case = TRUE)) {
  counts <- readRDS(counts_file)
} else {
  counts <- as.matrix(read.table(counts_file, header = TRUE, 
                                row.names = 1, sep = "\t", check.names = FALSE))
}

# Load metadata
metadata <- read.table(metadata_file, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, check.names = FALSE)

# Align samples
common_samples <- intersect(colnames(counts), metadata[[1]])
counts <- counts[, common_samples]
metadata <- metadata[match(common_samples, metadata[[1]]), ]

# Create standardized metadata
std_meta <- data.frame(
  specimenID = metadata[[1]]
)

# Try to find and standardize columns
for (col in colnames(metadata)) {
  if (grepl("sex|gender", col, ignore.case = TRUE)) {
    vals <- tolower(metadata[[col]])
    vals[vals %in% c("m", "male", "1")] <- "male"
    vals[vals %in% c("f", "female", "2")] <- "female"
    std_meta$sex <- vals
  }
  if (grepl("diagnosis|disease|condition|phenotype", col, ignore.case = TRUE)) {
    vals <- tolower(metadata[[col]])
    vals[vals %in% c("control", "ctrl", "normal", "0")] <- "control"
    vals[vals %in% c("ad", "alzheimer", "1")] <- "AD"
    std_meta$diagnosis <- vals
  }
  if (grepl("source|batch|tissue", col, ignore.case = TRUE) && col != colnames(metadata)[1]) {
    std_meta[[tolower(col)]] <- metadata[[col]]
  }
}

# Save
saveRDS(counts, sprintf("data/counts/%s_counts.rds", brain_area))
write.csv(std_meta, sprintf("data/metadata/%s_metadata.csv", brain_area), 
         row.names = FALSE)

cat(sprintf("%s: %d genes × %d samples\n", brain_area, nrow(counts), ncol(counts)))
RSCRIPT

# Run R script for each tissue
for tissue in "${TISSUES[@]}"; do
    echo -n "  Processing $tissue..."
    if Rscript /tmp/adapt_data_auto.R "$SOURCE_DIR" "$tissue" 2>/dev/null; then
        echo -e " ${GREEN}✓${NC}"
    else
        echo -e " ${RED}✗${NC}"
    fi
done

echo ""

# ==============================================================================
# STEP 4: UPDATE CONFIG
# ==============================================================================

echo -e "${BLUE}[4/5] Updating configuration...${NC}"

# Build brain_areas list
BRAIN_AREAS_STR=$(printf '%s\n' "${TISSUES[@]}" | 
                  sed 's/^/  - /' | 
                  sed 's/$//')

# Replace in config.yaml if it exists
if [ -f "config.yaml" ]; then
    # Create temporary config
    cat > config.yaml.tmp << YAML_CONFIG
# Configuration for Sex-Stratified DESeq2 Analysis

# Brain areas (auto-detected)
brain_areas:
${BRAIN_AREAS_STR}

# Input/output directories
counts_dir: "data/counts"
metadata_dir: "data/metadata"
outdir: "results"

# Statistical thresholds
padj_threshold: 0.05
lfc_threshold: 0.5

# Computational resources
cores: $CORES
YAML_CONFIG

    mv config.yaml.tmp config.yaml
    echo -e "  ${GREEN}✓${NC} Updated config.yaml"
else
    echo -e "  ${YELLOW}⚠${NC} config.yaml not found (create it manually)"
fi

echo ""

# ==============================================================================
# STEP 5: RUN SNAKEMAKE PIPELINE
# ==============================================================================

echo -e "${BLUE}[5/5] Running Snakemake pipeline...${NC}"
echo ""

# Check if Snakemake is installed
if ! command -v snakemake &> /dev/null; then
    echo -e "${RED}Error: Snakemake not found. Install with: pip install snakemake${NC}"
    exit 1
fi

# Run pipeline
snakemake \
    --cores "$CORES" \
    --use-conda \
    --keep-going \
    -p

# Check exit code
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}PIPELINE COMPLETED SUCCESSFULLY!${NC}"
    echo -e "${GREEN}========================================${NC}"
else
    echo ""
    echo -e "${RED}========================================${NC}"
    echo -e "${RED}PIPELINE FAILED - CHECK LOGS ABOVE${NC}"
    echo -e "${RED}========================================${NC}"
    exit 1
fi

# ==============================================================================
# POST-PIPELINE SUMMARY
# ==============================================================================

echo ""
echo -e "${BLUE}RESULTS SUMMARY:${NC}"
echo ""

for tissue in "${TISSUES[@]}"; do
    summary_file="results/$tissue/${tissue}_summary.csv"
    if [ -f "$summary_file" ]; then
        echo -e "${GREEN}✓ $tissue${NC}"
        head -3 "$summary_file" | tail -1
    else
        echo -e "${YELLOW}✗ $tissue (incomplete)${NC}"
    fi
done

echo ""
echo -e "${BLUE}Output directory:${NC} $(pwd)/results/"
echo ""
echo -e "${BLUE}Next steps:${NC}"
echo "  1. Check DEGs: head results/*/DEGs.csv"
echo "  2. View report: firefox results/final_report.html"
echo "  3. Examine logs: less results/logs/*.log"
echo ""
echo -e "${GREEN}Done!${NC}"