#!/bin/bash

# check_requirements.sh - Check system requirements for dbSNP to Parquet conversion

set -e

echo "🔍 Checking requirements for dbSNP to Parquet conversion..."

# Check Python version
echo "Checking Python version..."
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version 2>&1 | grep -oP '\d+\.\d+')
    if [[ $(echo "$PYTHON_VERSION >= 3.11" | bc -l) -eq 1 ]]; then
        echo "✅ Python $PYTHON_VERSION (>= 3.11 required)"
    else
        echo "❌ Python $PYTHON_VERSION found, but >= 3.11 is required"
        exit 1
    fi
else
    echo "❌ Python3 not found"
    exit 1
fi

# Check bcftools
echo "Checking bcftools..."
if command -v bcftools &> /dev/null; then
    BCFTOOLS_VERSION=$(bcftools --version | head -n1 | grep -oP '\d+\.\d+')
    echo "✅ bcftools $BCFTOOLS_VERSION"
else
    echo "❌ bcftools not found (required for VCF/BCF processing)"
    exit 1
fi

# Check if uv is available (preferred) or pip
echo "Checking package manager..."
if command -v uv &> /dev/null; then
    echo "✅ uv package manager found"
    PACKAGE_MANAGER="uv"
elif command -v pip &> /dev/null; then
    echo "✅ pip found (uv recommended for faster installs)"
    PACKAGE_MANAGER="pip"
else
    echo "❌ Neither uv nor pip found"
    exit 1
fi

# Check if virtual environment exists
echo "Checking virtual environment..."
if [ -d ".venv" ]; then
    echo "✅ Virtual environment found at .venv"
    if [ -f ".venv/pyvenv.cfg" ]; then
        VENV_PYTHON=$(grep "version" .venv/pyvenv.cfg | cut -d' ' -f3)
        echo "  Python version in venv: $VENV_PYTHON"
    fi
else
    echo "❌ Virtual environment not found at .venv"
    echo "  Run: $PACKAGE_MANAGER venv .venv"
    exit 1
fi

# Check Python dependencies
echo "Checking Python dependencies..."
source .venv/bin/activate

MISSING_DEPS=()

# Check each dependency from pyproject.toml
for dep in cyvcf2 snakemake pandas pyarrow polars; do
    if python3 -c "import $dep" 2>/dev/null; then
        VERSION=$(python3 -c "import $dep; print($dep.__version__)" 2>/dev/null || echo "unknown")
        echo "✅ $dep ($VERSION)"
    else
        echo "❌ $dep not found"
        MISSING_DEPS+=("$dep")
    fi
done

# Check snakemake executor plugin separately
if python3 -c "import snakemake_executor_plugin_cluster_generic" 2>/dev/null; then
    echo "✅ snakemake-executor-plugin-cluster-generic"
else
    echo "❌ snakemake-executor-plugin-cluster-generic not found"
    MISSING_DEPS+=("snakemake-executor-plugin-cluster-generic")
fi

if [ ${#MISSING_DEPS[@]} -gt 0 ]; then
    echo ""
    echo "❌ Missing dependencies: ${MISSING_DEPS[*]}"
    echo "  Run: $PACKAGE_MANAGER sync"
    exit 1
fi

# Check if download.sh has been executed
echo "Checking downloaded data files..."
DOWNLOAD_FILES=("GCF_000001405.40.gz" "GCF_000001405.40.gz.tbi")
MISSING_DOWNLOADS=()

for file in "${DOWNLOAD_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "✅ $file"
    else
        echo "❌ $file not found"
        MISSING_DOWNLOADS+=("$file")
    fi
done

if [ ${#MISSING_DOWNLOADS[@]} -gt 0 ]; then
    echo ""
    echo "❌ Downloaded data files missing: ${MISSING_DOWNLOADS[*]}"
    echo "  Please run: ./download.sh"
    echo "  This will download the required dbSNP VCF files from NCBI."
    exit 1
fi

# Check required files
echo "Checking required files..."
REQUIRED_FILES=("Snakefile" "convert_to_parquet.py" "contig_map.tsv")
for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "✅ $file"
    else
        echo "❌ $file not found"
        exit 1
    fi
done

# Check disk space (warn if less than 50GB available)
echo "Checking disk space..."
AVAILABLE_GB=$(df . | tail -1 | awk '{print int($4/1024/1024)}')
if [ $AVAILABLE_GB -lt 50 ]; then
    echo "⚠️  Warning: Only ${AVAILABLE_GB}GB available (dbSNP processing may require significant disk space)"
else
    echo "✅ Sufficient disk space (${AVAILABLE_GB}GB available)"
fi

echo ""
echo "✅ All requirements check passed!"
echo "Ready to run: snakemake --cores all"