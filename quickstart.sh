#!/bin/bash
# Variant-to-Function Reporter — Quick Start
# Run this once to set up your environment.

set -e
echo "============================================"
echo "  V2F Reporter — Environment Setup"
echo "============================================"

# Check Python version
PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
echo "Python: $PYTHON_VERSION"
if [[ "$PYTHON_VERSION" < "3.9" ]]; then
    echo "ERROR: Python 3.9+ required"
    exit 1
fi

# Install dependencies
echo ""
echo "Installing Python dependencies..."
python3 -m pip install -q -r requirements.txt 2>/dev/null || {
    echo "WARNING: Some packages failed to install. Core pipeline may still work."
}

# Check optional tools
echo ""
echo "Checking optional tools:"
command -v bcftools >/dev/null 2>&1 && echo "  ✓ bcftools $(bcftools --version | head -1)" || echo "  ○ bcftools not found (VCF normalization will be skipped)"
command -v docker >/dev/null 2>&1 && echo "  ✓ docker $(docker --version 2>/dev/null | head -1)" || echo "  ○ docker not found (set vep_mode: skip in config)"
command -v gsutil >/dev/null 2>&1 && echo "  ✓ gsutil available" || echo "  ○ gsutil not found (GCS features unavailable)"

# Create output directories
echo ""
echo "Creating output directories..."
mkdir -p data/{vcf,annotated,enriched,classified,phenotype} reports eval logs

# Check config
if [ -f config/pipeline_config.yaml ]; then
    echo ""
    echo "  ✓ Config found: config/pipeline_config.yaml"
else
    echo ""
    echo "  No config found. Run the setup wizard:"
    echo "    python3 setup_wizard.py"
fi

echo ""
echo "============================================"
echo "  Setup complete!"
echo "============================================"
echo ""
echo "  Next steps:"
echo "    1. python3 setup_wizard.py              # Configure pipeline"
echo "    2. python3 run_pipeline.py --dry-run     # Validate setup"
echo "    3. python3 run_pipeline.py               # Run pipeline"
echo "    4. python3 dashboard/app.py              # Start dashboard"
echo ""
