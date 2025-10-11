#!/bin/bash

echo "Language    Runtime"
echo "-------------------"

# Run Python tests
python3 test_phylo_python.py

# Set Python library path for Codon
export CODON_PYTHON=$(which python3)
codon run test_phylo_codon.py