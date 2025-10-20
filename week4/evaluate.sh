#!/bin/bash
# Evaluate script for Week 4 - runs all alignment tests and outputs timing table

cd "$(dirname "$0")"

echo "Method            Language    Runtime"
echo "--------------------------------------"

# Run Python tests
python src/algorithms_python.py 2>&1 | grep -E "^(global|local|semi|affine)-"

# Run Codon tests  
codon run src/algorithms_codon.py 2>&1 | grep -E "^(global|local|semi|affine)-"
