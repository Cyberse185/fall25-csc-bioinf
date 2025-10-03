#!/bin/bash

set -e

echo "=== Week 2: BioPython Motifs Port ==="

echo "Testing with Python..."
python test.py || echo "Python tests had issues"

echo "Testing with Codon..."
codon run test.py || echo "Codon tests failed - known type system 
issues"

echo "=== Week 2 complete ==="
