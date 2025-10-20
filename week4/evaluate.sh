cat > evaluate.sh << 'EOF'
#!/bin/bash
# Evaluate script for Week 4 - runs all alignment tests and outputs timing table

cd "$(dirname "$0")"

echo "Method            Language    Runtime"
echo "--------------------------------------"

# Run Python tests
python src/algorithms_python.py 2>/dev/null | grep -E "^(global|local|semi|affine)-" | awk '{print $1, $2, $3}'

# Run Codon tests
codon run src/algorithms_codon.py 2>/dev/null | grep -E "^(global|local|semi|affine)-" | awk '{print $1, $2, $3}'
EOF

chmod +x evaluate.sh
./evaluate.sh