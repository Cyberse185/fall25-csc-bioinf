#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# --- Python timing ---
PY_MS=$(python test_phylo_python.py --just-ms | tail -n1)

# --- Codon timing ---
export CODON_PYTHON=$(python -m find_libpython)
export CODONPATH="$SCRIPT_DIR"
CODON_MS=$(codon run test_phylo_codon.py --just-ms | tail -n1)

# --- Output ---
echo "Language    Runtime"
echo "-------------------"
printf "python      %sms\n" "$PY_MS"
printf "codon       %sms\n" "$CODON_MS"
