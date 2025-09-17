#!/bin/bash
set -euo pipefail

# Function to calculate N50
calculate_n50() {
    local contig_file="$1"
    if [[ ! -f "$contig_file" ]]; then
        echo "0"
        return
    fi
    
    python3 -c "
import sys
lengths = []
with open('$contig_file', 'r') as f:
    seq = ''
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if seq:
                lengths.append(len(seq))
                seq = ''
        else:
            seq += line
    if seq:
        lengths.append(len(seq))

if not lengths:
    print(0)
else:
    lengths.sort(reverse=True)
    total = sum(lengths)
    target = total / 2
    cumsum = 0
    for length in lengths:
        cumsum += length
        if cumsum >= target:
            print(length)
            break
    else:
        print(0)
"
}

# Print header
printf "%-10s\t%-10s\t%-10s\t%-10s\n" "Dataset" "Language" "Runtime" "N50"
printf "%.70s\n" "$(printf '%*s' 70 | tr ' ' '-')"

cd code

# Set stack size for data4
ulimit -s 65520

# Test all datasets including data4
for dataset in data1 data2 data3 data4; do
    dataset_path="../data/$dataset"
    
    if [[ ! -d "$dataset_path" ]]; then
        continue
    fi
    
    # Clean previous results
    rm -f "$dataset_path/contig.fasta"
    
    # Run Python version
    start_time=$(date +%s)
    python main.py "$dataset_path"
    end_time=$(date +%s)
    python_runtime=$((end_time - start_time))
    python_minutes=$((python_runtime / 60))
    python_seconds=$((python_runtime % 60))
    python_time=$(printf "%d:%02d" $python_minutes $python_seconds)
    
    python_n50=$(calculate_n50 "$dataset_path/contig.fasta")
    
    # Clean for Codon run
    rm -f "$dataset_path/contig.fasta"
    
    # Run Codon version
    start_time=$(date +%s)
    ~/.codon/bin/codon run -plugin seq -release main.codon "$dataset_path" > /dev/null 2>&1
    end_time=$(date +%s)
    codon_runtime=$((end_time - start_time))
    codon_minutes=$((codon_runtime / 60))
    codon_seconds=$((codon_runtime % 60))
    codon_time=$(printf "%d:%02d" $codon_minutes $codon_seconds)
    
    codon_n50=$(calculate_n50 "$dataset_path/contig.fasta")
    
    # Output results
    printf "%-10s\t%-10s\t%-10s\t%-10s\n" "$dataset" "python" "$python_time" "$python_n50"
    printf "%-10s\t%-10s\t%-10s\t%-10s\n" "$dataset" "codon" "$codon_time" "$codon_n50"
done
