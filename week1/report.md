# Week 1 Report

## Task A: Python Implementation and Reproduction

### Setup Process
- Cloned Zhongyu Chen's de Bruijn graph genome assembler repository
- Required resolving Chinese documentation using translation tools
- Set up Python environment with required dependencies (matplotlib)

### Reproduction Attempt
- Successfully ran Python version on data1-4
- Original README showed NGA50 values: 9118.8, 9129.2, 7859.2, 55757.8
- Assignment update clarified to use N50 instead of NGA50
- Direct comparison not possible as original used NGA50 (requires reference genomes)
- My Python results produced different N50 values, which is expected per assignment notes

### Runtime Observations
- data1-3: ~10-15 seconds each
- data4: ~13 minutes (required stack limit adjustment)

## Task B: Codon Conversion

### Conversion Process
- Converted main.py, dbg.py, utils.py to corresponding .codon files
- Main challenges encountered:
  - Python interop issues: Removed `from python import` statements
  - Type annotations: Added explicit typing where needed
  - File path handling: Replaced `os.path.join` with string concatenation
  - Syntax errors: Fixed line breaks in assignment statements

### Key Solutions
- Environment setup: Set `CODON_PYTHON` variable for Python bridge
- Stack limits: Found `ulimit -s 65520` works on macOS for data4
- Algorithm behavior: Codon and Python produce slightly different results due to internal data 
structure ordering differences

### Performance Results
- Codon generally 2x faster than Python on smaller datasets
- data4 showed unexpected behavior where Codon was slightly slower

## Task C: Automation and CI Integration

### Automation Script
- Created `evaluate.sh` that runs both versions and compares performance
- Calculates N50 from generated contig files
- Outputs formatted comparison table with runtime and N50 metrics
- Handles data4 stack limit requirements

### CI Integration
- Set up GitHub Actions with Codon installation and Python bridge
- Automated testing runs successfully on Ubuntu environment
- Script executes from correct directory structure

### Results Summary
All datasets run successfully with automated performance comparison demonstrating the 
conversion works correctly despite minor numerical differences in assembly results.
