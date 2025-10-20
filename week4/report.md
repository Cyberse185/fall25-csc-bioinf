# Week 4 Report: Sequence Alignment Algorithms

## Implementation Steps

1. Implemented four sequence alignment algorithms in Python:
   - Global alignment (Needleman-Wunsch)
   - Local alignment (Smith-Waterman)
   - Semi-global (fitting) alignment
   - Affine gap penalty global alignment

2. Created FASTA file parser for loading test sequences

3. Tested Python implementation on small sequences (q1.fa vs t1.fa)

4. Attempted testing on large mitochondrial DNA sequences (MT-human.fa vs MT-orang.fa)
   - Discovered Python runs out of memory on affine algorithm with large sequences

5. Ported all algorithms to Codon with proper type annotations

6. Verified Codon handles large sequences efficiently (5-6x faster than Python, no memory issues)

7. Fixed spec clarification: semi-global uses gap=-2, not gap_open=-5 (affine uses gap_open)

8. Created evaluate.sh script to run both Python and Codon tests

9. Documented implementation in report.md and ai.md

## Gotchas and Challenges

### Memory Management with Large Sequences
- MT sequences are ~16.5k bases each, creating 273 million cell DP tables
- Python integer overhead (~28 bytes each) = ~7.6 GB per table
- Affine algorithm needs 3 tables (M, D, I) = ~23 GB total
- Python ran out of system RAM attempting affine on MT sequences


### Gap Penalty Specification Clarification
- Student asked clarification, professor confirmed: affine uses gap_open=-5/gap_ext=-1
- Semi-global, global, and local all use simple gap=-2
- Had to update Python and Codon implementations

### Type Annotations for Codon
- Codon requires explicit type annotations: `def func(seq1: str, seq2: str) -> int:`
- Cannot use `float('-inf')` in Codon, must use explicit large negative constant like `-999999999`
- Used `from time import time` instead of Python's `time.time()`

### File Naming Conventions
- Test files use hyphens not underscores: `MT-human.fa` not `MT_human.fa`
- Had to update file paths in test_pairs list to match actual filenames

### Traceback vs Score-Only
- Initially implemented full alignment traceback (generating aligned strings with gaps)
- Removed traceback logic to improve performance since deliverable only needs scores and timing
- Final version returns only alignment score, not the aligned sequences

### Test Data Availability
- ksw2 repo only has q1/t1 and q2/t2 (compressed as .gz)
- Used only q1/t1 and MT sequences for testing
- No q3-q5 test pairs available in upstream repo

## Algorithm Details

### Global Alignment
- Needleman-Wunsch algorithm
- Aligns entire sequences end-to-end
- DP table: dp[i][j] = best score for first i chars of seq1, first j chars of seq2
- O(n*m) time and space

### Local Alignment
- Smith-Waterman algorithm
- Finds best local region match between sequences
- Allows cells to reset to 0 (start new alignment)
- Tracks maximum score anywhere in table
- O(n*m) time and space

### Semi-global Alignment
- Fitting alignment, aligns one sequence to substring of another
- First row initialized to 0 (query can start anywhere)
- Reports max score in last row (query fully aligned)
- O(n*m) time and space

### Affine Gap Penalty
- Gotoh algorithm for realistic gap costs
- Maintains 3 DP tables: M (match), D (deletion), I (insertion)
- Gap open = -5, gap extension = -1
- O(n*m) time, O(3nm) space

## Test Results

### Python Performance
Small sequences (q1 vs t1: 105bp × 143bp):
```
global-q1:       3ms
local-q1:        4ms
semi-global-q1:  3ms
affine-q1:       6ms
```

Large sequences (MT-human vs MT-orang: 16,569bp × 16,499bp):
```
global-mt_human:       80,542ms (80.5s)
local-mt_human:        92,207ms (92.2s)
semi-global-mt_human:  77,260ms (77.3s)
affine-mt_human:       FAILED - Out of memory
```

### Codon Performance
Small sequences:
```
global-q1:       0ms
local-q1:        0ms
semi-global-q1:  0ms
affine-q1:       1ms
```

Large sequences:
```
global-mt_human:       14,555ms (14.6s) - 5.5x faster
local-mt_human:        17,481ms (17.5s) - 5.3x faster
semi-global-mt_human:  16,225ms (16.2s) - 4.8x faster
affine-mt_human:       31,896ms (31.9s) - Works (Python crashed)
```

## Codon vs Python Comparison

**Performance improvements with Codon:**
- 5-6x speedup on large sequences
- Handles memory-intensive computations Python cannot
- Identical algorithm logic, just compiled instead of interpreted

**Code differences:**
- Codon needs explicit type annotations
- Codon uses integer constants instead of float('-inf')
- Codon time import differs from Python's time module
- Otherwise nearly identical code

## Time Spent

Rougly 7 hours. 