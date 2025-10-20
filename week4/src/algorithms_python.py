"""
Sequence alignment algorithms for bioinformatics.
Implements: global, local, semi-global, and affine gap penalty alignments.
"""
import time
import os
import sys

def parse_fasta(filename):
    """Parse a FASTA file and return the sequence."""
    seq = ""
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('>'):
                    seq += line
    except FileNotFoundError:
        print(f"Error: {filename} not found", file=sys.stderr)
        return None
    return seq


class SequenceAligner:
    def __init__(self, match=3, mismatch=-3, gap=-2, gap_open=-5, gap_ext=-1):
        """
        Initialize alignment parameters.
        
        Args:
            match: Score for matching bases
            mismatch: Score for mismatching bases
            gap: Score for gap (used in global/local/semi-global)
            gap_open: Score for opening a gap (affine)
            gap_ext: Score for extending a gap (affine)
        """
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.gap_open = gap_open
        self.gap_ext = gap_ext
    
    def _match_score(self, base1, base2):
        """Return score for comparing two bases."""
        return self.match if base1 == base2 else self.mismatch
    
    def global_alignment(self, seq1, seq2):
        """Needleman-Wunsch: Global alignment."""
        n, m = len(seq1), len(seq2)
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        
        for i in range(n + 1):
            dp[i][0] = i * self.gap
        for j in range(m + 1):
            dp[0][j] = j * self.gap
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = dp[i-1][j-1] + self._match_score(seq1[i-1], seq2[j-1])
                delete = dp[i-1][j] + self.gap
                insert = dp[i][j-1] + self.gap
                dp[i][j] = max(match, delete, insert)
        
        return dp[n][m]
    
    def local_alignment(self, seq1, seq2):
        """Smith-Waterman: Local alignment."""
        n, m = len(seq1), len(seq2)
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        
        max_score = 0
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = dp[i-1][j-1] + self._match_score(seq1[i-1], seq2[j-1])
                delete = dp[i-1][j] + self.gap
                insert = dp[i][j-1] + self.gap
                dp[i][j] = max(0, match, delete, insert)
                max_score = max(max_score, dp[i][j])
        
        return max_score
    
    def semi_global_alignment(self, seq1, seq2):
        """Fitting/Semi-global alignment."""
        n, m = len(seq1), len(seq2)
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        
        for j in range(m + 1):
            dp[0][j] = 0
        
        for i in range(1, n + 1):
            dp[i][0] = i * self.gap
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = dp[i-1][j-1] + self._match_score(seq1[i-1], seq2[j-1])
                delete = dp[i-1][j] + self.gap
                insert = dp[i][j-1] + self.gap
                dp[i][j] = max(match, delete, insert)
        
        return max(dp[n])
    
    def affine_gap_alignment(self, seq1, seq2):
        """Global alignment with affine gap penalty."""
        n, m = len(seq1), len(seq2)
        
        M = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
        D = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
        I = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
        
        M[0][0] = 0
        
        for i in range(1, n + 1):
            D[i][0] = self.gap_open + (i - 1) * self.gap_ext
        
        for j in range(1, m + 1):
            I[0][j] = self.gap_open + (j - 1) * self.gap_ext
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                M[i][j] = max(M[i-1][j-1], D[i-1][j-1], I[i-1][j-1]) + self._match_score(seq1[i-1], seq2[j-1])
                D[i][j] = max(M[i-1][j] + self.gap_open, D[i-1][j] + self.gap_ext)
                I[i][j] = max(M[i][j-1] + self.gap_open, I[i][j-1] + self.gap_ext)
        
        return max(M[n][m], D[n][m], I[n][m])


def run_tests():
    """Run all alignment tests and output results."""
    aligner = SequenceAligner(match=3, mismatch=-3, gap=-2, gap_open=-5, gap_ext=-1)
    
    test_pairs = [
        ("q1.fa", "t1.fa", "q1"),
        # ("q2.fa", "t2.fa", "q2"),  # Don't have this file
        ("MT-human.fa", "MT-orang.fa", "mt_human"),  # Enable MT tests
    ]
    
    methods = ["global", "local", "semi-global", "affine"]
    
    print(f"{'Method':<20} {'Language':<12} {'Runtime':<10}")
    print("-" * 42)
    
    for seq1_file, seq2_file, pair_name in test_pairs:
        if not os.path.exists(seq1_file):
            print(f"Skipping {pair_name}: {seq1_file} not found", file=sys.stderr)
            continue
        if not os.path.exists(seq2_file):
            print(f"Skipping {pair_name}: {seq2_file} not found", file=sys.stderr)
            continue
        
        seq1 = parse_fasta(seq1_file)
        seq2 = parse_fasta(seq2_file)
        
        if seq1 is None or seq2 is None:
            continue
        
        print(f"# Testing {pair_name}: seq1_len={len(seq1)}, seq2_len={len(seq2)}")
        
        for method in methods:
            # Skipped as Python runs out of memory causing error 143
            if method == "affine" and pair_name == "mt_human":
                print(f"affine-{pair_name:<13} {'python':<12} SKIPPED (out of memory)")
                continue
            
            start = time.time()
            
            if method == "global":
                score = aligner.global_alignment(seq1, seq2)
            elif method == "local":
                score = aligner.local_alignment(seq1, seq2)
            elif method == "semi-global":
                score = aligner.semi_global_alignment(seq1, seq2)
            else:  # affine
                score = aligner.affine_gap_alignment(seq1, seq2)
            
            elapsed = time.time() - start
            runtime_ms = int(elapsed * 1000)
            
            method_label = f"{method}-{pair_name}"
            print(f"{method_label:<20} {'python':<12} {runtime_ms}ms")


if __name__ == "__main__":
    run_tests()
