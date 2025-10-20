"""
Sequence alignment algorithms for bioinformatics (Codon version).
Implements: global, local, semi-global, and affine gap penalty alignments.
"""
from time import time

def parse_fasta(filename: str) -> str:
    """Parse a FASTA file and return the sequence."""
    seq = ""
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('>'):
                    seq += line
    except:
        print(f"Error: {filename} not found")
        return ""
    return seq


class SequenceAligner:
    match: int
    mismatch: int
    gap: int
    gap_open: int
    gap_ext: int
    
    def __init__(self, match: int = 3, mismatch: int = -3, gap: int = -2, 
                 gap_open: int = -5, gap_ext: int = -1):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.gap_open = gap_open
        self.gap_ext = gap_ext
    
    def _match_score(self, base1: str, base2: str) -> int:
        """Return score for comparing two bases."""
        if base1 == base2:
            return self.match
        else:
            return self.mismatch
    
    def global_alignment(self, seq1: str, seq2: str) -> int:
        """Needleman-Wunsch: Global alignment."""
        n = len(seq1)
        m = len(seq2)
        
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
                dp[i][j] = max(match, max(delete, insert))
        
        return dp[n][m]
    
    def local_alignment(self, seq1: str, seq2: str) -> int:
        """Smith-Waterman: Local alignment."""
        n = len(seq1)
        m = len(seq2)
        
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        
        max_score = 0
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = dp[i-1][j-1] + self._match_score(seq1[i-1], seq2[j-1])
                delete = dp[i-1][j] + self.gap
                insert = dp[i][j-1] + self.gap
                dp[i][j] = max(0, max(match, max(delete, insert)))
                if dp[i][j] > max_score:
                    max_score = dp[i][j]
        
        return max_score
    
    def semi_global_alignment(self, seq1: str, seq2: str) -> int:
        """Fitting/Semi-global alignment."""
        n = len(seq1)
        m = len(seq2)
        
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
                dp[i][j] = max(match, max(delete, insert))
        
        result = dp[n][0]
        for j in range(m + 1):
            if dp[n][j] > result:
                result = dp[n][j]
        return result
    
    def affine_gap_alignment(self, seq1: str, seq2: str) -> int:
        """Global alignment with affine gap penalty."""
        n = len(seq1)
        m = len(seq2)
        
        INF = -999999999
        
        M = [[INF] * (m + 1) for _ in range(n + 1)]
        D = [[INF] * (m + 1) for _ in range(n + 1)]
        I = [[INF] * (m + 1) for _ in range(n + 1)]
        
        M[0][0] = 0
        
        for i in range(1, n + 1):
            D[i][0] = self.gap_open + (i - 1) * self.gap_ext
        
        for j in range(1, m + 1):
            I[0][j] = self.gap_open + (j - 1) * self.gap_ext
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                m_val = max(M[i-1][j-1], max(D[i-1][j-1], I[i-1][j-1])) + self._match_score(seq1[i-1], seq2[j-1])
                d_val = max(M[i-1][j] + self.gap_open, D[i-1][j] + self.gap_ext)
                i_val = max(M[i][j-1] + self.gap_open, I[i][j-1] + self.gap_ext)
                
                M[i][j] = m_val
                D[i][j] = d_val
                I[i][j] = i_val
        
        return max(M[n][m], max(D[n][m], I[n][m]))


def run_tests():
    """Run all alignment tests and output results."""
    aligner = SequenceAligner(match=3, mismatch=-3, gap=-2, gap_open=-5, gap_ext=-1)
    
    test_pairs = [
        ("q1.fa", "t1.fa", "q1"),
        ("q2.fa", "t2.fa", "q2"),
        ("MT-human.fa", "MT-orang.fa", "mt_human"),
    ]
    
    methods = ["global", "local", "semi-global", "affine"]
    
    print(f"{'Method':<20} {'Language':<12} {'Runtime':<10}")
    print("-" * 42)
    
    for seq1_file, seq2_file, pair_name in test_pairs:
        seq1 = parse_fasta(seq1_file)
        if not seq1:
            continue
        
        seq2 = parse_fasta(seq2_file)
        if not seq2:
            continue
        
        print(f"# Testing {pair_name}: seq1_len={len(seq1)}, seq2_len={len(seq2)}")
        
        for method in methods:
            start = time()
            
            if method == "global":
                score = aligner.global_alignment(seq1, seq2)
            elif method == "local":
                score = aligner.local_alignment(seq1, seq2)
            elif method == "semi-global":
                score = aligner.semi_global_alignment(seq1, seq2)
            else:  # affine
                score = aligner.affine_gap_alignment(seq1, seq2)
            
            elapsed = time() - start
            runtime_ms = int(elapsed * 1000)
            
            method_label = f"{method}-{pair_name}"
            print(f"{method_label:<20} {'codon':<12} {runtime_ms}ms")


if __name__ == "__main__":
    run_tests()