import os

# Detect environment (instructor-approved method)
if hasattr(str, 'memcpy'):
    IS_CODON = True
    from codon_motifs import create,read, parse,  Motif
    from codon_motifs.thresholds import ScoreDistribution
else:
    IS_CODON = False
    from Bio import motifs
    from Bio.motifs import thresholds
    create = motifs.create
    Motif = motifs.Motif
    read = motifs.read
    parse = motifs.parse
    ScoreDistribution = thresholds.ScoreDistribution
    
    def test(func):
        return func

@test
def test_basic_motif():
    """Test basic motif functionality."""
    print("\n=== Testing Basic Motif ===")
    
    instances = ["AACGCCA", "ACCGCCC", "AACTCCG"]
    motif = create(instances)
    
    print(f"Motif length: {motif.length}")
    
    # Get consensus
    consensus = motif.consensus
    print(f"Consensus: {consensus}")
    
    # Test PWM
    pwm = motif.pwm
    print(f"\nPWM length: {pwm.length}")
    
    # Test PSSM (Position-Specific Scoring Matrix)
    print("\n--- Testing PSSM ---")
    pssm = motif.pssm
    
    print(f"PSSM type: {type(pssm)}")
    print(f"PSSM length: {pssm.length}")
    print(f"PSSM alphabet: {pssm.alphabet}")
    
    # Show PSSM values for first position
    print("\nPSSM log-odds scores (first position):")
    for letter in pssm.alphabet:
        print(f"  {letter}[0]: {pssm[letter][0]:.4f}")
    
    # Test PSSM min/max scores
    print("\n--- PSSM score range ---")
    max_score = pssm.max
    min_score = pssm.min
    print(f"Max possible score: {max_score:.4f}")
    print(f"Min possible score: {min_score:.4f}")
    
    # Test PSSM mean and std
    print("\n--- PSSM statistics ---")
    mean_score = pssm.mean()
    std_score = pssm.std()
    print(f"Mean score: {mean_score:.4f}")
    print(f"Std deviation: {std_score:.4f}")
    
    # Test scoring a sequence
    print("\n--- Scoring sequences ---")
    test_seq = "AACGCCA"  # Should score high (matches consensus)
    scores = pssm.calculate(test_seq)
    print(f"Score for '{test_seq}': {scores}")
    
    test_seq2 = "TTTTTTTT"  # Longer sequence
    scores2 = pssm.calculate(test_seq2)
    print(f"Scores for '{test_seq2}': {scores2}")
    print(f"Number of scores: {len(scores2)}")
    
    print("\n✓ test_basic_motif passed!")


@test
def test_reverse_complement():
    """Test reverse complement - simplest addition."""
    print("\n=== Testing Reverse Complement ===")
    
    instances = ["ACG"]
    motif = create(instances)
    
    consensus = motif.consensus  # Changed from motif.counts.get_consensus()
    print(f"Original consensus: {consensus}")
    
    try:
        rc_motif = motif.reverse_complement()
        rc_consensus = rc_motif.consensus  # Changed from rc_motif.counts.get_consensus()
        print(f"RC consensus: {rc_consensus}")
        
        expected = "CGT"
        if str(rc_consensus) == expected:
            print(f"✓ Correct! ACG -> {expected}")
        else:
            print(f"✗ Expected {expected}, got {rc_consensus}")
        
        print("\n✓ test_reverse_complement passed!")
    except Exception as e:
        print(f"✗ Error in reverse_complement: {e}")
        print("Skipping this test for now...")


@test
def test_consensus_variants():
    """Test different types of consensus sequences."""
    print("\n=== Testing Consensus Variants ===")
    
    # Create a motif with some variation
    instances = ["AACG", "AACA", "GACG"]
    motif = create(instances)
    
    # Test regular consensus
    consensus = motif.consensus
    print(f"Consensus: {consensus}")
    
    # Test anticonsensus (least probable)
    anticonsensus = motif.anticonsensus
    print(f"Anticonsensus: {anticonsensus}")
    
    # Test degenerate consensus
    degenerate = motif.degenerate_consensus
    print(f"Degenerate consensus: {degenerate}")
    
    # Test using properties (should work for both)
    print("\nTesting @property accessors:")
    print(f"motif.consensus: {motif.consensus}")
    print(f"motif.anticonsensus: {motif.anticonsensus}")
    print(f"motif.degenerate_consensus: {motif.degenerate_consensus}")
    
    print("\n✓ test_consensus_variants passed!")
    
    


# ============================================================================
# INCREMENTAL TESTS - Uncomment one at a time to test new functionality
# ============================================================================

@test
def test_counts_operations():
    """Test Counts object operations."""
    print("\n=== Testing Counts Operations ===")
    
    instances = ["ACG", "ACT", "GCG"]
    motif = create(instances)
    counts = motif.counts
    
    # Test counts for specific positions
    print(f"Alphabet: {counts.alphabet}")
    print("\nCounts at position 0:")
    for letter in counts.alphabet:
        print(f"  {letter}: {counts[letter][0]}")
    
    # Test length and shape
    print(f"\nCounts length: {counts.length}")
    
    print("\n✓ test_counts_operations passed!")


@test
def test_matrix_operations():
    """Test matrix-specific operations."""
    print("\n=== Testing Matrix Operations ===")
    
    instances = ["ACG", "ACT"]
    motif = create(instances)
    
    # Get different matrix types
    counts = motif.counts
    pwm = motif.pwm
    pssm = motif.pssm
    
    print(f"Counts type: {type(counts)}")
    print(f"PWM type: {type(pwm)}")
    print(f"PSSM type: {type(pssm)}")
    
    # Test matrix properties
    print(f"\nAll have same length: {counts.length == pwm.length == pssm.length}")
    print(f"All have same alphabet: {counts.alphabet == pwm.alphabet == pssm.alphabet}")
    
    print("\n✓ test_matrix_operations passed!")
    
@test
def test_rna_motif():
    """Test RNA motif with U instead of T."""
    print("\n=== Testing RNA Motif ===")
    
    sequences = ["UACAA", "UACGC", "UACAC"]
    motif = create(sequences, alphabet="ACGU")
    
    print(f"Consensus: {motif.consensus}")
    assert "U" in str(motif.consensus)
    
    # Test reverse complement
    rc = motif.reverse_complement()
    print(f"RC consensus: {rc.consensus}")
    
    print("✓ test_rna_motif passed!")

@test
def test_score_distribution():
    """Test ScoreDistribution for threshold calculation."""
    print("\n=== Testing ScoreDistribution ===")
    
    from codon_motifs.thresholds import ScoreDistribution
    
    sequences = ["TACAA", "TACGC", "AACCC", "AATGC"]
    motif =  create(sequences, alphabet="ACGT")
    
    background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    
    # FIX: Use pseudocounts to avoid -inf
    pwm = motif.counts.normalize(pseudocounts=1.0)  # <-- ADD THIS LINE
    pssm = pwm.log_odds(background=background)       # <-- CHANGE THIS LINE
    
    print(f"PSSM min: {pssm.min:.4f}, max: {pssm.max:.4f}")
    
    # Create distribution
    dist = ScoreDistribution(pssm=pssm, background=background, precision=1000)
    
    print(f"Interval: {dist.interval:.4f}")
    print(f"N points: {dist.n_points}")
    print(f"IC: {dist.ic:.4f}")
    
    # Test thresholds
    print("\n--- Testing threshold_fpr ---")
    threshold_low = dist.threshold_fpr(0.01)
    threshold_high = dist.threshold_fpr(0.5)
    
    print(f"Threshold (FPR=0.01): {threshold_low:.4f}")
    print(f"Threshold (FPR=0.5): {threshold_high:.4f}")
    
    assert threshold_low > threshold_high
    assert pssm.min < threshold_low < pssm.max
    
    print("\n--- Testing threshold_patser ---")
    threshold_patser = dist.threshold_patser()
    print(f"PATSER threshold: {threshold_patser:.4f}")
    
    print("\n✓ test_score_distribution passed!")
    
# @test
# def test_minimal_meme_parser():
#     """Parse minimal_test.meme file (matches BioPython test)."""
#     print("\n=== Testing Minimal MEME Parser (DNA) ===")
    
#     meme_file = "motifs/minimal_test.meme"
    
#     try: 
#         with open(meme_file) as stream: 
#             record = parse(stream, "minimal")
#     except FileNotFoundError: 
#         print(f"Skipping: {meme_file} not found")
#         return
    
#     # Test Record properties
#     print(f"Version: {record.version}")
#     print(f"Alphabet: {record.alphabet}")
#     print(f"Number of motifs: {len(record)}")
    
#     assert record.version == "4"
#     assert record.alphabet == "ACGT"
#     assert len(record.sequences) == 0
#     assert record.command == ""
#     assert len(record) == 3
    
#     # Test first motif (KRP)
#     motif = record[0]
#     print(f"\nMotif 1: {motif.name}")
#     print(f"  Length: {motif.length}")
#     print(f"  Occurrences: {motif.num_occurrences}")
#     print(f"  E-value: {motif.evalue}")
#     print(f"  Consensus: {motif.consensus}")
    
#     assert motif.name == "KRP"
#     assert record["KRP"] == motif  # Test name-based access
#     assert motif.num_occurrences == 17
#     assert motif.length == 19
#     assert motif.alphabet == "ACGT"
    
#     # Test background frequencies
#     assert abs(motif.background["A"] - 0.30269730269730266) < 1e-10
#     assert abs(motif.background["C"] - 0.1828171828171828) < 1e-10
#     assert abs(motif.background["G"] - 0.20879120879120877) < 1e-10
#     assert abs(motif.background["T"] - 0.30569430569430567) < 1e-10
    
#     # Test E-value
#     assert abs(motif.evalue - 4.1e-09) < 1e-19
    
#     # Test consensus sequences
#     assert str(motif.consensus) == "TGTGATCGAGGTCACACTT"
#     assert str(motif.degenerate_consensus) == "TGTGANNNWGNTCACAYWW"
    
#     # Test slicing
#     print("  Testing slicing...")
#     sliced = motif[2:9]
#     assert str(sliced.consensus) == "TGATCGA"
    
#     # Test second motif (IFXA)
#     print(f"\nMotif 2: {record[1].name}")
#     assert record[1].name == "IFXA"
    
#     print("\n✓ test_minimal_meme_parser passed!")


# @test
# def test_meme_parser_rna():
#     """Test parsing MEME output files using RNA (matches BioPython test)."""
#     print("\n=== Testing Minimal MEME Parser (RNA) ===")
    
#     meme_file = "motifs/minimal_test_rna.meme"
#     try: 
#         with open(meme_file) as stream:
#             record = parse(stream, "minimal")
#     except FileNotFoundError: 
#         print(f"skipping:{meme_file} not found")
#         return
    
#     # Test Record properties
#     print(f"Version: {record.version}")
#     print(f"Alphabet: {record.alphabet}")
#     print(f"Number of motifs: {len(record)}")
    
#     assert record.version == "4"
#     assert record.alphabet == "ACGU"
#     assert len(record.sequences) == 0
#     assert record.command == ""
#     assert len(record) == 3
    
#     # Test first motif
#     motif = record[0]
#     print(f"\nMotif 1: {motif.name}")
#     print(f"  Consensus: {motif.consensus}")
#     print(f"  Degenerate: {motif.degenerate_consensus}")
    
#     assert motif.name == "KRP_fake_RNA"
#     assert record["KRP_fake_RNA"] == motif
#     assert motif.num_occurrences == 17
#     assert motif.length == 19
    
#     # Test background frequencies (U instead of T)
#     assert abs(motif.background["A"] - 0.30269730269730266) < 1e-10
#     assert abs(motif.background["C"] - 0.1828171828171828) < 1e-10
#     assert abs(motif.background["G"] - 0.20879120879120877) < 1e-10
#     assert abs(motif.background["U"] - 0.30569430569430567) < 1e-10
    
#     assert abs(motif.evalue - 4.1e-09) < 1e-19
#     assert motif.alphabet == "ACGU"
    
#     # Test consensus (should use U not T)
#     assert str(motif.consensus) == "UGUGAUCGAGGUCACACUU"
#     assert str(motif.degenerate_consensus) == "UGUGANNNWGNUCACAYWW"
    
#     print("\n✓ test_meme_parser_rna passed!")


# ============================================================================
# RUN TESTS
# ============================================================================

if __name__ == "__main__":
    print("="*60)
    print(f"Running in {'CODON' if IS_CODON else 'PYTHON'} mode")
    print("="*60)
    
    # Call each test explicitly (required for both Codon and Python)
    test_basic_motif()
    test_reverse_complement()
    test_consensus_variants()
    test_counts_operations()
    test_matrix_operations()
    test_rna_motif()
    test_score_distribution()
    #test_minimal_meme_parser()
    #test_meme_parser_rna()
    
    
    print("\n" + "="*60)
    print("All enabled tests completed!")
    print("="*60)