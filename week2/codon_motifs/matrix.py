"""Support for various forms of sequence motif matrices."""

import math
from typing import Dict, List, Tuple, Optional
from .seq import Seq

class GenericPositionMatrix:
    """Base class for the support of position matrix operations."""
    
    _data: Dict[str, List[float]]
    length: int
    alphabet: List[str]
    
    
    def __init__(self, alphabet: List[str], values: Dict[str, List[float]]):
        """Initialize the class."""
        self._data: Dict[str, List[float]] = {}
        self.length: int = 0
        self.alphabet: List[str] = alphabet
        
        first = True
        for letter in alphabet:
            if first:
                self.length = len(values[letter])
                first = False
            elif self.length != len(values[letter]):
                raise ValueError("data has inconsistent lengths")
            # Cast any values to Python floats
            self._data[letter] = [float(v) for v in values[letter]]
        self.alphabet = alphabet
    
    def __getitem__(self, key: str) -> List[float]:
        """Get item like a dictionary."""
        return self._data[key]
    
    def __setitem__(self, key: str, value: List[float]):
        """Set item like a dictionary."""
        self._data[key] = value
    
    def __contains__(self, key: str) -> bool:
        """Check if key exists like a dictionary."""
        return key in self._data
    
    def __str__(self) -> str:
        """Return a string containing nucleotides and counts of the alphabet in the Matrix."""
        words = [f"{i:6d}" for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        for letter in self.alphabet:
            words = [f"{value:6.2f}" for value in self[letter]]
            line = f"{letter}: " + " ".join(words)
            lines.append(line)
        text = "\n".join(lines) + "\n"
        return text
    
    # ... keep all your other methods (get_consensus, get_anticonsensus, etc.) exactly as they are
    
    def get_consensus(self) -> Seq:
        """Return the consensus sequence."""
        sequence = ""
        for i in range(self.length):
            maximum = -math.inf
            sequence_letter = ""
            for letter in self.alphabet:
                count = self[letter][i]
                if count > maximum:
                    maximum = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq(sequence)
    
    @property
    def consensus(self) -> Seq:
        """Return the consensus sequence as a property."""
        return self.get_consensus()

    def get_anticonsensus(self) -> Seq:
        """Return the anticonsensus sequence."""
        sequence = ""
        for i in range(self.length):
            minimum = math.inf
            sequence_letter = ""
            for letter in self.alphabet:
                count = self[letter][i]
                if count < minimum:
                    minimum = count
                    sequence_letter = letter
            sequence += sequence_letter
        return Seq(sequence)
    @property
    def anticonsensus(self) -> Seq:
        """Return the anticonsensus sequence as a property."""
        return self.get_anticonsensus()

    
    def get_degenerate_consensus(self) -> Seq:
        """Return the degenerate consensus sequence."""
        degenerate_nucleotide = {
            "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
            "AC": "M", "AG": "R", "AT": "W", "AU": "W",
            "CG": "S", "CT": "Y", "CU": "Y", "GT": "K", "GU": "K",
            "ACG": "V", "ACT": "H", "ACU": "H", "AGT": "D", "AGU": "D",
            "CGT": "B", "CGU": "B", "ACGT": "N", "ACGU": "N",
        }
        sequence = ""
        for i in range(self.length):
            # Sort nucleotides by count at position i
            nucleotides = sorted(self.alphabet, 
                               key=lambda x: self[x][i], 
                               reverse=True)
            counts = [self[c][i] for c in nucleotides]
            
            # Follow the Cavener rules:
            if counts[0] > sum(counts[1:]) and counts[0] > 2 * counts[1]:
                key = nucleotides[0]
            elif 4 * sum(counts[:2]) > 3 * sum(counts):
                key = "".join(sorted(nucleotides[:2]))
            elif counts[3] == 0:
                key = "".join(sorted(nucleotides[:3]))
            else:
                key = "ACGT"
            
            nucleotide = degenerate_nucleotide.get(key, key)
            sequence += nucleotide
        return Seq(sequence)
    
    def calculate_consensus(self, 
                          substitution_matrix: Optional[Dict] = None,
                          plurality: Optional[float] = None,
                          identity: float = 0.0,
                          setcase: Optional[float] = None) -> str:
        """Return the consensus sequence (as a string) for the given parameters."""
        alphabet = self.alphabet
        
        # Determine undefined character
        alphabet_set = set(alphabet)
        if alphabet_set.union(set("ACGTUN-")) == set("ACGTUN-"):
            undefined = "N"
        else:
            undefined = "X"
        
        if substitution_matrix is None:
            if plurality is not None:
                raise ValueError(
                    "plurality must be None if substitution_matrix is None"
                )
            sequence = ""
            for i in range(self.length):
                maximum = 0.0
                total = 0.0
                consensus_letter = ""
                for letter in alphabet:
                    count = self[letter][i]
                    total += count
                    if count > maximum:
                        maximum = count
                        consensus_letter = letter
                
                if maximum < identity * total:
                    consensus_letter = undefined
                else:
                    if setcase is None:
                        setcase_threshold = total / 2.0
                    else:
                        setcase_threshold = setcase * total
                    if maximum <= setcase_threshold:
                        consensus_letter = consensus_letter.lower()
                sequence += consensus_letter
        else:
            raise NotImplementedError(
                "calculate_consensus currently only supports substitution_matrix=None"
            )
        return sequence
    
    def get_gc_content(self) -> float:
        """Compute the fraction GC content."""
        alphabet = self.alphabet
        gc_total = 0.0
        total = 0.0
        for i in range(self.length):
            for letter in alphabet:
                if letter in "CG":
                    gc_total += self[letter][i]
                total += self[letter][i]
        return gc_total / total
    
    @property
    def gc_content(self) -> float:
        """Return the GC content as a property."""
        return self.get_gc_content()
    
    def reverse_complement(self):
        """Compute reverse complement."""
        values: Dict[str, List[float]] = {}
        if self.alphabet == ["A", "C", "G", "U"]:
            values["A"] = self["U"][::-1]
            values["U"] = self["A"][::-1]
        else:
            values["A"] = self["T"][::-1]
            values["T"] = self["A"][::-1]
        values["G"] = self["C"][::-1]
        values["C"] = self["G"][::-1]
        alphabet = self.alphabet
        return self.__class__(alphabet, values)


class PositionWeightMatrix(GenericPositionMatrix):
    """Class for the support of weight calculations on the Position Matrix."""
    
    def __init__(self, alphabet: List[str], counts: Dict[str, List[float]]):
        """Initialize the class."""
        GenericPositionMatrix.__init__(self, alphabet, counts)
        
        # Normalize
        for i in range(self.length):
            total = sum(self[letter][i] for letter in alphabet)
            for letter in alphabet:
                self[letter][i] /= total
    
    def log_odds(self, background: Optional[Dict[str, float]] = None):  # Also remove -> "PositionSpecificScoringMatrix" here
        """Return the Position-Specific Scoring Matrix."""
        values: Dict[str, List[float]] = {}
        alphabet = self.alphabet
        
        if background is None:
            bg: Dict[str, float] = {letter: 1.0 for letter in self.alphabet}
        else:
            bg = dict(background)
        
        total = sum(bg.values())
        for letter in alphabet:
            bg[letter] /= total
            values[letter] = []
        
        for i in range(self.length):
            for letter in alphabet:
                b = bg[letter]
                if b > 0:
                    p = self[letter][i]
                    if p > 0:
                        logodds = math.log(p / b) / math.log(2.0)
                    else:
                        logodds = -math.inf
                else:
                    p = self[letter][i]
                    if p > 0:
                        logodds = math.inf
                    else:
                        logodds = math.nan
                values[letter].append(logodds)
        
        pssm = PositionSpecificScoringMatrix(alphabet, values)
        return pssm


class FrequencyPositionMatrix(GenericPositionMatrix):
    """Class for the support of frequency calculations on the Position Matrix."""
    
    def __init__(self, alphabet: List[str], values: Dict[str, List[float]]):
        """Initialize the class."""
        GenericPositionMatrix.__init__(self, alphabet, values)
        
    def normalize(self, pseudocounts = None):  # Remove -> "PositionWeightMatrix"
        """Create and return a position-weight matrix by normalizing the counts matrix."""
        counts: Dict[str, List[float]] = {}
        
        if pseudocounts is None:
            for letter in self.alphabet:
                counts[letter] = [0.0] * self.length
        elif isinstance(pseudocounts, dict):
            for letter in self.alphabet:
                counts[letter] = [float(pseudocounts[letter])] * self.length
        else:
            for letter in self.alphabet:
                counts[letter] = [float(pseudocounts)] * self.length
        
        for i in range(self.length):
            for letter in self.alphabet:
                counts[letter][i] += self[letter][i]
        
        return PositionWeightMatrix(self.alphabet, counts)

class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """Class for the support of Position Specific Scoring Matrix calculations."""
    
    def __init__(self, alphabet: List[str], values: Dict[str, List[float]]):
        """Initialize the class."""
        GenericPositionMatrix.__init__(self, alphabet, values)
    
    def calculate(self, sequence: str):
        """Return the PWM score for a given sequence for all positions."""
        if sorted(self.alphabet) != ["A", "C", "G", "T"]:
            raise ValueError(
                f"PSSM has wrong alphabet: {self.alphabet} - Use only with DNA motifs"
            )
        
        sequence = sequence.upper()
        n = len(sequence)
        m = self.length
        
        if n < m:
            return []
        
        scores: List[float] = []
        for pos in range(n - m + 1):
            score = 0.0
            for i in range(m):
                letter = sequence[pos + i]
                if letter in self:
                    score += self[letter][i]
                else:
                    score = math.nan
                    break
            scores.append(score)
        
        return scores
    
    def search(self, sequence: str, threshold: float = 0.0, 
               both: bool = True):
        """Find hits with PWM score above given threshold."""
        sequence = sequence.upper()
        seq_len = len(sequence)
        motif_l = self.length
        
        # Forward strand
        scores = self.calculate(sequence)
        
        for pos, score in enumerate(scores):
            if not math.isnan(score) and score >= threshold:
                yield (pos, score)
        
        # Reverse complement strand
        if both:
            rc = self.reverse_complement()
            rc_scores = rc.calculate(sequence)
            
            for pos, score in enumerate(rc_scores):
                if not math.isnan(score) and score >= threshold:
                    yield (pos - seq_len, score)
    
    def get_max(self) -> float:
        """Maximal possible score for this motif."""
        score = 0.0
        letters = self.alphabet
        for position in range(self.length):
            score += max(self[letter][position] for letter in letters)
        return score
    
    @property
    def max(self) -> float:
        """Return the maximal possible score as a property."""
        return self.get_max()
    
    def get_min(self) -> float:
        """Minimal possible score for this motif."""
        score = 0.0
        letters = self.alphabet
        for position in range(self.length):
            score += min(self[letter][position] for letter in letters)
        return score
    
    @property
    def min(self) -> float:
        """Return the minimal possible score as a property."""
        return self.get_min()
    
    def mean(self, background: Optional[Dict[str, float]] = None) -> float:
        """Return expected value of the score of a motif."""
        if background is None:
            bg: Dict[str, float] = {letter: 1.0 for letter in self.alphabet}
        else:
            bg = dict(background)
        
        total = sum(bg.values())
        for letter in self.alphabet:
            bg[letter] /= total
        
        sx = 0.0
        for i in range(self.length):
            for letter in self.alphabet:
                logodds = self[letter][i]
                if math.isnan(logodds):
                    continue
                if math.isinf(logodds) and logodds < 0:
                    continue
                b = bg[letter]
                p = b * math.pow(2.0, logodds)
                sx += p * logodds
        return sx
    def std(self, background: Optional[Dict[str, float]] = None) -> float:
        """Return standard deviation of the score of a motif."""
        if background is None:
            bg: Dict[str, float] = {letter: 1.0 for letter in self.alphabet}
        else:
            bg = dict(background)
    
        total = sum(bg.values())
        for letter in self.alphabet:
            bg[letter] /= total
    
        variance = 0.0
        for i in range(self.length):
            sx = 0.0
            sxx = 0.0
            for letter in self.alphabet:
                logodds = self[letter][i]
                if math.isnan(logodds):
                    continue
                if math.isinf(logodds) and logodds < 0:
                    continue
                b = bg[letter]
                p = b * math.pow(2.0, logodds)
                sx += p * logodds
                sxx += p * logodds * logodds
            sxx -= sx * sx
            variance += sxx
        variance = max(variance, 0.0)  # to avoid roundoff problems
        return math.sqrt(variance)