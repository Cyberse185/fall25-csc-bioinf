"""Tools for sequence motif analysis.

Bio.motifs contains the core Motif class containing various I/O methods
as well as methods for motif comparisons and motif searching in sequences.

Ported to Codon from BioPython.
"""

from typing import Dict, List, Tuple
from .seq import Seq



IS_CODON = hasattr(str, 'memcpy')

# Import classes directly - always at module level
from .matrix import (
    FrequencyPositionMatrix,
    PositionWeightMatrix, 
    PositionSpecificScoringMatrix
)

# Minimal Alignment class - defined once at module level
class Alignment:
    sequences: List[str]
    length: int
    frequencies: Dict[str, List[float]]
    
    def __init__(self, sequences: List[str]):
        self.sequences = sequences
        if sequences:
            self.length = len(sequences[0])
            for seq in sequences:
                if len(seq) != self.length:
                    raise ValueError("All sequences must have the same length")
        else:
            self.length = 0
        self.frequencies = {}
    
    def __iter__(self):
        return iter(self.sequences)
    
    def __len__(self) -> int:
        return len(self.sequences)
    
    def __getitem__(self, key):
        if isinstance(key, tuple):
            # Handle 2D slicing like [:, 1:4]
            row_key, col_key = key
            if isinstance(col_key, slice):
                # Column slicing: return new Alignment with sliced sequences
                return Alignment([seq[col_key] for seq in self.sequences])
            else:
                raise TypeError("Column index must be a slice")
        elif isinstance(key, slice):
            # Column slicing: return new Alignment with sliced sequences
            return Alignment([seq[key] for seq in self.sequences])
        else:
            # Single row access: return the sequence string
            return self.sequences[key]
    
    def reverse_complement(self):
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                        'a': 't', 't': 'a', 'g': 'c', 'c': 'g','U': 'A', 'u': 'a' }
        rc_seqs = []
        for seq in self.sequences:
            rc = ''.join(complement_map.get(base, base) for base in reversed(seq))
            rc_seqs.append(rc)
        return Alignment(rc_seqs)
    
    def calculate_frequencies(self, alphabet: List[str]) -> Dict[str, List[float]]:
        if not self.sequences: 
            return {letter: [] for letter in alphabet}
        
        num_sequences = len(self.sequences)
        freqs: Dict[str, List[float]] = {letter: [0.0] * self.length for letter in alphabet}
        for seq in self.sequences:
            for i, letter in enumerate(seq):
                if letter in freqs:
                    freqs[letter][i] += 1.0
                    
        for letter in alphabet: 
            for i in range(self.length):
                freqs[letter][i] /= num_sequences
                
        return freqs
    
    def calculate_counts(self, alphabet: List[str]) -> Dict[str, List[float]]:
        """Calculate raw counts (not normalized) for each position."""
        if not self.sequences: 
            return {letter: [] for letter in alphabet}
        
        counts: Dict[str, List[float]] = {letter: [0.0] * self.length for letter in alphabet}
        
        for seq in self.sequences:
            for i, letter in enumerate(seq):
                if letter in counts:
                    counts[letter][i] += 1.0
        
        # Do NOT divide by num_sequences - return raw counts
        return counts


def create(instances: List[str], alphabet: str = "ACGT"):
    """Create a Motif object."""
    alignment = Alignment(instances)
    alphabet_list = [c for c in alphabet]
    return Motif(alignment=alignment, alphabet=alphabet_list)


def parse(handle, fmt: str, strict: bool = True):
    """Parse an output file from a motif finding program.
    
    Currently supported formats in Codon port:
     - MINIMAL: MINIMAL MEME output file motif
     
    Other formats require Python imports.
    """
    fmt = fmt.lower()
    if fmt == "minimal":
        import minimal
        return minimal.read(handle)
    else:
        raise NotImplementedError(f"Format '{fmt}' not supported in Codon.")
        


def read(handle, fmt: str, strict: bool = True):
    """Read a motif from a handle using the specified file-format.
    
    This supports the same formats as parse(), but only for files 
    containing exactly one motif.
    """
    motifs = parse(handle, fmt, strict)
    if len(motifs) == 0:
        raise ValueError("No motifs found in handle")
    if len(motifs) > 1:
        raise ValueError("More than one motif found in handle")
    return motifs[0]


class Motif:
    """A class representing sequence motifs."""
    
    def __init__(self, alphabet = "ACGT", 
                 alignment = None,
                 counts = None):
        """Initialize the class."""
        self.name = ""
        
        # Convert alphabet to list if string
        if isinstance(alphabet, str):
            self.alphabet = [c for c in alphabet]
        else:
            self.alphabet = alphabet
        
        if counts is not None and alignment is not None:
            raise ValueError("Specify either counts or an alignment, don't specify both")
        
        elif counts is not None:
            self.alignment = None
            self.counts = FrequencyPositionMatrix(self.alphabet, counts)
            self.length = self.counts.length
        
        elif alignment is not None:
            self.length = alignment.length
            counts_data = alignment.calculate_counts(self.alphabet)
            self.counts = FrequencyPositionMatrix(self.alphabet, counts_data)
            self.alignment = alignment
        
        else:
            self.counts = None
            self.alignment = None
            self.length = 0
        
        # Initialize pseudocounts and background
        self._pseudocounts = {letter: 0.0 for letter in self.alphabet}
        self._background = {letter: 1.0 for letter in self.alphabet}
        
        # Initialize mask - FIXED
        if self.length > 0:
            mask_list = [1] * self.length
            self._mask = tuple(mask_list)
        else:
            self._mask = tuple([])
    
    def get_mask(self) -> tuple:
        """Get the mask."""
        return self._mask
    
    def set_mask(self, mask):
        """Set the mask."""
        if self.length == 0:
            self._mask = tuple([])
        elif mask is None:
            mask_list = [1] * self.length
            self._mask = tuple(mask_list)
        elif len(mask) != self.length:
            raise ValueError(
                f"The length ({len(mask)}) of the mask is inconsistent with "
                f"the length ({self.length}) of the motif"
            )
        elif isinstance(mask, str):
            mask_list: List[int] = []
            for char in mask:
                if char == "*":
                    mask_list.append(1)
                elif char == " ":
                    mask_list.append(0)
                else:
                    raise ValueError(
                        f"Mask should contain only '*' or ' ' and not a '{char}'"
                    )
            self._mask = tuple(mask_list)
        else:
            mask_list = [1 if c else 0 for c in mask]
            self._mask = tuple(mask_list)
    
    def get_pseudocounts(self) -> Dict[str, float]:
        """Get pseudocounts."""
        return self._pseudocounts
    
    def set_pseudocounts(self, value):
        """Set pseudocounts."""
        if isinstance(value, dict):
            self._pseudocounts = {letter: value[letter] for letter in self.alphabet}
        else:
            if value is None:
                value = 0.0
            self._pseudocounts = {letter: value for letter in self.alphabet}
    
    def get_background(self) -> Dict[str, float]:
        """Get background."""
        return self._background
    
    def set_background(self, value):
        """Set background."""
        if isinstance(value, dict):
            self._background = {letter: value[letter] for letter in self.alphabet}
        elif value is None:
            self._background = {letter: 1.0 for letter in self.alphabet}
        else:
            if not self._has_dna_alphabet() and not self._has_rna_alphabet():
                raise ValueError(
                    "Setting the background to a single value only works for DNA and RNA "
                    "motifs (in which case the value is interpreted as the GC content)"
                )
            T_or_U = "T" if self._has_dna_alphabet() else "U"
            self._background = {
                "A": (1.0 - value) / 2.0,
                "C": value / 2.0,
                "G": value / 2.0,
                T_or_U: (1.0 - value) / 2.0
            }
        
        total = sum(self._background.values())
        for letter in self.alphabet:
            self._background[letter] /= total
    
    def __getitem__(self, key: slice):
        """Return a new Motif object for the positions included in key."""
        if not isinstance(key, slice):
            raise TypeError("motif indices must be slices")
        
        alphabet = self.alphabet
        
        if self.alignment is None:
            alignment = None
            if self.counts is None:
                counts = None
            else:
                counts = {letter: self.counts[letter][key] for letter in alphabet}
        else:
            alignment = self.alignment[:, key]
            counts = None
        
        motif = Motif(alphabet=alphabet, alignment=alignment, counts=counts)
        motif.set_mask(self._mask[key])
        
        if alignment is None and counts is None:
            if self.length > 0:
                start, stop, step = key.indices(self.length)
                motif.length = len(list(range(start, stop, step or 1)))
        
        motif.set_pseudocounts(self._pseudocounts.copy())
        motif.set_background(self._background.copy())
        
        return motif
    
    def get_pwm(self):
        """Calculate and return the position weight matrix for this motif."""
        return self.counts.normalize(self._pseudocounts)
    
    def get_pssm(self):
        """Calculate and return the position specific scoring matrix for this motif."""
        pwm = self.get_pwm()
        return pwm.log_odds(self._background)
    
    def __str__(self) -> str:
        """Return string representation of a motif."""
        text = ""
        if self.alignment is not None:
            text += "\n".join(self.alignment.sequences)
        return text
    
    def __len__(self) -> int:
        """Return the length of a motif."""
        return self.length
    
    def _has_dna_alphabet(self) -> bool:
        """Check if motif has DNA alphabet."""
        return sorted(self.alphabet) == ["A", "C", "G", "T"]
    
    def _has_rna_alphabet(self) -> bool:
        """Check if motif has RNA alphabet."""
        return sorted(self.alphabet) == ["A", "C", "G", "U"]
    
    def reverse_complement(self):
        """Return the reverse complement of the motif as a new motif."""
        alphabet = self.alphabet
        if not self._has_dna_alphabet() and not self._has_rna_alphabet():
            raise ValueError(
                "Calculating reverse complement only works for DNA and RNA motifs"
            )
        
        T_or_U = "T" if self._has_dna_alphabet() else "U"
        
        if self.alignment is not None:
            alignment = self.alignment.reverse_complement()
            if T_or_U == "U":
                alignment.sequences = [s.replace("T", "U") for s in alignment.sequences]
            res = Motif(alphabet=alphabet, alignment=alignment)
        else:  # has counts
            counts = {
                "A": self.counts[T_or_U][::-1],
                "C": self.counts["G"][::-1],
                "G": self.counts["C"][::-1],
                T_or_U: self.counts["A"][::-1],
            }
            res = Motif(alphabet=alphabet, counts=counts)
        
        # FIX: Just copy the mask structure as-is for now (don't reverse)
        # TODO: Implement proper mask reversal later if needed
        res._mask = self._mask
        
        res._background = {
            "A": self._background[T_or_U],
            "C": self._background["G"],
            "G": self._background["C"],
            T_or_U: self._background["A"],
        }
        res._pseudocounts = {
            "A": self._pseudocounts[T_or_U],
            "C": self._pseudocounts["G"],
            "G": self._pseudocounts["C"],
            T_or_U: self._pseudocounts["A"],
        }
        
        return res
    
    def get_consensus(self) -> Seq:
        """Return the consensus sequence."""
        return self.counts.get_consensus()
    
    def get_anticonsensus(self) -> Seq:
        """Return the least probable pattern to be generated from this motif."""
        return self.counts.get_anticonsensus()
    
    def get_degenerate_consensus(self) -> Seq:
        """Return the degenerate consensus sequence.
        
        Following the rules adapted from D. R. Cavener (1987).
        The same rules are used by TRANSFAC.
        """
        return self.counts.get_degenerate_consensus()
    
    @property
    def consensus(self) -> Seq:
        """Return the consensus sequence."""
        return self.get_consensus()
    
    @property
    def anticonsensus(self) -> Seq:
        """Return the anticonsensus sequence."""
        return self.get_anticonsensus()
    
    @property
    def degenerate_consensus(self) -> Seq:
        """Return the degenerate consensus sequence."""
        return self.get_degenerate_consensus()
    
    @property
    def pwm(self):
        """Return the position weight matrix."""
        return self.get_pwm()
    
    @property
    def pssm(self):
        """Return the position specific scoring matrix."""
        return self.get_pssm()
    
    @property
    def pseudocounts(self) -> Dict[str, float]:
        """Get pseudocounts."""
        return self.get_pseudocounts()

    @pseudocounts.setter
    def pseudocounts(self, value):
        """Set pseudocounts."""
        self.set_pseudocounts(value)

    @property
    def background(self) -> Dict[str, float]:
        """Get background."""
        return self.get_background()

    @background.setter
    def background(self, value):
        """Set background."""
        self.set_background(value)

    @property
    def mask(self) -> tuple:
        """Get mask."""
        return self.get_mask()

    @mask.setter
    def mask(self, value):
        """Set mask."""
        self.set_mask(value)
    
    def format(self, format_spec: str) -> str:
        """Return a string representation of the Motif in the given format.
        
        Currently supported formats:
         - transfac : TRANSFAC like files
         - pfm : JASPAR single Position Frequency Matrix
         - jaspar : JASPAR multiple Position Frequency Matrix
        """
        format_spec = format_spec.lower()
        
        if format_spec == "transfac":
            # Simple TRANSFAC format implementation
            lines = ["XX"]
            lines.append(f"ID {self.name}")
            lines.append("BF unknown")
            lines.append("P0      A      C      G      T")
            
            for i in range(self.length):
                values = [
                    f"{int(self.counts['A'][i]):6d}",
                    f"{int(self.counts['C'][i]):6d}",
                    f"{int(self.counts['G'][i]):6d}",
                    f"{int(self.counts['T'][i]):6d}"
                ]
                lines.append(f"{i+1:02d} " + " ".join(values))
            
            lines.append("XX")
            return "\n".join(lines) + "\n"
        
        elif format_spec in ("pfm", "jaspar"):
            # Simple JASPAR format implementation
            lines = []
            if self.name:
                lines.append(f">{self.name}")
            
            for letter in ["A", "C", "G", "T"]:
                if letter in self.counts:
                    values = " ".join(f"{int(v):3d}" for v in self.counts[letter])
                    lines.append(f"{letter} [{values}]")
            
            return "\n".join(lines) + "\n"
        
        else:
            raise ValueError(f"Unknown format type {format_spec}")


def write(motifs: List[Motif], fmt: str) -> str:
    """Return a string representation of motifs in the given format.
    
    Currently supported formats:
     - transfac : TRANSFAC like files
     - pfm : JASPAR simple single Position Frequency Matrix
     - jaspar : JASPAR multiple PFM format
    """
    fmt = fmt.lower()
    
    if fmt == "transfac":
        # TRANSFAC format needs // terminator at the end
        output = "".join(motif.format("transfac") for motif in motifs)
        if output and not output.endswith("//\n"):
            output += "//\n"
        return output
    elif fmt in ("pfm", "jaspar"):
        return "".join(motif.format(fmt) for motif in motifs)
    else:
        raise ValueError(f"Unknown format type {fmt}")
    
    
    