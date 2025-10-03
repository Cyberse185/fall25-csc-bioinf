"""Module for the support of MEME minimal motif format.

Ported to Codon from BioPython.
"""

from typing import Dict, List, Optional, Tuple
from __init__ import Motif


def read(handle):
    """Parse the text output of the MEME program into a Record object.
    
    This function won't retrieve instances, as there are none in minimal meme format.
    """
 
    
    
    motif_number = 0
    record = Record()
    _read_version(record, handle)
    _read_alphabet(record, handle)
    _read_background(record, handle)
    
    while True:
        # Look for next motif
        found_motif = False
        name = ""
        for line in handle:
            if line.startswith("MOTIF"):
                found_motif = True
                name = line.split()[1]
                break
        
        if not found_motif:
            return record
        
        motif_number += 1
        
        # Read the letter-probability matrix line directly
        for line in handle:
            line_stripped = line.strip()
            if line_stripped.startswith("letter-probability matrix:"):
                # Parse statistics from this line
                length, num_occurrences, evalue = _parse_motif_statistics(line_stripped)
                break
        else:
            raise ValueError("Could not find letter-probability matrix line")
        
        # Now read the frequency matrix
        counts = _read_lpm(record, handle, length, num_occurrences)
        
        # Create motif using locally imported module
        alphabet_list = [c for c in record.alphabet]

        motif =Motif(alphabet=alphabet_list, counts=counts)
        motif.background = record.background
        motif.length = motif.counts.length
        motif.num_occurrences = num_occurrences
        motif.evalue = evalue
        motif.name = name
        record.append(motif)
        assert len(record) == motif_number
    
    return record


class Record(list):
    """Class for holding the results of a minimal MEME run."""
    version: str
    datafile: str
    command: str
    alphabet: str
    background: Dict[str, float]  
    
    
    def __init__(self):
        """Initialize record class values."""
        super().__init__()
        self.version = ""
        self.datafile = ""
        self.command = ""
        self.alphabet = ""
        self.background = {}
        self.sequences = []
    
    def __getitem__(self, key):
        """Return the motif of index key."""
        if isinstance(key, str):
            for motif in self:
                if motif.name == key:
                    return motif
        else:
            return list.__getitem__(self, key)


# Everything below is private

def _read_version(record: Record, handle):
    """Read MEME version (PRIVATE)."""
    for line in handle:
        if line.startswith("MEME version"):
            break
    else:
        raise ValueError(
            "Improper input file. File should contain a line starting MEME version."
        )
    
    line = line.strip()
    ls = line.split()
    record.version = ls[2]


def _read_alphabet(record: Record, handle):
    """Read alphabet (PRIVATE)."""
    for line in handle:
        if line.startswith("ALPHABET"):
            break
    else:
        raise ValueError(
            "Unexpected end of stream: Expected to find line starting with 'ALPHABET'"
        )
    
    if not line.startswith("ALPHABET= "):
        raise ValueError(f"Line does not start with 'ALPHABET':\n{line}")
    
    line = line.strip().replace("ALPHABET= ", "")
    
    if line == "ACGT":
        al = "ACGT"
    elif line == "ACGU":
        al = "ACGU"
    else:
        raise ValueError("Only parsing of DNA and RNA motifs is implemented")
    
    record.alphabet = al


def _read_background(record: Record, handle):
    """Read background letter frequencies (PRIVATE)."""
    for line in handle:
        if line.startswith("Background letter frequencies"):
            background_freqs = []
            for line in handle:
                line = line.rstrip()
                if line:
                    background_freqs.extend(
                        [
                            float(freq)
                            for i, freq in enumerate(line.split(" "))
                            if i % 2 == 1
                        ]
                    )
                else:
                    break
            if not background_freqs:
                raise ValueError(
                    "Unexpected end of stream: Expected to find line starting background frequencies."
                )
            break
    else:
        raise ValueError(
            "Improper input file. File should contain a line starting background frequencies."
        )
    record.background = dict(zip(record.alphabet, background_freqs))


def _parse_motif_statistics(line: str) -> Tuple[Optional[int], int, float]:
    """Parse motif statistics from letter-probability matrix line (PRIVATE).
    
    Returns:
        Tuple of (length, num_occurrences, evalue)
    
    Example line:
         letter-probability matrix: alength= 4 w= 4 nsites= 17 E= 4.1e-009
    """
    # Parse nsites (defaults to 20)
    num_occurrences = 20
    if "nsites=" in line:
        try:
            nsites_part = line.split("nsites=")[1].split()[0]
            num_occurrences = int(nsites_part)
        except (IndexError, ValueError):
            num_occurrences = 20
    
    # Parse length (can be None, inferred later)
    length: Optional[int] = None
    if "w=" in line:
        try:
            w_part = line.split("w=")[1].split()[0]
            length = int(w_part)
        except (IndexError, ValueError):
            length = None
    
    # Parse E-value (defaults to 0.0)
    evalue = 0.0
    if "E=" in line:
        try:
            e_part = line.split("E=")[1].split()[0]
            evalue = float(e_part)
        except (IndexError, ValueError):
            evalue = 0.0
    
    return length, num_occurrences, evalue


def _read_lpm(record: Record, handle, length: Optional[int], 
              num_occurrences: int) -> Dict[str, List[float]]:
    """Read letter probability matrix (PRIVATE)."""
    counts: List[List[float]] = [[], [], [], []]
    
    for line in handle:
        line = line.strip()
        
        # Skip empty lines
        if not line:
            continue
        
        freqs = line.split()
        
        # Check if this line has exactly 4 frequency values
        if len(freqs) != 4:
            break
        
        # Try to parse as floats - if they're not numbers, break
        try:
            f0 = float(freqs[0])
            f1 = float(freqs[1])
            f2 = float(freqs[2])
            f3 = float(freqs[3])
        except ValueError:
            break
        
        counts[0].append(round(f0 * num_occurrences))
        counts[1].append(round(f1 * num_occurrences))
        counts[2].append(round(f2 * num_occurrences))
        counts[3].append(round(f3 * num_occurrences))
        
        # If we've read the expected length, stop
        if length is not None and len(counts[0]) == length:
            break
    
    # Validate that we got some data
    if not counts[0]:
        raise ValueError("No valid frequency data found in letter-probability matrix")
    
    # Create counts dictionary
    c: Dict[str, List[float]] = {}
    for i, letter in enumerate(record.alphabet):
        if i < len(counts):
            c[letter] = counts[i]
        else:
            c[letter] = []
    
    return c


def _read_motif_name(handle) -> str:
    """Read motif name (PRIVATE)."""
    for line in handle:
        if "sorted by position p-value" in line:
            break
    else:
        raise ValueError("Unexpected end of stream: Failed to find motif name")
    
    line = line.strip()
    words = line.split()
    name = " ".join(words[0:2])
    return name