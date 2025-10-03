"""Shared Seq class for codon_motifs package.

This module contains the Seq class that is shared across all motif modules
to avoid type conflicts in Codon.
"""

class Seq:
    """Minimal sequence class."""
    sequence: str
    
    def __init__(self, data: str):
        self.sequence = data
    
    def __str__(self) -> str:
        return self.sequence
    
    def __repr__(self) -> str:
        return f"Seq('{self.sequence}')"
    
    def __len__(self) -> int:
        return len(self.sequence)
    
    def upper(self) -> "Seq":
        return Seq(self.sequence.upper())
    
    def replace(self, old: str, new: str) -> "Seq":
        return Seq(self.sequence.replace(old, new))