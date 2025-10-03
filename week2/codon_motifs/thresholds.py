"""Approximate calculation of appropriate thresholds for motif finding. Ported to Codon from BioPython."""

from typing import Dict, List, Optional, Tuple

# Detect Codon vs Python using str.memcpy (Codon-specific attribute)
IS_CODON = hasattr(str, 'memcpy')

if IS_CODON:
    from .matrix import PositionSpecificScoringMatrix
else:
    from .matrix import PositionSpecificScoringMatrix


class ScoreDistribution:
    """Class representing approximate score distribution for a given motif.

    Utilizes a dynamic programming approach to calculate the distribution of
    scores with a predefined precision. Provides a number of methods for calculating
    thresholds for motif occurrences.
    """
    
    min_score: float
    interval: float
    n_points: int
    ic: float
    step: float
    mo_density: List[float]
    bg_density: List[float]
    
    def __init__(self, motif=None, precision: int = 1000, 
                 pssm: Optional[PositionSpecificScoringMatrix] = None,
                 background: Optional[Dict[str, float]] = None):
        """Initialize the class."""
        if pssm is None:
            # Original motif-based initialization (not implemented in our port)
            if motif is None:
                raise ValueError("Either motif or pssm must be provided")
            # This branch would require the full Motif class with ic() method
            raise NotImplementedError(
                "Initialization from motif object not yet ported. Use pssm parameter instead."
            )
        else:
            # PSSM-based initialization
            self.min_score = min(0.0, pssm.min)
            self.interval = max(0.0, pssm.max) - self.min_score
            self.n_points = precision * pssm.length
            self.ic = pssm.mean(background)
        
        self.step = self.interval / float(self.n_points - 1)
        
        # Initialize density arrays
        self.mo_density = [0.0] * self.n_points
        self.mo_density[self._add(0, -self._index_diff(self.min_score))] = 1.0
        
        self.bg_density = [0.0] * self.n_points
        self.bg_density[self._add(0, -self._index_diff(self.min_score))] = 1.0
        
        if pssm is not None:
            if background is None:
                background = {letter: 1.0 for letter in pssm.alphabet}
                total = sum(background.values())
                for letter in pssm.alphabet:
                    background[letter] /= total
            
            # Process each position in the PSSM
            for position in range(pssm.length):
                mo_new = [0.0] * self.n_points
                bg_new = [0.0] * self.n_points
                
                # Get scores for this position
                for letter in pssm.alphabet:
                    score = pssm[letter][position]
                    bg = background[letter]
                    mo = pow(2.0, score) * bg
                    d = self._index_diff(score)
                    
                    for i in range(self.n_points):
                        mo_new[self._add(i, d)] += self.mo_density[i] * mo
                        bg_new[self._add(i, d)] += self.bg_density[i] * bg
                
                self.mo_density = mo_new
                self.bg_density = bg_new
    
    def _index_diff(self, x: float, y: float = 0.0) -> int:
        """Calculate index difference for a score difference."""
        return int((x - y + 0.5 * self.step) / self.step)
    
    def _add(self, i: int, j: int) -> int:
        """Add indices with bounds checking."""
        return max(0, min(self.n_points - 1, i + j))
    
    def modify(self, scores: Dict[str, float], 
               mo_probs: Dict[str, float], 
               bg_probs: Dict[str, float]):
        """Modify motifs and background density."""
        mo_new = [0.0] * self.n_points
        bg_new = [0.0] * self.n_points
        
        for k, v in scores.items():
            d = self._index_diff(v)
            for i in range(self.n_points):
                mo_new[self._add(i, d)] += self.mo_density[i] * mo_probs[k]
                bg_new[self._add(i, d)] += self.bg_density[i] * bg_probs[k]
        
        self.mo_density = mo_new
        self.bg_density = bg_new
    
    def threshold_fpr(self, fpr: float) -> float:
        """Approximate the log-odds threshold which makes the type I error (false positive rate).
        
        Args:
            fpr: Desired false positive rate
            
        Returns:
            Threshold score
        """
        i = self.n_points
        prob = 0.0
        while prob < fpr:
            i -= 1
            if i < 0:
                break
            prob += self.bg_density[i]
        return self.min_score + float(i) * self.step
    
    def threshold_fnr(self, fnr: float) -> float:
        """Approximate the log-odds threshold which makes the type II error (false negative rate).
        
        Args:
            fnr: Desired false negative rate
            
        Returns:
            Threshold score
        """
        i = -1
        prob = 0.0
        while prob < fnr:
            i += 1
            if i >= self.n_points:
                break
            prob += self.mo_density[i]
        return self.min_score + float(i) * self.step
    
    def threshold_balanced(self, rate_proportion: float = 1.0, 
                          return_rate: bool = False) -> Tuple[float, float]:
        """Approximate log-odds threshold making FNR equal to FPR times rate_proportion.
        
        Args:
            rate_proportion: Proportion between FPR and FNR
            return_rate: If True, return both threshold and FPR
            
        Returns:
            If return_rate is False: threshold score
            If return_rate is True: (threshold score, fpr)
        """
        i = self.n_points
        fpr = 0.0
        fnr = 1.0
        
        while fpr * rate_proportion < fnr:
            i -= 1
            if i < 0:
                break
            fpr += self.bg_density[i]
            fnr -= self.mo_density[i]
        
        threshold = self.min_score + float(i) * self.step
        
        if return_rate:
            return (threshold, fpr)
        else:
            return (threshold, 0.0)  # Return tuple for consistency
    
    def threshold_patser(self) -> float:
        """Threshold selection mimicking the behaviour of patser (Hertz, Stormo 1999) software.

        It selects such a threshold that the log(fpr)=-ic(M)
        
        Note: the actual patser software uses natural logarithms instead of log_2, so the numbers
        are not directly comparable.
        
        Returns:
            Threshold score
        """
        return self.threshold_fpr(fpr=pow(2.0, -self.ic))