"""
utils.py — Helper functions: reverse complement, sequence stats.
"""

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(complement)[::-1]

def gc_content(sequence: str) -> float:
    """Return GC content as a percentage (0–100)."""
    if not sequence:
        return 0.0
    gc = sum(1 for base in sequence.upper() if base in "GC")
    return (gc / len(sequence)) * 100