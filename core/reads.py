"""
reads.py — Read ingestion and k-mer decomposition for De Bruijn graph assembly.

Handles loading DNA reads from lists or FASTA files, decomposing them into
overlapping k-mers, and producing (prefix, suffix) edge tuples for graph construction.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class ReadSet:
    """Container for a collection of DNA reads and their derived k-mer edges.

    Attributes:
        reads: The raw DNA sequences loaded from input.
        k: The k-mer length used for decomposition.
        edges: List of (prefix, suffix) tuples where prefix = kmer[:-1]
               and suffix = kmer[1:], representing edges in the De Bruijn graph.
    """
    reads: list[str]
    k: int
    edges: list[tuple[str, str]] = field(default_factory=list)

    @property
    def total_kmers(self) -> int:
        """Total number of k-mer edges extracted from all reads."""
        return len(self.edges)


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

_VALID_DNA_RE = re.compile(r"^[ACGTNacgtn]+$")


def _validate_dna_sequence(sequence: str, index: int) -> None:
    """Raise ValueError if *sequence* is empty or contains non-DNA characters.

    Raises:
        ValueError: If the sequence is empty or contains invalid characters.
    """
    if not sequence:
        raise ValueError(f"Read at index {index} is empty.")
    if not _VALID_DNA_RE.match(sequence):
        raise ValueError(
            f"Read at index {index} contains invalid characters: '{sequence}'. "
            "Only A, C, G, T, N (case-insensitive) are allowed."
        )


def _validate_k(k: int, min_read_length: int) -> None:
    """Ensure k is a sensible value relative to the shortest read.

    Raises:
        ValueError: If k is less than 2 or longer than the shortest read.
    """
    if k < 2:
        raise ValueError(f"k must be >= 2, got {k}.")
    if k > min_read_length:
        raise ValueError(
            f"k ({k}) is longer than the shortest read ({min_read_length} bp). "
            "Reduce k or provide longer reads."
        )


# ---------------------------------------------------------------------------
# FASTA parser
# ---------------------------------------------------------------------------

def load_fasta(filepath: str | Path) -> list[str]:
    """Parse a FASTA file and return a list of uppercase DNA sequences.

    Multi-line sequences are concatenated.  Comment lines (starting with ';')
    and blank lines are ignored.

    Args:
        filepath: Path to the FASTA file.

    Returns:
        List of DNA sequences as uppercase strings.

    Raises:
        FileNotFoundError: If *filepath* does not exist.
        ValueError: If the file is empty or contains no valid sequences.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    sequences: list[str] = []
    current_sequence_parts: list[str] = []

    with path.open("r", encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if current_sequence_parts:
                    sequences.append("".join(current_sequence_parts).upper())
                    current_sequence_parts = []
            else:
                current_sequence_parts.append(line)

        # Flush the last record
        if current_sequence_parts:
            sequences.append("".join(current_sequence_parts).upper())

    if not sequences:
        raise ValueError(f"No sequences found in FASTA file: {path}")

    return sequences


# ---------------------------------------------------------------------------
# K-mer decomposition
# ---------------------------------------------------------------------------

def decompose_read(read: str, k: int) -> Iterator[tuple[str, str]]:
    """Yield (prefix, suffix) edge tuples from a single DNA read.

    For a read of length L and k-mer length k, this yields L - k + 1 edges.
    Each k-mer is split as:
        prefix = kmer[:-1]  (length k-1)
        suffix = kmer[1:]   (length k-1)

    *read* is expected to be uppercase; call ``read.upper()`` beforehand if needed.

    Example:
        >>> list(decompose_read("ACGTTGCA", 4))
        [('ACG', 'CGT'), ('CGT', 'GTT'), ('GTT', 'TTG'), ('TTG', 'TGC'), ('TGC', 'GCA')]
    """
    for start in range(len(read) - k + 1):
        kmer = read[start : start + k]
        yield kmer[:-1], kmer[1:]


def decompose_reads(reads: list[str], k: int) -> ReadSet:
    """Decompose a list of DNA reads into a :class:`ReadSet` with k-mer edges.

    Validates every read and the chosen k before decomposition.

    Args:
        reads: List of DNA sequences (strings).  Case-insensitive.
        k: The k-mer length.  Typical values: 4 (testing), 21–31 (real genomes).

    Returns:
        A :class:`ReadSet` containing the original reads, the k value, and all
        (prefix, suffix) edge tuples.

    Raises:
        ValueError: If *reads* is empty, any read is invalid, or k is out of range.

    Example:
        >>> rs = decompose_reads(["ACGTTGCA", "TTGCATGC"], k=4)
        >>> rs.total_kmers
        10
    """
    if not reads:
        raise ValueError("reads list is empty — nothing to assemble.")

    normalised: list[str] = []
    for index, read in enumerate(reads):
        _validate_dna_sequence(read, index)
        normalised.append(read.upper())

    _validate_k(k, min(len(r) for r in normalised))

    all_edges: list[tuple[str, str]] = []
    for read in normalised:
        all_edges.extend(decompose_read(read, k))

    return ReadSet(reads=normalised, k=k, edges=all_edges)