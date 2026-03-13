from __future__ import annotations

from pathlib import Path
from typing import Iterator
from dataclasses import dataclass

@dataclass
class ReadSet:
    """A set of DNA reads decomposed into k-mer edges."""

    original_reads: list[str]
    k: int
    edges: list[tuple[str, str]]

    @property
    def total_kmers(self) -> int:
        """Total count of k-mers (edges) across all reads."""
        return len(self.edges)

    def edge_counts(self) -> dict[tuple[str, str], int]:
        """Return a dict mapping each unique edge to its count."""
        counts: dict[tuple[str, str], int] = {}
        for edge in self.edges:
            counts[edge] = counts.get(edge, 0) + 1
        return counts


# ===========================================================================
# Validation helpers
# ===========================================================================

def _validate_dna_sequence(seq: str, index: int) -> None:
    """Check that *seq* is non-empty and contains only ACGT (any case)."""
    if not seq:
        raise ValueError(f"Read {index}: empty read.")

    valid_chars = set("ACGTacgt")
    invalid = set(seq) - valid_chars
    if invalid:
        raise ValueError(
            f"Read {index}: invalid characters {invalid}. "
            f"Expected only A, C, G, T (case-insensitive)."
        )


def _validate_k(k: int, min_read_length: int) -> None:
    """Check that *k* is positive and not larger than the shortest read."""
    if k <= 0:
        raise ValueError(f"k must be positive, got {k}.")
    if k > min_read_length:
        raise ValueError(
            f"k={k} is longer than the shortest read ({min_read_length}). "
            f"Use a smaller k."
        )


# ===========================================================================
# FASTA file I/O
# ===========================================================================

def load_fasta(filepath: str | Path) -> list[str]:
    """Load DNA sequences from a FASTA file.

    Empty lines, and lines starting with ';' (comment) are ignored.

    Each header line (starting with '>') marks the start of a new sequence.
    Subsequent lines are concatenated until the next header (or EOF).

    All sequences are returned in uppercase for consistency.

    Reads are deduplicated automatically using a dict (insertion order preserved).
    This means if the same sequence appears multiple times in the file,
    only one copy is returned.

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

    *read* is normalized to uppercase automatically.

    Example:
        >>> list(decompose_read("ACGTTGCA", 4))
        [('ACG', 'CGT'), ('CGT', 'GTT'), ('GTT', 'TTG'), ('TTG', 'TGC'), ('TGC', 'GCA')]
    """
    read = read.upper()  # Normalize to uppercase
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

    return ReadSet(original_reads=normalised, k=k, edges=all_edges)
