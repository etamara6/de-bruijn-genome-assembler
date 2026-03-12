"""
assembler.py — Contig reconstruction from Eulerian paths in a De Bruijn graph.

Ties together read ingestion, graph construction, and path traversal to produce
assembled contigs from raw DNA reads.  Also provides FASTA output and a CLI demo.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path


from reads import ReadSet, decompose_reads, load_fasta
from debruijn import DeBruijnGraph


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class AssemblyResult:
    """Complete output of a genome assembly run.

    Attributes:
        contigs: Assembled DNA sequences, ordered longest-first.
        total_kmers: Total k-mer edges used to build the graph.
        graph_nodes: Number of unique (k-1)-mer nodes in the graph.
        graph_edges: Total directed edges in the graph.
        elapsed_ms: Wall-clock time for the assembly in milliseconds.
    """
    contigs: list[str]
    total_kmers: int
    graph_nodes: int
    graph_edges: int
    elapsed_ms: float

    @property
    def longest_contig(self) -> str:
        """Return the longest assembled contig, or an empty string if none."""
        return max(self.contigs, key=len, default="")

    @property
    def n50(self) -> int:
        """N50 length: the contig length at which 50 % of total assembly is covered."""
        if not self.contigs:
            return 0
        sorted_lengths = sorted((len(c) for c in self.contigs), reverse=True)
        total_length = sum(sorted_lengths)
        cumulative = 0
        for contig_length in sorted_lengths:
            cumulative += contig_length
            if cumulative >= total_length / 2:
                return contig_length
        return 0


# ---------------------------------------------------------------------------
# Path → contig conversion
# ---------------------------------------------------------------------------

def path_to_contig(path: list[str]) -> str:
    """Reconstruct a DNA string from an ordered list of (k-1)-mer nodes.

    Consecutive nodes share a (k-2)-mer overlap, so the assembled sequence is
    the first node followed by the last character of every subsequent node.

    Args:
        path: Ordered list of (k-1)-mer strings from an Eulerian traversal.

    Returns:
        The reconstructed DNA contig as a single uppercase string.

    Raises:
        ValueError: If *path* is empty.

    Example:
        >>> path_to_contig(["ACG", "CGT", "GTT", "TTG", "TGC", "GCA"])
        'ACGTTGCA'
    """
    if not path:
        raise ValueError("Cannot reconstruct a contig from an empty path.")
    return path[0] + "".join(node[-1] for node in path[1:])


def write_fasta(contigs: list[str], filepath: str | Path, prefix: str = "contig") -> None:
    """Write assembled contigs to a FASTA file.

    Args:
        contigs: List of DNA sequences to write, ordered as-is.
        filepath: Destination file path.
        prefix: Header prefix for each record (default ``"contig"``).
    """
    path = Path(filepath)
    with path.open("w") as fasta_out:
        for index, contig in enumerate(contigs, start=1):
            fasta_out.write(f">{prefix}_{index} length={len(contig)}\n")
            # Wrap at 60 characters per line (FASTA convention)
            for chunk_start in range(0, len(contig), 60):
                fasta_out.write(contig[chunk_start : chunk_start + 60] + "\n")


# ---------------------------------------------------------------------------
# Component assembly helper
# ---------------------------------------------------------------------------

def _assemble_component(graph: "DeBruijnGraph") -> list[str]:
    """Assemble contigs from a single weakly-connected graph component.

    Handles three cases:
    1. Standard Eulerian path (one start, one end node) → one contig.
    2. Eulerian circuit (all degrees balanced) → one contig.
    3. Partially unbalanced graph (multiple start/end nodes) → greedy traversal
       from each semi-balanced start node to maximise contig recovery.

    Args:
        graph: A weakly-connected De Bruijn graph component.

    Returns:
        List of assembled contig strings from this component.
    """
    from debruijn import DeBruijnGraph

    stats = graph.get_stats()

    if stats.has_eulerian_path or stats.has_eulerian_circuit:
        try:
            path = graph.find_eulerian_path()
            return [path_to_contig(path)]
        except RuntimeError:
            pass

    # Fallback: find all nodes where out > in (potential starts) and
    # greedily traverse from each, sharing the mutable adjacency copy.
    adjacency_copy: dict[str, list[str]] = {
        node: list(neighbours) for node, neighbours in graph.adjacency.items()
    }

    all_nodes = set(graph.adjacency.keys())
    start_nodes = [
        node for node in all_nodes
        if graph.out_degree[node] - graph.in_degree[node] == 1
    ]
    if not start_nodes:
        # circuit fallback: pick any node with outgoing edges
        start_nodes = [n for n in all_nodes if graph.out_degree[n] > 0]

    contigs: list[str] = []
    visited_all_edges = False

    for start_node in start_nodes:
        if not adjacency_copy.get(start_node):
            continue
        path = DeBruijnGraph._hierholzer(adjacency_copy, start_node)
        if len(path) >= 2:
            contigs.append(path_to_contig(path))

    return contigs




class GenomeAssembler:
    """End-to-end genome assembler using a De Bruijn graph and Eulerian traversal.

    Usage::

        assembler = GenomeAssembler(reads=["ACGTTGCA", "TTGCATGC"], k=4)
        result = assembler.assemble()
        print(result.longest_contig)

    Attributes:
        reads: The input DNA sequences.
        k: The k-mer length.
    """

    def __init__(self, reads: list[str], k: int = 4) -> None:
        """Initialise the assembler with reads and k-mer length.

        Args:
            reads: List of DNA sequences.  Case-insensitive.
            k: K-mer length (>= 2).  Use k=4 for testing, k=21–31 for real data.

        Raises:
            ValueError: If *reads* is empty or *k* is invalid.
        """
        if not reads:
            raise ValueError("reads list must not be empty.")
        self.reads: list[str] = reads
        self.k: int = k

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def assemble(self) -> AssemblyResult:
        """Run the full assembly pipeline and return an :class:`AssemblyResult`.

        Pipeline steps:
            1. Decompose reads into (prefix, suffix) k-mer edges.
            2. Build a De Bruijn graph.
            3. Find weakly-connected components.
            4. For each component, attempt Eulerian traversal(s) to produce contigs.
            5. Sort contigs longest-first and return results.

        Returns:
            An :class:`AssemblyResult` with contigs and assembly statistics.

        Raises:
            ValueError: If reads are empty or k is invalid (propagated from
                        :func:`~reads.decompose_reads`).
            RuntimeError: If Eulerian traversal fails for any component.
        """
        start_time = time.perf_counter()

        # Step 1 — k-mer decomposition
        read_set: ReadSet = decompose_reads(self.reads, self.k)

        # Step 2 — graph construction
        graph = DeBruijnGraph(read_set.edges)

        # Step 3 — find connected components and assemble each
        components = graph.find_connected_components()
        contigs: list[str] = []

        for component_nodes in components:
            try:
                component_graph = graph.subgraph(component_nodes)
            except ValueError:
                continue

            component_contigs = _assemble_component(component_graph)
            contigs.extend(component_contigs)

        # Sort longest-first
        contigs.sort(key=len, reverse=True)

        elapsed_ms = (time.perf_counter() - start_time) * 1000.0

        return AssemblyResult(
            contigs=contigs,
            total_kmers=read_set.total_kmers,
            graph_nodes=len(graph.nodes),
            graph_edges=graph.edge_count,
            elapsed_ms=elapsed_ms,
        )

    # ------------------------------------------------------------------
    # Convenience: assemble from FASTA
    # ------------------------------------------------------------------

    @classmethod
    def from_fasta(cls, filepath: str | Path, k: int = 31) -> "GenomeAssembler":
        """Construct a :class:`GenomeAssembler` by loading reads from a FASTA file.

        Args:
            filepath: Path to the FASTA input file.
            k: K-mer length (default 31 for real genome data).

        Returns:
            A ready-to-use :class:`GenomeAssembler` instance.

        Raises:
            FileNotFoundError: If the FASTA file does not exist.
            ValueError: If no sequences are found or k is invalid.
        """
        reads = load_fasta(filepath)
        return cls(reads=reads, k=k)


# ---------------------------------------------------------------------------
# CLI demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # These reads tile the genome ACGTTGCATGCAAC with sufficient overlap
    # so that multi-edge coverage produces a single Eulerian path.
    reads = ["ACGTTGCATGC", "TGCATGCAAC"]
    k = 4

    print("=" * 50)
    print("De Bruijn Graph Genome Assembler — Demo")
    print("=" * 50)
    print(f"Input reads : {reads}")
    print(f"k           : {k}")
    print()

    assembler = GenomeAssembler(reads=reads, k=k)
    result = assembler.assemble()

    print(f"Contigs assembled : {len(result.contigs)}")
    print(f"Longest contig    : {result.longest_contig}")
    print(f"Total k-mers      : {result.total_kmers}")
    print(f"Graph nodes       : {result.graph_nodes}")
    print(f"Graph edges       : {result.graph_edges}")
    print(f"N50               : {result.n50}")
    print(f"Time              : {result.elapsed_ms:.2f} ms")

    if len(result.contigs) > 1:
        print()
        print("All contigs:")
        for idx, contig in enumerate(result.contigs, start=1):
            print(f"  [{idx}] {contig}  (len={len(contig)})")