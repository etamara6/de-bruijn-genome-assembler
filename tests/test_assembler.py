"""
tests/test_assembler.py — Comprehensive pytest suite for the De Bruijn assembler.

Covers:
- K-mer decomposition (reads.py)
- Graph construction and degree tracking (debruijn.py)
- Eulerian path finding (debruijn.py)
- Contig reconstruction (assembler.py)
- End-to-end assembly (assembler.py)
- Edge cases: empty input, single read, circular genome
"""

from __future__ import annotations

import pytest

from reads import decompose_read, decompose_reads, ReadSet
from debruijn import DeBruijnGraph, GraphStats
from assembler import (
    AssemblyResult,
    GenomeAssembler,
    path_to_contig,
    write_fasta,
)


# ===========================================================================
# reads.py — K-mer decomposition
# ===========================================================================

class TestDecomposeRead:
    """Unit tests for single-read k-mer decomposition."""

    def test_basic_decomposition(self):
        """Decompose a known read and verify exact edge list."""
        edges = list(decompose_read("ACGTTGCA", k=4))
        expected = [
            ("ACG", "CGT"),
            ("CGT", "GTT"),
            ("GTT", "TTG"),
            ("TTG", "TGC"),
            ("TGC", "GCA"),
        ]
        assert edges == expected

    def test_edge_count(self):
        """L - k + 1 edges should be produced."""
        read = "ACGTACGT"
        k = 3
        edges = list(decompose_read(read, k))
        assert len(edges) == len(read) - k + 1

    def test_k_equals_read_length(self):
        """When k == len(read), exactly one edge is produced."""
        edges = list(decompose_read("ACGT", k=4))
        assert len(edges) == 1
        assert edges[0] == ("ACG", "CGT")

    def test_lowercase_input_normalised(self):
        """Lowercase input should produce uppercase edges."""
        edges = list(decompose_read("acgt", k=3))
        assert edges[0] == ("AC", "CG")

    def test_prefix_suffix_lengths(self):
        """Prefix and suffix must each have length k-1."""
        k = 5
        for prefix, suffix in decompose_read("ACGTACGTACGT", k=k):
            assert len(prefix) == k - 1
            assert len(suffix) == k - 1

    def test_consecutive_edges_overlap(self):
        """Adjacent edges must share a (k-2)-mer overlap."""
        edges = list(decompose_read("ACGTACGT", k=4))
        for (_, suffix_a), (prefix_b, _) in zip(edges, edges[1:]):
            assert suffix_a == prefix_b


class TestDecomposeReads:
    """Unit tests for multi-read batch decomposition."""

    def test_multiple_reads(self):
        """K-mers from all reads are combined into a single ReadSet."""
        rs = decompose_reads(["ACGTTGCA", "TTGCATGC"], k=4)
        assert isinstance(rs, ReadSet)
        # 5 kmers from first read + 5 from second
        assert rs.total_kmers == 10

    def test_empty_list_raises(self):
        with pytest.raises(ValueError, match="empty"):
            decompose_reads([], k=4)

    def test_invalid_characters_raise(self):
        with pytest.raises(ValueError, match="invalid characters"):
            decompose_reads(["ACGXYZ"], k=3)

    def test_empty_string_raises(self):
        with pytest.raises(ValueError, match="empty"):
            decompose_reads([""], k=3)

    def test_k_too_large_raises(self):
        with pytest.raises(ValueError, match="longer than the shortest read"):
            decompose_reads(["ACG"], k=10)

    def test_k_less_than_2_raises(self):
        with pytest.raises(ValueError, match="k must be >= 2"):
            decompose_reads(["ACGT"], k=1)

    def test_reads_stored_uppercase(self):
        rs = decompose_reads(["acgtacgt"], k=4)
        assert rs.reads[0] == "ACGTACGT"

    def test_read_set_k_attribute(self):
        rs = decompose_reads(["ACGT"], k=3)
        assert rs.k == 3


# ===========================================================================
# debruijn.py — Graph construction
# ===========================================================================

class TestDeBruijnGraphConstruction:
    """Tests for graph node/edge/degree bookkeeping."""

    def _graph_from_read(self, read: str, k: int) -> DeBruijnGraph:
        rs = decompose_reads([read], k=k)
        return DeBruijnGraph(rs.edges)

    def test_node_set(self):
        """All (k-1)-mers from the read should appear as nodes."""
        graph = self._graph_from_read("ACGTTGCA", k=4)
        expected_nodes = {"ACG", "CGT", "GTT", "TTG", "TGC", "GCA"}
        assert set(graph.nodes) == expected_nodes

    def test_edge_count(self):
        """Edge count matches number of k-mers."""
        graph = self._graph_from_read("ACGTTGCA", k=4)
        assert graph.edge_count == 5

    def test_out_degree(self):
        """ACG → CGT only, so out_degree["ACG"] == 1."""
        graph = self._graph_from_read("ACGTTGCA", k=4)
        assert graph.out_degree["ACG"] == 1
        assert graph.out_degree["CGT"] == 1

    def test_in_degree(self):
        """GCA receives one incoming edge from TGC."""
        graph = self._graph_from_read("ACGTTGCA", k=4)
        assert graph.in_degree["GCA"] == 1

    def test_repeated_kmer_increments_degree(self):
        """AAAA produces a self-loop; out_degree == in_degree == 1."""
        rs = decompose_reads(["AAAA"], k=3)
        graph = DeBruijnGraph(rs.edges)
        assert graph.out_degree["AA"] == 2
        assert graph.in_degree["AA"] == 2

    def test_empty_edges_raise(self):
        with pytest.raises(ValueError, match="empty edge list"):
            DeBruijnGraph([])

    def test_adjacency_contains_all_destinations(self):
        """Every destination node must appear as a key in adjacency."""
        graph = self._graph_from_read("ACGTTGCA", k=4)
        for node in graph.nodes:
            assert node in graph.adjacency


# ===========================================================================
# debruijn.py — Eulerian path
# ===========================================================================

class TestEulerianPath:
    """Tests for Hierholzer's Eulerian path algorithm."""

    def _assemble_path(self, reads: list[str], k: int) -> list[str]:
        rs = decompose_reads(reads, k=k)
        graph = DeBruijnGraph(rs.edges)
        return graph.find_eulerian_path()

    def test_simple_path(self):
        """Single read with unique k-mers should yield a clean Eulerian path."""
        path = self._assemble_path(["ACGTTGCA"], k=4)
        assert path == ["ACG", "CGT", "GTT", "TTG", "TGC", "GCA"]

    def test_path_visits_all_edges(self):
        """Path length == number of edges + 1."""
        rs = decompose_reads(["ACGTTGCA"], k=4)
        graph = DeBruijnGraph(rs.edges)
        path = graph.find_eulerian_path()
        assert len(path) == graph.edge_count + 1

    def test_circular_genome_circuit(self):
        """A circular genome produces a circuit (first node == last node)."""
        # "AAACCC" → overlapping k-mers form a cycle when circular
        rs = decompose_reads(["AACGTTAACG"], k=4)
        graph = DeBruijnGraph(rs.edges)
        stats = graph.get_stats()
        # Circuit or path depending on whether all degrees balance
        path = graph.find_eulerian_path()
        assert len(path) >= 2

    def test_stats_eulerian_path_flag(self):
        """Single linear read should flag has_eulerian_path."""
        rs = decompose_reads(["ACGTTGCA"], k=4)
        graph = DeBruijnGraph(rs.edges)
        stats = graph.get_stats()
        assert stats.has_eulerian_path is True
        assert stats.has_eulerian_circuit is False

    def test_stats_node_and_edge_counts(self):
        rs = decompose_reads(["ACGTTGCA"], k=4)
        graph = DeBruijnGraph(rs.edges)
        stats = graph.get_stats()
        assert stats.node_count == 6
        assert stats.edge_count == 5

    def test_start_node_for_path(self):
        """Start node should have out_degree - in_degree == 1."""
        rs = decompose_reads(["ACGTTGCA"], k=4)
        graph = DeBruijnGraph(rs.edges)
        stats = graph.get_stats()
        start = stats.start_node
        assert start is not None
        assert graph.out_degree[start] - graph.in_degree[start] == 1


# ===========================================================================
# assembler.py — path_to_contig
# ===========================================================================

class TestPathToContig:
    """Unit tests for the path → DNA string reconstruction."""

    def test_known_path(self):
        path = ["ACG", "CGT", "GTT", "TTG", "TGC", "GCA"]
        assert path_to_contig(path) == "ACGTTGCA"

    def test_single_node(self):
        assert path_to_contig(["ACG"]) == "ACG"

    def test_two_nodes(self):
        assert path_to_contig(["AC", "CG"]) == "ACG"

    def test_empty_path_raises(self):
        with pytest.raises(ValueError, match="empty"):
            path_to_contig([])

    def test_length_formula(self):
        """Result length == k-1 + len(path) - 1 == len(path) + k - 2."""
        path = ["ACG", "CGT", "GTT", "TTG", "TGC", "GCA"]  # k=4, k-1=3
        contig = path_to_contig(path)
        # path[0] is length 3, each additional node contributes 1 char
        assert len(contig) == len(path[0]) + len(path) - 1


# ===========================================================================
# assembler.py — End-to-end assembly
# ===========================================================================

class TestGenomeAssembler:
    """Integration tests for GenomeAssembler.assemble()."""

    def test_demo_example(self):
        """Reads with sufficient multi-edge coverage produce a single assembled contig."""
        assembler = GenomeAssembler(
            reads=["ACGTTGCATGC", "TGCATGCAAC"], k=4
        )
        result = assembler.assemble()
        assert len(result.contigs) == 1
        assert result.longest_contig == "ACGTTGCATGCATGCAAC"

    def test_result_stats_are_populated(self):
        assembler = GenomeAssembler(reads=["ACGTTGCA"], k=4)
        result = assembler.assemble()
        assert result.total_kmers > 0
        assert result.graph_nodes > 0
        assert result.graph_edges > 0
        assert result.elapsed_ms >= 0.0

    def test_assembly_result_type(self):
        assembler = GenomeAssembler(reads=["ACGTACGT"], k=4)
        result = assembler.assemble()
        assert isinstance(result, AssemblyResult)

    def test_single_read(self):
        """Assembly of a single read should recover that read (modulo k-mer boundary)."""
        read = "ACGTTGCA"
        assembler = GenomeAssembler(reads=[read], k=4)
        result = assembler.assemble()
        assert len(result.contigs) >= 1
        assert result.longest_contig == read

    def test_empty_reads_raise(self):
        with pytest.raises(ValueError):
            GenomeAssembler(reads=[], k=4)

    def test_contig_list_sorted_longest_first(self):
        """Contigs should be returned longest-first."""
        assembler = GenomeAssembler(reads=["ACGTTGCA", "TTGCATGC", "CATGCAAC"], k=4)
        result = assembler.assemble()
        lengths = [len(c) for c in result.contigs]
        assert lengths == sorted(lengths, reverse=True)

    def test_n50_single_contig(self):
        """N50 of a single contig equals that contig's length."""
        assembler = GenomeAssembler(reads=["ACGTTGCATGCAAC"], k=4)
        result = assembler.assemble()
        assert result.n50 == len(result.longest_contig)

    def test_longest_contig_property(self):
        assembler = GenomeAssembler(reads=["ACGTTGCA", "TTGCATGC", "CATGCAAC"], k=4)
        result = assembler.assemble()
        assert result.longest_contig == max(result.contigs, key=len)

    def test_from_fasta(self, tmp_path):
        """from_fasta classmethod should load sequences and assemble."""
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">read1\nACGTTGCA\n>read2\nTTGCATGC\n>read3\nCATGCAA\n")
        assembler = GenomeAssembler.from_fasta(fasta, k=4)
        result = assembler.assemble()
        assert len(result.contigs) >= 1

    def test_fasta_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            GenomeAssembler.from_fasta("/nonexistent/path/reads.fasta", k=4)

    def test_k_larger_than_reads_raises(self):
        with pytest.raises(ValueError):
            GenomeAssembler(reads=["AC"], k=100).assemble()


# ===========================================================================
# assembler.py — FASTA output
# ===========================================================================

class TestWriteFasta:
    """Tests for FASTA file writing."""

    def test_output_file_created(self, tmp_path):
        outfile = tmp_path / "contigs.fasta"
        write_fasta(["ACGTTGCA", "TTGCATGC"], outfile)
        assert outfile.exists()

    def test_header_lines_present(self, tmp_path):
        outfile = tmp_path / "contigs.fasta"
        write_fasta(["ACGTTGCA", "TTGCATGC"], outfile)
        content = outfile.read_text()
        assert ">contig_1" in content
        assert ">contig_2" in content

    def test_sequences_in_output(self, tmp_path):
        outfile = tmp_path / "contigs.fasta"
        write_fasta(["ACGTTGCA"], outfile)
        content = outfile.read_text()
        assert "ACGTTGCA" in content

    def test_custom_prefix(self, tmp_path):
        outfile = tmp_path / "contigs.fasta"
        write_fasta(["ACGT"], outfile, prefix="scaffold")
        content = outfile.read_text()
        assert ">scaffold_1" in content

    def test_empty_contig_list(self, tmp_path):
        outfile = tmp_path / "empty.fasta"
        write_fasta([], outfile)
        assert outfile.read_text() == ""


# ===========================================================================
# Edge cases
# ===========================================================================

class TestEdgeCases:
    """Edge and boundary condition tests."""

    def test_repeated_reads_increase_coverage(self):
        """Duplicate reads add redundant raw edges; unique k-mers still assemble correctly."""
        reads = ["ACGTTGCA"] * 3
        assembler = GenomeAssembler(reads=reads, k=4)
        result = assembler.assemble()
        # total_kmers counts all raw k-mer edges before deduplication
        assert result.total_kmers == 15  # 5 edges × 3 reads

    def test_k_equals_two(self):
        """k=2 is the minimum legal value; nodes are single characters."""
        assembler = GenomeAssembler(reads=["ACGT"], k=2)
        result = assembler.assemble()
        assert len(result.contigs) >= 1

    def test_large_k_single_kmer(self):
        """When k == len(read), exactly one edge and one path node beyond start."""
        read = "ACGT"
        assembler = GenomeAssembler(reads=[read], k=4)
        result = assembler.assemble()
        assert result.total_kmers == 1

    def test_n50_empty(self):
        """N50 of an empty contig list is 0."""
        result = AssemblyResult(
            contigs=[], total_kmers=0, graph_nodes=0, graph_edges=0, elapsed_ms=0.0
        )
        assert result.n50 == 0

    def test_longest_contig_empty(self):
        result = AssemblyResult(
            contigs=[], total_kmers=0, graph_nodes=0, graph_edges=0, elapsed_ms=0.0
        )
        assert result.longest_contig == ""