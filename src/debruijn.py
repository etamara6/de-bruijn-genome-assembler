"""
debruijn.py — De Bruijn graph construction and Eulerian path/circuit traversal.

Builds a directed graph from (prefix, suffix) k-mer edge tuples and implements
Hierholzer's O(E) algorithm for finding Eulerian paths and circuits.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from reads import ReadSet


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class GraphStats:
    """Summary statistics for a De Bruijn graph.

    Attributes:
        node_count: Number of unique (k-1)-mer nodes.
        edge_count: Total number of directed edges (k-mer occurrences).
        has_eulerian_path: True if an Eulerian path (not circuit) exists.
        has_eulerian_circuit: True if an Eulerian circuit exists.
        start_node: Designated start node for traversal, or None if the graph
                    is disconnected / has no valid Eulerian traversal.
    """
    node_count: int
    edge_count: int
    has_eulerian_path: bool
    has_eulerian_circuit: bool
    start_node: Optional[str]


# ---------------------------------------------------------------------------
# De Bruijn graph
# ---------------------------------------------------------------------------

class DeBruijnGraph:
    """Directed De Bruijn graph built from (prefix, suffix) k-mer edges.

    Nodes are (k-1)-mers.  Each edge represents a k-mer whose left half is the
    source node and whose right half is the destination node.

    Attributes:
        adjacency: Mapping from each node to its list of neighbour nodes.
        in_degree: Mapping from each node to its in-degree count.
        out_degree: Mapping from each node to its out-degree count.

    Example:
        >>> from reads import decompose_reads
        >>> rs = decompose_reads(["ACGTTGCA"], k=4)
        >>> g = DeBruijnGraph(rs.edges)
        >>> path = g.find_eulerian_path()
        >>> path
        ['ACG', 'CGT', 'GTT', 'TTG', 'TGC', 'GCA']
    """

    def __init__(self, edges: List[Tuple[str, str]]) -> None:
        """Build the graph from a list of (prefix, suffix) edge tuples.

        Args:
            edges: Each tuple (u, v) adds a directed edge u → v.

        Raises:
            ValueError: If *edges* is empty.
        """
        if not edges:
            raise ValueError(
                "Cannot build a De Bruijn graph from an empty edge list. "
                "Check that your reads are long enough for the chosen k."
            )

        self.adjacency: Dict[str, List[str]] = defaultdict(list)
        self.in_degree: Dict[str, int] = defaultdict(int)
        self.out_degree: Dict[str, int] = defaultdict(int)

        self._build(edges)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def _build(self, edges: List[Tuple[str, str]]) -> None:
        """Populate adjacency list and degree tables from *edges*.

        Args:
            edges: List of (source, destination) node pairs.
        """
        for source, destination in edges:
            self.adjacency[source].append(destination)
            self.out_degree[source] += 1
            self.in_degree[destination] += 1
            # Ensure every destination node exists as a key in adjacency
            if destination not in self.adjacency:
                self.adjacency[destination] = []

    # ------------------------------------------------------------------
    # Degree inspection
    # ------------------------------------------------------------------

    @property
    def nodes(self) -> List[str]:
        """Sorted list of all node labels in the graph."""
        return sorted(self.adjacency.keys())

    @property
    def edge_count(self) -> int:
        """Total number of directed edges."""
        return sum(self.out_degree.values())

    # ------------------------------------------------------------------
    # Eulerian traversal helpers
    # ------------------------------------------------------------------

    def _classify_nodes(
        self,
    ) -> Tuple[List[str], List[str], List[str]]:
        """Classify nodes by their degree imbalance for Eulerian analysis.

        Returns:
            A 3-tuple of:
                - start_candidates: nodes where out - in == 1
                - end_candidates:   nodes where in  - out == 1
                - balanced:         nodes where in  == out
        """
        start_candidates: List[str] = []
        end_candidates: List[str] = []
        balanced: List[str] = []

        all_nodes = set(self.adjacency.keys()) | set(self.in_degree.keys())
        for node in all_nodes:
            diff = self.out_degree[node] - self.in_degree[node]
            if diff == 1:
                start_candidates.append(node)
            elif diff == -1:
                end_candidates.append(node)
            elif diff == 0:
                balanced.append(node)
            else:
                # diff outside {-1, 0, 1} means no Eulerian path/circuit exists
                pass

        return start_candidates, end_candidates, balanced

    def get_stats(self) -> GraphStats:
        """Return summary statistics and Eulerian property flags for this graph.

        Returns:
            A :class:`GraphStats` dataclass instance.
        """
        start_candidates, end_candidates, _ = self._classify_nodes()

        has_circuit = len(start_candidates) == 0 and len(end_candidates) == 0
        has_path = len(start_candidates) == 1 and len(end_candidates) == 1

        if has_path:
            start_node: Optional[str] = start_candidates[0]
        elif has_circuit:
            # Any node with outgoing edges can start the circuit
            start_node = next(
                (n for n in self.adjacency if self.out_degree[n] > 0), None
            )
        else:
            start_node = None

        return GraphStats(
            node_count=len(self.adjacency),
            edge_count=self.edge_count,
            has_eulerian_path=has_path,
            has_eulerian_circuit=has_circuit,
            start_node=start_node,
        )

    # ------------------------------------------------------------------
    # Hierholzer's algorithm
    # ------------------------------------------------------------------

    def find_eulerian_path(self) -> List[str]:
        """Find an Eulerian path or circuit using Hierholzer's algorithm.

        Hierholzer's algorithm runs in O(V + E) time by maintaining a stack
        and appending nodes to the result only when they have no remaining
        outgoing edges.

        Returns:
            An ordered list of nodes representing the Eulerian path/circuit.
            The path visits every edge exactly once.

        Raises:
            RuntimeError: If the graph has no valid Eulerian path or circuit
                          (i.e., more than one node has degree imbalance != {-1, 0, 1}).
        """
        stats = self.get_stats()

        if not stats.has_eulerian_path and not stats.has_eulerian_circuit:
            raise RuntimeError(
                "Graph does not have an Eulerian path or circuit. "
                "The input reads may be insufficient to cover the genome uniformly, "
                "or k is too large for the given read set."
            )

        start_node = stats.start_node
        if start_node is None:
            raise RuntimeError(
                "Could not determine a valid start node for Eulerian traversal."
            )

        # Work on a mutable copy so the original adjacency is preserved
        adjacency_copy: Dict[str, List[str]] = {
            node: list(neighbours) for node, neighbours in self.adjacency.items()
        }

        path: List[str] = self._hierholzer(adjacency_copy, start_node)
        return path

    @staticmethod
    def _hierholzer(
        adjacency: Dict[str, List[str]],
        start_node: str,
    ) -> List[str]:
        """Core iterative implementation of Hierholzer's algorithm.

        Args:
            adjacency: Mutable adjacency list (will be consumed during traversal).
            start_node: The node from which traversal begins.

        Returns:
            Ordered list of nodes forming the Eulerian path/circuit.
        """
        stack: List[str] = [start_node]
        result: List[str] = []

        while stack:
            current_node = stack[-1]
            neighbours = adjacency.get(current_node, [])
            if neighbours:
                # Follow the next available edge
                next_node = neighbours.pop(0)
                stack.append(next_node)
            else:
                # Dead end — this node belongs in the path
                result.append(stack.pop())

        result.reverse()
        return result

    # ------------------------------------------------------------------
    # Connected components for multi-contig assembly
    # ------------------------------------------------------------------

    def find_connected_components(self) -> List[List[str]]:
        """Return weakly-connected components as lists of node labels.

        Treats the directed graph as undirected for connectivity analysis.
        Each component is a candidate for independent Eulerian traversal.

        Returns:
            List of components; each component is a list of node labels.
        """
        all_nodes = set(self.adjacency.keys())
        visited: set[str] = set()
        components: List[List[str]] = []

        # Build an undirected neighbour map for BFS
        undirected: Dict[str, set[str]] = defaultdict(set)
        for node, neighbours in self.adjacency.items():
            for neighbour in neighbours:
                undirected[node].add(neighbour)
                undirected[neighbour].add(node)

        for node in all_nodes:
            if node in visited:
                continue
            # BFS to collect this component
            component: List[str] = []
            queue: List[str] = [node]
            visited.add(node)
            while queue:
                current = queue.pop(0)
                component.append(current)
                for neighbour in undirected[current]:
                    if neighbour not in visited:
                        visited.add(neighbour)
                        queue.append(neighbour)
            components.append(component)

        return components

    def subgraph(self, nodes: List[str]) -> "DeBruijnGraph":
        """Create a new :class:`DeBruijnGraph` restricted to *nodes*.

        Args:
            nodes: The subset of node labels to include.

        Returns:
            A new graph containing only edges where both endpoints are in *nodes*.

        Raises:
            ValueError: If the induced subgraph has no edges.
        """
        node_set = set(nodes)
        sub_edges: List[Tuple[str, str]] = [
            (source, dest)
            for source, neighbours in self.adjacency.items()
            if source in node_set
            for dest in neighbours
            if dest in node_set
        ]
        return DeBruijnGraph(sub_edges)