/**
 * algorithms.js — Client-side De Bruijn genome assembler (JavaScript port).
 *
 * Mirrors the Python pipeline in reads.py → debruijn.py → assembler.py.
 * Used as an offline fallback when the Flask API is unreachable, and
 * independently tested via Jest (algorithms.test.js).
 *
 * Pipeline:
 *   decomposeReads(reads, k)  →  edges[]
 *   buildGraph(edges)         →  { adjacency, inDegree, outDegree }
 *   findEulerianPath(graph)   →  path[]
 *   pathToContig(path)        →  string
 *   assemble(reads, k)        →  AssemblyResult
 */

// ─── Constants ────────────────────────────────────────────────────────────────

const VALID_DNA_RE = /^[ACGTNacgtn]+$/;

// ─── Validation ───────────────────────────────────────────────────────────────

/**
 * Validate a list of DNA read strings.
 *
 * @param {string[]} reads
 * @param {number}   k
 * @returns {string[]} Array of error messages (empty = valid).
 */
export function validateReads(reads, k) {
  const errors = [];

  if (!Array.isArray(reads) || reads.length === 0) {
    errors.push("reads must be a non-empty array.");
    return errors;
  }

  if (!Number.isInteger(k) || k < 2) {
    errors.push(`k must be an integer >= 2, got ${k}.`);
  }

  let minLen = Infinity;

  reads.forEach((read, i) => {
    if (typeof read !== "string" || read.trim().length === 0) {
      errors.push(`Read ${i + 1}: empty or not a string.`);
      return;
    }
    const upper = read.trim().toUpperCase();
    if (!VALID_DNA_RE.test(upper)) {
      errors.push(
        `Read ${i + 1}: invalid characters in "${upper}". Only A, C, G, T, N are allowed.`
      );
    }
    if (upper.length < minLen) minLen = upper.length;
  });

  if (errors.length === 0 && Number.isInteger(k) && k > minLen) {
    errors.push(
      `k (${k}) is longer than the shortest read (${minLen} bp). Reduce k or provide longer reads.`
    );
  }

  return errors;
}

// ─── K-mer decomposition ──────────────────────────────────────────────────────

/**
 * Decompose a single DNA read into overlapping k-mer edge tuples.
 *
 * Each k-mer kmer[i : i+k] produces an edge:
 *   prefix = kmer.slice(0, k-1)
 *   suffix = kmer.slice(1, k)
 *
 * @param {string} read  - Uppercase DNA string.
 * @param {number} k     - K-mer length (>= 2).
 * @returns {{ prefix: string, suffix: string }[]}
 *
 * @example
 * decomposeRead("ACGTTGCA", 4)
 * // → [{ prefix:"ACG", suffix:"CGT" }, { prefix:"CGT", suffix:"GTT" }, ...]
 */
export function decomposeRead(read, k) {
  const upper = read.toUpperCase();
  const edges = [];
  for (let i = 0; i <= upper.length - k; i++) {
    const kmer = upper.slice(i, i + k);
    edges.push({ prefix: kmer.slice(0, k - 1), suffix: kmer.slice(1) });
  }
  return edges;
}

/**
 * Decompose a list of reads into all (prefix, suffix) edge tuples.
 *
 * @param {string[]} reads
 * @param {number}   k
 * @returns {{ prefix: string, suffix: string }[]}
 */
export function decomposeReads(reads, k) {
  const edges = [];
  for (const read of reads) {
    const trimmed = read.trim();
    if (trimmed.length >= k) {
      edges.push(...decomposeRead(trimmed, k));
    }
  }
  return edges;
}

// ─── De Bruijn graph ──────────────────────────────────────────────────────────

/**
 * @typedef {Object} DeBruijnGraph
 * @property {Object.<string, string[]>} adjacency  - Node → list of neighbour nodes.
 * @property {Object.<string, number>}   inDegree   - Node → in-degree count.
 * @property {Object.<string, number>}   outDegree  - Node → out-degree count.
 * @property {number}                    edgeCount  - Total directed edges.
 */

/**
 * Build a De Bruijn graph from a list of (prefix, suffix) edge tuples.
 *
 * @param {{ prefix: string, suffix: string }[]} edges
 * @returns {DeBruijnGraph}
 */
export function buildGraph(edges) {
  /** @type {Object.<string, string[]>} */
  const adjacency = {};
  /** @type {Object.<string, number>} */
  const inDegree  = {};
  /** @type {Object.<string, number>} */
  const outDegree = {};

  const ensureNode = (node) => {
    if (!(node in adjacency)) adjacency[node] = [];
    if (!(node in inDegree))   inDegree[node]  = 0;
    if (!(node in outDegree))  outDegree[node] = 0;
  };

  for (const { prefix, suffix } of edges) {
    ensureNode(prefix);
    ensureNode(suffix);
    adjacency[prefix].push(suffix);
    outDegree[prefix] += 1;
    inDegree[suffix]  += 1;
  }

  const edgeCount = Object.values(outDegree).reduce((sum, d) => sum + d, 0);

  return { adjacency, inDegree, outDegree, edgeCount };
}

// ─── Eulerian path detection ──────────────────────────────────────────────────

/**
 * @typedef {Object} GraphStats
 * @property {boolean}       hasEulerianPath    - True if exactly one start and one end node.
 * @property {boolean}       hasEulerianCircuit - True if all degrees are balanced.
 * @property {string|null}   startNode          - Node to begin traversal from.
 * @property {number}        nodeCount
 * @property {number}        edgeCount
 */

/**
 * Analyse the graph and determine Eulerian traversal properties.
 *
 * @param {DeBruijnGraph} graph
 * @returns {GraphStats}
 */
export function getGraphStats(graph) {
  const { adjacency, inDegree, outDegree, edgeCount } = graph;
  const nodes = Object.keys(adjacency);

  const startCandidates = [];
  const endCandidates   = [];

  for (const node of nodes) {
    const diff = (outDegree[node] || 0) - (inDegree[node] || 0);
    if      (diff ===  1) startCandidates.push(node);
    else if (diff === -1) endCandidates.push(node);
  }

  const hasEulerianPath    = startCandidates.length === 1 && endCandidates.length === 1;
  const hasEulerianCircuit = startCandidates.length === 0 && endCandidates.length === 0;

  let startNode = null;
  if (hasEulerianPath) {
    startNode = startCandidates[0];
  } else if (hasEulerianCircuit) {
    startNode = nodes.find((n) => (outDegree[n] || 0) > 0) ?? null;
  }

  return {
    hasEulerianPath,
    hasEulerianCircuit,
    startNode,
    nodeCount: nodes.length,
    edgeCount,
  };
}

// ─── Hierholzer's algorithm ───────────────────────────────────────────────────

/**
 * Find an Eulerian path or circuit using Hierholzer's iterative algorithm.
 *
 * Runs in O(V + E). Works on a mutable copy of the adjacency list so the
 * original graph is not modified.
 *
 * @param {DeBruijnGraph} graph
 * @returns {string[]} Ordered list of nodes forming the Eulerian path/circuit.
 * @throws {Error} If no valid Eulerian traversal exists.
 */
export function findEulerianPath(graph) {
  const stats = getGraphStats(graph);

  if (!stats.hasEulerianPath && !stats.hasEulerianCircuit) {
    throw new Error(
      "Graph does not have an Eulerian path or circuit. " +
      "The reads may not provide sufficient coverage for the chosen k."
    );
  }

  if (stats.startNode === null) {
    throw new Error("Could not determine a valid start node for Eulerian traversal.");
  }

  // Deep-copy adjacency so we can consume edges during traversal
  const adj = {};
  for (const [node, neighbours] of Object.entries(graph.adjacency)) {
    adj[node] = [...neighbours];
  }

  return _hierholzer(adj, stats.startNode);
}

/**
 * Core iterative Hierholzer traversal.
 *
 * @param {Object.<string, string[]>} adj   - Mutable adjacency list (consumed in place).
 * @param {string}                    start - Start node.
 * @returns {string[]}
 */
function _hierholzer(adj, start) {
  const stack  = [start];
  const result = [];

  while (stack.length > 0) {
    const current = stack[stack.length - 1];
    const neighbours = adj[current];
    if (neighbours && neighbours.length > 0) {
      stack.push(neighbours.shift()); // follow next available edge
    } else {
      result.push(stack.pop());       // dead end — belongs in path
    }
  }

  result.reverse();
  return result;
}

// ─── Contig reconstruction ────────────────────────────────────────────────────

/**
 * Reconstruct a DNA string from an ordered list of (k-1)-mer nodes.
 *
 * Consecutive nodes share a (k-2)-mer overlap, so:
 *   contig = path[0] + path[1].slice(-1) + path[2].slice(-1) + ...
 *
 * @param {string[]} path - Ordered list of (k-1)-mer nodes.
 * @returns {string}
 * @throws {Error} If path is empty.
 *
 * @example
 * pathToContig(["ACG","CGT","GTT","TTG","TGC","GCA"])
 * // → "ACGTTGCA"
 */
export function pathToContig(path) {
  if (!path || path.length === 0) {
    throw new Error("Cannot reconstruct a contig from an empty path.");
  }
  return path[0] + path.slice(1).map((node) => node[node.length - 1]).join("");
}

// ─── Connected components ─────────────────────────────────────────────────────

/**
 * Find weakly-connected components of the graph using BFS on an undirected view.
 *
 * @param {DeBruijnGraph} graph
 * @returns {string[][]} Array of components; each component is an array of node labels.
 */
export function findConnectedComponents(graph) {
  const { adjacency } = graph;
  const allNodes = Object.keys(adjacency);
  const visited  = new Set();
  const components = [];

  // Build undirected neighbour map
  /** @type {Object.<string, Set<string>>} */
  const undirected = {};
  for (const [node, neighbours] of Object.entries(adjacency)) {
    if (!undirected[node]) undirected[node] = new Set();
    for (const nb of neighbours) {
      undirected[node].add(nb);
      if (!undirected[nb]) undirected[nb] = new Set();
      undirected[nb].add(node);
    }
  }

  for (const startNode of allNodes) {
    if (visited.has(startNode)) continue;
    const component = [];
    const queue = [startNode];
    visited.add(startNode);

    while (queue.length > 0) {
      const current = queue.shift();
      component.push(current);
      for (const nb of (undirected[current] || [])) {
        if (!visited.has(nb)) {
          visited.add(nb);
          queue.push(nb);
        }
      }
    }
    components.push(component);
  }

  return components;
}

/**
 * Extract a subgraph containing only the edges where both endpoints are in nodeSet.
 *
 * @param {DeBruijnGraph} graph
 * @param {string[]}      nodes
 * @returns {{ prefix: string, suffix: string }[]} Edge list for the subgraph.
 */
function subgraphEdges(graph, nodes) {
  const nodeSet = new Set(nodes);
  const edges   = [];
  for (const [source, neighbours] of Object.entries(graph.adjacency)) {
    if (!nodeSet.has(source)) continue;
    for (const dest of neighbours) {
      if (nodeSet.has(dest)) edges.push({ prefix: source, suffix: dest });
    }
  }
  return edges;
}

// ─── Component-level assembly ─────────────────────────────────────────────────

/**
 * Assemble contigs from a single weakly-connected component graph.
 *
 * Handles standard Eulerian path, Eulerian circuit, and partially-unbalanced
 * graphs (greedy multi-start traversal) — matching the Python assembler behaviour.
 *
 * @param {DeBruijnGraph} componentGraph
 * @returns {string[]} Assembled contig strings from this component.
 */
function assembleComponent(componentGraph) {
  const stats = getGraphStats(componentGraph);

  if (stats.hasEulerianPath || stats.hasEulerianCircuit) {
    try {
      const path   = findEulerianPath(componentGraph);
      const contig = pathToContig(path);
      return [contig];
    } catch (_) {
      // fall through to greedy
    }
  }

  // Greedy fallback: traverse from each semi-balanced start node
  const { adjacency, outDegree, inDegree } = componentGraph;
  const allNodes = Object.keys(adjacency);

  // Mutable copy shared across starts
  const adj = {};
  for (const [node, nbs] of Object.entries(adjacency)) {
    adj[node] = [...nbs];
  }

  const startNodes = allNodes.filter(
    (n) => (outDegree[n] || 0) - (inDegree[n] || 0) === 1
  );
  const fallbackStarts = startNodes.length > 0
    ? startNodes
    : allNodes.filter((n) => (outDegree[n] || 0) > 0);

  const contigs = [];
  for (const start of fallbackStarts) {
    if (!adj[start] || adj[start].length === 0) continue;
    const path = _hierholzer(adj, start);
    if (path.length >= 2) contigs.push(pathToContig(path));
  }

  return contigs;
}

// ─── AssemblyResult ───────────────────────────────────────────────────────────

/**
 * @typedef {Object} AssemblyResult
 * @property {string[]} contigs      - Assembled contigs, longest-first.
 * @property {number}   totalKmers   - Total k-mer edges (including duplicates).
 * @property {number}   graphNodes   - Unique (k-1)-mer nodes in the graph.
 * @property {number}   graphEdges   - Total directed edges in the graph.
 * @property {number}   elapsedMs    - Wall-clock assembly time in milliseconds.
 * @property {string}   longestContig
 * @property {number}   n50
 */

/**
 * Compute the N50 statistic for a list of contigs.
 *
 * @param {string[]} contigs
 * @returns {number}
 */
export function computeN50(contigs) {
  if (!contigs || contigs.length === 0) return 0;
  const lengths = contigs.map((c) => c.length).sort((a, b) => b - a);
  const total   = lengths.reduce((s, l) => s + l, 0);
  let cumulative = 0;
  for (const len of lengths) {
    cumulative += len;
    if (cumulative >= total / 2) return len;
  }
  return 0;
}

// ─── Main entry point ─────────────────────────────────────────────────────────

/**
 * Assemble contigs from raw DNA reads — full client-side pipeline.
 *
 * Drop-in offline equivalent of POST /api/assemble. Returns an object with
 * the same shape as the Flask API response so App.js can use either source.
 *
 * @param {string[]} reads - Raw DNA sequences (case-insensitive).
 * @param {number}   k     - K-mer length (default 4).
 * @returns {AssemblyResult}
 * @throws {Error} If reads are invalid or k is out of range.
 *
 * @example
 * const result = assemble(["ACGTTGCATGC", "TGCATGCAAC"], 4);
 * console.log(result.longestContig); // "ACGTTGCATGCATGCAAC"
 */
export function assemble(reads, k = 4) {
  const t0 = performance.now();

  // Validate
  const errors = validateReads(reads, k);
  if (errors.length > 0) throw new Error(errors.join("\n"));

  const normalised = reads.map((r) => r.trim().toUpperCase());

  // K-mer decomposition
  const edges = decomposeReads(normalised, k);
  if (edges.length === 0) {
    throw new Error("No k-mer edges produced. Check that reads are longer than k.");
  }

  // Build graph
  const graph = buildGraph(edges);

  // Find weakly-connected components and assemble each
  const components = findConnectedComponents(graph);
  const contigs = [];

  for (const componentNodes of components) {
    const subEdges = subgraphEdges(graph, componentNodes);
    if (subEdges.length === 0) continue;

    const componentGraph = buildGraph(subEdges);
    const componentContigs = assembleComponent(componentGraph);
    contigs.push(...componentContigs);
  }

  // Sort longest-first (matches Python assembler output order)
  contigs.sort((a, b) => b.length - a.length);

  const elapsedMs = performance.now() - t0;

  return {
    contigs,
    totalKmers:    edges.length,
    graphNodes:    graph.edgeCount > 0 ? Object.keys(graph.adjacency).length : 0,
    graphEdges:    graph.edgeCount,
    elapsedMs:     parseFloat(elapsedMs.toFixed(4)),
    longestContig: contigs[0] ?? "",
    n50:           computeN50(contigs),
  };
}