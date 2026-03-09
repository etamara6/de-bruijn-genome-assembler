# 🧬 De Bruijn Genome Assembler

[![Python](https://img.shields.io/badge/python-3.8%2B-3776ab?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Status](https://img.shields.io/badge/status-in%20development-orange)]()

> A genome assembler that reconstructs DNA sequences from short k-mer reads using De Bruijn graphs and Eulerian path traversal — implemented from scratch in Python.

---

## 🧠 What This Project Demonstrates

| Skill | How it shows up |
|---|---|
| Graph algorithm implementation | De Bruijn graph construction, Eulerian path/circuit detection, Hierholzer's algorithm |
| Bioinformatics knowledge | k-mer decomposition, read overlap, contig assembly, coverage depth |
| Software engineering | Clean OOP design, type hints, dataclasses, separation of concerns |
| Testing discipline | pytest suite covering graph construction, path finding, and assembly edge cases |
| CS fundamentals | Directed graphs, adjacency lists, DFS, Eulerian vs Hamiltonian paths |

---

## 🔬 How It Works — Algorithm Overview

### 1. K-mer Decomposition

A DNA read (e.g. `ACGTTGCA`) is broken into overlapping substrings of length `k`:

```
read:   ACGTTGCA   (k=4)
kmers:  ACGT → CGTT → GTTG → TTGC → TGCA
```

### 2. De Bruijn Graph Construction

Each k-mer `ACGT` becomes a directed edge from node `ACG` → `CGT` (prefix → suffix of length k-1).

```
        ACGT         CGTT         GTTG
ACG  ──────────▶  CGT  ──────────▶  GTT  ──────────▶  TTG  ...
```

All reads are decomposed and their edges merged into a single graph. Repeated k-mers increase edge multiplicity (coverage).

### 3. Eulerian Path / Circuit

The assembled genome corresponds to an **Eulerian path** through the graph — a path that visits every edge exactly once.

**Conditions checked:**
- Eulerian circuit: every node has `in-degree == out-degree`
- Eulerian path: exactly one node has `out − in = 1` (start), one has `in − out = 1` (end)

**Algorithm:** Hierholzer's algorithm — O(E) time complexity.

```
Start node → follow edges greedily → backtrack when stuck → splice in sub-tours
```

### 4. Contig Output

If the graph is not fully connected (real sequencing has gaps), multiple contigs are output — one per Eulerian path component.

---

## 🏗️ Planned Architecture

```
┌─────────────────────────┐      ┌──────────────────────────┐      ┌─────────────────────┐
│   Input Layer           │      │   Graph Engine           │      │   Output Layer      │
│   reads.py              │─────▶│   debruijn.py            │─────▶│   assembler.py      │
│   FASTA / raw reads     │      │   build_graph()          │      │   contigs → FASTA   │
│   k-mer size config     │      │   find_eulerian_path()   │      │   stats report      │
└─────────────────��───────┘      └──────────────────────────┘      └─────────────────────┘
```

---

## 📥 Example (Planned)

**Input reads:**
```
ACGTTGCA
TTGCATGC
CATGCAAC
```

**k = 4, assembled contig:**
```
ACGTTGCATGCAAC
```

---

## 🗺️ Roadmap

### Core Engine
- [ ] K-mer decomposition from raw reads (`reads.py`)
- [ ] De Bruijn graph construction with adjacency list (`debruijn.py`)
- [ ] Eulerian path detection (in/out-degree analysis)
- [ ] Hierholzer's algorithm for Eulerian path traversal
- [ ] Contig extraction and FASTA output (`assembler.py`)

### Robustness
- [ ] Handle sequencing errors (low-coverage edge pruning)
- [ ] Support for paired-end reads
- [ ] Bubble detection and resolution (heterozygous variants)
- [ ] Tip removal (dead-end branches from read errors)

### Visualisation & Interface
- [ ] DOT/Graphviz export of De Bruijn graph
- [ ] Interactive graph visualisation (React + D3.js, matching Sequence Alignment Studio style)
- [ ] REST API (`POST /api/assemble`) with Flask backend
- [ ] CLI interface: `python assemble.py --reads input.fa --k 31`

### Testing & Quality
- [ ] pytest suite: graph construction, Eulerian path, edge cases (empty, single read, circular)
- [ ] CI via GitHub Actions on every push
- [ ] Coverage reporting

### Performance
- [ ] NumPy/SciPy sparse matrix graph representation for large genomes
- [ ] Benchmarking against real sequencing datasets (E. coli, phage genomes)

---

## 📁 Planned Project Structure

```
de-bruijn-genome-assembler/
├── .github/workflows/ci.yml     # CI — pytest + coverage
├── src/
│   ├── reads.py                 # Read ingestion, k-mer decomposition
│   ├── debruijn.py              # De Bruijn graph class + Eulerian path
│   ├── assembler.py             # Contig assembly + FASTA output
│   └── utils.py                 # Helpers: reverse complement, stats
├── tests/
│   └── test_assembler.py        # pytest suite
├── data/
│   └── example_reads.fa         # Sample input for demo
├── requirements.txt
└── README.md
```

---

## 🔗 Related Project

This assembler is the natural successor to [Sequence Alignment Studio](https://github.com/etamara6/sequence-alignment-tool) — once reads are assembled into contigs, alignment algorithms (NW/SW) compare them to reference genomes.

---

## 📄 License

MIT © [etamara6](https://github.com/etamara6)
