# 🧬 De Bruijn Genome Assembler

[![Python](https://img.shields.io/badge/python-3.8%2B-3776ab?logo=python&logoColor=white)](https://www.python.org/)
[![React](https://img.shields.io/badge/react-18-61dafb?logo=react&logoColor=white)](https://react.dev/)
[![Flask](https://img.shields.io/badge/flask-2.3%2B-000000?logo=flask&logoColor=white)](https://flask.palletsprojects.com/)
[![Tests](https://img.shields.io/badge/tests-53%20passed-brightgreen)](tests/test_assembler.py)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

> A genome assembler that reconstructs DNA sequences from short k-mer reads using De Bruijn graphs and Eulerian path traversal — implemented from scratch in Python with a React frontend and Flask REST API.

---

## 🧠 What This Project Demonstrates

| Skill | How it shows up |
|---|---|
| Graph algorithm implementation | De Bruijn graph construction, Eulerian path/circuit detection, Hierholzer's O(E) algorithm |
| Bioinformatics knowledge | k-mer decomposition, read overlap, contig assembly, N50 coverage depth |
| Software engineering | Clean OOP design, type hints, dataclasses, separation of concerns |
| Testing discipline | 53-test pytest suite covering graph construction, path finding, and assembly edge cases |
| Full-stack development | React frontend + Flask REST API + JavaScript algorithm port |
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

All reads are decomposed and their edges merged into a single graph. Repeated k-mers increase edge multiplicity (coverage depth).

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

If the graph is not fully connected (real sequencing has gaps), multiple contigs are output — one per weakly-connected component, assembled independently.

---

## 🏗️ Architecture

```
┌─────────────────────────┐      ┌──────────────────────────┐      ┌─────────────────────┐
│   Input Layer           │      │   Graph Engine           │      │   Output Layer      │
│   reads.py              │─────▶│   debruijn.py            │─────▶│   assembler.py      │
│   FASTA / raw reads     │      │   DeBruijnGraph()        │      │   contigs → FASTA   │
│   k-mer size config     │      │   find_eulerian_path()   │      │   AssemblyResult    │
└─────────────────────────┘      └──────────────────────────┘      └─────────────────────┘
                                            │
                                            ▼
                               ┌──────────────────────────┐
                               │   REST API + Frontend    │
                               │   api.py  (Flask)        │
                               │   POST /api/assemble     │
                               │   src/App.js  (React)    │
                               └──────────────────────────┘
```

---

## 📁 Project Structure

```
de-bruijn-genome-assembler/
├── .github/workflows/ci.yml     # CI — pytest + coverage on Python 3.9/3.10/3.11
├── core/
│   ├── reads.py                 # Read ingestion, k-mer decomposition
│   ├── debruijn.py              # De Bruijn graph class + Eulerian path (Hierholzer's)
│   ├── assembler.py             # Contig assembly, N50, FASTA output
│   └── utils.py                 # Helpers: reverse complement, GC content
├── src/
│   ├── App.js                   # React frontend — input, visualisation, results
│   ├── algorithms.js            # Full JavaScript port of the Python pipeline
│   └── index.js                 # React entry point
├── public/
│   └── index.html               # Animated splash screen, floating base particles
├── tests/
│   └── test_assembler.py        # 53-test pytest suite
├── data/
│   └── example_reads.fa         # Demo input reads
├── api.py                       # Flask REST API
├── conftest.py                  # PYTHONPATH config for pytest
├── pytest.ini                   # Test discovery config
├── requirements.txt
└── README.md
```

---

## 📥 Example

**Input reads:**
```
ACGTTGCATGC
TGCATGCAAC
```

**k = 4, assembled contig:**
```
ACGTTGCATGCATGCAAC  (18 bp, 1 contig, 100% coverage)
```

---

## 🚀 Running Locally

**Terminal 1 — Flask API:**
```bash
pip install -r requirements.txt
PYTHONPATH=core python api.py        # Linux/Mac
$env:PYTHONPATH="core"; python api.py  # Windows PowerShell
```

**Terminal 2 — React frontend:**
```bash
npm install
npm start
```

Browser opens at `http://localhost:3000`. The React app proxies `/api/*` to Flask on port 5000 automatically.

**Run tests:**
```bash
PYTHONPATH=core python -m pytest tests/ -v        # Linux/Mac
$env:PYTHONPATH="core"; python -m pytest tests/ -v  # Windows PowerShell
```

---

## ✅ Test Coverage

53 tests across 8 classes:

| Class | What it covers |
|---|---|
| `TestDecomposeRead` | Single-read k-mer decomposition, edge counts, lowercase normalisation |
| `TestDecomposeReads` | Multi-read batching, validation errors, ReadSet structure |
| `TestDeBruijnGraphConstruction` | Node/edge/degree bookkeeping, adjacency correctness |
| `TestEulerianPath` | Hierholzer's algorithm, path length, circuit detection, start node selection |
| `TestPathToContig` | DNA string reconstruction from node paths |
| `TestGenomeAssembler` | End-to-end assembly, N50, multi-contig sorting, FASTA loading |
| `TestWriteFasta` | FASTA file output, headers, custom prefixes |
| `TestEdgeCases` | Repeated reads, k=2, single k-mer, empty result handling |

---

## 🗺️ Roadmap

- [ ] DOT/Graphviz export of De Bruijn graph
- [ ] Bubble detection and resolution (heterozygous variants)
- [ ] Tip removal (dead-end branches from sequencing errors)
- [ ] Paired-end read support
- [ ] Benchmarking against real datasets (E. coli, phage genomes)
- [ ] NumPy sparse matrix representation for large genomes

---

## 🔗 Related Project

This assembler is the natural successor to [Sequence Alignment Studio](https://github.com/etamara6/sequence-alignment-tool) — once reads are assembled into contigs, alignment algorithms (NW/SW) compare them to reference genomes.

---

## 📄 License

MIT © [etamara6](https://github.com/etamara6)

---

## 📄 License

MIT © [etamara6](https://github.com/etamara6)
