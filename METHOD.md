# Method and Implementation

[Project overview](README.md) · [Full paper](paper/main.pdf) · [Run the analysis](REPRODUCE.md)

## Pipeline

The saved analysis graph contains **37,421 proteins, 92,790 edges, and 5,262
connected components across 50 species**. These counts can be checked directly
against the tracked input and output tables under `golden/`. Graph construction
uses k-mer overlap and Jaccard similarity in place of pairwise protein alignment.

```text
Protein FASTAs from RefSeq
        |
        v
Step 1: FASTA parsing and k-mer extraction
        Build per-protein k-mer sets and an inverted index.
        |
        v
Step 2: Candidate edge generation
        Keep protein pairs with enough shared k-mers and top-M similarity.
        |
        v
Step 3: Graph pruning
        Apply species-pair Jaccard quantiles and top-X rows per source endpoint.
        |
        v
Step 4: Connected component analysis
        Treat proteins as nodes and retained similarities as edges.
        |
        v
Step 5: Species-pair robust normalization
        Score each edge relative to its species-pair background.
        |
        v
Step 6: Node and component features
        Compute betweenness, clustering, high-z concentration, and graph context.
        |
        v
Step 7: Candidate ranking
        Write ranked HGT candidates and supporting feature tables.
```

Key signal used by the ranking:

- **Edge surprise**: robust z-score of an edge within its species-pair background.
- **Node context**: whether a protein participates in unusually strong cross-species edges.
- **Graph position**: centrality and local topology inside the component.
- **Component concentration**: whether high-z edges cluster in a meaningful connected region.

The top-X pruning step groups rows by the stored source (`u`) endpoint. It is
not a strict undirected degree cap. The saved graph is the comparison input
for the tracked analysis; a fresh RefSeq download may produce a different graph.

## Code map

| Path | Notes |
|---|---|
| `src/graft/graph_construction/refseq_fetch_proteins.py` | Builds a manifest from NCBI assembly metadata and optionally downloads protein FASTAs. |
| `src/graft/graph_construction/orchestrator.py` | Builds candidate edges and applies graph pruning in one CLI. |
| `src/graft/pipeline.py` | Main analysis entry point for feature extraction and HGT candidate ranking. |
| `src/graft/stages/` | Modular stages for I/O, graph operations, pair statistics, node features, component features, and ranking. |
| `tools/reporting/` | Post-run scripts used for the report figures and explanation files. |
| `golden/reference_inputs/` | Small and canonical graph inputs used for quick runs and tests. |

The editable distribution is `graft-hgt`; the importable package is `graft`.
Installation also provides the `graft`, `graft-graph`, and `graft-fetch-refseq`
command-line entry points.

## Inputs and outputs

- [Saved input and regression baselines](golden/README.md)
- [Output tables, component figures, and plotting commands](artifacts/README.md)
- [Saved-graph and fresh-data workflows](REPRODUCE.md)
- [Regression test modes](tests/README.md)
