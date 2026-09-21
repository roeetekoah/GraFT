<table>
  <tr>
    <td width="220" align="center">
      <img src="assets/logo-readme-lockup-card.png" alt="GraFT logo" width="190">
    </td>
    <td>
      <h1>GraFT: Graph-based Finder of Transfers</h1>
      <p><strong>Alignment-Free Detection of Horizontal Gene Transfer Candidates via Protein Similarity Graphs and Anomaly Scoring</strong></p>
    </td>
  </tr>
</table>

**GraFT** stands for **Gra**ph-based **F**inder of **T**ransfers. This repository
contains the GraFT pipeline: an alignment-free k-mer Jaccard workflow for
detecting horizontal gene transfer (HGT) candidates across bacterial proteomes
without whole-genome alignment or phylogenetic reconstruction.

GraFT asks which proteins are unusually similar across species relative to each
species pair's own background, and whether those similarities form informative
graph structures. Its output is a ranked set of candidates for biological
follow-up; graph anomalies alone do not establish a transfer event.

**Explore:** [read the paper](paper/main.pdf), [run the saved graph](#getting-started),
or [inspect the component figures](#results-and-figures).

<p align="center">
  <img src="https://img.shields.io/badge/Python-%3E%3D3.9-3776AB?logo=python&logoColor=white" alt="python >= 3.9">
  <img src="https://img.shields.io/badge/NetworkX-graph%20analysis-0A66C2" alt="networkx">
  <img src="https://img.shields.io/badge/Pandas-tabular%20outputs-150458?logo=pandas" alt="pandas">
  <img src="https://img.shields.io/badge/Matplotlib-report%20plots-11557C" alt="matplotlib">
  <img src="https://img.shields.io/badge/NCBI%20RefSeq-proteomes-4CAF50" alt="NCBI RefSeq">
</p>

The active manuscript is self-contained under [`paper/`](paper/), including
its figures in [`paper/figures/`](paper/figures/).  
Read the full paper at [`paper/main.pdf`](paper/main.pdf).

---

## Table of Contents

- [Overview](#overview)
- [What This Repository Contains](#what-this-repository-contains)
- [Method](#method)
- [Project Structure](#project-structure)
- [Getting Started](#getting-started)
- [GraFT Outputs](#graft-outputs)
- [Results and Figures](#results-and-figures)
- [Reproducibility](#reproducibility)
- [Testing](#testing)
- [Data Policy](#data-policy)
- [Citation and License](#citation-and-license)

---

## Overview

Horizontal gene transfer can leave a local signature: a protein may look
suspiciously similar to proteins in phylogenetically distant species, even when
the surrounding genomes do not share the same history. This repository models
that problem as a sparse cross-species protein similarity graph.

GraFT treats proteins as graph nodes. Edges connect proteins with high
alignment-free similarity, measured by shared k-mers and Jaccard score. The
analysis then normalizes each edge against the background distribution for its
species pair and ranks proteins/components by graph anomaly features.

The practical goal is not to prove HGT alone. The goal is to produce a compact,
ranked set of HGT candidates that are worth biological validation.

---

## What This Repository Contains

This repository contains **GraFT** and the supporting reproducibility artifacts:

| Area | Purpose |
|---|---|
| `graft.graph_construction` | Builds candidate protein-pair edges from RefSeq FASTA files and prunes them into a graph input. |
| `graft.pipeline` | Computes edge, protein, and component features, then ranks HGT candidates. |
| `tools/reporting` | Produces summary tables, anomaly plots, component explanations, and candidate explanations. |
| `golden` | Canonical inputs and regression baselines for deterministic checks. |
| `paper` | Active manuscript source, generated PDF, and manuscript figures. |

Two entry paths are supported:

| Path | Start From | Best For |
|---|---|---|
| Shortcut | `golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv` | Fast reproduction of the analysis pipeline. |
| Fresh-data run | `data/assembly_summary_refseq.txt` + `config/species.txt` | Building a new graph from RefSeq metadata and downloaded proteomes; not an exact replay of the historical graph. |

---

## Method

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

---

## Project Structure

```text
.
|-- config/
|   `-- species.txt
|-- data/
|   `-- out_refseq/                 # local RefSeq-derived outputs when present
|-- paper/
|   `-- figures/
|-- golden/
|   |-- reference_inputs/
|   |-- bw_pipeline/
|   |-- no_bw_pipeline/
|   `-- hackathon_report_refs/
|-- src/
|   `-- graft/
|       |-- graph_construction/
|       |-- pipeline.py
|       `-- stages/
|-- tests/
|-- tools/
|   |-- reporting/
|   |-- REPRODUCE.md
|   `-- reproduce.py
|-- REPRODUCE.md
|-- pyproject.toml
`-- README.md
```

| Path | Notes |
|---|---|
| `src/graft/graph_construction/refseq_fetch_proteins.py` | Builds a manifest from NCBI assembly metadata and optionally downloads protein FASTAs. |
| `src/graft/graph_construction/orchestrator.py` | Builds candidate edges and applies graph pruning in one CLI. |
| `src/graft/pipeline.py` | Main analysis entry point for feature extraction and HGT candidate ranking. |
| `src/graft/stages/` | Modular stages for I/O, graph operations, pair statistics, node features, component features, and ranking. |
| `tools/reporting/` | Post-run scripts used for the report figures and explanation files. |
| `golden/reference_inputs/` | Small and canonical graph inputs used for quick runs and tests. |

---

## Getting Started

### Prerequisites

- Python 3.9 or newer
- `pip`
- `curl.exe` or another downloader if you want to fetch the full RefSeq assembly summary

Install the package from the repository root:

```powershell
python -m pip install -e .
```

This installs the editable distribution `graft-hgt`. The importable Python
package is `graft`, and the installed command-line entry points are `graft`,
`graft-graph`, and `graft-fetch-refseq`.

### Quickstart: Run GraFT on the Canonical Graph

This uses the preincluded pruned graph input and skips the expensive graph
construction stage.

With betweenness centrality:

```powershell
python -m graft.pipeline `
  --in_edges golden\reference_inputs\edges_PRUNED_JACCARD_92790.tsv `
  --out_dir tmp_run_bw
```

Without betweenness, which is faster:

```powershell
python -m graft.pipeline `
  --in_edges golden\reference_inputs\edges_PRUNED_JACCARD_92790.tsv `
  --out_dir tmp_run_no_bw `
  --no_betweenness
```

### Generate Basic Reports

```powershell
python tools\reporting\top_anomaly_edges.py `
  --edges tmp_run_bw\edge_features.tsv `
  --protein_features tmp_run_bw\protein_features.tsv `
  --top_n 25 `
  --out_dir tmp_run_bw\results

python tools\reporting\summarize_global_stats.py `
  --component_features tmp_run_bw\component_features.tsv `
  --protein_features tmp_run_bw\protein_features.tsv `
  --edge_features tmp_run_bw\edge_features.tsv `
  --hgt_candidates tmp_run_bw\hgt_candidates.tsv `
  --out_prefix tmp_run_bw\results\global_stats
```

---

## GraFT Outputs

Each pipeline run writes TSV deliverables to the chosen `--out_dir`:

| File | Description |
|---|---|
| `hgt_candidates.tsv` | Top-ranked HGT candidate proteins with final scores. |
| `all_scores.tsv` | HGT-likeness score for every scored protein. |
| `edge_features.tsv` | Per-edge Jaccard/shared-k-mer values and robust species-pair z-scores. |
| `protein_features.tsv` | Per-protein graph features, component membership, and local anomaly context. |
| `component_features.tsv` | Per-component size, concentration, and high-z edge summaries. |

Reporting scripts can additionally create:

| File | Description |
|---|---|
| `results/top_anomaly_edges.tsv` | Highest-z edges in ranked tabular form. |
| `results/top_anomaly_edges.png` | Plot of top anomalous edges. |
| `results/global_stats.tsv` | Compact global summary for reproducibility checks. |
| `results/global_stats.json` | Machine-readable version of the global summary. |
| `explanations/*.txt` | Human-readable component and candidate explanations. |

---

## Results and Figures

The canonical analysis highlights graph components where cross-species protein
similarity is concentrated in unusually strong edges. Components 5, 8, and 32
are shown below as compact previews. Click any component to open the
high-resolution render with legible species colors, candidate-node sizing, and
highlighted anomalous edges.

| Component 5 | Component 8 | Component 32 |
|---|---|---|
| [![Component 5 high-resolution graph](artifacts/updated_plots/topk/component_5.png)](artifacts/updated_plots/topk/component_5.png) | [![Component 8 high-resolution graph](artifacts/updated_plots/topk/component_8.png)](artifacts/updated_plots/topk/component_8.png) | [![Component 32 high-resolution graph](artifacts/updated_plots/topk/component_32.png)](artifacts/updated_plots/topk/component_32.png) |

Recreate these plots:

```powershell
python tools\reporting\plot_components.py `
  --edges golden\bw_pipeline\rerun_pruned\edge_features.tsv `
  --protein_features golden\bw_pipeline\rerun_pruned\protein_features.tsv `
  --hgt_candidates golden\bw_pipeline\rerun_pruned\hgt_candidates.tsv `
  --component_ids 5,8,32 `
  --z_min_highlight -999 `
  --max_highlight_edges 60 `
  --highlight_policy_label "Top-60 edges by z_robust" `
  --node_size_mode score `
  --out_dir artifacts\updated_plots\topk
```

---

## Reproducibility

For a concise reviewer-facing recipe, see [`REPRODUCE.md`](REPRODUCE.md).
Start with the tracked graph to reproduce the saved analysis. A fresh-data run
can select different RefSeq assemblies and produce a different graph.

<details>
<summary>Fresh-data workflow from RefSeq assembly metadata</summary>

Install once:

```powershell
python -m pip install -e .
```

Download the NCBI RefSeq assembly summary if it is not already available:

```powershell
curl.exe -L -o data\assembly_summary_refseq.txt https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
```

Build a manifest and download protein FASTAs:

```powershell
python -m graft.graph_construction.refseq_fetch_proteins `
  --assembly_summary data\assembly_summary_refseq.txt `
  --species_list config\species.txt `
  --out_dir tmp\e2e\graph_construction `
  --max_assemblies_per_species 2 `
  --require_latest `
  --download
```

Construct candidates and pruned graph edges:

```powershell
python -m graft.graph_construction.orchestrator construct-edges `
  --manifest tmp\e2e\graph_construction\manifest.tsv `
  --downloads_dir tmp\e2e\graph_construction\downloads `
  --out_candidates tmp\e2e\candidates.tsv `
  --out_edges tmp\e2e\edges_pruned.tsv `
  --k 6 --min_len 50 --max_postings 100 --min_shared 6 --top_m 10 --q 0.9 --top_x 20
```

Run the GraFT analysis:

```powershell
python -m graft.pipeline `
  --in_edges tmp\e2e\edges_pruned.tsv `
  --out_dir tmp\e2e\pipeline_bw

python -m graft.pipeline `
  --in_edges tmp\e2e\edges_pruned.tsv `
  --out_dir tmp\e2e\pipeline_no_bw `
  --no_betweenness
```

Generate reports:

```powershell
python tools\reporting\top_anomaly_edges.py `
  --edges tmp\e2e\pipeline_bw\edge_features.tsv `
  --protein_features tmp\e2e\pipeline_bw\protein_features.tsv `
  --top_n 25 `
  --out_dir tmp\e2e\pipeline_bw\results

python tools\reporting\summarize_global_stats.py `
  --component_features tmp\e2e\pipeline_bw\component_features.tsv `
  --protein_features tmp\e2e\pipeline_bw\protein_features.tsv `
  --edge_features tmp\e2e\pipeline_bw\edge_features.tsv `
  --hgt_candidates tmp\e2e\pipeline_bw\hgt_candidates.tsv `
  --out_prefix tmp\e2e\pipeline_bw\results\global_stats
```

</details>

<details>
<summary>Optional explanation commands</summary>

```powershell
python tools\reporting\explain_component.py --component_id 5 --edges tmp\e2e\pipeline_bw\edge_features.tsv --protein_features tmp\e2e\pipeline_bw\protein_features.tsv --component_features tmp\e2e\pipeline_bw\component_features.tsv --hgt_candidates tmp\e2e\pipeline_bw\hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools\reporting\explain_component.py --component_id 8 --edges tmp\e2e\pipeline_bw\edge_features.tsv --protein_features tmp\e2e\pipeline_bw\protein_features.tsv --component_features tmp\e2e\pipeline_bw\component_features.tsv --hgt_candidates tmp\e2e\pipeline_bw\hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools\reporting\explain_component.py --component_id 32 --edges tmp\e2e\pipeline_bw\edge_features.tsv --protein_features tmp\e2e\pipeline_bw\protein_features.tsv --component_features tmp\e2e\pipeline_bw\component_features.tsv --hgt_candidates tmp\e2e\pipeline_bw\hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools\reporting\explain_top_candidates.py --edges tmp\e2e\pipeline_bw\edge_features.tsv --protein_features tmp\e2e\pipeline_bw\protein_features.tsv --component_features tmp\e2e\pipeline_bw\component_features.tsv --hgt_candidates tmp\e2e\pipeline_bw\hgt_candidates.tsv --top_n 20 --top_k_neighbors 12
```

</details>

---

## Testing

Run the graph regression suite from the repository root:

```powershell
python tests\test_regression_baselines.py --mode graph
```

Other modes:

| Mode | Coverage |
|---|---|
| `fast` | Small graph-construction fixtures and reporting checks. |
| `graph` | Default no-betweenness regression and fast checks. |
| `graph_bw` | Includes betweenness-on regression. |
| `full` | Includes the heavier full k-mer regression when local data exists. |

---

## Data Policy

Large raw data is intentionally not tracked in git. In particular,
`data/assembly_summary_refseq.txt` is large, so fetch it from NCBI when needed:

```powershell
curl.exe -L -o data\assembly_summary_refseq.txt https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
```

Canonical regression inputs and outputs live in [`golden/`](golden/) so the
analysis can be checked without requiring every raw RefSeq download.

---

## Citation and License

Citation metadata is available in [`CITATION.cff`](CITATION.cff).

This project is licensed under the GNU General Public License v3.0. See
[`LICENSE`](LICENSE) for details.
