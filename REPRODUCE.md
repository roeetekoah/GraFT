# Reproduce the Saved Analysis or Run Fresh Data

This document uses only the two module entrypoints used in the codebase:
- `graft.graph_construction.*` (including the graph-construction orchestrator)
- `graft.pipeline`

No `reproduce.py` wrapper is required. For comparison with the saved results,
start with the tracked graph below. The fresh-data workflow is a separate route
and is not an exact replay of the historical graph construction.

## Environment

From repository root (PowerShell):

```powershell
python -m pip install -e .
```

The `graph_bw` regression suite was checked on Windows with Python 3.13.13,
NetworkX 3.6.1, Matplotlib 3.10.9, and pandas 2.3.3 (clean installation)
and 3.0.2. This verifies the saved-graph analysis and small construction
fixtures; it does not cover a fresh full-FASTA reconstruction.

## A) Start With the Tracked Analysis Graph

This route needs no RefSeq download. From the repository root:

```powershell
python -m graft.pipeline `
  --in_edges golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv `
  --out_dir tmp/e2e_review/run/pipeline_no_bw `
  --no_betweenness

python -m graft.pipeline `
  --in_edges golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv `
  --out_dir tmp/e2e_review/run/pipeline_bw
```

Both commands use the same 92,790-edge input. Their rankings differ because
betweenness is a scoring feature; compare each run only with its matching
`golden/no_bw_pipeline` or `golden/bw_pipeline` baseline. Continue to the reporting
commands below, or use `python tests/test_regression_baselines.py --mode graph_bw`
to run the existing regression checks.

## B) Optional Fresh-Data Workflow

This route requires RefSeq metadata and FASTA downloads. A current assembly
summary may select different assemblies from the historical run. The current
orchestrator applies species-pair quantiles and top-X rows grouped by the stored
`u` endpoint; this is not a strict undirected degree cap. Use route A when
comparing with the saved analysis of the 92,790-edge graph.

Canonical manifest is tracked in the repo at:
- `data/out_refseq/manifest.tsv`

Download the assembly summary first if needed (see the command at the end).

1. Build a manifest and download its protein FASTAs:

```powershell
python -m graft.graph_construction.refseq_fetch_proteins `
  --assembly_summary data/assembly_summary_refseq.txt `
  --species_list config/species.txt `
  --out_dir tmp/e2e_review/graph_construction `
  --max_assemblies_per_species 2 `
  --require_latest `
  --download
```

2. Construct graph inputs with the graph-construction orchestrator:

```powershell
python -m graft.graph_construction.orchestrator construct-edges `
  --manifest tmp/e2e_review/graph_construction/manifest.tsv `
  --downloads_dir tmp/e2e_review/graph_construction/downloads `
  --out_candidates tmp/e2e_review/run/candidates.tsv `
  --out_edges tmp/e2e_review/run/edges_pruned.tsv `
  --k 6 --min_len 50 --max_postings 100 --min_shared 6 --top_m 10 --q 0.9 --top_x 20
```

3. Run the analysis pipeline (no betweenness, then betweenness):

```powershell
python -m graft.pipeline `
  --in_edges tmp/e2e_review/run/edges_pruned.tsv `
  --out_dir tmp/e2e_review/run/pipeline_no_bw `
  --no_betweenness

python -m graft.pipeline `
  --in_edges tmp/e2e_review/run/edges_pruned.tsv `
  --out_dir tmp/e2e_review/run/pipeline_bw
```

## C) Reporting / Explanations

The component IDs below refer to the tracked graph in route A. A newly
constructed graph may have different components and identifiers.

No-betweenness reports:

```powershell
python tools/reporting/top_anomaly_edges.py `
  --edges tmp/e2e_review/run/pipeline_no_bw/edge_features.tsv `
  --top_n 25 `
  --out_dir tmp/e2e_review/run/pipeline_no_bw/results

python tools/reporting/summarize_global_stats.py `
  --component_features tmp/e2e_review/run/pipeline_no_bw/component_features.tsv `
  --protein_features tmp/e2e_review/run/pipeline_no_bw/protein_features.tsv `
  --edge_features tmp/e2e_review/run/pipeline_no_bw/edge_features.tsv `
  --hgt_candidates tmp/e2e_review/run/pipeline_no_bw/hgt_candidates.tsv `
  --out_prefix tmp/e2e_review/run/pipeline_no_bw/results/global_stats
```

Betweenness-on reports + explanations:

```powershell
python tools/reporting/top_anomaly_edges.py `
  --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv `
  --top_n 25 `
  --out_dir tmp/e2e_review/run/pipeline_bw/results

python tools/reporting/summarize_global_stats.py `
  --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv `
  --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv `
  --edge_features tmp/e2e_review/run/pipeline_bw/edge_features.tsv `
  --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv `
  --out_prefix tmp/e2e_review/run/pipeline_bw/results/global_stats
```

```powershell
python tools/reporting/explain_component.py --component_id 5 --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools/reporting/explain_component.py --component_id 8 --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools/reporting/explain_component.py --component_id 32 --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools/reporting/explain_top_candidates.py --edges tmp/e2e_review/run/pipeline_bw/edge_features.tsv --protein_features tmp/e2e_review/run/pipeline_bw/protein_features.tsv --component_features tmp/e2e_review/run/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e_review/run/pipeline_bw/hgt_candidates.tsv --top_n 20 --top_k_neighbors 12
```

## D) Reference Artifacts and Verification Scope

Tracked input and comparison targets for route A:
- `golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv`
- `golden/no_bw_pipeline/rerun_pruned/*.tsv`
- `golden/no_bw_pipeline/results/*.tsv`
- `golden/bw_pipeline/rerun_pruned/*.tsv`
- `golden/bw_pipeline/results/*.tsv`
- `golden/bw_pipeline/explanations/*.txt`

The existing suite checks table values and selected explanation text. Compare
each analysis mode against its corresponding reference artifacts. These checks
verify implementation consistency, not biological confirmation of HGT.

Note:
- `data/assembly_summary_refseq.txt` is large and not tracked by git.
- Download from NCBI when needed:

```powershell
curl.exe -L -o data/assembly_summary_refseq.txt https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
```
