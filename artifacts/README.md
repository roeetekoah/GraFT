# Outputs and Component Figures

[Project overview](../README.md) · [Reproduction guide](../REPRODUCE.md) · [Method](../METHOD.md)

## Pipeline outputs

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

## Component figures

The tracked previews show components 5, 8, and 32 of the saved analysis graph.
Node color indicates species, node size reflects candidate score, and black
edges show the top 25 edges by robust species-pair z-score. Other edges remain
visible as context. These are candidate patterns for biological investigation,
not confirmed transfer events.

- [Component 5](updated_plots/topk/component_5.png)
- [Component 8](updated_plots/topk/component_8.png)
- [Component 32](updated_plots/topk/component_32.png)

The README uses equal-sized previews under `updated_plots/previews/`. Each
complete figure is scaled to fit the same white canvas, without stretching or
cropping, and links to its unchanged full-resolution image under `topk/`.

Recreate the views from the repository root:

```powershell
python tools\reporting\plot_components.py `
  --edges golden\bw_pipeline\rerun_pruned\edge_features.tsv `
  --protein_features golden\bw_pipeline\rerun_pruned\protein_features.tsv `
  --hgt_candidates golden\bw_pipeline\rerun_pruned\hgt_candidates.tsv `
  --component_ids 5,8,32 `
  --z_min_highlight -999 `
  --max_highlight_edges 25 `
  --highlight_policy_label "Top-25 edges by z_robust" `
  --node_size_mode score `
  --out_dir artifacts\updated_plots\topk
```

The commands for summary tables and candidate/component explanations are in
the [reproduction guide](../REPRODUCE.md#c-reporting--explanations).
Component identifiers refer to the saved graph; fresh graphs may have different
components and identifiers.

## Historical locations

Canonical comparison inputs and outputs live in [`golden/`](../golden/).
The ignored `legacy_runs/` and `historical_results/` directories hold local,
historical outputs when present; they are not required for the saved-graph route.
