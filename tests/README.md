# Regression Tests

Run from the repository root after installation:

```sh
python tests/test_regression_baselines.py --mode graph
```

| Mode | Coverage |
|---|---|
| `fast` | Small graph-construction fixtures and reporting checks. |
| `graph` (default) | No-betweenness pipeline regression and fast checks. |
| `graph_bw` | Also includes betweenness-on regression. |
| `full` | Also includes the heavier full k-mer regression, when local data exists. |

Compare each analysis mode with its matching baseline under [`golden/`](../golden/).
Regression checks verify implementation consistency, not biological confirmation
of HGT. See the [reproduction guide](../REPRODUCE.md) for the saved-graph workflow
and the recorded verification environment.
