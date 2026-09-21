import argparse
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
REPORTING_DIR = REPO_ROOT / "tools" / "reporting"
STAGES_DIR = SRC_DIR / "graft" / "stages"
GRAPH_CONSTRUCTION_DIR = SRC_DIR / "graft" / "graph_construction"
GOLDEN_DIR = REPO_ROOT / "golden"
RUN_GRAPH_PIPELINE_REGRESSION = False
RUN_GRAPH_PIPELINE_BW_REGRESSION = False
RUN_KMER_FULL_REGRESSION = False


def assert_pipeline_table_equal(actual, expected, name):
    """Keep identifiers and scoring exact; allow one serialized diagnostic unit."""
    diagnostic = "top5_mean_jac"
    if name == "protein_features.tsv":
        # Different floating-point summation methods can round a mean to
        # adjacent six-decimal values. No other field gets this allowance.
        pd.testing.assert_frame_equal(
            actual.drop(columns=[diagnostic]),
            expected.drop(columns=[diagnostic]),
            check_dtype=False,
            check_exact=True,
        )
        pd.testing.assert_series_equal(
            actual[diagnostic], expected[diagnostic], check_dtype=False,
            check_exact=False, rtol=0, atol=1e-6 + 1e-12,
        )
    else:
        pd.testing.assert_frame_equal(
            actual, expected, check_dtype=False, check_exact=True,
        )


def load_module(module_path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_python(args: List[str], cwd: Path, timeout_sec: int = 3600, env: Optional[Dict[str, str]] = None):
    cmd = [sys.executable, *args]
    completed = subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env,
        text=True,
        capture_output=True,
        timeout=timeout_sec,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"Command failed ({completed.returncode}): {' '.join(cmd)}\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )


def run_python_capture(args: List[str], cwd: Path, timeout_sec: int = 3600) -> str:
    cmd = [sys.executable, *args]
    completed = subprocess.run(
        cmd,
        cwd=str(cwd),
        text=True,
        capture_output=True,
        timeout=timeout_sec,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"Command failed ({completed.returncode}): {' '.join(cmd)}\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    return completed.stdout.replace("\r\n", "\n")


class RegressionBaselines(unittest.TestCase):
    def _assert_graph_regression(self, include_betweenness: bool, golden_pipeline_dir: Path):
        with tempfile.TemporaryDirectory() as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)
            out_dir = tmp_dir / "rerun_pruned"
            out_dir.mkdir(parents=True, exist_ok=True)

            cmd = [
                "-m",
                "graft.pipeline",
                "--in_edges",
                str(GOLDEN_DIR / "reference_inputs" / "edges_PRUNED_JACCARD_92790.tsv"),
                "--out_dir",
                str(out_dir),
            ]
            if not include_betweenness:
                cmd.append("--no_betweenness")

            run_python(
                cmd,
                cwd=SRC_DIR,
                timeout_sec=7200,
            )

            for name in [
                "edge_features.tsv",
                "protein_features.tsv",
                "component_features.tsv",
                "hgt_candidates.tsv",
                "all_scores.tsv",
            ]:
                actual = pd.read_csv(out_dir / name, sep="\t")
                expected = pd.read_csv(golden_pipeline_dir / "rerun_pruned" / name, sep="\t")
                assert_pipeline_table_equal(actual, expected, name)

            results_dir = tmp_dir / "results"
            results_dir.mkdir(parents=True, exist_ok=True)

            run_python(
                [
                    str(REPORTING_DIR / "top_anomaly_edges.py"),
                    "--edges",
                    str(out_dir / "edge_features.tsv"),
                    "--top_n",
                    "25",
                    "--out_dir",
                    str(results_dir),
                ],
                cwd=REPO_ROOT,
                timeout_sec=900,
            )

            run_python(
                [
                    str(REPORTING_DIR / "summarize_global_stats.py"),
                    "--component_features",
                    str(out_dir / "component_features.tsv"),
                    "--protein_features",
                    str(out_dir / "protein_features.tsv"),
                    "--edge_features",
                    str(out_dir / "edge_features.tsv"),
                    "--hgt_candidates",
                    str(out_dir / "hgt_candidates.tsv"),
                    "--out_prefix",
                    str(results_dir / "global_stats"),
                ],
                cwd=SRC_DIR,
                timeout_sec=900,
            )

            top_actual = pd.read_csv(results_dir / "top_anomaly_edges.tsv", sep="\t")
            top_expected = pd.read_csv(golden_pipeline_dir / "results" / "top_anomaly_edges.tsv", sep="\t")
            pd.testing.assert_frame_equal(top_actual, top_expected, check_dtype=False)

            stats_actual = pd.read_csv(results_dir / "global_stats.tsv", sep="\t")
            stats_expected = pd.read_csv(golden_pipeline_dir / "results" / "global_stats.tsv", sep="\t")
            pd.testing.assert_frame_equal(stats_actual, stats_expected, check_dtype=False)

    def test_fasta_parsing_smoke_against_repo_data(self):
        fasta_parsing = load_module(GRAPH_CONSTRUCTION_DIR / "fasta_parsing.py", "fasta_parsing")
        manifest = GOLDEN_DIR / "reference_inputs" / "manifest_tiny_set.tsv"
        species_map = fasta_parsing.load_manifest_species_map(manifest)

        self.assertIn("GCF_000005845.2_ASM584v2", species_map)
        self.assertEqual(species_map["GCF_000005845.2_ASM584v2"], "Escherichia coli")

        faa_path = GOLDEN_DIR / "reference_inputs" / "downloads_tiny" / "GCF_000005845.2_ASM584v2_protein.faa.gz"
        first_header, first_seq = next(fasta_parsing.iter_fasta_gz(faa_path))
        self.assertTrue(first_header)
        self.assertGreater(len(first_seq), 0)
        self.assertTrue(fasta_parsing.protein_id_from_header(first_header))

    def test_kmer_candidates_tiny_regression_matches_baseline(self):
        with tempfile.TemporaryDirectory() as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)
            out_candidates = tmp_dir / "candidates_tiny.tsv"

            run_python(
                [
                    "-m",
                    "graft.graph_construction.orchestrator",
                    "build-candidates",
                    "--manifest",
                    str(GOLDEN_DIR / "reference_inputs" / "manifest_tiny_set.tsv"),
                    "--downloads_dir",
                    str(GOLDEN_DIR / "reference_inputs" / "downloads_tiny"),
                    "--out_candidates",
                    str(out_candidates),
                    "--k",
                    "6",
                    "--min_len",
                    "50",
                    "--max_postings",
                    "100",
                    "--min_shared",
                    "6",
                    "--top_m",
                    "10",
                ],
                cwd=SRC_DIR,
                timeout_sec=900,
            )

            actual = pd.read_csv(out_candidates, sep="\t")
            expected = pd.read_csv(GOLDEN_DIR / "reference_inputs" / "candidates_tiny_set.tsv", sep="\t")
            pd.testing.assert_frame_equal(actual, expected, check_dtype=False)

    def test_pruning_logic_matches_edges_pruned_baseline(self):
        graph_pruning = load_module(GRAPH_CONSTRUCTION_DIR / "graph_pruning.py", "graph_pruning")

        candidates_path = GOLDEN_DIR / "reference_inputs" / "candidates_tiny_set.tsv"
        baseline_path = GOLDEN_DIR / "reference_inputs" / "edges_pruned_tiny.tsv"

        df = pd.read_csv(candidates_path, sep="\t")
        filtered_df = graph_pruning.keep_q_percentile_edges(df, q=0.9)

        # Baseline preserves species_pair string representation.
        filtered_df = filtered_df.copy()
        filtered_df["species_pair"] = filtered_df.apply(
            lambda r: str(tuple(sorted([r["species_u"], r["species_v"]]))),
            axis=1,
        )
        top_edges = graph_pruning.keep_top_X_edges_per_node(filtered_df, X=20)

        expected = pd.read_csv(baseline_path, sep="\t")
        cols = ["u", "v", "shared_kmers", "jaccard", "species_u", "species_v", "species_pair"]

        top_edges = top_edges[cols].sort_values(cols).reset_index(drop=True)
        expected = expected[cols].sort_values(cols).reset_index(drop=True)

        self.assertEqual(len(top_edges), len(expected))
        pd.testing.assert_frame_equal(top_edges, expected, check_like=False, check_dtype=False)

    def test_quantile_parameter_and_unordered_species_pairs(self):
        pruning = load_module(GRAPH_CONSTRUCTION_DIR / "graph_pruning.py", "graph_pruning")
        edges = pd.DataFrame({
            "u": ["a", "b", "c", "d", "e", "f"],
            "v": ["v1", "v2", "v3", "v4", "v5", "v6"],
            "species_u": ["A", "B", "A", "B", "C", "C"],
            "species_v": ["B", "A", "B", "A", "D", "D"],
            "jaccard": [0.1, 0.2, 0.3, 0.4, 0.8, 0.9],
        })
        original = edges.copy(deep=True)
        self.assertEqual(pruning.keep_q_percentile_edges(edges, q=0.5).u.tolist(), ["c", "d", "f"])
        self.assertEqual(pruning.keep_q_percentile_edges(edges, q=0.9).u.tolist(), ["d", "f"])
        pd.testing.assert_frame_equal(edges, original)

    def test_endpoint_pruning_preserves_columns_order_and_ties(self):
        pruning = load_module(GRAPH_CONSTRUCTION_DIR / "graph_pruning.py", "graph_pruning")
        edges = pd.DataFrame({
            "u": ["b", "a", "a", "a", "b"],
            "v": ["x", "y", "z", "w", "q"],
            "jaccard": [0.2, 0.8, 0.8, 0.1, 0.9],
            "annotation": ["b-low", "first-tie", "second-tie", "a-low", "b-high"],
        }, index=[9, 9, 4, 7, 4])
        original = edges.copy(deep=True)
        expected = edges.iloc[[1, 4]].reset_index(drop=True)
        pd.testing.assert_frame_equal(pruning.keep_top_X_edges_per_node(edges, X=1), expected)
        pd.testing.assert_frame_equal(pruning.keep_top_X_edges_per_node(edges.iloc[:0], X=1), edges.iloc[:0].reset_index(drop=True))
        pd.testing.assert_frame_equal(edges, original)

    def test_numeric_allowance_is_limited_to_one_diagnostic_unit(self):
        expected = pd.DataFrame({"u": ["p"], "top5_mean_jac": [0.095102], "max_z": [3.0]})
        actual = expected.copy()
        actual.loc[0, "top5_mean_jac"] = 0.095103
        assert_pipeline_table_equal(actual, expected, "protein_features.tsv")
        actual.loc[0, "top5_mean_jac"] = 0.095104
        with self.assertRaises(AssertionError):
            assert_pipeline_table_equal(actual, expected, "protein_features.tsv")
        actual = expected.copy()
        actual.loc[0, "max_z"] += 1e-6
        with self.assertRaises(AssertionError):
            assert_pipeline_table_equal(actual, expected, "protein_features.tsv")
        scores = pd.DataFrame({"u": ["p"], "score": [3.0]})
        changed = scores.copy()
        changed.loc[0, "score"] += 1e-6
        with self.assertRaises(AssertionError):
            assert_pipeline_table_equal(changed, scores, "all_scores.tsv")

    def test_graph_pipeline_matches_pruned_and_results_baselines(self):
        if not RUN_GRAPH_PIPELINE_REGRESSION:
            self.skipTest("Graph pipeline regression disabled for this run.")
        self._assert_graph_regression(
            include_betweenness=False,
            golden_pipeline_dir=GOLDEN_DIR / "no_bw_pipeline",
        )

    def test_graph_pipeline_with_betweenness_matches_baseline(self):
        if not RUN_GRAPH_PIPELINE_BW_REGRESSION:
            self.skipTest("Graph pipeline (with betweenness) regression disabled for this run.")
        self._assert_graph_regression(
            include_betweenness=True,
            golden_pipeline_dir=GOLDEN_DIR / "bw_pipeline",
        )

    def test_kmer_candidates_full_regression_matches_baseline(self):
        if not RUN_KMER_FULL_REGRESSION:
            self.skipTest("Full k-mer regression disabled for this run.")
        if not (REPO_ROOT / "data" / "out_refseq" / "downloads").exists():
            self.skipTest("Full downloads set not available locally.")
        with tempfile.TemporaryDirectory() as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)

            out_candidates = tmp_dir / "candidates.tsv"
            run_python(
                [
                    "-m",
                    "graft.graph_construction.orchestrator",
                    "build-candidates",
                    "--manifest",
                    str(REPO_ROOT / "data" / "out_refseq" / "manifest.tsv"),
                    "--downloads_dir",
                    str(REPO_ROOT / "data" / "out_refseq" / "downloads"),
                    "--out_candidates",
                    str(out_candidates),
                    "--k",
                    "6",
                    "--min_len",
                    "50",
                    "--max_postings",
                    "100",
                    "--min_shared",
                    "6",
                    "--top_m",
                    "10",
                ],
                cwd=SRC_DIR,
                timeout_sec=14400,
            )

            actual = pd.read_csv(out_candidates, sep="\t")
            expected = pd.read_csv(REPO_ROOT / "data" / "out_refseq" / "candidates.tsv", sep="\t")
            pd.testing.assert_frame_equal(actual, expected, check_dtype=False)

    def test_bw_explanations_match_golden(self):
        if not RUN_GRAPH_PIPELINE_BW_REGRESSION:
            self.skipTest("Explanation regression disabled for this run.")

        explain_dir = GOLDEN_DIR / "bw_pipeline" / "explanations"
        rerun_dir = GOLDEN_DIR / "bw_pipeline" / "rerun_pruned"

        component_runs = [
            (5, "comp_5_explained.txt"),
            (8, "comp_8_explained.txt"),
            (32, "comp_32_explained.txt"),
        ]
        for component_id, out_name in component_runs:
            stdout = run_python_capture(
                [
                    str(REPORTING_DIR / "explain_component.py"),
                    "--component_id",
                    str(component_id),
                    "--edges",
                    str(rerun_dir / "edge_features.tsv"),
                    "--protein_features",
                    str(rerun_dir / "protein_features.tsv"),
                    "--component_features",
                    str(rerun_dir / "component_features.tsv"),
                    "--hgt_candidates",
                    str(rerun_dir / "hgt_candidates.tsv"),
                    "--top_nodes",
                    "20",
                    "--top_edges",
                    "25",
                ],
                cwd=REPO_ROOT,
                timeout_sec=900,
            )
            expected = (explain_dir / out_name).read_text(encoding="utf-8").replace("\r\n", "\n")
            self.assertEqual(stdout, expected, f"Mismatch in {out_name}")

        top_stdout = run_python_capture(
            [
                str(REPORTING_DIR / "explain_top_candidates.py"),
                "--edges",
                str(rerun_dir / "edge_features.tsv"),
                "--protein_features",
                str(rerun_dir / "protein_features.tsv"),
                "--component_features",
                str(rerun_dir / "component_features.tsv"),
                "--hgt_candidates",
                str(rerun_dir / "hgt_candidates.tsv"),
                "--top_n",
                "20",
                "--top_k_neighbors",
                "12",
            ],
            cwd=REPO_ROOT,
            timeout_sec=900,
        )
        top_expected = (explain_dir / "top_candidates_explained.txt").read_text(encoding="utf-8").replace("\r\n", "\n")
        self.assertEqual(top_stdout, top_expected, "Mismatch in top_candidates_explained.txt")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument(
        "--mode",
        choices=["fast", "graph", "graph_bw", "full"],
        default="graph",
        help="fast: quick checks only, graph: no-betweenness graph regression, graph_bw: include both no-bw and bw graph regressions, full: all regressions",
    )
    known, remaining = parser.parse_known_args()

    RUN_GRAPH_PIPELINE_REGRESSION = known.mode in {"graph", "graph_bw", "full"}
    RUN_GRAPH_PIPELINE_BW_REGRESSION = known.mode in {"graph_bw", "full"}
    RUN_KMER_FULL_REGRESSION = known.mode == "full"

    unittest.main(argv=[sys.argv[0], *remaining])
