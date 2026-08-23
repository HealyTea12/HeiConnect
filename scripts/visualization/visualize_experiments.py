#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "matplotlib",
#     "pandas",
#     "seaborn",
# ]
# ///

import argparse
import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/heiconnect-matplotlib")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


EXACT_ALGORITHMS = {
    "block_tree_ilp_nored",
    "block_tree_ilp_red",
    "eilp",
    "scilp_nored",
    "scilp_red",
}
HEURISTIC_ALGORITHMS = {
    "block_tree_nored",
    "gwc",
    "sc_double_nored",
}

MARKER_OPTIONS = ["o", "s", "^", "D", "v", "<", ">", "p", "P", "X", "h", "H", "8", "*"]
LINE_STYLES = ["-", "--", "-.", ":"]
MARKER_CACHE_PATH = Path(__file__).with_name("algorithm_markers.json")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create runtime, memory, and solution quality scaling visualizations from a results CSV."
    )
    parser.add_argument("results_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    return parser.parse_args()


def _safe_filename(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def _load_marker_cache(path):
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def _save_marker_cache(path, markers):
    path.write_text(json.dumps(markers, indent=2, sort_keys=True))


def _next_marker(used_markers):
    used = set(used_markers)
    for marker in MARKER_OPTIONS:
        if marker not in used:
            return marker
    return MARKER_OPTIONS[len(used) % len(MARKER_OPTIONS)]


def _get_marker_map(algorithms, cache_path):
    marker_map = _load_marker_cache(cache_path)
    changed = False

    for algorithm in sorted(set(algorithms)):
        if algorithm not in marker_map:
            marker_map[algorithm] = _next_marker(marker_map.values())
            changed = True

    if changed:
        _save_marker_cache(cache_path, marker_map)

    return marker_map


def save_figure(path):
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()
    print(f"Saved {path}")


def plot_scaling(data, output_dir, value_column, ylabel, filename_prefix):
    plot_data = data.dropna(
        subset=["dataset", "link_configuration", "algorithm", "nodes", value_column]
    ).copy()
    plot_data["nodes"] = pd.to_numeric(plot_data["nodes"], errors="coerce")
    plot_data[value_column] = pd.to_numeric(plot_data[value_column], errors="coerce")
    plot_data = plot_data.dropna(subset=["nodes", value_column])
    if plot_data.empty:
        return

    marker_map = _get_marker_map(plot_data["algorithm"].unique(), MARKER_CACHE_PATH)
    algorithms = sorted(plot_data["algorithm"].unique())
    palette = {
        algorithm: sns.color_palette("colorblind", len(algorithms))[index]
        for index, algorithm in enumerate(algorithms)
    }

    for (dataset, link_distribution), group in plot_data.groupby(
        ["dataset", "link_configuration"]
    ):
        plt.figure(figsize=(10.5, 6.5))
        axis = plt.gca()
        group_algorithms = sorted(group["algorithm"].unique())
        for idx, algorithm in enumerate(group_algorithms):
            alg_data = group[group["algorithm"] == algorithm].copy()
            points = alg_data.groupby("nodes", as_index=False)[value_column].median()
            points = points.sort_values("nodes")
            style = LINE_STYLES[idx % len(LINE_STYLES)]
            axis.plot(
                points["nodes"],
                points[value_column],
                marker=marker_map.get(algorithm, "o"),
                linestyle=style,
                linewidth=2.1,
                markersize=8,
                markeredgewidth=1.0,
                markeredgecolor=palette[algorithm],
                markerfacecolor="white",
                label=algorithm,
                color=palette[algorithm],
                alpha=0.95,
            )
            axis.scatter(
                points["nodes"],
                points[value_column],
                color=palette[algorithm],
                marker=marker_map.get(algorithm, "o"),
                s=24,
                alpha=0.45,
            )

        print(
            f"{dataset}/{link_distribution} algorithms: "
            + ", ".join(group_algorithms)
        )
        if (group[value_column] > 0).all():
            axis.set_yscale("log")
        axis.set(
            xlabel="Number of nodes",
            ylabel=ylabel,
            title=f"{dataset} — {link_distribution}",
        )
        axis.legend(title="Algorithm", bbox_to_anchor=(1.02, 1), loc="upper left")
        axis.grid(True, alpha=0.2)

        filename = f"{filename_prefix}_{_safe_filename(dataset)}_{_safe_filename(link_distribution)}.png"
        save_figure(output_dir / filename)


def prepare_solution_quality(data):
    required_columns = [
        "dataset",
        "link_configuration",
        "instance",
        "algorithm",
        "nodes",
        "solution_cost",
    ]
    if not set(required_columns).issubset(data.columns):
        print("Missing columns required for solution quality plots")
        return pd.DataFrame()

    quality_data = data.dropna(subset=required_columns).copy()
    quality_data["nodes"] = pd.to_numeric(quality_data["nodes"], errors="coerce")
    quality_data["solution_cost"] = pd.to_numeric(
        quality_data["solution_cost"], errors="coerce"
    )
    quality_data = quality_data.dropna(subset=["nodes", "solution_cost"])

    known_algorithms = EXACT_ALGORITHMS | HEURISTIC_ALGORITHMS
    unknown_algorithms = sorted(set(quality_data["algorithm"]) - known_algorithms)
    if unknown_algorithms:
        print("Ignoring unclassified algorithms: " + ", ".join(unknown_algorithms))
    quality_data = quality_data[quality_data["algorithm"].isin(known_algorithms)]

    instance_columns = ["dataset", "link_configuration", "instance"]
    best_observed_costs = (
        quality_data.groupby(instance_columns, as_index=False)["solution_cost"]
        .min()
        .rename(columns={"solution_cost": "best_observed_cost"})
    )
    exact_results = quality_data[quality_data["algorithm"].isin(EXACT_ALGORITHMS)]
    reference_costs = (
        exact_results.groupby(instance_columns, as_index=False)["solution_cost"]
        .min()
        .rename(columns={"solution_cost": "exact_reference_cost"})
    )
    quality_data = quality_data.merge(best_observed_costs, on=instance_columns, how="left")
    quality_data = quality_data.merge(reference_costs, on=instance_columns, how="left")
    quality_data["reference_cost"] = quality_data["exact_reference_cost"].fillna(
        quality_data["best_observed_cost"]
    )
    quality_data["reference_is_exact"] = quality_data["exact_reference_cost"].notna()
    quality_data = quality_data[quality_data["reference_cost"] > 0].copy()
    quality_data["relative_gap_percent"] = 100 * (
        quality_data["solution_cost"] / quality_data["reference_cost"] - 1
    )
    return quality_data


def plot_solution_quality_gap(quality_data, output_dir):
    if quality_data.empty:
        return

    marker_map = _get_marker_map(quality_data["algorithm"].unique(), MARKER_CACHE_PATH)
    algorithms = sorted(quality_data["algorithm"].unique())
    palette = {
        algorithm: sns.color_palette("colorblind", len(algorithms))[index]
        for index, algorithm in enumerate(algorithms)
    }

    for (dataset, link_distribution), group in quality_data.groupby(
        ["dataset", "link_configuration"]
    ):
        plt.figure(figsize=(10.5, 6.5))
        axis = plt.gca()
        group_algorithms = sorted(group["algorithm"].unique())
        for idx, algorithm in enumerate(group_algorithms):
            alg_data = group[group["algorithm"] == algorithm]
            points = (
                alg_data.groupby("nodes")["relative_gap_percent"]
                .agg(
                    median="median",
                    q1=lambda values: values.quantile(0.25),
                    q3=lambda values: values.quantile(0.75),
                )
                .reset_index()
                .sort_values("nodes")
            )
            style = LINE_STYLES[idx % len(LINE_STYLES)]
            is_exact = algorithm in EXACT_ALGORITHMS
            axis.fill_between(
                points["nodes"].to_numpy(),
                points["q1"].to_numpy(),
                points["q3"].to_numpy(),
                color=palette[algorithm],
                alpha=0.08 if is_exact else 0.16,
            )
            axis.plot(
                points["nodes"],
                points["median"],
                marker=marker_map.get(algorithm, "o"),
                linestyle=style,
                linewidth=1.4 if is_exact else 2.1,
                markersize=7,
                markerfacecolor="white",
                label=algorithm,
                color=palette[algorithm],
                alpha=0.55 if is_exact else 0.95,
            )

        axis.axhline(0, color="black", linewidth=1.0, alpha=0.5)
        axis.set(
            xlabel="Number of nodes",
            ylabel="Gap to exact/best-known cost (%)",
            title=f"{dataset} — {link_distribution} — relative solution quality",
        )
        axis.legend(title="Algorithm", bbox_to_anchor=(1.02, 1), loc="upper left")
        axis.grid(True, alpha=0.2)

        filename = f"solution_quality_gap_{_safe_filename(dataset)}_{_safe_filename(link_distribution)}.png"
        save_figure(output_dir / filename)


def plot_solution_quality_ecdf(quality_data, output_dir):
    if quality_data.empty:
        return

    algorithms = sorted(quality_data["algorithm"].unique())
    palette = {
        algorithm: sns.color_palette("colorblind", len(algorithms))[index]
        for index, algorithm in enumerate(algorithms)
    }

    for (dataset, link_distribution), group in quality_data.groupby(
        ["dataset", "link_configuration"]
    ):
        plt.figure(figsize=(10.5, 6.5))
        axis = plt.gca()
        group_algorithms = sorted(group["algorithm"].unique())
        for idx, algorithm in enumerate(group_algorithms):
            gaps = (
                group[group["algorithm"] == algorithm]["relative_gap_percent"]
                .sort_values()
                .to_numpy()
            )
            fractions = [(index + 1) / len(gaps) for index in range(len(gaps))]
            is_exact = algorithm in EXACT_ALGORITHMS
            axis.step(
                gaps,
                fractions,
                where="post",
                linestyle=LINE_STYLES[idx % len(LINE_STYLES)],
                linewidth=1.4 if is_exact else 2.1,
                label=f"{algorithm} (n={len(gaps)})",
                color=palette[algorithm],
                alpha=0.55 if is_exact else 0.95,
            )

        axis.set(
            xlabel="Gap to exact/best-known cost (%)",
            ylabel="Fraction of instances",
            title=f"{dataset} — {link_distribution} — solution quality ECDF",
            ylim=(0, 1.02),
        )
        axis.set_xlim(left=0)
        axis.legend(title="Algorithm", bbox_to_anchor=(1.02, 1), loc="upper left")
        axis.grid(True, alpha=0.2)

        filename = f"solution_quality_ecdf_{_safe_filename(dataset)}_{_safe_filename(link_distribution)}.png"
        save_figure(output_dir / filename)


def main():
    arguments = parse_arguments()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook", palette="colorblind")

    data = pd.read_csv(arguments.results_csv)
    print(f"Loaded {len(data)} runs from {arguments.results_csv}")

    plot_scaling(data, arguments.output_dir, "runtime_seconds", "Runtime (seconds)", "runtime")
    plot_scaling(data, arguments.output_dir, "peak_memory_mb", "Peak memory (MiB)", "memory")
    quality_data = prepare_solution_quality(data)
    plot_solution_quality_gap(quality_data, arguments.output_dir)
    plot_solution_quality_ecdf(quality_data, arguments.output_dir)


if __name__ == "__main__":
    main()
