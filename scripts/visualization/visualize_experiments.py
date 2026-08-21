#!/usr/bin/env python3

import argparse
import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/heiconnect-matplotlib")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


MARKER_OPTIONS = ["o", "s", "^", "D", "v", "<", ">", "p", "P", "X", "h", "H", "8", "*"]
LINE_STYLES = ["-", "--", "-.", ":"]
MARKER_CACHE_PATH = Path(__file__).with_name("algorithm_markers.json")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create compact runtime scaling visualizations from a results CSV."
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


def plot_runtime_scaling(data, output_dir):
    plot_data = data.dropna(
        subset=["dataset", "link_configuration", "algorithm", "nodes", "runtime_seconds"]
    )
    plot_data["nodes"] = pd.to_numeric(plot_data["nodes"], errors="coerce")
    plot_data["runtime_seconds"] = pd.to_numeric(
        plot_data["runtime_seconds"], errors="coerce"
    )
    plot_data = plot_data.dropna(subset=["nodes", "runtime_seconds"])
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
            points = alg_data.groupby("nodes", as_index=False)["runtime_seconds"].median()
            points = points.sort_values("nodes")
            style = LINE_STYLES[idx % len(LINE_STYLES)]
            axis.plot(
                points["nodes"],
                points["runtime_seconds"],
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
                points["runtime_seconds"],
                color=palette[algorithm],
                marker=marker_map.get(algorithm, "o"),
                s=24,
                alpha=0.45,
            )

        print(
            f"{dataset}/{link_distribution} algorithms: "
            + ", ".join(group_algorithms)
        )
        if (group["runtime_seconds"] > 0).all():
            axis.set_yscale("log")
        axis.set(
            xlabel="Number of nodes",
            ylabel="Runtime (seconds)",
            title=f"{dataset} — {link_distribution}",
        )
        axis.legend(title="Algorithm", bbox_to_anchor=(1.02, 1), loc="upper left")
        axis.grid(True, alpha=0.2)

        filename = f"runtime_{_safe_filename(dataset)}_{_safe_filename(link_distribution)}.png"
        save_figure(output_dir / filename)


def main():
    arguments = parse_arguments()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook", palette="colorblind")

    data = pd.read_csv(arguments.results_csv)
    print(f"Loaded {len(data)} runs from {arguments.results_csv}")

    plot_runtime_scaling(data, arguments.output_dir)


if __name__ == "__main__":
    main()
