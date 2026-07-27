#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", "/tmp/heiconnect-matplotlib")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ALGORITHM_COLORS = {
    "block_tree_nored": "#4E79A7",
    "block_tree_red": "#F28E2B",
    "double_nored": "#59A14F",
    "double_red": "#E15759",
    "gwc": "#B07AA1",
}

ALGORITHM_MARKERS = {
    "block_tree_nored": "o",
    "block_tree_red": "s",
    "double_nored": "^",
    "double_red": "D",
    "gwc": "X",
}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create Seaborn visualizations from a collected results CSV."
    )
    parser.add_argument("results_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    return parser.parse_args()


def save_figure(path):
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()
    print(f"Saved {path}")


def plot_runtime_scaling(data, output_dir):
    plot_data = data.dropna(subset=["nodes", "runtime_seconds"])
    if plot_data.empty:
        return

    if plot_data["dataset"].nunique() == 1:
        plt.figure(figsize=(11, 6.5))
        axis = sns.lineplot(
            data=plot_data,
            x="nodes",
            y="runtime_seconds",
            hue="algorithm",
            style="link_configuration",
            markers=True,
            dashes=False,
            estimator="median",
            errorbar=None,
        )
        if (plot_data["runtime_seconds"] > 0).all():
            axis.set_yscale("log")
        axis.set(
            xlabel="Number of nodes",
            ylabel="Runtime (seconds)",
            title=f"Runtime scaling on {plot_data['dataset'].iloc[0]}",
        )
        axis.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
        save_figure(output_dir / "runtime_scaling.png")
        return

    grid = sns.relplot(
        data=plot_data,
        x="nodes",
        y="runtime_seconds",
        hue="algorithm",
        style="link_configuration",
        col="dataset",
        col_wrap=2,
        kind="line",
        markers=True,
        dashes=False,
        estimator="median",
        errorbar=None,
        facet_kws={"sharex": False, "sharey": False},
        height=4.2,
        aspect=1.25,
    )
    grid.set_axis_labels("Number of nodes", "Runtime (seconds)")
    grid.set_titles("{col_name}")
    grid.figure.suptitle("Runtime scaling by dataset", y=1.03, fontsize=15)
    for axis in grid.axes.flat:
        if (plot_data["runtime_seconds"] > 0).all():
            axis.set_yscale("log")
        axis.grid(True, alpha=0.2)
    grid.figure.savefig(
        output_dir / "runtime_scaling.png", dpi=220, bbox_inches="tight"
    )
    plt.close(grid.figure)


def plot_runtime_distribution(data, output_dir):
    plot_data = data.dropna(subset=["runtime_seconds"])
    if plot_data.empty:
        return

    plt.figure(figsize=(12, 6))
    sns.boxplot(
        data=plot_data,
        x="algorithm",
        y="runtime_seconds",
        hue="link_configuration",
        showfliers=False,
    )
    sns.stripplot(
        data=plot_data,
        x="algorithm",
        y="runtime_seconds",
        hue="link_configuration",
        dodge=True,
        alpha=0.45,
        size=4,
        legend=False,
    )
    if (plot_data["runtime_seconds"] > 0).all():
        plt.yscale("log")
    plt.xlabel("Algorithm configuration")
    plt.ylabel("Runtime (seconds, log scale)")
    plt.title("Runtime distribution across instances")
    plt.xticks(rotation=25, ha="right")
    plt.legend(
        title="Link configuration", bbox_to_anchor=(1.02, 1), loc="upper left"
    )
    save_figure(output_dir / "runtime_distribution.png")


def plot_solution_quality(data, output_dir):
    plot_data = data.dropna(subset=["solution_cost"]).copy()
    if plot_data.empty:
        return

    groups = ["dataset", "instance", "link_configuration"]
    best_cost = plot_data.groupby(groups)["solution_cost"].transform("min")
    plot_data["cost_over_best"] = plot_data["solution_cost"] / best_cost

    plt.figure(figsize=(11, 6))
    sns.barplot(
        data=plot_data,
        x="algorithm",
        y="cost_over_best",
        hue="link_configuration",
        estimator="mean",
        errorbar=("ci", 95),
    )
    plt.axhline(1.0, color="black", linestyle="--", linewidth=1, alpha=0.7)
    plt.xlabel("Algorithm configuration")
    plt.ylabel("Solution cost / best cost")
    plt.title("Relative solution quality (lower is better)")
    plt.xticks(rotation=25, ha="right")
    plt.legend(
        title="Link configuration", bbox_to_anchor=(1.02, 1), loc="upper left"
    )
    save_figure(output_dir / "relative_solution_quality.png")


def plot_runtime_memory(data, output_dir):
    plot_data = data.dropna(
        subset=["runtime_seconds", "peak_memory_mb", "nodes"]
    )
    if plot_data.empty:
        return

    plt.figure(figsize=(10, 7))
    sns.scatterplot(
        data=plot_data,
        x="runtime_seconds",
        y="peak_memory_mb",
        hue="algorithm",
        style="link_configuration",
        size="nodes",
        sizes=(50, 350),
        alpha=0.8,
    )
    if (plot_data["runtime_seconds"] > 0).all():
        plt.xscale("log")
    if (plot_data["peak_memory_mb"] > 0).all():
        plt.yscale("log")
    plt.xlabel("Runtime (seconds, log scale)")
    plt.ylabel("Peak memory (MiB, log scale)")
    plt.title("Runtime–memory trade-off")
    plt.grid(True, alpha=0.2)
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(output_dir / "runtime_memory_tradeoff.png")


def plot_algorithm_scaling(data, output_dir):
    plot_data = data.dropna(subset=["nodes", "runtime_seconds"])
    for algorithm, algorithm_data in plot_data.groupby("algorithm"):
        if algorithm_data["nodes"].nunique() < 2:
            continue

        plt.figure(figsize=(10.5, 6.5))
        axis = sns.lineplot(
            data=algorithm_data,
            x="nodes",
            y="runtime_seconds",
            hue="dataset",
            style="link_configuration",
            markers=True,
            dashes=True,
            estimator="median",
            errorbar=None,
            linewidth=2,
        )
        if (algorithm_data["nodes"] > 0).all():
            axis.set_xscale("log")
        if (algorithm_data["runtime_seconds"] > 0).all():
            axis.set_yscale("log")
        axis.set(
            xlabel="Number of nodes (log scale)",
            ylabel="Runtime (seconds, log scale)",
            title=f"Runtime scaling — {algorithm}",
        )
        axis.grid(True, which="both", alpha=0.2)
        axis.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", algorithm)
        save_figure(output_dir / f"scaling_{safe_name}.png")


def plot_dataset_scaling(data, output_dir):
    plot_data = data.dropna(subset=["nodes", "runtime_seconds"])
    for dataset, dataset_data in plot_data.groupby("dataset"):
        if dataset_data["nodes"].nunique() < 2:
            continue

        link_configurations = sorted(dataset_data["link_configuration"].unique())
        algorithms = sorted(dataset_data["algorithm"].unique())
        colors = {
            algorithm: ALGORITHM_COLORS.get(
                algorithm, sns.color_palette("husl", len(algorithms))[index]
            )
            for index, algorithm in enumerate(algorithms)
        }
        markers = {
            algorithm: ALGORITHM_MARKERS.get(algorithm, "o")
            for algorithm in algorithms
        }
        figure, axes = plt.subplots(
            1,
            len(link_configurations),
            figsize=(7 * len(link_configurations), 6),
            sharex=True,
            sharey=True,
            squeeze=False,
        )
        for index, link_configuration in enumerate(link_configurations):
            axis = axes[0, index]
            link_data = dataset_data[
                dataset_data["link_configuration"] == link_configuration
            ]
            sns.lineplot(
                data=link_data,
                x="nodes",
                y="runtime_seconds",
                hue="algorithm",
                style="algorithm",
                palette=colors,
                markers=markers,
                dashes=False,
                estimator="median",
                errorbar=None,
                linewidth=2.2,
                markersize=8,
                ax=axis,
            )
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set(
                xlabel="Number of nodes (log scale)",
                ylabel="Runtime (seconds, log scale)" if index == 0 else "",
                title=link_configuration.replace("_", " "),
            )
            axis.grid(True, which="both", alpha=0.2)
            if index == 0:
                axis.legend(
                    title="Algorithm",
                    bbox_to_anchor=(0.5, -0.2),
                    loc="upper center",
                    ncol=3,
                )
            else:
                axis.get_legend().remove()

        figure.suptitle(f"Algorithm runtime comparison — {dataset}", fontsize=15)
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", dataset)
        save_figure(output_dir / f"algorithm_comparison_{safe_name}.png")


def main():
    arguments = parse_arguments()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook", palette="colorblind")

    data = pd.read_csv(arguments.results_csv)
    print(f"Loaded {len(data)} runs from {arguments.results_csv}")

    plot_runtime_scaling(data, arguments.output_dir)
    plot_runtime_distribution(data, arguments.output_dir)
    plot_solution_quality(data, arguments.output_dir)
    plot_runtime_memory(data, arguments.output_dir)
    plot_algorithm_scaling(data, arguments.output_dir)
    plot_dataset_scaling(data, arguments.output_dir)


if __name__ == "__main__":
    main()
