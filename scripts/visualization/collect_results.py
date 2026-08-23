#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "pandas",
# ]
# ///

import argparse
from pathlib import Path
import re

import pandas as pd


RESULT_PATTERN = re.compile(r"res-(.+)\.txt")
PAGE_SIZE_BYTES = 4096

def parse_path_metadata(parts):
    if len(parts) >= 4:
        return parts[-4], parts[-3], parts[-2]

    if len(parts) == 3:
        return parts[0], "unknown", parts[1]

    if len(parts) == 2:
        return "unknown", "unknown", parts[0]

    return "unknown", "unknown", "unknown"


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Collect experiment result files into one tidy CSV."
    )
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("output_csv", type=Path)
    return parser.parse_args()


def parse_value(value):
    if value == "true":
        return True
    if value == "false":
        return False
    try:
        return float(value)
    except ValueError:
        return value


def parse_result_file(path):
    runs = []
    current = {}
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key == "run.instance" and current:
            runs.append(current)
            current = {}
        current[key] = parse_value(value)
    if current:
        runs.append(current)
    return runs


def load_results(results_dir):
    rows = []
    for result_file in sorted(results_dir.rglob("res-*.txt")):
        match = RESULT_PATTERN.fullmatch(result_file.name)
        if not match:
            continue

        relative_path = result_file.relative_to(results_dir)
        algorithm_configuration, dataset, instance = parse_path_metadata(
            relative_path.parts
        )
        link_configuration = match.group(1)

        for run in parse_result_file(result_file):
            source_instance = Path(str(run.get("run.instance", "")))
            if dataset == "unknown" and source_instance.parent.name:
                dataset = source_instance.parent.name
            node_match = re.search(r"\d+", instance)
            node_count = run.get("input_graph.n")
            if node_count is None and node_match:
                node_count = float(node_match.group())

            memory_pages = run.get("run.peak_memory_pages")
            rows.append(
                {
                    "algorithm": algorithm_configuration,
                    "link_configuration": link_configuration,
                    "dataset": dataset,
                    "instance": instance,
                    "nodes": node_count,
                    "edges": run.get("input_graph.m"),
                    "links": run.get("input_graph.num_links"),
                    "runtime_seconds": run.get("run.total_time_seconds"),
                    "solution_cost": run.get("solution.cost"),
                    "solution_size": run.get("solution.size"),
                    "peak_memory_mb": (
                        memory_pages * PAGE_SIZE_BYTES / (1024 * 1024)
                        if isinstance(memory_pages, float)
                        else None
                    ),
                }
            )

    data = pd.DataFrame(rows)
    if data.empty:
        raise ValueError(f"no named result files found below {results_dir}")
    return data


def main():
    arguments = parse_arguments()
    data = load_results(arguments.results_dir)
    arguments.output_csv.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(arguments.output_csv, index=False)
    print(f"Collected {len(data)} runs into {arguments.output_csv}")


if __name__ == "__main__":
    main()
