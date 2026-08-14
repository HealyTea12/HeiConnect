#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "tqdm",
# ]
# ///

import argparse
from collections import defaultdict
import os
from pathlib import Path
import re
import resource
import signal
import subprocess
import sys
import tempfile
import tomllib

from tqdm import tqdm


VALID_NAME = re.compile(r"[A-Za-z0-9_.-]+")
VALID_DISTRIBUTIONS = {"constant", "float_uniform", "integer_uniform"}


def positive_integer(value):
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ValueError("must be a positive integer")
    return value


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate links and run named configurations on every graph instance."
    )
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "-c", "--configurations", required=True, type=Path, metavar="FILE"
    )
    parser.add_argument(
        "--no-repeat",
        action="store_true",
        help="skip runs that already have a non-empty named result file",
    )
    parser.add_argument(
        "--mode",
        choices=("exhaustive", "adaptive"),
        default="exhaustive",
        help="run every instance or discover each algorithm's solvability frontier",
    )
    parser.add_argument(
        "--adaptive-probes",
        type=int,
        default=8,
        metavar="COUNT",
        help="number of widely spaced size levels to consider in adaptive mode",
    )
    return parser.parse_args()


def validate_name(name, kind):
    if not VALID_NAME.fullmatch(name):
        raise ValueError(
            f"{kind} name {name!r} may only contain letters, numbers, '.', '_', and '-'"
        )


def load_configurations(path):
    with path.open("rb") as file:
        document = tomllib.load(file)

    configurations = document.get("configurations")
    if not isinstance(configurations, dict) or not configurations:
        raise ValueError("the TOML file must contain [configurations.<name>] tables")

    for name, configuration in configurations.items():
        validate_name(name, "configuration")
        if not isinstance(configuration, dict):
            raise ValueError(f"configuration {name!r} must be a TOML table")
        if not isinstance(configuration.get("algorithm"), str):
            raise ValueError(f"configuration {name!r} requires an algorithm string")
        if configuration.get("timeout") is not None:
            positive_integer(configuration["timeout"])
        if configuration.get("max_memory_mb") is not None:
            positive_integer(configuration["max_memory_mb"])
        if not isinstance(configuration.get("params", {}), dict):
            raise ValueError(f"configuration {name!r} params must be a TOML table")

    link_configurations = document.get("link_configurations")
    if not isinstance(link_configurations, dict) or not link_configurations:
        raise ValueError(
            "the TOML file must contain [link_configurations.<name>] tables"
        )

    for name, configuration in link_configurations.items():
        validate_name(name, "link configuration")
        if not isinstance(configuration, dict):
            raise ValueError(f"link configuration {name!r} must be a TOML table")
        if configuration.get("distribution") not in VALID_DISTRIBUTIONS:
            raise ValueError(
                f"link configuration {name!r} distribution must be constant, "
                "float_uniform, or integer_uniform"
            )
        seed = configuration.get("seed")
        if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
            raise ValueError(
                f"link configuration {name!r} requires a non-negative integer seed"
            )

    dataset_selection = document.get("dataset_selection", {})
    if not isinstance(dataset_selection, dict):
        raise ValueError("dataset_selection must be a TOML table")
    if dataset_selection.get("default_max_nodes") is not None:
        positive_integer(dataset_selection["default_max_nodes"])
    included_datasets = dataset_selection.get("include", [])
    if not isinstance(included_datasets, list) or not all(
        isinstance(name, str) for name in included_datasets
    ):
        raise ValueError("dataset_selection.include must be an array of strings")
    dataset_limits = dataset_selection.get("datasets", {})
    if not isinstance(dataset_limits, dict):
        raise ValueError("dataset_selection.datasets must be a TOML table")
    for name, selection in dataset_limits.items():
        validate_name(name, "dataset")
        if not isinstance(selection, dict):
            raise ValueError(f"dataset selection {name!r} must be a TOML table")
        positive_integer(selection.get("max_nodes"))

    return configurations, link_configurations, dataset_selection


def read_node_count(graph_file):
    with graph_file.open() as file:
        for line in file:
            stripped = line.strip()
            if stripped and not stripped.startswith("%"):
                return int(stripped.split()[0])
    raise ValueError(f"could not read a node count from {graph_file}")


def find_instances(input_dir, dataset_selection):
    instances = []
    skipped = 0
    included_datasets = set(dataset_selection.get("include", []))
    dataset_limits = dataset_selection.get("datasets", {})
    default_max_nodes = dataset_selection.get("default_max_nodes")
    for graph_file in sorted(input_dir.rglob("*.xml")):
        metis_file = graph_file.with_suffix(".graph")
        if not metis_file.is_file():
            print(
                f"Skipping {graph_file}: matching .graph file not found",
                file=sys.stderr,
            )
            skipped += 1
            continue

        relative_path = graph_file.relative_to(input_dir)
        dataset_name = (
            relative_path.parts[0] if len(relative_path.parts) > 1 else input_dir.name
        )
        if included_datasets and dataset_name not in included_datasets:
            skipped += 1
            continue

        node_count = read_node_count(metis_file)
        max_nodes = dataset_limits.get(dataset_name, {}).get(
            "max_nodes", default_max_nodes
        )
        if max_nodes is not None and node_count > max_nodes:
            skipped += 1
            continue

        instances.append((graph_file, metis_file, dataset_name, node_count))
    return instances, skipped


def evenly_spaced_indices(size_levels, count):
    if count == 1:
        return [0]
    if len(size_levels) <= count:
        return list(range(len(size_levels)))

    smallest = size_levels[0]
    size_range = size_levels[-1] - smallest
    targets = [smallest + index * size_range / (count - 1) for index in range(count)]
    return sorted(
        {
            min(
                range(len(size_levels)),
                key=lambda index: abs(size_levels[index] - target),
            )
            for target in targets
        }
    )


def group_instances_by_dataset_and_size(instances):
    grouped = defaultdict(lambda: defaultdict(list))
    for graph_file, metis_file, dataset_name, node_count in instances:
        grouped[dataset_name][node_count].append((graph_file, metis_file))
    return grouped


def run_adaptive_series(size_levels, probe_count, run_level):
    outcomes = {}
    consecutive_timeout_levels = 0

    for index in evenly_spaced_indices(size_levels, probe_count):
        outcome = run_level(size_levels[index])
        outcomes[index] = outcome
        if outcome == "timeout":
            consecutive_timeout_levels += 1
            if consecutive_timeout_levels == 3:
                break
        else:
            consecutive_timeout_levels = 0

    solved_indices = [index for index, outcome in outcomes.items() if outcome == "solved"]
    last_solved = max(solved_indices, default=-1)
    timeout_indices = [
        index
        for index, outcome in outcomes.items()
        if index > last_solved and outcome == "timeout"
    ]
    boundary = min(timeout_indices, default=max(outcomes))

    for index in range(boundary + 1):
        if index not in outcomes:
            outcomes[index] = run_level(size_levels[index])

    return {size_levels[index]: outcome for index, outcome in outcomes.items()}


def serialize_parameter(value):
    if isinstance(value, bool):
        return str(value).lower()
    return str(value)


def write_configuration(
    path, name, configuration, link_configurations, dataset_selection
):
    lines = [
        f"name={name}",
        f"algorithm={configuration['algorithm']}",
        f"timeout={configuration.get('timeout', '')}",
        f"max_memory_mb={configuration.get('max_memory_mb', '')}",
    ]
    for key, value in configuration.get("params", {}).items():
        lines.append(f"param.{key}={serialize_parameter(value)}")
    for link_name, link_configuration in link_configurations.items():
        for key, value in link_configuration.items():
            lines.append(
                f"link_configuration.{link_name}.{key}={serialize_parameter(value)}"
            )
    for dataset_name in dataset_selection.get("include", []):
        lines.append(f"dataset_selection.include={dataset_name}")
    if dataset_selection.get("default_max_nodes") is not None:
        lines.append(
            "dataset_selection.default_max_nodes="
            f"{dataset_selection['default_max_nodes']}"
        )
    for dataset_name, selection in dataset_selection.get("datasets", {}).items():
        lines.append(
            f"dataset_selection.{dataset_name}.max_nodes={selection['max_nodes']}"
        )
    path.write_text("\n".join(lines) + "\n")


def memory_limit(max_memory_mb):
    if max_memory_mb is None:
        return None

    limit = positive_integer(max_memory_mb) * 1024 * 1024

    def set_limit():
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))

    return set_limit


def run_command(command, timeout=None, max_memory_mb=None):
    process = subprocess.Popen(
        command,
        preexec_fn=memory_limit(max_memory_mb),
        start_new_session=True,
    )
    try:
        return process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        raise


def link_generator_command(executable, graph_file, output_file, configuration):
    command = [
        executable,
        "--type",
        "links",
        "--generator",
        "complete",
        "--input_graph",
        graph_file,
        "--output",
        output_file,
        "--distribution",
        configuration["distribution"],
        "--seed",
        str(configuration["seed"]),
    ]
    distribution = configuration["distribution"]
    if distribution == "constant":
        command.extend(
            ["--constant_weight", str(configuration.get("constant_weight", 1.0))]
        )
    elif distribution == "float_uniform":
        command.extend(
            [
                "--float_uniform_lower",
                str(configuration.get("float_uniform_lower", 0.0)),
                "--float_uniform_upper",
                str(configuration.get("float_uniform_upper", 1.0)),
            ]
        )
    else:
        command.extend(
            [
                "--integer_uniform_lower",
                str(configuration.get("integer_uniform_lower", 1)),
                "--integer_uniform_upper",
                str(configuration.get("integer_uniform_upper", 9)),
            ]
        )
    return command


def run_experiment(
    executable,
    configuration_name,
    configuration,
    graph_file,
    links_file,
    instance_output,
    result_file_name,
):
    command = [
        executable,
        "-a",
        configuration["algorithm"],
        "-f",
        graph_file,
        links_file,
        "-o",
        instance_output,
        "--result_file_name",
        result_file_name,
    ]
    for key, value in configuration.get("params", {}).items():
        command.extend(["-p", f"{key}={serialize_parameter(value)}"])

    result_file = instance_output / result_file_name
    previous_size = result_file.stat().st_size if result_file.exists() else 0
    try:
        return_code = run_command(
            command,
            configuration.get("timeout"),
            configuration.get("max_memory_mb"),
        )
    except subprocess.TimeoutExpired:
        print(
            f"[{configuration_name}] Timed out after "
            f"{configuration['timeout']}s: {graph_file}",
            file=sys.stderr,
        )
        return "timeout"

    result_size = result_file.stat().st_size if result_file.exists() else 0
    if return_code == 0 and result_size > previous_size:
        return "completed"

    reason = f"exit code {return_code}" if return_code else "no result was written"
    print(
        f"[{configuration_name}] Failed ({reason}): {graph_file}",
        file=sys.stderr,
    )
    return "failed"


def main():
    arguments = parse_arguments()
    positive_integer(arguments.adaptive_probes)
    input_dir = arguments.input_dir.resolve()
    output_dir = arguments.output_dir.resolve()
    project_dir = Path(__file__).resolve().parents[2]
    experiments_executable = Path(
        os.environ.get(
            "EXPERIMENTS_BINARY",
            project_dir / "build" / "experiments" / "experiments",
        )
    )
    link_generator_executable = Path(
        os.environ.get(
            "DATASET_GENERATOR_BINARY",
            project_dir / "build" / "experiments" / "generate_datasets",
        )
    )

    if not input_dir.is_dir():
        raise FileNotFoundError(f"input directory not found: {input_dir}")
    for executable in (experiments_executable, link_generator_executable):
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise FileNotFoundError(
                f"required executable not found or not executable: {executable}"
            )

    configurations, link_configurations, dataset_selection = load_configurations(
        arguments.configurations
    )
    instances, filtered_instances = find_instances(input_dir, dataset_selection)
    if not instances:
        raise ValueError("no .xml files with matching .graph files were found")
    total_runs = len(instances) * len(link_configurations) * len(configurations)
    print(f"Selected instances: {len(instances)}")
    print(f"Filtered or invalid instances: {filtered_instances}")
    print(f"Total experiment runs: {total_runs}")
    output_dir.mkdir(parents=True, exist_ok=True)

    for name, configuration in configurations.items():
        configuration_dir = output_dir / name
        configuration_dir.mkdir(parents=True, exist_ok=True)
        write_configuration(
            configuration_dir / "configuration.txt",
            name,
            configuration,
            link_configurations,
            dataset_selection,
        )

    completed = 0
    failed = 0
    timed_out = 0
    skipped = 0
    processed = 0
    progress = tqdm(total=total_runs, desc="Running experiments", unit="run")
    adaptive_reports = defaultdict(list)

    def run_graph(graph_file, metis_file, selected_configurations):
        nonlocal completed, failed, timed_out, skipped, processed
        relative_path = graph_file.relative_to(input_dir).with_suffix("")
        outcomes = {name: [] for name in selected_configurations}
        for link_name, link_configuration in link_configurations.items():
            result_file_name = f"res-{link_name}.txt"
            pending = []
            for configuration_name in selected_configurations:
                configuration = configurations[configuration_name]
                instance_output = output_dir / configuration_name / relative_path
                result_file = instance_output / result_file_name
                if (
                    arguments.no_repeat
                    and result_file.is_file()
                    and result_file.stat().st_size > 0
                ):
                    print(
                        f"[{configuration_name}/{link_name}] "
                        f"Skipping completed {relative_path}"
                    )
                    skipped += 1
                    processed += 1
                    outcomes[configuration_name].append("completed")
                    progress.set_postfix_str(
                        f"skipped {configuration_name}/{link_name}/{relative_path}"
                    )
                    progress.update()
                else:
                    pending.append(
                        (configuration_name, configuration, instance_output)
                    )

            if not pending:
                continue

            with tempfile.TemporaryDirectory(prefix="heiconnect-links-") as temp_dir:
                links_file = Path(temp_dir) / f"{graph_file.stem}-{link_name}.links"
                print(f"[{link_name}] Generating links for {relative_path}")
                generator_result = run_command(
                    link_generator_command(
                        link_generator_executable,
                        metis_file,
                        links_file,
                        link_configuration,
                    )
                )
                if generator_result != 0:
                    print(
                        f"[{link_name}] Link generation failed for {graph_file}",
                        file=sys.stderr,
                    )
                    failed += len(pending)
                    processed += len(pending)
                    for configuration_name, _, _ in pending:
                        outcomes[configuration_name].append("failed")
                        progress.set_postfix_str(
                            f"failed {configuration_name}/{link_name}/{relative_path}"
                        )
                        progress.update()
                    continue

                for configuration_name, configuration, instance_output in pending:
                    instance_output.mkdir(parents=True, exist_ok=True)
                    print(
                        f"[{configuration_name}/{link_name}] "
                        f"Running {relative_path}"
                    )
                    outcome = run_experiment(
                        experiments_executable,
                        configuration_name,
                        configuration,
                        graph_file,
                        links_file,
                        instance_output,
                        result_file_name,
                    )
                    outcomes[configuration_name].append(outcome)
                    processed += 1
                    if outcome == "completed":
                        completed += 1
                        status = "completed"
                    elif outcome == "timeout":
                        timed_out += 1
                        status = "timed out"
                    else:
                        failed += 1
                        status = "failed"
                    progress.set_postfix_str(
                        f"{status} {configuration_name}/{link_name}/{relative_path}"
                    )
                    progress.update()
        return outcomes

    if arguments.mode == "exhaustive":
        for graph_file, metis_file, _, _ in instances:
            run_graph(graph_file, metis_file, configurations)
    else:
        grouped_instances = group_instances_by_dataset_and_size(instances)
        for dataset_name, instances_by_size in sorted(grouped_instances.items()):
            size_levels = sorted(instances_by_size)
            for configuration_name in configurations:
                print(f"[{configuration_name}] Discovering frontier for {dataset_name}")

                def run_level(node_count):
                    level_outcomes = []
                    print(
                        f"[{configuration_name}] Testing {dataset_name} "
                        f"at {node_count} nodes"
                    )
                    for graph_file, metis_file in instances_by_size[node_count]:
                        outcomes = run_graph(
                            graph_file, metis_file, [configuration_name]
                        )
                        level_outcomes.extend(outcomes[configuration_name])
                    if level_outcomes and all(
                        outcome == "timeout" for outcome in level_outcomes
                    ):
                        return "timeout"
                    if "completed" in level_outcomes:
                        return "solved"
                    return "failed"

                outcomes = run_adaptive_series(
                    size_levels, arguments.adaptive_probes, run_level
                )
                tested_sizes = ", ".join(
                    f"{size}={outcome}" for size, outcome in sorted(outcomes.items())
                )
                print(f"[{configuration_name}/{dataset_name}] {tested_sizes}")
                adaptive_reports[configuration_name].append(
                    f"dataset={dataset_name}"
                )
                for size in size_levels:
                    outcome = outcomes.get(size, "censored")
                    adaptive_reports[configuration_name].append(
                        f"size={size} status={outcome}"
                    )

        for configuration_name, report in adaptive_reports.items():
            report_path = output_dir / configuration_name / "adaptive-frontier.txt"
            report_path.write_text("\n".join(report) + "\n")

    censored = total_runs - processed
    if censored:
        progress.update(censored)
    progress.close()
    print(f"Completed: {completed}")
    print(f"Skipped: {skipped}")
    print(f"Timed out: {timed_out}")
    print(f"Failed: {failed}")
    print(f"Censored: {censored}")
    return 1 if failed or (arguments.mode == "exhaustive" and timed_out) else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, tomllib.TOMLDecodeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
