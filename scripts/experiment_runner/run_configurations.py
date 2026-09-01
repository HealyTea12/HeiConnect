#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "tqdm",
# ]
# ///

import argparse
from collections import defaultdict
from itertools import product
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
        help=(
            "run every instance or use bounded probes to discover each algorithm's "
            "solvability frontier"
        ),
    )
    parser.add_argument(
        "--adaptive-probes",
        type=int,
        default=6,
        metavar="COUNT",
        help="maximum frontier probes before adaptive mode fills lower levels",
    )
    return parser.parse_args()


def validate_name(name, kind):
    if not VALID_NAME.fullmatch(name):
        raise ValueError(
            f"{kind} name {name!r} may only contain letters, numbers, '.', '_', and '-'"
        )


def configuration_name_component(value):
    component = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).lower()).strip("_")
    if not component:
        raise ValueError(f"cannot use {value!r} in an expanded configuration name")
    return component


def expand_configuration_matrices(matrices):
    if not isinstance(matrices, dict):
        raise ValueError("configuration_matrices must be a TOML table")

    expanded = {}
    for matrix_name, matrix in matrices.items():
        validate_name(matrix_name, "configuration matrix")
        if not isinstance(matrix, dict):
            raise ValueError(f"configuration matrix {matrix_name!r} must be a table")

        allowed_keys = {"algorithms", "timeout", "max_memory_mb", "params"}
        unknown_keys = set(matrix) - allowed_keys
        if unknown_keys:
            raise ValueError(
                f"configuration matrix {matrix_name!r} has unknown keys: "
                f"{', '.join(sorted(unknown_keys))}"
            )

        algorithms = matrix.get("algorithms")
        if not isinstance(algorithms, list) or not algorithms or not all(
            isinstance(algorithm, str) for algorithm in algorithms
        ):
            raise ValueError(
                f"configuration matrix {matrix_name!r} algorithms must be "
                "a non-empty array of strings"
            )
        params = matrix.get("params", {})
        if not isinstance(params, dict):
            raise ValueError(
                f"configuration matrix {matrix_name!r} params must be a table"
            )

        fixed_params = {}
        dimensions = []
        for key, value in params.items():
            validate_name(key, "matrix parameter")
            if isinstance(value, list):
                if not value:
                    raise ValueError(
                        f"configuration matrix {matrix_name!r} parameter "
                        f"{key!r} must not be an empty array"
                    )
                if not all(
                    isinstance(item, (str, int, float, bool)) for item in value
                ):
                    raise ValueError(
                        f"configuration matrix {matrix_name!r} parameter "
                        f"{key!r} values must be strings, numbers, or booleans"
                    )
                dimensions.append((key, value))
            else:
                if not isinstance(value, (str, int, float, bool)):
                    raise ValueError(
                        f"configuration matrix {matrix_name!r} parameter "
                        f"{key!r} must be a string, number, boolean, or array"
                    )
                fixed_params[key] = value

        dimension_values = [values for _, values in dimensions]
        for algorithm in algorithms:
            algorithm_name = configuration_name_component(algorithm)
            combinations = product(*dimension_values) if dimensions else [()]
            for values in combinations:
                configuration_name_parts = [matrix_name, algorithm_name]
                configuration_params = dict(fixed_params)
                for (key, _), value in zip(dimensions, values):
                    configuration_params[key] = value
                    configuration_name_parts.extend(
                        [key, configuration_name_component(value)]
                    )

                configuration_name = "_".join(configuration_name_parts)
                if configuration_name in expanded:
                    raise ValueError(
                        f"expanded configuration name {configuration_name!r} "
                        "is not unique"
                    )
                configuration = {
                    "algorithm": algorithm,
                    "params": configuration_params,
                }
                for key in ("timeout", "max_memory_mb"):
                    if key in matrix:
                        configuration[key] = matrix[key]
                expanded[configuration_name] = configuration
    return expanded


def expand_link_configurations(link_configurations):
    expanded = {}
    for name, configuration in link_configurations.items():
        has_seed = "seed" in configuration
        has_seeds = "seeds" in configuration
        if has_seed == has_seeds:
            raise ValueError(
                f"link configuration {name!r} requires exactly one of seed or seeds"
            )

        seeds = [configuration["seed"]] if has_seed else configuration["seeds"]
        if not isinstance(seeds, list) or not seeds:
            raise ValueError(
                f"link configuration {name!r} seeds must be a non-empty array"
            )
        for seed in seeds:
            if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
                raise ValueError(
                    f"link configuration {name!r} seeds must be non-negative integers"
                )
        if len(seeds) != len(set(seeds)):
            raise ValueError(f"link configuration {name!r} seeds must be unique")

        for seed in seeds:
            expanded_name = name if has_seed else f"{name}_seed_{seed}"
            if expanded_name in expanded:
                raise ValueError(
                    f"expanded link configuration name {expanded_name!r} is not unique"
                )
            expanded_configuration = {
                key: value
                for key, value in configuration.items()
                if key not in {"seed", "seeds"}
            }
            expanded_configuration["seed"] = seed
            expanded_configuration["_family"] = name
            expanded[expanded_name] = expanded_configuration
    return expanded


def load_configurations(path):
    with path.open("rb") as file:
        document = tomllib.load(file)

    configurations = document.get("configurations", {})
    if not isinstance(configurations, dict):
        raise ValueError("configurations must be a TOML table")
    configurations = dict(configurations)
    matrix_configurations = expand_configuration_matrices(
        document.get("configuration_matrices", {})
    )
    duplicate_names = set(configurations) & set(matrix_configurations)
    if duplicate_names:
        raise ValueError(
            "configuration names are not unique: "
            + ", ".join(sorted(duplicate_names))
        )
    configurations.update(matrix_configurations)
    if not configurations:
        raise ValueError(
            "the TOML file must contain configurations or configuration_matrices"
        )

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
    link_configurations = expand_link_configurations(link_configurations)

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


def exponential_probe_indices(level_count, probe_count):
    indices = []
    index = 0
    step = 1
    while len(indices) < probe_count:
        indices.append(index)
        if index == level_count - 1:
            break
        if len(indices) == probe_count - 1:
            index = level_count - 1
        else:
            index = min(level_count - 1, index + step)
            step *= 2
    return indices


def group_instances_by_dataset_and_size(instances):
    grouped = defaultdict(lambda: defaultdict(list))
    for graph_file, metis_file, dataset_name, node_count in instances:
        grouped[dataset_name][node_count].append((graph_file, metis_file))
    return grouped


def group_link_configurations(link_configurations):
    grouped = defaultdict(list)
    for name, configuration in link_configurations.items():
        grouped[configuration["_family"]].append(name)
    return grouped


def run_adaptive_series(size_levels, probe_count, run_level):
    outcomes = {}
    last_solved_index = -1
    terminal_index = None
    terminal_outcome = None

    for index in exponential_probe_indices(len(size_levels), probe_count):
        outcome = run_level(size_levels[index])
        outcomes[index] = outcome
        if outcome != "solved":
            terminal_index = index
            terminal_outcome = outcome
            break
        last_solved_index = index

    fill_through = last_solved_index
    if terminal_outcome == "timeout":
        fill_through = terminal_index - 1

    for index in range(fill_through + 1):
        if index in outcomes:
            continue
        outcome = run_level(size_levels[index])
        outcomes[index] = outcome
        if outcome != "solved":
            break

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
            if key.startswith("_"):
                continue
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
        with result_file.open("rb") as file:
            file.seek(previous_size)
            if b"run.total_time_seconds=" in file.read():
                return "completed"

    reason = (
        f"exit code {return_code}"
        if return_code
        else "no complete result was written"
    )
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

    def run_graph(
        graph_file,
        metis_file,
        selected_configurations,
        selected_link_configurations=None,
        link_cache_dir=None,
    ):
        nonlocal completed, failed, timed_out, skipped, processed
        relative_path = graph_file.relative_to(input_dir).with_suffix("")
        outcomes = {name: [] for name in selected_configurations}
        selected_link_configurations = (
            link_configurations
            if selected_link_configurations is None
            else selected_link_configurations
        )
        for link_name in selected_link_configurations:
            link_configuration = link_configurations[link_name]
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

            temporary_links = None
            if link_cache_dir is None:
                temporary_links = tempfile.TemporaryDirectory(
                    prefix="heiconnect-links-"
                )
                links_dir = Path(temporary_links.name)
            else:
                links_dir = link_cache_dir / relative_path.parent
                links_dir.mkdir(parents=True, exist_ok=True)
            links_file = links_dir / f"{graph_file.stem}-{link_name}.links"
            if not links_file.is_file() or links_file.stat().st_size == 0:
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
                    links_file.unlink(missing_ok=True)
                    failed += len(pending)
                    processed += len(pending)
                    for configuration_name, _, _ in pending:
                        outcomes[configuration_name].append("failed")
                        progress.set_postfix_str(
                            f"failed {configuration_name}/{link_name}/{relative_path}"
                        )
                        progress.update()
                    if temporary_links is not None:
                        temporary_links.cleanup()
                    continue
            else:
                print(f"[{link_name}] Reusing links for {relative_path}")

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
            if temporary_links is not None:
                temporary_links.cleanup()
        return outcomes

    if arguments.mode == "exhaustive":
        for graph_file, metis_file, _, _ in instances:
            run_graph(graph_file, metis_file, configurations)
    else:
        grouped_instances = group_instances_by_dataset_and_size(instances)
        grouped_link_configurations = group_link_configurations(link_configurations)
        for dataset_name, instances_by_size in sorted(grouped_instances.items()):
            size_levels = sorted(instances_by_size)
            for link_family, link_names in grouped_link_configurations.items():
                with tempfile.TemporaryDirectory(
                    prefix="heiconnect-adaptive-links-"
                ) as cache_dir:
                    link_cache_dir = Path(cache_dir)
                    for configuration_name in configurations:
                        print(
                            f"[{configuration_name}] Discovering frontier for "
                            f"{dataset_name}/{link_family}"
                        )

                        def run_level(node_count):
                            level_completed = False
                            print(
                                f"[{configuration_name}] Testing {dataset_name}/"
                                f"{link_family} at {node_count} nodes"
                            )
                            for graph_file, metis_file in instances_by_size[
                                node_count
                            ]:
                                for link_name in link_names:
                                    outcomes = run_graph(
                                        graph_file,
                                        metis_file,
                                        [configuration_name],
                                        [link_name],
                                        link_cache_dir,
                                    )
                                    outcome = outcomes[configuration_name][0]
                                    if outcome == "timeout":
                                        return "timeout"
                                    if outcome == "failed":
                                        return "failed"
                                    level_completed = True
                            return "solved" if level_completed else "failed"

                        outcomes = run_adaptive_series(
                            size_levels, arguments.adaptive_probes, run_level
                        )
                        tested_sizes = ", ".join(
                            f"{size}={outcome}"
                            for size, outcome in sorted(outcomes.items())
                        )
                        print(
                            f"[{configuration_name}/{dataset_name}/{link_family}] "
                            f"{tested_sizes}"
                        )
                        adaptive_reports[configuration_name].append(
                            f"dataset={dataset_name} "
                            f"link_configuration={link_family}"
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
