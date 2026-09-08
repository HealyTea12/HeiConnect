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
import json
import os
from pathlib import Path
import re
import resource
import signal
import subprocess
import sys
import tempfile
import time
import tomllib

from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dataset_generation.generate_graph_datasets import (
    cactus_cycles,
    graph_generator_command,
)


VALID_NAME = re.compile(r"[A-Za-z0-9_.-]+")
VALID_DISTRIBUTIONS = {"constant", "float_uniform", "integer_uniform"}


def compact_label(value, max_length=64):
    text = str(value)
    if len(text) <= max_length:
        return text
    side_length = (max_length - 1) // 2
    return f"{text[:side_length]}…{text[-side_length:]}"


def format_duration(seconds):
    seconds = int(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes:02d}m {seconds:02d}s"
    if minutes:
        return f"{minutes}m {seconds:02d}s"
    return f"{seconds}s"


def positive_integer(value):
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ValueError("must be a positive integer")
    return value


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run existing datasets and search synthetic solvability frontiers."
    )
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "-c", "--configurations", required=True, type=Path, metavar="FILE"
    )
    parser.add_argument(
        "--no-repeat",
        action="store_true",
        help="skip runs that already contain a completed named result",
    )
    parser.add_argument(
        "--mode",
        choices=("exhaustive", "adaptive"),
        default="exhaustive",
        help=(
            "scheduling for existing files; synthetic families always search their frontier"
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
    if "dataset_selection" in document:
        input_path = dataset_selection.get("input_dir")
        if not isinstance(input_path, str) or not input_path.strip():
            raise ValueError("dataset_selection.input_dir must be a non-empty path string")
        input_dir = Path(input_path).expanduser()
        if not input_dir.is_absolute():
            input_dir = path.resolve().parent / input_dir
        dataset_selection["input_dir"] = str(input_dir.resolve())
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
        if selection.get("min_nodes") is not None:
            positive_integer(selection["min_nodes"])
        positive_integer(selection.get("max_nodes"))
        if selection.get("min_nodes", 1) > selection["max_nodes"]:
            raise ValueError(
                f"dataset selection {name!r} min_nodes must not exceed max_nodes"
            )

    synthetic = validate_synthetic_datasets(document.get("synthetic_datasets", {}))
    if "dataset_selection" not in document:
        dataset_selection = None
    if dataset_selection is None and not synthetic:
        raise ValueError("configure dataset_selection or synthetic_datasets to select instances")
    return configurations, link_configurations, dataset_selection, synthetic


def validate_synthetic_datasets(datasets):
    if not isinstance(datasets, dict):
        raise ValueError("synthetic_datasets must be a TOML table")
    validated = {}
    for name, settings in datasets.items():
        validate_name(name, "synthetic dataset")
        if name in {".", ".."}:
            raise ValueError("synthetic dataset names must not be '.' or '..'")
        if not isinstance(settings, dict):
            raise ValueError(f"synthetic dataset {name!r} must be a table")
        settings = dict(settings)
        allowed = {
            "generator", "min_nodes", "max_nodes", "samples", "resolution",
            "seeds", "cycle_length", "cycle_mass", "generation_timeout",
            "generation_max_memory_mb",
            "min_cycle_size", "max_cycle_size",
        }
        if set(settings) - allowed:
            raise ValueError(f"unknown synthetic settings for {name!r}: {set(settings) - allowed}")
        generator = settings.get("generator")
        if not isinstance(generator, str) or generator not in {"cycle", "star", "tree", "cactus", "cactus_variable"}:
            raise ValueError(f"unsupported synthetic generator for {name!r}")
        for key, default in {
            "min_nodes": 10, "samples": 10, "resolution": 1,
            "generation_timeout": 300, "generation_max_memory_mb": 4096,
        }.items():
            settings.setdefault(key, default)
            positive_integer(settings[key])
        minimum_nodes = {"cycle": 3, "star": 2, "tree": 2, "cactus": 1, "cactus_variable": 1}
        if settings["min_nodes"] < minimum_nodes[generator]:
            raise ValueError(f"min_nodes is too small for {generator}")
        if "max_nodes" in settings:
            positive_integer(settings["max_nodes"])
            if settings["max_nodes"] < settings["min_nodes"]:
                raise ValueError("synthetic max_nodes must be at least min_nodes")
        seeds = settings.setdefault("seeds", [42])
        if not isinstance(seeds, list) or not seeds or any(
            not isinstance(seed, int) or isinstance(seed, bool)
            or not 0 <= seed <= 2**32 - 1 for seed in seeds
        ):
            raise ValueError("synthetic seeds must be a non-empty array of unsigned 32-bit integers")
        if len(set(seeds)) != len(seeds):
            raise ValueError("synthetic seeds must be unique")
        if generator in {"cycle", "star"}:
            if len(seeds) > 1:
                raise ValueError(f"{generator} has only one graph per size; use at most one graph seed")
            settings["seeds"] = [42]
        if generator == "cactus":
            positive_integer(settings.get("cycle_length"))
            mass = settings.get("cycle_mass")
            if not isinstance(mass, (int, float)) or isinstance(mass, bool):
                raise ValueError("cactus cycle_mass must be a number")
            cactus_cycles(settings["min_nodes"], settings["cycle_length"], mass)
        elif "cycle_length" in settings or "cycle_mass" in settings:
            raise ValueError("cycle_length and cycle_mass apply only to cactus graphs")
        if generator == "cactus_variable":
            positive_integer(settings.get("min_cycle_size"))
            positive_integer(settings.get("max_cycle_size"))
            if not 2 <= settings["min_cycle_size"] <= settings["max_cycle_size"]:
                raise ValueError("cycle sizes must satisfy 2 <= min_cycle_size <= max_cycle_size")
        elif "min_cycle_size" in settings or "max_cycle_size" in settings:
            raise ValueError("min_cycle_size and max_cycle_size apply only to cactus_variable graphs")
        validated[name] = settings
    return validated


def evenly_spaced_sizes(lower, upper, count):
    count = min(count, upper - lower + 1)
    if count == 1:
        return [upper]
    return [lower + (upper - lower) * index // (count - 1) for index in range(count)]


def run_synthetic_series(settings, run_level):
    outcomes = {}
    lower = None
    upper = None
    size = settings["min_nodes"]
    cap = settings.get("max_nodes")
    status = "capped"

    def probe(size):
        if size not in outcomes:
            outcomes[size] = run_level(size)
        return outcomes[size]

    while True:
        outcome = probe(size)
        if outcome != "solved":
            status = outcome
            if outcome == "timeout":
                upper = size
            break
        lower = size
        if size == cap:
            break
        size = size * 2 if cap is None else min(size * 2, cap)

    if lower is not None and upper is not None:
        while upper - lower > settings["resolution"]:
            size = (lower + upper) // 2
            outcome = probe(size)
            if outcome == "solved":
                lower = size
            elif outcome == "timeout":
                upper = size
            else:
                status = outcome
                break
        if status == "timeout":
            status = "bracketed"

    samples = []
    if lower is not None:
        samples = evenly_spaced_sizes(settings["min_nodes"], lower, settings["samples"])
        for size in reversed(samples):
            probe(size)
        if any(outcomes[size] != "solved" for size in samples):
            status = "non_monotone"
    smallest_timeout = min(
        (size for size, outcome in outcomes.items() if outcome == "timeout"),
        default=None,
    )
    return {
        "status": status, "largest_solved": lower, "smallest_timeout": smallest_timeout,
        "samples": samples, "outcomes": outcomes,
    }


def generate_synthetic_graph(executable, output_dir, node_count, seed, settings):
    try:
        return_code = run_command(
            graph_generator_command(executable, output_dir, node_count, seed, settings),
            settings["generation_timeout"],
            settings["generation_max_memory_mb"],
        )
    except subprocess.TimeoutExpired:
        return None
    graphs = list(output_dir.glob("*.xml"))
    if return_code != 0 or len(graphs) != 1:
        return None
    graph_file = graphs[0]
    metis_file = graph_file.with_suffix(".graph")
    if not metis_file.is_file() or read_node_count(metis_file) != node_count:
        return None
    return graph_file, metis_file


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
        limits = dataset_limits.get(dataset_name, {})
        min_nodes = limits.get("min_nodes")
        max_nodes = limits.get(
            "max_nodes", default_max_nodes
        )
        if min_nodes is not None and node_count < min_nodes:
            skipped += 1
            continue
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
    path, name, configuration, link_configurations, dataset_selection, synthetic=None
):
    dataset_selection = dataset_selection or {}
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
    if "input_dir" in dataset_selection:
        lines.append(f"dataset_selection.input_dir={dataset_selection['input_dir']}")
    for dataset_name in dataset_selection.get("include", []):
        lines.append(f"dataset_selection.include={dataset_name}")
    if dataset_selection.get("default_max_nodes") is not None:
        lines.append(
            "dataset_selection.default_max_nodes="
            f"{dataset_selection['default_max_nodes']}"
        )
    for dataset_name, selection in dataset_selection.get("datasets", {}).items():
        if selection.get("min_nodes") is not None:
            lines.append(
                f"dataset_selection.{dataset_name}.min_nodes="
                f"{selection['min_nodes']}"
            )
        lines.append(
            f"dataset_selection.{dataset_name}.max_nodes={selection['max_nodes']}"
        )
    for dataset_name, settings in (synthetic or {}).items():
        for key, value in settings.items():
            lines.append(f"synthetic_dataset.{dataset_name}.{key}={serialize_parameter(value)}")
    path.write_text("\n".join(lines) + "\n")


def has_complete_result(path):
    if not path.is_file():
        return False
    with path.open() as file:
        return any(line.startswith("run.total_time_seconds=") for line in file)


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
    start_time = time.monotonic()
    arguments = parse_arguments()
    positive_integer(arguments.adaptive_probes)
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

    configurations, link_configurations, dataset_selection, synthetic = load_configurations(
        arguments.configurations
    )
    input_dir = Path(dataset_selection["input_dir"]) if dataset_selection is not None else None
    if input_dir is not None and not input_dir.is_dir():
        raise FileNotFoundError(f"input directory not found: {input_dir}")
    if synthetic and any(configuration.get("timeout") is None for configuration in configurations.values()):
        raise ValueError("every algorithm configuration needs a timeout for synthetic frontier search")
    for executable in (experiments_executable, link_generator_executable):
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise FileNotFoundError(
                f"required executable not found or not executable: {executable}"
            )

    instances, filtered_instances = (
        find_instances(input_dir, dataset_selection)
        if dataset_selection is not None else ([], 0)
    )
    existing_names = {instance[2] for instance in instances}
    if existing_names & synthetic.keys():
        raise ValueError(
            "existing and synthetic dataset names overlap; use distinct names: "
            + ", ".join(sorted(existing_names & synthetic.keys()))
        )
    if not instances and not synthetic:
        raise ValueError("no .xml files with matching .graph files were found")
    total_runs = len(instances) * len(link_configurations) * len(configurations)
    instances_by_dataset = defaultdict(int)
    for _, _, dataset_name, _ in instances:
        instances_by_dataset[dataset_name] += 1
    dataset_summary = ", ".join(
        f"{name}: {count:,}" for name, count in sorted(instances_by_dataset.items())
    )
    print("\nExperiment plan")
    print(f"  Mode             {arguments.mode}")
    print(f"  Configurations   {len(configurations):,}")
    print(f"  Link variants    {len(link_configurations):,}")
    print(f"  Instances        {len(instances):,} ({dataset_summary})")
    print(f"  Filtered         {filtered_instances:,}")
    print(f"  Matrix size      {total_runs:,} runs")
    if synthetic:
        print(f"  Synthetic        {', '.join(synthetic)} (run count discovered during search)")
    print(f"  Output           {output_dir}\n")
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
            synthetic,
        )

    completed = 0
    failed = 0
    timed_out = 0
    skipped = 0
    processed = 0
    progress = tqdm(
        total=None if synthetic else total_runs,
        desc="Experiments",
        unit="run",
        dynamic_ncols=True,
        bar_format=None if synthetic else (
            "{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} "
            "[{elapsed}<{remaining}, {rate_fmt}] {postfix}"
        ),
    )
    adaptive_reports = defaultdict(list)

    def run_graph(
        graph_file,
        metis_file,
        selected_configurations,
        selected_link_configurations=None,
        link_cache_dir=None,
        relative_path=None,
        generation_settings=None,
    ):
        nonlocal completed, failed, timed_out, skipped, processed
        if relative_path is None:
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
                    and has_complete_result(result_file)
                ):
                    skipped += 1
                    processed += 1
                    outcomes[configuration_name].append("completed")
                    progress.set_postfix_str(
                        compact_label(
                            f"skip · {configuration_name} · {link_name} · {relative_path}"
                        ),
                        refresh=False,
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
                progress.set_postfix_str(
                    compact_label(f"links · {link_name} · {relative_path}")
                )
                generation_settings = generation_settings or {}
                try:
                    generator_result = run_command(
                        link_generator_command(
                            link_generator_executable,
                            metis_file,
                            links_file,
                            link_configuration,
                        ),
                        generation_settings.get("generation_timeout"),
                        generation_settings.get("generation_max_memory_mb"),
                    )
                except subprocess.TimeoutExpired:
                    generator_result = -1
                if generator_result != 0:
                    print(
                        f"[{link_name}] Link generation failed for {graph_file}",
                        file=sys.stderr,
                    )
                    links_file.unlink(missing_ok=True)
                    failed += len(pending)
                    processed += len(pending)
                    for configuration_name, _, _ in pending:
                        outcomes[configuration_name].append(
                            "generation_failed" if generation_settings else "failed"
                        )
                        progress.set_postfix_str(
                            compact_label(
                                f"failed · {configuration_name} · {link_name} · {relative_path}"
                            ),
                            refresh=False,
                        )
                        progress.update()
                    if temporary_links is not None:
                        temporary_links.cleanup()
                    continue
            for configuration_name, configuration, instance_output in pending:
                instance_output.mkdir(parents=True, exist_ok=True)
                progress.set_postfix_str(
                    compact_label(
                        f"run · {configuration_name} · {link_name} · {relative_path}"
                    )
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
                    compact_label(
                        f"{status} · {configuration_name} · {link_name} · {relative_path}"
                    ),
                    refresh=False,
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
                        progress.set_postfix_str(
                            compact_label(
                                f"frontier · {configuration_name} · {dataset_name} · {link_family}"
                            )
                        )

                        def run_level(node_count):
                            level_completed = False
                            progress.set_postfix_str(
                                compact_label(
                                    f"probe · {configuration_name} · {dataset_name} · "
                                    f"{link_family} · n={node_count}"
                                )
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
    existing_timed_out = timed_out
    if censored and not synthetic:
        progress.update(censored)

    synthetic_reports = defaultdict(list)
    for dataset_name, settings in synthetic.items():
        for link_family, link_names in group_link_configurations(link_configurations).items():
            for configuration_name in configurations:
                def run_level(node_count):
                    nonlocal failed
                    progress.set_postfix_str(compact_label(
                        f"synthetic · {configuration_name} · {dataset_name} · n={node_count}"
                    ))
                    for seed in settings["seeds"]:
                        with tempfile.TemporaryDirectory(prefix="heiconnect-graph-") as directory:
                            graph = generate_synthetic_graph(
                                link_generator_executable, Path(directory),
                                node_count, seed, settings,
                            )
                            if graph is None:
                                failed += 1
                                print(f"Graph generation failed: {dataset_name}, n={node_count}, seed={seed}", file=sys.stderr)
                                return "generation_failed"
                            graph_file, metis_file = graph
                            for link_name in link_names:
                                outcomes = run_graph(
                                    graph_file, metis_file, [configuration_name], [link_name],
                                    relative_path=Path(dataset_name) / graph_file.stem,
                                    generation_settings=settings,
                                )
                                outcome = outcomes[configuration_name][0]
                                if outcome != "completed":
                                    return outcome
                    return "solved"

                report = run_synthetic_series(settings, run_level)
                report.update(dataset=dataset_name, link_configuration=link_family)
                synthetic_reports[configuration_name].append(report)
                report_path = output_dir / configuration_name / "synthetic-frontier.json"
                report_path.write_text(json.dumps(synthetic_reports[configuration_name], indent=2) + "\n")
                tqdm.write(
                    f"[{configuration_name}] {dataset_name}/{link_family}: "
                    f"{report['status']}, largest solved={report['largest_solved']}, "
                    f"smallest timeout={report['smallest_timeout']}"
                )
    progress.close()
    print("\nExperiment summary")
    print(f"  Completed        {completed:,}")
    print(f"  Skipped          {skipped:,}")
    print(f"  Timed out        {timed_out:,}")
    print(f"  Failed           {failed:,}")
    print(f"  Censored         {censored:,}")
    print(f"  Elapsed          {format_duration(time.monotonic() - start_time)}")
    print(f"  Results          {output_dir}")
    return 1 if failed or (arguments.mode == "exhaustive" and existing_timed_out) else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, tomllib.TOMLDecodeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
