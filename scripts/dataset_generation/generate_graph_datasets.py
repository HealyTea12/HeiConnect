#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "tqdm",
# ]
# ///

from pathlib import Path
import subprocess

from tqdm import tqdm


CONFIG = {
    "stars": {
        "generator": "star",
        "ranges": [(10, 400, 10), (500, 2000, 100), (2000, 10000, 400)],
        "seeds": [42],
    },
    "cycles": {
        "generator": "cycle",
        "ranges": [(10, 500, 10), (500, 2000, 100), (2000, 3000, 200)],
        "seeds": [42],
    },
    "trees": {
        "generator": "tree",
        "ranges": [(10, 400, 10), (500, 2000, 100), (2000, 10000, 400)],
        "seeds": [42, 43, 44, 45, 46],
    },
    "cacti": {
        "generator": "cactus",
        "ranges": [(10, 400, 10), (500, 2000, 100), (2000, 10000, 400)],
        "seeds": [42, 43, 44, 45, 46],
        "cycle_length": 4,
        "cycle_mass": 0.75,
    },
}


def sizes_from_ranges(ranges):
    sizes = set()
    for start, stop, step in ranges:
        sizes.update(range(start, stop + 1, step))
    return sorted(sizes)


def cactus_cycles(node_count, cycle_length, cycle_mass):
    if cycle_length < 3:
        raise ValueError("cactus cycle_length must be at least 3")
    if not 0.0 <= cycle_mass <= 1.0:
        raise ValueError("cactus cycle_mass must be between 0 and 1")
    return int((node_count - 1) * cycle_mass) // (cycle_length - 1)


def graph_generator_command(executable, output_dir, node_count, seed, settings):
    command = [
        executable,
        "--type",
        "graph",
        "--generator",
        settings["generator"],
        "--nodes",
        str(node_count),
        "--seed",
        str(seed),
        "--output",
        output_dir,
    ]
    if settings["generator"] == "cactus":
        cycle_length = settings["cycle_length"]
        cycles = cactus_cycles(node_count, cycle_length, settings["cycle_mass"])
        command.extend(
            ["--cycles", str(cycles), "--cycle-length", str(cycle_length)]
        )
    elif settings["generator"] == "cactus_cycles":
        command.extend(["--cycles", str(settings["cycles"])])
    elif settings["generator"] == "cactus_variable":
        command.extend([
            "--min-cycle-size", str(settings["min_cycle_size"]),
            "--max-cycle-size", str(settings["max_cycle_size"]),
        ])
    return command


def main():
    project_root = Path(__file__).resolve().parents[2]
    executable = project_root / "build" / "experiments" / "generate_datasets"
    datasets_dir = project_root / "datasets"

    if not executable.is_file():
        raise FileNotFoundError(
            f"{executable} does not exist. Build the generate_datasets target first."
        )

    total = sum(
        len(sizes_from_ranges(settings["ranges"])) * len(settings["seeds"])
        for settings in CONFIG.values()
    )

    with tqdm(total=total, desc="Generating graph datasets", unit="graph") as progress:
        for folder, settings in CONFIG.items():
            output_dir = datasets_dir / folder
            for size in sizes_from_ranges(settings["ranges"]):
                for seed in settings["seeds"]:
                    progress.set_postfix(type=folder, nodes=size, seed=seed)
                    subprocess.run(
                        graph_generator_command(
                            executable,
                            output_dir,
                            size,
                            seed,
                            settings,
                        ),
                        check=True,
                    )
                    progress.update()


if __name__ == "__main__":
    main()
