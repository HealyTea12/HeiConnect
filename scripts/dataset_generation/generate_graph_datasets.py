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
}


def sizes_from_ranges(ranges):
    sizes = set()
    for start, stop, step in ranges:
        sizes.update(range(start, stop + 1, step))
    return sorted(sizes)


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
                        [
                            executable,
                            "--type",
                            "graph",
                            "--generator",
                            settings["generator"],
                            "--nodes",
                            str(size),
                            "--seed",
                            str(seed),
                            "--output",
                            output_dir,
                        ],
                        check=True,
                    )
                    progress.update()


if __name__ == "__main__":
    main()
