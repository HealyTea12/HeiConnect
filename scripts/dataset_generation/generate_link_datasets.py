#!/usr/bin/env python3

from pathlib import Path
import subprocess


CONFIG = {
    "folders": ["stars", "cycles"],
    "distributions": {
        "float_uniform_0_1": [
            "--distribution",
            "float_uniform",
            "--float_uniform_lower",
            "0",
            "--float_uniform_upper",
            "1",
        ],
        "int_uniform_1_10": [
            "--distribution",
            "integer_uniform",
            "--integer_uniform_lower",
            "1",
            "--integer_uniform_upper",
            "10",
        ],
    },
}


def main():
    project_root = Path(__file__).resolve().parents[2]
    executable = project_root / "build" / "experiments" / "generate_datasets"
    datasets_dir = project_root / "datasets"

    if not executable.is_file():
        raise FileNotFoundError(
            f"{executable} does not exist. Build the generate_datasets target first."
        )

    for folder in CONFIG["folders"]:
        for graph_file in sorted((datasets_dir / folder).glob("*.graph")):
            for suffix, arguments in CONFIG["distributions"].items():
                output_file = graph_file.with_name(f"{graph_file.stem}-{suffix}.links")
                subprocess.run(
                    [
                        executable,
                        "--type",
                        "links",
                        "--generator",
                        "complete",
                        "--input_graph",
                        graph_file,
                        "--output",
                        output_file,
                        *arguments,
                    ],
                    check=True,
                )


if __name__ == "__main__":
    main()
