#!/usr/bin/env python3

from pathlib import Path
import subprocess


CONFIG = {
    "folders": ["stars", "cycles", "trees", "cacti"],
    "distributions": {
        "float_uniform_0_1": {
            "arguments": [
                "--distribution",
                "float_uniform",
                "--float_uniform_lower",
                "0",
                "--float_uniform_upper",
                "1",
            ],
            "seeds": [42, 43, 44, 45, 46],
        },
        "int_uniform_1_5": {
            "arguments": [
                "--distribution",
                "integer_uniform",
                "--integer_uniform_lower",
                "1",
                "--integer_uniform_upper",
                "5",
            ],
            "seeds": [42, 43, 44, 45, 46],
        },
    },
}


def link_generator_command(
    executable, graph_file, output_file, distribution_arguments, seed
):
    return [
        executable,
        "--type",
        "links",
        "--generator",
        "complete",
        "--input_graph",
        graph_file,
        "--output",
        output_file,
        *distribution_arguments,
        "--seed",
        str(seed),
    ]


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
            for suffix, settings in CONFIG["distributions"].items():
                for seed in settings["seeds"]:
                    output_file = graph_file.with_name(
                        f"{graph_file.stem}-{suffix}_seed_{seed}.links"
                    )
                    subprocess.run(
                        link_generator_command(
                            executable,
                            graph_file,
                            output_file,
                            settings["arguments"],
                            seed,
                        ),
                        check=True,
                    )


if __name__ == "__main__":
    main()
