import re
from sys import argv

import matplotlib.pyplot as plt

if __name__ == "__main__":
    state = 0
    data = {}
    for i in range(1, len(argv)):
        filename = argv[i]
        data[filename] = {"cycle_sizes": [], "times_ms": []}
        with open(filename, "r") as f:
            for line in f.readlines():
                if state == 0:
                    if line.find(r"{") != -1:
                        state = 1
                        continue
                    if line.find(": ") == -1:
                        continue
                    cycle_size = int(line.split(":")[0].split(" ")[1])
                    time_ms = line.split(": ")[1].split(" ")[0]
                    if cycle_size > 300:
                        continue
                    data[filename]["cycle_sizes"].append(cycle_size)
                    data[filename]["times_ms"].append(float(time_ms))
                if state == 1:
                    if line.find(r"}"):
                        state = 0

        plt.scatter(
            data[filename]["cycle_sizes"],
            data[filename]["times_ms"],
            label=filename,
        )
    plt.xlabel("Cycle Size")
    plt.ylabel("Time (ms)")
    plt.title("Graph Augmentation Time for Cycle Sizes")
    plt.legend()
    output_file = "plot_results.png"
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"Plot saved to {output_file}")
