# Experiment configurations

For a single server-ready thesis configuration, run:

```sh
uv run scripts/experiment_runner/run_configurations.py results/thesis-server \
  -c scripts/experiment_runner/configurations.thesis.server.toml --no-repeat
```

This generates all synthetic inputs and runs 1,800 comparisons. Build the current
experiment/generator binaries first; see [THESIS_EXPERIMENTS.md](THESIS_EXPERIMENTS.md)
for build commands and the rationale behind the design.

For the bounded thesis pilot and main comparison, see
[THESIS_EXPERIMENTS.md](THESIS_EXPERIMENTS.md). These use cycles, stars, trees,
and variable cacti with cycle-size bounds 2 and 16.
The claim-driven main setup reuses 15 configurations across representation,
reduction, and Supernova comparisons. Calibrate it on the experiment machine;
`configurations.thesis.real_world.toml` provides a separate finalist comparison.

Lazy Block Tree ILP results include `solver.status` and `solution.optimal`.
Only an optimal solver status (or a fully reduced trivial instance) sets optimality
to true. False means optimality was not established; consult the status and whether
a solution cost exists. The CSV collector includes `solver_status` and
`solution_optimal`; missing fields in older or interrupted runs remain unknown.

Synthetic datasets can set `sizes = [20, 40, 80, 160, 320]` to run a fixed grid
instead of frontier search. Do not combine `sizes` with `min_nodes`, `max_nodes`,
`samples`, or `resolution`. Every selected repetition is attempted even after a
timeout. Outcomes and execution order are recorded in `synthetic-fixed.json`.

With explicit sizes, `seed_mode = "paired"` pairs tree/cactus graph seeds with
link seeds by array position within each distribution. The counts must match.
Stars and cycles use their single graph with all link seeds. The default
`seed_mode = "cross"` retains the full cross-product. Fixed-grid execution shares
generated graphs and links across algorithms and shuffles algorithm order with
scheduling seed 42. Without `sizes`, the frontier workflow below applies.

For synthetic graphs, start with `configurations.synthetic.example.toml`:

```sh
uv run scripts/experiment_runner/run_configurations.py results/synthetic \
  -c scripts/experiment_runner/configurations.synthetic.example.toml
```

The configuration tables select which instances to run and where existing datasets
are stored. The only positional argument is the output directory. Graphs and complete
links are generated in temporary directories as sizes are requested. Results keep
the usual `<output>/<algorithm configuration>/<dataset>/<instance>/res-<link>.txt`
layout, which is supported by `scripts/visualization/collect_results.py`.

Each `[synthetic_datasets.<name>]` selects a `generator`: `cycle`, `star`, `tree`,
`cactus`, `cactus_variable`, or `cactus_cycles`. Parameters are:

| Setting | Default | Meaning |
| --- | --- | --- |
| `min_nodes` | 10 | Smallest size to test and lower endpoint for plotting samples |
| `samples` | 10 | Number of evenly spaced integer sizes from the minimum through the largest solved probe |
| `resolution` | 1 | Stop bisection when the solved/timeout bracket is this many nodes wide |
| `max_nodes` | none | Optional search cap; success here only establishes a lower bound |
| `seeds` | `[42]` | Graph realizations; multiple seeds are allowed for trees and cacti only |
| `cycle_length` | required for cactus | Length of each cactus cycle, at least 3 |
| `cycle_mass` | required for cactus | Fraction in `[0, 1]` used to compute `floor(floor((n-1)*mass)/(length-1))` cycles |
| `cycles` | required for cactus_cycles | Exact cycle count including bridges as length-two cycles: 1 through `n-1`, or 0 for a single-node graph |
| `min_cycle_size` | required for cactus_variable | Inclusive lower bound, at least 2; size 2 adds a bridge |
| `max_cycle_size` | required for cactus_variable | Inclusive upper bound, at least `min_cycle_size` |
| `generation_timeout` | 300 | Seconds allowed for each graph or link generation command |
| `generation_max_memory_mb` | 4096 | Address-space limit for each generation command |

`cactus_variable` starts with one vertex. At each step it chooses an existing
attachment vertex uniformly, samples a cycle size uniformly from the configured
inclusive bounds, and adds that cycle using the attachment vertex plus new
vertices. A sampled size of 2 adds a bridge to one new vertex. A cycle of size `k`
consumes `k - 1` new vertices. If the sampled cycle
would exceed the remaining budget, it adds the largest cycle that fits, even when
that is smaller than `min_cycle_size`. If only one new vertex remains, it adds a
bridge instead. This produces exactly the requested number of vertices.
Cycle edges have weight 1; bridges have weight 2. Setting both bounds to 2 produces
a random tree. With a minimum of at least 3, only the final addition can be a bridge.

This generator takes no `cycle_mass` or fixed `cycle_length`. Cycle-size samples
are uniform before truncation; this construction does not sample uniformly from
all possible cactus graphs. Seeds control the attachment choices, cycle sizes,
and final vertex-label permutation. It is also available directly:

```sh
build/experiments/generate_datasets --type graph --generator cactus_variable \
  --nodes 100 --min-cycle-size 3 --max-cycle-size 8 --seed 42 --output datasets/variable_cacti
```

`cactus_cycles` generates exactly `n` vertices and `cycles` blocks, counting each
bridge as a cycle of length two without storing parallel edges. It starts each
block with one attachment vertex and one new vertex. Each remaining vertex chooses
a block uniformly to extend; blocks that stay at size two become bridges.
The blocks are shuffled, each attaches at a uniformly chosen existing
vertex, and vertex labels are permuted. This takes linear time and space and does
not sample uniformly from all cactus graphs. Cycle edges have weight 1 and bridges
have weight 2. Setting `cycles = n-1` produces a tree, and `cycles = 1` produces
a single cycle for `n >= 3` (a bridge for `n = 2`). Zero is valid only for `n = 1`.
The same inputs and seed reproduce the same graph. For example:

```sh
build/experiments/generate_datasets --type graph --generator cactus_cycles \
  --nodes 100 --cycles 12 --seed 42 --output datasets/counted_cacti
```

This writes `cactus_cycles_n100_q12_seed42.graph` and `.xml`. The GraphML file
records the generator, node and edge counts, seed, total `cycles` including bridges,
`proper_cycles` of length at least three, and `bridges`. Thus
`cycles = proper_cycles + bridges`. In a
synthetic dataset table, use `generator = "cactus_cycles"` and `cycles = 12` with
`min_nodes >= 13` (or explicit `sizes` all at least 13). The cycle count stays
fixed as the node count grows.

Define shared algorithm limits once per configuration file:

```toml
[budget]
timeout = 60
max_memory_mb = 4096
```

`timeout` is in seconds and `max_memory_mb` is in MiB; both must be positive
integers. These defaults apply to explicit configurations and configuration
matrices. Individual configurations or matrices can override either limit.
Existing files with only per-algorithm limits remain supported. Every algorithm
needs an effective `timeout` for synthetic experiments. Generation limits are
configured separately with `generation_timeout` and `generation_max_memory_mb`.
Generation failures are reported separately from algorithm timeouts and do not
establish an algorithm frontier.

For each algorithm configuration, dataset, and link distribution family, the
runner doubles the node count until a timeout, then bisects between the last
solved size and the timeout. A size is solved only if every graph seed and every
link seed completes. The first unsuccessful repetition ends a probe. Stars and
cycles have one graph per size; their link weights can still vary across seeds.
Use one link seed for constant weights to avoid identical repetitions.

The runner then fills the requested samples in descending order, reusing probes
already measured. Integer spacing differs by at most one node. A single sample
selects the frontier; short intervals include each available integer size once.
Search probes remain as additional measurements, so total plotted sizes may
exceed `samples`. Different algorithms can have different sample grids; solution
quality comparisons can only use their shared instances.

`synthetic-frontier.json` in each algorithm directory records the largest solved
probe, the timeout bound, selected samples, and every measured outcome. `bracketed`
means the requested resolution was reached; `capped` means the cap was solved.
Failures leave the search incomplete. As with any binary search, the method assumes
solvability decreases with size. Random instances and timing noise can violate
that assumption: a backfill failure is marked `non_monotone`, and the observations
are retained. The report describes an empirical boundary for the configured seeds
and resource limits, rather than a guaranteed largest solvable graph.

`--no-repeat` reuses result files containing a completion marker. Failed and timed
out probes are retried when resuming. Use a fresh output directory when changing
generation parameters, seeds, algorithm parameters, or resource limits.

Existing files retain their previous workflow and `[dataset_selection]` settings:

```sh
uv run scripts/experiment_runner/run_configurations.py results/existing \
  -c scripts/experiment_runner/configurations.all.toml --mode adaptive
```

`--mode exhaustive` (the default) runs all selected existing files; `--mode adaptive`
uses the previous bounded probes and backfill, with `--adaptive-probes` controlling
the probe budget. Synthetic families without explicit sizes use the new frontier search.
Include both table types in the same configuration file to run both sources:

```toml
[dataset_selection]
input_dir = "../../datasets"
include = ["real_world"]

[dataset_selection.datasets.real_world]
max_nodes = 10000

[synthetic_datasets.trees]
generator = "tree"
seeds = [42, 43, 44]
samples = 10
```

Omit `[dataset_selection]` to run only synthetic instances, or omit
`[synthetic_datasets.<name>]` tables to run only existing instances.
`dataset_selection.input_dir` is required whenever `[dataset_selection]` is present.
Relative input paths are resolved from the configuration file's directory;
absolute paths and `~` are also supported. A selection containing only `input_dir`
selects all existing files below that directory. Existing and synthetic
dataset names must be distinct in a combined run to prevent result collisions.
For older configurations, move the former positional dataset path into
`dataset_selection.input_dir` and remove it from the command line.
