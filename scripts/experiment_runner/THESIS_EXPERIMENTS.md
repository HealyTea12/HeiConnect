# Thesis experiment protocol

The initial revised pilot completed 474 of 480 runs, with six timeouts. It exposed
a bug in element domination: equivalent constraints could all be removed. The
reducer now retains one representative, and all four exact configurations agree
on the previously failing tree. Rebuild before running the server configuration
and use fresh result directories; pre-fix reduced results are not reliable.
See [THESIS_REPRODUCTION.md](THESIS_REPRODUCTION.md) for the diagnosis and validation.

## Run this on the server

Use **`scripts/experiment_runner/configurations.thesis.server.toml`** for the
complete synthetic comparison. It is self-contained: no pre-generated datasets,
configuration merging, or pilot-dependent algorithm selection is needed. It fixes
the main design below at 1,800 executions. The pilot is a separate calibration and
validation run; it is not a prerequisite built into the server configuration.

From the repository root, with the project's C++ dependencies and a working
Gurobi license installed on the server:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DHEICONNECT_BUILD_TESTS=OFF -DHEICONNECT_BUILD_EXPERIMENTS=ON
cmake --build build --target experiments generate_datasets -j 2

uv run scripts/experiment_runner/run_configurations.py results/thesis-server \
  -c scripts/experiment_runner/configurations.thesis.server.toml

uv run scripts/visualization/collect_results.py results/thesis-server results/thesis-server.csv
```

Run this in a persistent server session or scheduled job so disconnecting SSH does
not interrupt the experiment. Use the same runner command to resume: completed
results are reused; failed and timed-out attempts are retried. Do not repeatedly
resume a finished run just to remove its timeouts. Record the revision, machine,
and actual solver thread count alongside the results. The configuration limits
each execution to 60 seconds and 4 GiB of address space. Generation is separately
limited to 300 seconds and 4 GiB; its overhead is outside the execution budget.
Algorithms run sequentially, but Gurobi can use multiple threads internally.
The configuration does not pin its thread count; keep the server's CPU allocation
fixed and report that setting when interpreting runtimes.

The real-world comparison remains separate because it needs a dataset directory
on the server. Its ready-to-edit configuration and input semantics are below.

This design follows `../report/chapters/plan/exp-results.txt` (relative to the
repository root). Its purpose is to measure the planned claims, including negative
results, rather than screen every available option. Run the preflight, pilot,
and measurements on the experiment machine. Local Gurobi failures do not determine
which algorithms belong in the thesis comparison.

## Comparisons and reuse

The main configuration contains 15 configurations on one shared set of instances.
Reuse their measurements across the following sections; do not rerun each section
independently. Configuration names below refer to the TOML/result directory names.

| Thesis question | Configurations | Main measurements |
| --- | --- | --- |
| Is SetCoverCyc/CycContext more efficient than GWC? | `sc_double`, `sc_no_reductions`, `gwc` | End-to-end runtime and peak memory; also cost to expose different solutions/tie-breaking |
| How does greedy-cheapest compare with MSTConnect? | `sc_cheapest`, `mst_connect` | Runtime, memory, and solution cost |
| Do reductions accelerate SCILP? | `sc_ilp_no_reductions`, `sc_ilp` | Runtime, peak memory, completion rates, equal optimal cost where verified |
| Do reductions accelerate BTSC with ILP? | `bt_ilp_no_reductions`, `bt_ilp` | Runtime, peak memory, completion rates, solution cost |
| Do reductions accelerate Lazy Supernova? | `lazy_block_tree_ilp_no_reductions`, `lazy_block_tree_ilp` | Runtime, peak memory, completion rates, equal optimal cost where verified |
| Do reductions improve greedy quality? | `sc_no_reductions` vs `sc_greedy`; `bt_no_reductions` vs `bt_greedy` | Paired cost differences, wins/ties/losses; runtime as a secondary measure |
| How does Supernova compare with other approaches? | `sc_greedy`, `sc_double`, `gwc`, `mst_connect`, `mst_connect_local_search`, `bt_greedy`, `bt_ilp` | Cost versus runtime, memory, completion rates |
| How do the exact approaches compare? | `sc_ilp`, `lazy_block_tree_ilp` | Runtime and memory on common solved instances, plus completion rates on all instances |

`sc_double` selects `reduction_type = "double_csr_cyc"`, which constructs the
CSR/Cyc representation and uses `CycContext`. `sc_no_reductions` provides the
ordinary CSR control. All three representation-comparison candidates run without
our reductions, trimming, or local search; GWC uses its bounds variant
(`sampling = 0`). This comparison measures the implemented pipelines, not the
context in isolation. Use stage timings additionally if isolating construction
versus solving is important.

`sc_cheapest` also has reductions, trimming, and local search disabled to make the
MSTConnect comparison about the base approaches. All five reduction pairs differ
only in `reductions`; representation and all other options are held fixed.
Trimming and set-cover local search remain disabled throughout. Oracle and the
trimming cross-product have been removed because the plan does not need them.

Lazy Supernova maps to `Lazy Block Tree ILP` and **supports reductions**: its runner
invokes `FullSingleDomReducer` when `reductions = true`. Treat BTSC with an ILP
backend as its own method; solving its intermediate ILP does not by itself make
it an exact baseline for the complete augmentation problem.

## Instance design and budget

Both stages use cycles, stars, trees, and **variable cacti only**, with
`min_cycle_size = 2` and `max_cycle_size = 16`. Size 2 adds a bridge. Keep these
bounds fixed; conclusions describe this generator, not every cactus topology.
Links are complete, with float-uniform [0,1] and integer-uniform [1,5] weights.

| Factor | Pilot | Main / server |
| --- | --- | --- |
| Node counts in every family | 20, 80 | 20, 80, 320 |
| Graph seeds for trees and variable cacti | 11, 12 | 101, 102, 103, 104, 105 |
| Link seeds in each distribution | 21, 22 | 201, 202, 203, 204, 205 |
| Independent weighted instances per family/distribution/size | 2 | 5 |
| Weighted instances | 32 | 120 |
| Algorithm configurations | 15 | 15 |
| Algorithm executions | 480 | 1,800 |
| Timeout per execution | 30 s | 60 s |
| Memory limit per execution (adjust on target machine) | 4 GiB | 4 GiB |
| Sum of execution timeout budgets | 4 h | 30 h |

The timeout sums exclude generation and overhead and are not expected runtimes.
The main experiment has 720 ILP executions (12 h of timeout budget); the remaining
1,080 heuristic executions should be budgeted using pilot measurements. Three
geometrically spaced sizes are a compact initial comparison, not a detailed
scaling curve. Add sizes only for a specific unresolved scaling question.

Use the pilot to choose meaningful small/medium/large sizes and resource limits
on the target machine. The largest size should put useful pressure on the slower
methods; the smallest should allow exact references. Keep the five reduction
pairs and the named baselines even if one loses the pilot. Do not screen out the
very controls needed to support the thesis claims. Two pilot repetitions are for
calibration, not inference. Five main repetitions are a starting budget, not a
precision guarantee; increase the planned count before the main run if needed.

Cycles and stars have one topology per size and use all link seeds. Trees and
cacti pair graph/link seeds by array position: (101,201), ..., (105,205) in the
main experiment, separately for each distribution. This samples joint topology/
weight variation without the full cross-product; it does not separate those two
variance sources. Every algorithm receives the same weighted instances. Synthetic
algorithm order is shuffled per instance using scheduling seed 42 and recorded.
Keep grids shared across algorithms and distributions within each family.

## Preflight and exact-solver interpretation

On the experiment machine, `configurations.thesis.preflight.toml` checks all 15
configurations on 8 tiny weighted instances (120 executions), with a 5 s timeout.
This checks the build and solver environment and supplies no thesis measurements.
The pilot also belongs on that machine so its timings can inform the main budget.

Lazy Supernova is the intended exact baseline. Its result files now contain
`solver.status` and `solution.optimal`; the CSV collector preserves these as
`solver_status` and `solution_optimal`. Optimality is true only for `GRB_OPTIMAL`,
including the trivial case where reductions leave at most one vertex. A returned
`SUBOPTIMAL` solution is explicitly marked false, meaning optimality was not
established. An unsuccessful solve reports its status and false without a solution
cost. External timeouts kill the process, so they may have no result/status record;
missing optimality fields (including old results) are unknown, not true or false.
Use only verified optimal results for exact-reference costs. The shared ILP helper
also accepts `GRB_SUBOPTIMAL` for SCILP; check its solver-stage status rather than
inferring an optimum from the generic completion marker alone.

## Real-world experiment

`configurations.thesis.real_world.toml` selects the `real_world` dataset directory
with a provisional 1,000-node cap, both distributions, and three link seeds.
It includes six provisional finalists: GWC, MSTConnect, MSTConnect+LS, Cyc greedy,
CSR greedy with reductions, and Lazy Supernova with reductions. After the synthetic
results, replace the three new-method slots with the fastest/useful-quality/exact
finalists; keep the GWC and MSTConnect baselines. This costs **36 executions per
selected graph** with six configurations. Do not repeat the synthetic reduction
and representation matrix on real-world data.

Set `dataset_selection.input_dir` to the target machine's dataset location and
review the selected graphs and size cap before running. The runner generates
complete weighted links and ignores existing `.links` files: this design compares
real-world-derived base graphs under the same synthetic link distributions.
It is not a comparison using native real-world candidate links. Use a separate
workflow if native links are required. Existing-dataset execution uses configuration
order rather than the synthetic shuffle; account for run order on the target host.

## Running and analysis

Build on the experiment machine and record the revision/build options, hardware,
solver version, solver thread settings, memory limits, and distribution bounds.
Run one experiment process at a time on an otherwise idle machine. All paths below
are relative to the repository root. Use fresh result directories for this revised
design: some configuration names were reused with changed parameters.

```sh
uv run scripts/experiment_runner/run_configurations.py results/thesis-v2-preflight \
  -c scripts/experiment_runner/configurations.thesis.preflight.toml

uv run scripts/experiment_runner/run_configurations.py results/thesis-v2-pilot \
  -c scripts/experiment_runner/configurations.thesis.pilot.toml

# After calibrating and freezing the main setup on the experiment machine:
uv run scripts/experiment_runner/run_configurations.py results/thesis-v2-main \
  -c scripts/experiment_runner/configurations.thesis.main.toml

# After selecting real-world finalists and input graphs:
uv run scripts/experiment_runner/run_configurations.py results/thesis-v2-real-world \
  -c scripts/experiment_runner/configurations.thesis.real_world.toml

uv run scripts/visualization/collect_results.py results/thesis-v2-main results/thesis-v2-main.csv
```

`runner.no_repeat = true` skips completed results but retries timeouts and failures. Use a fresh
output directory whenever seeds, sizes, parameters, limits, or selected methods
change. Fixed synthetic outcomes, including failures and timeouts, are recorded in
each algorithm's `synthetic-fixed.json`. The grid continues after timeouts. The
CSV contains available result records; use the JSON outcomes as the denominator
for synthetic completion/failure rates. Expected synthetic algorithm timeouts do
not make the runner exit unsuccessfully; inspect the reports and counts.

For speed comparisons, use `run.total_time_seconds`/CSV `runtime_seconds`, which
includes reductions, construction, and solving. Lazy's internal `time_total`
starts after reductions and would misrepresent their total cost. Use CSV
`peak_memory_mb` for completed runs; verify the collector's 4,096-byte page-size
assumption against the experiment host. Timed-out runs do not currently have a
complete peak-memory measurement; do not substitute zero.

Report each family/distribution separately. Match results by dataset, instance,
and link configuration. For reductions, plot paired runtime/memory ratios and
cost differences; cost can increase, decrease, or tie for heuristics. For overall
comparison, plot runtime and solution quality by size and include completion rates.
Calculate speedups on shared solved instances and show how many were excluded.
Use exact-reference gaps only for verified optima; label other references as best
known. Do not count a timeout as a runtime observation or a completion as optimal
without checking. Distinguish repeated timings from independent instance seeds.
Keep pilot measurements out of the final comparison.
