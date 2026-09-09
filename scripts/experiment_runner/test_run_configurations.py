import importlib.util
from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock


TQDM = types.ModuleType("tqdm")
TQDM.tqdm = object
sys.modules.setdefault("tqdm", TQDM)
MODULE_PATH = Path(__file__).with_name("run_configurations.py")
SPEC = importlib.util.spec_from_file_location("run_configurations", MODULE_PATH)
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


class OutputFormattingTest(unittest.TestCase):
    def test_compacts_long_labels(self):
        label = RUNNER.compact_label("abcdefghijklmnopqrstuvwxyz", 9)

        self.assertEqual(len(label), 9)
        self.assertEqual(label, "abcd…wxyz")

    def test_formats_durations(self):
        self.assertEqual(RUNNER.format_duration(7), "7s")
        self.assertEqual(RUNNER.format_duration(67), "1m 07s")
        self.assertEqual(RUNNER.format_duration(3667), "1h 01m 07s")


class BudgetTest(unittest.TestCase):
    def load(self, budget):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.toml"
            path.write_text(budget + '''
[dataset_selection]
input_dir = "inputs"
[link_configurations.unit]
distribution = "constant"
seed = 42
[configurations.inherited]
algorithm = "test"
[configurations.override]
algorithm = "test"
timeout = 120
[configuration_matrices.inherited]
algorithms = ["test"]
[configuration_matrices.override]
algorithms = ["test"]
max_memory_mb = 2048
''')
            return RUNNER.load_configurations(path)[0]

    def test_shared_limits_and_individual_overrides(self):
        configs = self.load("[budget]\ntimeout = 60\nmax_memory_mb = 4096\n")
        for name in ("inherited", "inherited_test"):
            self.assertEqual(configs[name]["timeout"], 60)
            self.assertEqual(configs[name]["max_memory_mb"], 4096)
        self.assertEqual(configs["override"]["timeout"], 120)
        self.assertEqual(configs["override"]["max_memory_mb"], 4096)
        self.assertEqual(configs["override_test"]["timeout"], 60)
        self.assertEqual(configs["override_test"]["max_memory_mb"], 2048)

    def test_optional_and_partial_budget(self):
        for budget in ("", "[budget]\n", "[budget]\ntimeout = 60\n"):
            with self.subTest(budget=budget):
                configs = self.load(budget)
                self.assertNotIn("max_memory_mb", configs["inherited"])
                self.assertEqual(configs["override"]["timeout"], 120)
                self.assertEqual(configs["override_test"]["max_memory_mb"], 2048)
                if "timeout" in budget:
                    self.assertEqual(configs["inherited"]["timeout"], 60)
                else:
                    self.assertNotIn("timeout", configs["inherited"])

    def test_rejects_invalid_budget(self):
        invalid = ['budget = 60\n', '[budget]\ntimeot = 60\n']
        for key in ("timeout", "max_memory_mb"):
            for value in ('0', '-1', 'true', '1.5', '"60"', '[]', '{}'):
                invalid.append(f"[budget]\n{key} = {value}\n")
        for budget in invalid:
            with self.subTest(budget=budget):
                with self.assertRaisesRegex(ValueError, "budget"):
                    self.load(budget)


class ResumeSettingsTest(unittest.TestCase):
    def setUp(self):
        self.previous = {
            "runner": {"no_repeat": True, "scheduling_seed": 42},
            "configurations": {"a": {
                "algorithm": "test", "timeout": 60, "max_memory_mb": 4096,
                "params": {"time_limit_seconds": 50},
            }},
            "synthetic_datasets": {"trees": {
                "seeds": [42], "generation_timeout": 300, "generation_max_memory_mb": 4096,
            }},
        }

    def test_only_external_budget_increases_are_allowed(self):
        for section, key in (
            ("configurations", "timeout"), ("configurations", "max_memory_mb"),
            ("synthetic_datasets", "generation_timeout"),
            ("synthetic_datasets", "generation_max_memory_mb"),
        ):
            for factor in (0.5, 1, 2):
                with self.subTest(section=section, key=key, factor=factor):
                    current = deepcopy(self.previous)
                    settings = next(iter(current[section].values()))
                    settings[key] = int(settings[key] * factor)
                    original = deepcopy(self.previous)
                    self.assertEqual(RUNNER.compatible_resume_settings(self.previous, current), factor >= 1)
                    self.assertEqual(self.previous, original)

    def test_removing_a_limit_is_an_increase_but_adding_one_is_a_decrease(self):
        unlimited = deepcopy(self.previous)
        del unlimited["configurations"]["a"]["max_memory_mb"]
        self.assertTrue(RUNNER.compatible_resume_settings(self.previous, unlimited))
        self.assertFalse(RUNNER.compatible_resume_settings(unlimited, self.previous))

    def test_other_changes_are_rejected_even_with_a_larger_budget(self):
        for change in ("internal_limit", "seeds", "new_algorithm", "removed_algorithm", "scheduling"):
            with self.subTest(change=change):
                current = deepcopy(self.previous)
                current["configurations"]["a"]["timeout"] = 120
                if change == "internal_limit":
                    current["configurations"]["a"]["params"]["time_limit_seconds"] = 100
                elif change == "seeds":
                    current["synthetic_datasets"]["trees"]["seeds"] = [43]
                elif change == "new_algorithm":
                    current["configurations"]["b"] = deepcopy(current["configurations"]["a"])
                elif change == "removed_algorithm":
                    del current["configurations"]["a"]
                else:
                    current["runner"]["scheduling_seed"] = 43
                self.assertFalse(RUNNER.compatible_resume_settings(self.previous, current))
        current = deepcopy(self.previous)
        current["runner"]["no_repeat"] = False
        self.assertTrue(RUNNER.compatible_resume_settings(self.previous, current))


class DatasetSelectionTest(unittest.TestCase):
    def test_configured_input_paths_and_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.toml"
            base = '''
[configurations.test]
algorithm = "test"
[link_configurations.unit]
distribution = "constant"
seed = 42
[dataset_selection]
'''
            for value, expected in (
                ("inputs", root / "inputs"),
                (str(root / "absolute"), root / "absolute"),
                ("~/datasets", Path.home() / "datasets"),
            ):
                with self.subTest(value=value):
                    config.write_text(base + f"input_dir = {json.dumps(value)}\n")
                    _, _, selection, _ = RUNNER.load_configurations(config)
                    self.assertEqual(selection["input_dir"], str(expected.resolve()))
            for setting in ('', 'input_dir = ""', 'input_dir = " "', 'input_dir = 42'):
                with self.subTest(setting=setting):
                    config.write_text(base + setting)
                    with self.assertRaisesRegex(ValueError, "dataset_selection.input_dir"):
                        RUNNER.load_configurations(config)

    def test_selects_instances_between_minimum_and_maximum_node_count(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_dir = Path(temp_dir)
            cycles_dir = input_dir / "cycles"
            cycles_dir.mkdir()
            for node_count in (999, 1000, 1001):
                (cycles_dir / f"cycle_{node_count}.xml").write_text("<graphml/>")
                (cycles_dir / f"cycle_{node_count}.graph").write_text(
                    f"{node_count} 0\n"
                )

            instances, skipped = RUNNER.find_instances(
                input_dir,
                {
                    "include": ["cycles"],
                    "datasets": {
                        "cycles": {"min_nodes": 1000, "max_nodes": 1000}
                    },
                },
            )

        self.assertEqual([instance[3] for instance in instances], [1000])
        self.assertEqual(skipped, 2)


class ArgumentsTest(unittest.TestCase):
    def test_only_output_directory_is_positional(self):
        with mock.patch.object(sys, "argv", [
            "run_configurations.py", "results", "-c", "config.toml",
        ]):
            arguments = RUNNER.parse_arguments()
        self.assertEqual(arguments.output_dir, Path("results"))
        self.assertEqual(arguments.configurations, Path("config.toml"))
        self.assertEqual(set(vars(arguments)), {"output_dir", "configurations"})

    def test_runner_settings_are_configured_and_validate_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / "config.toml"
            config.write_text('''[runner]
mode = "adaptive"
adaptive_probes = 3
no_repeat = true
scheduling_seed = 7
experiments_binary = "bin/experiments"
''')
            settings = RUNNER.load_runner_settings(config)
            self.assertEqual(settings["mode"], "adaptive")
            self.assertEqual(settings["adaptive_probes"], 3)
            self.assertTrue(settings["no_repeat"])
            self.assertEqual(settings["scheduling_seed"], 7)
            self.assertEqual(settings["experiments_binary"], str(config.parent / "bin/experiments"))
            for invalid in ('mode = "other"', 'no_repeat = "yes"', 'adaptive_probes = 0',
                            'scheduling_seed = true', 'experiments_binary = 4', 'typo = 1'):
                with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                    config.write_text('[runner]\n' + invalid)
                    RUNNER.load_runner_settings(config)


class AdaptiveSchedulingTest(unittest.TestCase):
    def test_probes_grow_exponentially_and_include_largest_level(self):
        indices = RUNNER.exponential_probe_indices(10, 4)

        self.assertEqual(indices, [0, 1, 3, 9])

    def test_fills_every_size_below_the_frontier(self):
        sizes = list(range(100, 1100, 100))
        calls = []

        def run_level(size):
            calls.append(size)
            return "solved" if size <= 600 else "timeout"

        outcomes = RUNNER.run_adaptive_series(sizes, 4, run_level)

        self.assertTrue(all(outcomes[size] == "solved" for size in sizes[:6]))
        self.assertEqual(outcomes[700], "timeout")
        self.assertEqual(outcomes[1000], "timeout")
        self.assertNotIn(800, calls)
        self.assertNotIn(900, calls)

    def test_all_levels_are_filled_when_the_largest_probe_solves(self):
        sizes = list(range(100, 1100, 100))
        calls = []

        def run_level(size):
            calls.append(size)
            return "solved"

        outcomes = RUNNER.run_adaptive_series(sizes, 4, run_level)

        self.assertEqual(set(calls), set(sizes))
        self.assertEqual(set(outcomes), set(sizes))


class LinkConfigurationTest(unittest.TestCase):
    def test_expands_seed_repetitions_into_named_results(self):
        configurations = {
            "uniform": {
                "distribution": "float_uniform",
                "float_uniform_lower": 0.0,
                "float_uniform_upper": 1.0,
                "seeds": [42, 43],
            }
        }

        expanded = RUNNER.expand_link_configurations(configurations)

        self.assertEqual(list(expanded), ["uniform_seed_42", "uniform_seed_43"])
        self.assertEqual(expanded["uniform_seed_42"]["seed"], 42)
        self.assertEqual(expanded["uniform_seed_43"]["_family"], "uniform")
        self.assertEqual(
            RUNNER.group_link_configurations(expanded),
            {"uniform": ["uniform_seed_42", "uniform_seed_43"]},
        )

    def test_keeps_legacy_single_seed_name(self):
        configurations = {
            "unit": {"distribution": "constant", "seed": 42}
        }

        expanded = RUNNER.expand_link_configurations(configurations)

        self.assertEqual(list(expanded), ["unit"])

    def test_all_algorithms_configuration_loads(self):
        configuration_path = MODULE_PATH.with_name("configurations.all.toml")

        configurations, links, datasets, synthetic = RUNNER.load_configurations(
            configuration_path
        )

        self.assertEqual(len(configurations), 76)
        self.assertEqual(synthetic, {})
        self.assertEqual(len(links), 10)
        self.assertIn("cacti", datasets["include"])
        algorithms = {
            configuration["algorithm"] for configuration in configurations.values()
        }
        self.assertTrue(
            {
                "Set Cover CSR",
                "Block Tree Set Cover CSR",
                "Lazy Block Tree ILP",
                "WheelCon",
                "mst-connect",
                "mst-connect-ls",
            }.issubset(algorithms)
        )
        self.assertNotIn("gwc", algorithms)
        self.assertNotIn("eilp", algorithms)

        set_cover_configurations = [
            configuration
            for configuration in configurations.values()
            if configuration["algorithm"]
            in {"Set Cover CSR", "Block Tree Set Cover CSR"}
        ]
        self.assertEqual(len(set_cover_configurations), 72)
        combinations = {
            (
                configuration["algorithm"],
                configuration["params"]["reduction_type"],
                configuration["params"]["solver"],
                configuration["params"]["reductions"],
                configuration["params"]["trimming"],
            )
            for configuration in set_cover_configurations
        }
        self.assertEqual(len(combinations), 72)
        self.assertEqual(
            {
                configuration["params"]["reduction_type"]
                for configuration in set_cover_configurations
            },
            {"csr", "oracle", "double_csr_cyc"},
        )
        self.assertEqual(
            {
                configuration["params"]["solver"]
                for configuration in set_cover_configurations
            },
            {"greedy", "greedy_cheapest", "ilp"},
        )
        self.assertEqual(
            {
                configuration["params"]["reductions"]
                for configuration in set_cover_configurations
            },
            {False, True},
        )
        self.assertEqual(
            {
                configuration["params"]["trimming"]
                for configuration in set_cover_configurations
            },
            {False, True},
        )


class ConfigurationMatrixTest(unittest.TestCase):
    def test_expands_cartesian_product_with_fixed_parameters(self):
        matrices = {
            "cover": {
                "algorithms": ["Set Cover CSR", "Block Tree Set Cover CSR"],
                "timeout": 120,
                "params": {
                    "solver": ["greedy", "ilp"],
                    "trimming": [False, True],
                    "local_search": False,
                },
            }
        }

        configurations = RUNNER.expand_configuration_matrices(matrices)

        self.assertEqual(len(configurations), 8)
        self.assertEqual(
            {configuration["algorithm"] for configuration in configurations.values()},
            {"Set Cover CSR", "Block Tree Set Cover CSR"},
        )
        for configuration in configurations.values():
            self.assertEqual(configuration["timeout"], 120)
            self.assertFalse(configuration["params"]["local_search"])


class ExperimentResultTest(unittest.TestCase):
    def test_resume_requires_a_complete_result(self):
        with tempfile.TemporaryDirectory() as directory:
            result = Path(directory) / "res.txt"
            self.assertFalse(RUNNER.has_complete_result(result))
            result.write_text("run.peak_memory_pages=100\n")
            self.assertFalse(RUNNER.has_complete_result(result))
            result.write_text("run.total_time_seconds=0.5\n")
            self.assertTrue(RUNNER.has_complete_result(result))

    def test_partial_child_output_is_not_treated_as_completed(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            result_file = output_dir / "res.txt"

            def write_partial_result(*_):
                result_file.write_text("run.peak_memory_pages=100\n")
                return 0

            with mock.patch.object(
                RUNNER, "run_command", side_effect=write_partial_result
            ), mock.patch.object(RUNNER.sys, "stderr"):
                outcome = RUNNER.run_experiment(
                    Path("experiments"),
                    "test",
                    {"algorithm": "Set Cover CSR"},
                    Path("graph.xml"),
                    Path("links.links"),
                    output_dir,
                    "res.txt",
                )

            self.assertEqual(outcome, "failed")


class SyntheticSchedulingTest(unittest.TestCase):
    def test_exact_cycle_count_requires_a_feasible_integer(self):
        for cycles in (1, 4, 9):
            settings = self.settings(generator="cactus_cycles", cycles=cycles)
            self.assertEqual(settings["cycles"], cycles)
        for cycles in (None, -1, True, 1.5, 0, 10):
            with self.subTest(cycles=cycles), self.assertRaises(ValueError):
                self.settings(generator="cactus_cycles", cycles=cycles)
        with self.assertRaises(ValueError):
            self.settings(generator="tree", cycles=1)
        settings = RUNNER.validate_synthetic_datasets({
            "test": {"generator": "cactus_cycles", "sizes": [21, 30], "cycles": 10}
        })["test"]
        self.assertEqual(settings["cycles"], 10)
        with self.assertRaises(ValueError):
            RUNNER.validate_synthetic_datasets({
                "test": {"generator": "cactus_cycles", "sizes": [2, 30], "cycles": 2}
            })

    def test_zero_count_allows_only_single_node_datasets(self):
        for sizes in ([1], [1, 2]):
            document = {"test": {"generator": "cactus_cycles", "sizes": sizes, "cycles": 0}}
            if sizes == [1]:
                self.assertEqual(RUNNER.validate_synthetic_datasets(document)["test"]["cycles"], 0)
            else:
                with self.assertRaises(ValueError):
                    RUNNER.validate_synthetic_datasets(document)
        self.settings(generator="cactus_cycles", min_nodes=1, max_nodes=1, cycles=0)
        with self.assertRaises(ValueError):
            self.settings(generator="cactus_cycles", min_nodes=1, cycles=0)

    def test_variable_cactus_requires_valid_bounds(self):
        settings = self.settings(generator="cactus_variable", min_cycle_size=2, max_cycle_size=8)
        self.assertEqual(settings["min_cycle_size"], 2)
        self.assertEqual(settings["max_cycle_size"], 8)
        self.settings(generator="cactus_variable", min_cycle_size=2, max_cycle_size=2)
        for bounds in (
            {}, {"min_cycle_size": 3}, {"min_cycle_size": 1, "max_cycle_size": 5},
            {"min_cycle_size": 5, "max_cycle_size": 3},
            {"min_cycle_size": 3.5, "max_cycle_size": 8},
            {"min_cycle_size": 3, "max_cycle_size": 8, "cycle_mass": 0.5},
        ):
            with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                self.settings(generator="cactus_variable", **bounds)
        with self.assertRaises(ValueError):
            self.settings(generator="tree", min_cycle_size=3, max_cycle_size=8)

    def settings(self, **overrides):
        return RUNNER.validate_synthetic_datasets({
            "test": {"generator": "tree", "min_nodes": 10, "samples": 5, **overrides}
        })["test"]

    def test_searches_beyond_initial_size_then_bisects_and_samples(self):
        calls = []

        def run(size):
            calls.append(size)
            return "solved" if size <= 73 else "timeout"

        report = RUNNER.run_synthetic_series(self.settings(), run)
        self.assertEqual(calls[:4], [10, 20, 40, 80])
        self.assertEqual(report["largest_solved"], 73)
        self.assertEqual(report["smallest_timeout"], 74)
        self.assertEqual(report["status"], "bracketed")
        self.assertEqual(report["samples"], [10, 25, 41, 57, 73])
        self.assertEqual(len(calls), len(set(calls)))
        self.assertEqual(calls[-3:], [57, 41, 25])

    def test_timeout_at_minimum_has_no_solved_frontier(self):
        run = mock.Mock(return_value="timeout")
        report = RUNNER.run_synthetic_series(self.settings(), run)
        run.assert_called_once_with(10)
        self.assertIsNone(report["largest_solved"])
        self.assertEqual(report["samples"], [])

    def test_cap_is_reported_as_a_lower_bound(self):
        report = RUNNER.run_synthetic_series(
            self.settings(max_nodes=35), lambda size: "solved"
        )
        self.assertEqual(report["status"], "capped")
        self.assertEqual(report["largest_solved"], 35)
        self.assertIsNone(report["smallest_timeout"])
        self.assertLessEqual(max(report["outcomes"]), 35)

    def test_resolution_bounds_the_bracket(self):
        report = RUNNER.run_synthetic_series(
            self.settings(resolution=8),
            lambda size: "solved" if size <= 73 else "timeout",
        )
        self.assertLessEqual(report["smallest_timeout"] - report["largest_solved"], 8)

    def test_generation_failure_is_not_a_timeout_boundary(self):
        report = RUNNER.run_synthetic_series(
            self.settings(), lambda size: "solved" if size <= 40 else "generation_failed"
        )
        self.assertEqual(report["status"], "generation_failed")
        self.assertEqual(report["largest_solved"], 40)
        self.assertIsNone(report["smallest_timeout"])

    def test_failure_during_bisection_leaves_search_incomplete(self):
        def run(size):
            if size == 60:
                return "failed"
            return "solved" if size <= 73 else "timeout"

        report = RUNNER.run_synthetic_series(self.settings(), run)
        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["largest_solved"], 40)
        self.assertEqual(report["smallest_timeout"], 80)

    def test_backfill_detects_non_monotone_outcomes_and_keeps_sampling(self):
        report = RUNNER.run_synthetic_series(
            self.settings(max_nodes=50),
            lambda size: "timeout" if size == 30 else "solved",
        )
        self.assertEqual(report["status"], "non_monotone")
        self.assertEqual(report["smallest_timeout"], 30)
        self.assertTrue(set(report["samples"]).issubset(report["outcomes"]))

    def test_small_ranges_and_single_sample_do_not_duplicate_sizes(self):
        self.assertEqual(RUNNER.evenly_spaced_sizes(10, 12, 10), [10, 11, 12])
        self.assertEqual(RUNNER.evenly_spaced_sizes(10, 12, 1), [12])
        self.assertEqual(RUNNER.evenly_spaced_sizes(10, 10, 10), [10])

    def test_seed_rules_and_invalid_settings(self):
        for generator in ("cycle", "star"):
            self.assertEqual(self.settings(generator=generator)["seeds"], [42])
            with self.assertRaises(ValueError):
                self.settings(generator=generator, seeds=[42, 43])
        self.assertEqual(self.settings(seeds=[42, 43])["seeds"], [42, 43])
        for settings in (
            {"seeds": []}, {"seeds": [42, 42]}, {"seeds": [True]},
            {"seeds": [2**32]}, {"samples": 0}, {"resolution": 0},
            {"max_nodes": 5}, {"generator": "cycle", "min_nodes": 2},
            {"generator": "cactus", "cycle_length": 2, "cycle_mass": 0.5},
            {"cycle_mass": 0.5}, {"typo": 1},
        ):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                self.settings(**settings)


class ScalabilitySchedulingTest(unittest.TestCase):
    def settings(self, **changes):
        return RUNNER.validate_synthetic_datasets({
            "trees": dict(generator="tree", mode="scalability", min_nodes=4) | changes,
        })["trees"]

    def test_mixed_failures_continue_and_algorithms_stop_independently(self):
        calls = []

        def run(size, names):
            calls.append((size, list(names)))
            return {name: [
                {"status": "timeout"},
                {"status": "completed" if size < {"a": 8, "b": 32}[name] else "failed"},
            ] for name in names}

        checkpoint = mock.Mock()
        reports = RUNNER.run_scalability_series(self.settings(), ["a", "b"], run, checkpoint)
        self.assertEqual(calls, [(4, ["a", "b"]), (8, ["a", "b"]), (16, ["b"]), (32, ["b"])])
        self.assertEqual(checkpoint.call_count, 4)
        for name, last, stop in (("a", 4, 8), ("b", 16, 32)):
            self.assertEqual(reports[name]["status"], "all_failed")
            self.assertEqual(reports[name]["largest_successful_size"], last)
            self.assertEqual(reports[name]["stop_size"], stop)
            self.assertEqual(reports[name]["samples"], [])
            self.assertEqual(reports[name]["levels"][0]["succeeded"], 1)

    def test_fractional_growth_and_cap(self):
        run = mock.Mock(side_effect=lambda size, names: {name: [{"status": "completed"}] for name in names})
        report = RUNNER.run_scalability_series(
            self.settings(growth_factor=1.5, max_nodes=11), ["a"], run, mock.Mock(),
        )["a"]
        self.assertEqual([call.args[0] for call in run.call_args_list], [4, 6, 9, 11])
        self.assertEqual(report["status"], "capped")
        self.assertIsNone(report["stop_size"])

    def test_optional_intermediate_sizes_are_measured_once_even_if_they_fail(self):
        run = mock.Mock(side_effect=lambda size, names: {
            name: [{"status": "completed" if size in (4, 8, 16) else "timeout"}]
            for name in names
        })
        report = RUNNER.run_scalability_series(
            self.settings(fill_intermediate=True, samples=5), ["a"], run, mock.Mock(),
        )["a"]
        self.assertEqual([call.args[0] for call in run.call_args_list], [4, 8, 16, 32, 7, 10, 13])
        self.assertEqual(report["samples"], [4, 7, 10, 13, 16])
        self.assertEqual([level["phase"] for level in report["levels"]], ["growth"] * 4 + ["intermediate"] * 3)

    def test_generation_failure_and_failure_at_start(self):
        for status, expected in (("generation_failed", "generation_failed"), ("timeout", "all_failed")):
            run = mock.Mock(return_value={"a": [{"status": status}]})
            report = RUNNER.run_scalability_series(
                self.settings(fill_intermediate=True), ["a"], run, mock.Mock(),
            )["a"]
            self.assertEqual(run.call_count, 1)
            self.assertEqual(report["status"], expected)
            self.assertIsNone(report["largest_successful_size"])
            self.assertEqual(report["stop_size"], 4)
            self.assertEqual(report["samples"], [])

    def test_invalid_settings_and_paired_seeds(self):
        for changes in (
            {"growth_factor": 1}, {"growth_factor": True}, {"growth_factor": float("inf")},
            {"growth_factor": float("nan")}, {"growth_factor": "2"}, {"fill_intermediate": 1},
            {"resolution": 2}, {"sizes": [4, 8]}, {"mode": "unknown"},
            {"mode": "frontier", "growth_factor": 2},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.settings(**changes)
        settings = self.settings(seed_mode="paired", seeds=[101, 102])
        links = RUNNER.expand_link_configurations({"unit": dict(distribution="constant", seeds=[201, 202])})
        self.assertEqual(RUNNER.synthetic_seed_groups(settings, links), {
            101: ["unit_seed_201"], 102: ["unit_seed_202"],
        })


class SyntheticIntegrationTest(unittest.TestCase):
    def test_increased_budget_resumes_successes_retries_failures_and_grows_further(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.toml"
            config.write_text(f'''[runner]
no_repeat = true
experiments_binary = {json.dumps(sys.executable)}
dataset_generator_binary = {json.dumps(sys.executable)}
[budget]
timeout = 1
max_memory_mb = 128
[synthetic_datasets.trees]
generator = "tree"
mode = "scalability"
min_nodes = 4
seed_mode = "paired"
seeds = [101, 102]
[link_configurations.uniform]
distribution = "float_uniform"
seeds = [201, 202]
[configurations.a]
algorithm = "a"
''')
            arguments = types.SimpleNamespace(output_dir=root / "results", configurations=config)

            def generate(executable, output, size, seed, settings):
                graph = output / f"tree_{size}_seed_{seed}.xml"
                graph.write_text("<graphml/>")
                return graph, graph.with_suffix(".graph")

            def experiment(executable, name, settings, graph, links, output, result_name):
                size = int(graph.stem.split("_")[1])
                threshold = 16 if settings["timeout"] == 2 and settings["max_memory_mb"] == 256 else 8
                if graph.stem.endswith("101") or size >= threshold:
                    return "timeout"
                (output / result_name).write_text("run.total_time_seconds=0.01\n")
                return "completed"

            with mock.patch.object(RUNNER, "parse_arguments", return_value=arguments), \
                 mock.patch.object(RUNNER, "generate_synthetic_graph", side_effect=generate), \
                 mock.patch.object(RUNNER, "run_command", return_value=0), \
                 mock.patch.object(RUNNER, "run_experiment", side_effect=experiment) as run, \
                 mock.patch.object(RUNNER, "tqdm"), mock.patch("builtins.print"):
                self.assertEqual(RUNNER.main(), 0)
                self.assertEqual(run.call_count, 4)
                result = root / "results/a/trees/tree_4_seed_102/res-uniform_seed_202.txt"
                original_result = result.read_bytes()
                config.write_text(config.read_text().replace("timeout = 1", "timeout = 2")
                                  .replace("max_memory_mb = 128", "max_memory_mb = 256"))
                run.reset_mock()
                self.assertEqual(RUNNER.main(), 0)
                self.assertEqual(run.call_count, 5)
                self.assertNotIn("tree_4_seed_102", [call.args[3].stem for call in run.call_args_list])
                self.assertEqual(result.read_bytes(), original_result)
                report = json.loads((root / "results/a/synthetic-scalability.json").read_text())[0]
                self.assertEqual(report["stop_size"], 16)
                self.assertEqual(report["largest_successful_size"], 8)
                manifest_path = root / "results/run-settings.json"
                saved_manifest = manifest_path.read_text()
                settings = json.loads(saved_manifest)["configurations"]["a"]
                self.assertEqual((settings["timeout"], settings["max_memory_mb"]), (2, 256))
                self.assertEqual((root / "results/configuration.source.toml").read_text(), config.read_text())
                metadata = (root / "results/a/configuration.txt").read_text()
                self.assertIn("timeout=2\n", metadata)
                self.assertIn("max_memory_mb=256\n", metadata)
                config.write_text(config.read_text().replace("timeout = 2", "timeout = 1"))
                run.reset_mock()
                with self.assertRaisesRegex(ValueError, "only increased external"):
                    RUNNER.main()
                run.assert_not_called()
                self.assertEqual(manifest_path.read_text(), saved_manifest)

    def test_scalability_attempts_remaining_seeds_after_generation_failures(self):
        for stage in ("graph", "links"):
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                config = root / "config.toml"
                config.write_text(f'''[runner]
experiments_binary = {json.dumps(sys.executable)}
dataset_generator_binary = {json.dumps(sys.executable)}
[synthetic_datasets.trees]
generator = "tree"
mode = "scalability"
min_nodes = 4
seed_mode = "paired"
seeds = [101, 102]
[link_configurations.uniform]
distribution = "float_uniform"
seeds = [201, 202]
[configurations.a]
algorithm = "a"
timeout = 1
''')

                def generate(executable, output, size, seed, settings):
                    if stage == "graph" and (size == 8 or seed == 101):
                        return None
                    graph = output / f"tree_{size}_seed_{seed}.xml"
                    graph.write_text("<graphml/>")
                    return graph, graph.with_suffix(".graph")

                def generate_links(command, *args):
                    graph = command[command.index("--input_graph") + 1]
                    return int(stage == "links" and ("tree_8_" in graph.stem or graph.stem.endswith("101")))

                arguments = types.SimpleNamespace(output_dir=root / "results", configurations=config)
                with mock.patch.object(RUNNER, "parse_arguments", return_value=arguments), \
                     mock.patch.object(RUNNER, "generate_synthetic_graph", side_effect=generate) as generation, \
                     mock.patch.object(RUNNER, "run_command", side_effect=generate_links), \
                     mock.patch.object(RUNNER, "run_experiment", return_value="completed") as run, \
                     mock.patch.object(RUNNER, "tqdm"), mock.patch("builtins.print"):
                    self.assertEqual(RUNNER.main(), 1)
                    self.assertEqual(generation.call_count, 4)
                    self.assertEqual(run.call_count, 1)
                    report = json.loads((root / "results/a/synthetic-scalability.json").read_text())[0]
                    self.assertEqual(report["status"], "generation_failed")
                    self.assertEqual(report["stop_size"], 8)
                    self.assertEqual([(level["succeeded"], level["total"]) for level in report["levels"]], [(1, 2), (0, 2)])

    def test_scalability_shares_inputs_records_all_attempts_and_resumes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.toml"
            config.write_text(f'''[runner]
no_repeat = true
experiments_binary = {json.dumps(sys.executable)}
dataset_generator_binary = {json.dumps(sys.executable)}
[synthetic_datasets.trees]
generator = "tree"
mode = "scalability"
min_nodes = 4
seed_mode = "paired"
seeds = [101, 102]
[link_configurations.uniform]
distribution = "float_uniform"
seeds = [201, 202]
[link_configurations.integer]
distribution = "integer_uniform"
seeds = [201, 202]
[configurations.a]
algorithm = "a"
timeout = 1
[configurations.b]
algorithm = "b"
timeout = 1
''')
            arguments = types.SimpleNamespace(output_dir=root / "results", configurations=config)

            def generate(executable, output, size, seed, settings):
                graph = output / f"tree_{size}_seed_{seed}.xml"
                graph.write_text("<graphml/>")
                return graph, graph.with_suffix(".graph")

            def experiment(executable, name, settings, graph, links, output, result_name):
                size = int(graph.stem.split("_")[1])
                threshold = 8 if name == "a" or "integer" in result_name else 16
                if graph.stem.endswith("101"):
                    return "timeout"
                if size >= threshold:
                    return "failed"
                (output / result_name).write_text("run.total_time_seconds=0.01\n")
                return "completed"

            with mock.patch.object(RUNNER, "parse_arguments", return_value=arguments), \
                 mock.patch.object(RUNNER, "generate_synthetic_graph", side_effect=generate) as generation, \
                 mock.patch.object(RUNNER, "run_command", return_value=0) as links, \
                 mock.patch.object(RUNNER, "run_experiment", side_effect=experiment) as run, \
                 mock.patch.object(RUNNER, "tqdm"), mock.patch("builtins.print"):
                self.assertEqual(RUNNER.main(), 1)
                self.assertEqual(run.call_count, 18)
                self.assertEqual(generation.call_count, 10)
                self.assertEqual(links.call_count, 10)
                for name in ("a", "b"):
                    reports = json.loads((root / "results" / name / "synthetic-scalability.json").read_text())
                    self.assertEqual(len(reports), 2)
                    for report in reports:
                        stop = 16 if name == "b" and report["link_configuration"] == "uniform" else 8
                        self.assertEqual(report["status"], "all_failed")
                        self.assertEqual(report["stop_size"], stop)
                        self.assertEqual(report["levels"][0]["succeeded"], 1)
                        for level in report["levels"]:
                            self.assertEqual(level["total"], 2)
                            self.assertEqual({(row["graph_seed"], row["link_seed"]) for row in level["instances"]},
                                             {(101, 201), (102, 202)})
                self.assertEqual((root / "results/configuration.source.toml").read_text(), config.read_text())
                manifest = json.loads((root / "results/run-settings.json").read_text())
                self.assertEqual(manifest["synthetic_datasets"]["trees"]["growth_factor"], 2)
                run.reset_mock()
                self.assertEqual(RUNNER.main(), 1)
                self.assertEqual(run.call_count, 13)
                config.write_text(config.read_text().replace("min_nodes = 4", "min_nodes = 5"))
                run.reset_mock()
                with self.assertRaisesRegex(ValueError, "different experiment settings"):
                    RUNNER.main()
                run.assert_not_called()

    def test_fixed_sizes_pair_seeds_and_continue_after_timeout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.toml"
            config.write_text('''
[synthetic_datasets.trees]
generator = "tree"
sizes = [4, 8]
seeds = [101, 102]
seed_mode = "paired"
[synthetic_datasets.stars]
generator = "star"
sizes = [4, 8]
seed_mode = "paired"
[link_configurations.uniform]
distribution = "float_uniform"
seeds = [201, 202]
[configurations.a]
algorithm = "a"
timeout = 1
[configurations.b]
algorithm = "b"
timeout = 1
''')
            config.write_text(
                "[runner]\nno_repeat = true\n"
                + f"experiments_binary = {json.dumps(sys.executable)}\n"
                + f"dataset_generator_binary = {json.dumps(sys.executable)}\n"
                + config.read_text()
            )
            arguments = types.SimpleNamespace(
                output_dir=root / "results", configurations=config,
                mode="exhaustive", adaptive_probes=6, no_repeat=True,
            )

            def generate(executable, output, size, seed, settings):
                graph = output / f'{settings["generator"]}_{size}_seed_{seed}.xml'
                graph.write_text("<graphml/>")
                return graph, graph.with_suffix(".graph")

            def experiment(executable, name, settings, graph, links, output, result_name):
                if name == "a":
                    return "timeout"
                (output / result_name).write_text("run.total_time_seconds=0.01\n")
                return "completed"

            with mock.patch.object(RUNNER, "parse_arguments", return_value=arguments), \
                 mock.patch.object(RUNNER, "generate_synthetic_graph", side_effect=generate) as generation, \
                 mock.patch.object(RUNNER, "run_command", return_value=0) as links, \
                 mock.patch.object(RUNNER, "run_experiment", side_effect=experiment) as run, \
                 mock.patch.object(RUNNER, "run_synthetic_series") as frontier, \
                 mock.patch.object(RUNNER, "tqdm"), mock.patch("builtins.print"):
                self.assertEqual(RUNNER.main(), 0)
                self.assertEqual(run.call_count, 16)
                self.assertEqual(generation.call_count, 6)
                self.assertEqual(links.call_count, 8)
                frontier.assert_not_called()
                for name, status in (("a", "timeout"), ("b", "completed")):
                    report = json.loads((root / "results" / name / "synthetic-fixed.json").read_text())
                    self.assertEqual(len(report), 8)
                    self.assertEqual({row["status"] for row in report}, {status})
                    self.assertEqual({row["size"] for row in report}, {4, 8})
                    self.assertEqual(
                        {(row["graph_seed"], row["link_seed"]) for row in report if row["dataset"] == "trees"},
                        {(101, 201), (102, 202)},
                    )
                run.reset_mock()
                self.assertEqual(RUNNER.main(), 0)
                self.assertEqual(run.call_count, 8)
                self.assertTrue(all(call.args[1] == "a" for call in run.call_args_list))

    def test_thesis_budgets_and_variable_cacti(self):
        for stage, expected in (("preflight", 120), ("pilot", 480), ("main", 1800), ("server", 1800)):
            configurations, links, _, datasets = RUNNER.load_configurations(
                MODULE_PATH.with_name(f"configurations.thesis.{stage}.toml")
            )
            count = sum(
                len(settings["sizes"]) * sum(map(len, RUNNER.synthetic_seed_groups(settings, links).values()))
                for settings in datasets.values()
            ) * len(configurations)
            self.assertEqual(count, expected)
            self.assertEqual(set(datasets), {"cycles", "stars", "trees", "variable_cacti"})
            self.assertEqual(datasets["variable_cacti"]["generator"], "cactus_variable")
            self.assertEqual(datasets["variable_cacti"]["min_cycle_size"], 2)
            self.assertEqual(datasets["variable_cacti"]["max_cycle_size"], 16)

    def test_thesis_claims_have_matched_controls(self):
        for stage in ("preflight", "pilot", "main", "server"):
            configurations, _, _, _ = RUNNER.load_configurations(
                MODULE_PATH.with_name(f"configurations.thesis.{stage}.toml")
            )
            self.assertEqual(len(configurations), 15)
            self.assertEqual(configurations["gwc"]["algorithm"], "gwc")
            self.assertEqual(configurations["gwc"]["params"]["sampling"], 0)
            for without, with_reductions in (
                ("sc_no_reductions", "sc_greedy"),
                ("bt_no_reductions", "bt_greedy"),
                ("sc_ilp_no_reductions", "sc_ilp"),
                ("bt_ilp_no_reductions", "bt_ilp"),
                ("lazy_block_tree_ilp_no_reductions", "lazy_block_tree_ilp"),
            ):
                control = configurations[without]
                reduced = configurations[with_reductions]
                self.assertFalse(control["params"]["reductions"])
                self.assertTrue(reduced["params"]["reductions"])
                self.assertEqual(
                    control | {"params": control["params"] | {"reductions": True}},
                    reduced,
                )
            csr = configurations["sc_no_reductions"]
            self.assertEqual(
                configurations["sc_double"],
                csr | {"params": csr["params"] | {"reduction_type": "double_csr_cyc"}},
            )
            self.assertEqual(
                configurations["sc_cheapest"],
                csr | {"params": csr["params"] | {"solver": "greedy_cheapest"}},
            )

    def test_thesis_real_world_selection_keeps_baselines(self):
        configurations, links, datasets, synthetic = RUNNER.load_configurations(
            MODULE_PATH.with_name("configurations.thesis.real_world.toml")
        )
        self.assertEqual(synthetic, {})
        self.assertEqual(datasets["include"], ["real_world"])
        self.assertEqual(len(configurations) * len(links), 36)
        self.assertTrue({"gwc", "mst_connect", "mst_connect_local_search"} <= configurations.keys())

    def test_invalid_fixed_designs(self):
        for changes in (
            {"sizes": []}, {"sizes": [4, 4]}, {"sizes": [8, 4]},
            {"sizes": [True]}, {"sizes": [1]}, {"min_nodes": 4},
            {"seed_mode": "zip"},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                RUNNER.validate_synthetic_datasets({
                    "trees": dict(generator="tree", sizes=[4, 8]) | changes,
                })
        settings = RUNNER.validate_synthetic_datasets({
            "trees": dict(generator="tree", sizes=[4], seed_mode="paired", seeds=[1, 2]),
        })["trees"]
        with self.assertRaisesRegex(ValueError, "equal graph and link seed counts"):
            RUNNER.synthetic_seed_groups(settings, RUNNER.expand_link_configurations({
                "uniform": dict(distribution="float_uniform", seed=3),
            }))

    def test_tables_select_existing_synthetic_or_both(self):
        for existing, synthetic in ((True, False), (False, True), (True, True)):
            with self.subTest(existing=existing, synthetic=synthetic), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                input_dir = root / "inputs"
                if existing:
                    input_dir.mkdir()
                config = root / "config.toml"
                config.write_text('''
[link_configurations.unit]
distribution = "constant"
seed = 42
[configurations.test]
algorithm = "test"
timeout = 1
''' + ('\n[dataset_selection]\ninput_dir = "inputs"\n' if existing else "") + ('''
[synthetic_datasets.trees]
generator = "tree"
max_nodes = 10
''' if synthetic else ""))
                config.write_text(
                    "[runner]\nno_repeat = true\n"
                    + f"experiments_binary = {json.dumps(sys.executable)}\n"
                    + f"dataset_generator_binary = {json.dumps(sys.executable)}\n"
                    + config.read_text()
                )
                arguments = types.SimpleNamespace(
                    output_dir=root / "results",
                    configurations=config, mode="exhaustive",
                    adaptive_probes=6, no_repeat=False,
                )
                graph = input_dir / "real_world" / "instance.xml"
                report = {"status": "capped", "largest_solved": 10, "smallest_timeout": None}
                with mock.patch.object(RUNNER, "parse_arguments", return_value=arguments), \
                     mock.patch.object(RUNNER, "find_instances", return_value=(
                         [(graph, graph.with_suffix(".graph"), "real_world", 10)], 0,
                     )) as find, \
                     mock.patch.object(RUNNER, "run_synthetic_series", return_value=report) as search, \
                     mock.patch.object(RUNNER, "run_command", return_value=0), \
                     mock.patch.object(RUNNER, "run_experiment", return_value="completed") as experiment, \
                     mock.patch.object(RUNNER, "tqdm"), \
                     mock.patch("builtins.print"):
                    self.assertEqual(RUNNER.main(), 0)
                    self.assertEqual(find.call_count, int(existing))
                    if existing:
                        self.assertEqual(find.call_args.args[0], input_dir)
                        metadata = (root / "results" / "test" / "configuration.txt").read_text()
                        self.assertIn(f"dataset_selection.input_dir={input_dir}\n", metadata)
                    self.assertEqual(experiment.call_count, int(existing))
                    self.assertEqual(search.call_count, int(synthetic))

    def test_separate_frontiers_require_all_graph_and_link_seeds(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.toml"
            config.write_text('''
[synthetic_datasets.trees]
generator = "tree"
min_nodes = 4
max_nodes = 16
samples = 3
seeds = [42, 43]
[link_configurations.uniform]
distribution = "float_uniform"
seeds = [42, 43]
[configurations.small]
algorithm = "small"
timeout = 1
[configurations.large]
algorithm = "large"
timeout = 1
''')
            config.write_text(
                "[runner]\nno_repeat = true\n"
                + f"experiments_binary = {json.dumps(sys.executable)}\n"
                + f"dataset_generator_binary = {json.dumps(sys.executable)}\n"
                + config.read_text()
            )
            arguments = types.SimpleNamespace(
                output_dir=root / "results",
                configurations=config, mode="exhaustive",
                adaptive_probes=6, no_repeat=True,
            )
            calls = []

            def generate(executable, output, size, seed, settings):
                graph = output / f"tree_{size}_seed_{seed}.xml"
                graph.write_text("<graphml/>")
                graph.with_suffix(".graph").write_text(f"{size} 0\n")
                return graph, graph.with_suffix(".graph")

            def experiment(executable, name, settings, graph, links, output, result_name):
                size = RUNNER.read_node_count(graph.with_suffix(".graph"))
                seed = int(graph.stem.split("_")[-1])
                calls.append((name, size, seed, result_name))
                if size > {"small": 7, "large": 11}[name] and seed == 43 and "43" in result_name:
                    return "timeout"
                (output / result_name).write_text("run.total_time_seconds=0.01\n")
                return "completed"

            with mock.patch.object(RUNNER, "parse_arguments", return_value=arguments), \
                 mock.patch.object(RUNNER, "generate_synthetic_graph", side_effect=generate), \
                 mock.patch.object(RUNNER, "run_command", return_value=0), \
                 mock.patch.object(RUNNER, "run_experiment", side_effect=experiment), \
                 mock.patch.object(RUNNER, "tqdm"), \
                 mock.patch("builtins.print"):
                self.assertEqual(RUNNER.main(), 0)
                for name, frontier in (("small", 7), ("large", 11)):
                    report = json.loads((root / "results" / name / "synthetic-frontier.json").read_text())[0]
                    self.assertEqual(report["largest_solved"], frontier)
                    self.assertEqual(report["smallest_timeout"], frontier + 1)
                    self.assertEqual(report["status"], "bracketed")
                    self.assertEqual(
                        {(seed, link) for algorithm, size, seed, link in calls if algorithm == name and size == frontier},
                        {(seed, f"res-uniform_seed_{link}.txt") for seed in (42, 43) for link in (42, 43)},
                    )
                calls.clear()
                self.assertEqual(RUNNER.main(), 0)
                self.assertTrue(calls)
                self.assertTrue(all(seed == 43 and "43" in link for _, _, seed, link in calls))


if __name__ == "__main__":
    unittest.main()
