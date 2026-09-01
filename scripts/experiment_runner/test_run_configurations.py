import importlib.util
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

        configurations, links, datasets = RUNNER.load_configurations(
            configuration_path
        )

        self.assertEqual(len(configurations), 76)
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


if __name__ == "__main__":
    unittest.main()
