import importlib.util
from pathlib import Path
import sys
import types
import unittest


TQDM = types.ModuleType("tqdm")
TQDM.tqdm = object
sys.modules.setdefault("tqdm", TQDM)
MODULE_PATH = Path(__file__).with_name("run_configurations.py")
SPEC = importlib.util.spec_from_file_location("run_configurations", MODULE_PATH)
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


class AdaptiveSchedulingTest(unittest.TestCase):
    def test_probes_are_spaced_by_node_count(self):
        sizes = [100, 200, 300, 400, 500, 1000, 2000, 3000]

        indices = RUNNER.evenly_spaced_indices(sizes, 4)

        self.assertEqual([sizes[index] for index in indices], [100, 1000, 2000, 3000])

    def test_fills_sizes_through_first_timeout_after_last_solve(self):
        sizes = list(range(100, 1100, 100))
        calls = []

        def run_level(size):
            calls.append(size)
            return "solved" if size <= 300 else "timeout"

        outcomes = RUNNER.run_adaptive_series(sizes, 6, run_level)

        self.assertEqual(set(outcomes), {100, 200, 300, 400, 500, 600, 800})
        self.assertEqual(outcomes[300], "solved")
        self.assertEqual(outcomes[400], "timeout")
        self.assertNotIn(700, calls)
        self.assertNotIn(1000, calls)

    def test_larger_solved_probe_moves_the_frontier(self):
        sizes = list(range(100, 1100, 100))

        def run_level(size):
            if size in {100, 300, 800}:
                return "solved"
            return "timeout"

        outcomes = RUNNER.run_adaptive_series(sizes, 6, run_level)

        self.assertEqual(set(outcomes), set(sizes))


if __name__ == "__main__":
    unittest.main()
