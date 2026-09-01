import importlib.util
from pathlib import Path
import sys
import types
import unittest


TQDM = types.ModuleType("tqdm")
TQDM.tqdm = object
sys.modules.setdefault("tqdm", TQDM)
MODULE_PATH = Path(__file__).with_name("generate_graph_datasets.py")
SPEC = importlib.util.spec_from_file_location("generate_graph_datasets", MODULE_PATH)
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)
LINK_MODULE_PATH = Path(__file__).with_name("generate_link_datasets.py")
LINK_SPEC = importlib.util.spec_from_file_location(
    "generate_link_datasets", LINK_MODULE_PATH
)
LINK_GENERATOR = importlib.util.module_from_spec(LINK_SPEC)
LINK_SPEC.loader.exec_module(LINK_GENERATOR)


class CactusGenerationTest(unittest.TestCase):
    def test_cycle_mass_is_converted_to_a_valid_cycle_count(self):
        self.assertEqual(GENERATOR.cactus_cycles(10, 4, 0.75), 2)
        self.assertEqual(GENERATOR.cactus_cycles(100, 4, 1.0), 33)

    def test_cactus_parameters_are_forwarded_to_the_generator(self):
        command = GENERATOR.graph_generator_command(
            Path("generate_datasets"),
            Path("datasets/cacti"),
            10,
            43,
            {
                "generator": "cactus",
                "cycle_length": 4,
                "cycle_mass": 0.75,
            },
        )

        self.assertIn("--cycles", command)
        self.assertEqual(command[command.index("--cycles") + 1], "2")
        self.assertEqual(command[command.index("--cycle-length") + 1], "4")
        self.assertEqual(command[command.index("--seed") + 1], "43")

    def test_invalid_cactus_settings_are_rejected(self):
        with self.assertRaises(ValueError):
            GENERATOR.cactus_cycles(10, 2, 0.75)
        with self.assertRaises(ValueError):
            GENERATOR.cactus_cycles(10, 4, 1.1)


class LinkGenerationTest(unittest.TestCase):
    def test_cacti_and_repeated_seeds_are_configured(self):
        self.assertIn("cacti", LINK_GENERATOR.CONFIG["folders"])
        self.assertEqual(
            LINK_GENERATOR.CONFIG["distributions"]["int_uniform_1_5"]["seeds"],
            [42, 43, 44, 45, 46],
        )

    def test_seed_is_forwarded_to_the_link_generator(self):
        command = LINK_GENERATOR.link_generator_command(
            Path("generate_datasets"),
            Path("cactus.graph"),
            Path("cactus-seed_43.links"),
            ["--distribution", "float_uniform"],
            43,
        )

        self.assertEqual(command[command.index("--seed") + 1], "43")


if __name__ == "__main__":
    unittest.main()
