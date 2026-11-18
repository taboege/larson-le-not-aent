"""Exact rank deductions, real-valued possibilities, and the complete example."""

from contextlib import redirect_stdout
from fractions import Fraction
from io import StringIO
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

from draw_im_steps import compact_steps
from draw_matroid import load_matroid
from im_projection import RankBounds, from_rank_table, main, mask, projection


ROOT = Path(__file__).resolve().parents[1]


class RankBoundsTests(unittest.TestCase):
    def test_exact_fractional_input(self):
        # A modular polymatroid with singleton ranks 1/3 and 2/3.
        table = [0, Fraction(1, 3), Fraction(2, 3), 1]
        ranks = from_rank_table(table)
        self.assertEqual([ranks.rank_bounds(s) for s in range(4)],
                         [(r, r) for r in table])

    def test_floating_point_input_is_rejected(self):
        with self.assertRaisesRegex(AssertionError, "exact rational"):
            RankBounds(1, range(2), {0: 0, 1: 1.0})

    def test_submodularity_forces_ranks_and_parallel_closure(self):
        # Points 1 and 2 are parallel, and point 3 is independent of them.
        ranks = RankBounds(3, range(8), {0: 0, 1: 1, 2: 1, 3: 1, 4: 1, 7: 2})
        expected = [0, 1, 1, 1, 1, 2, 2, 2]
        self.assertEqual([ranks.rank_bounds(s) for s in range(8)],
                         [(r, r) for r in expected])
        compressed = ranks.compress()
        self.assertEqual(compressed.closure(mask([1])), mask([1, 2]))
        self.assertEqual(compressed.closure(mask([1, 3])), mask([1, 2, 3]))
        self.assertEqual([compressed.rank_bounds(s) for s in range(8)],
                         [(r, r) for r in expected])

    def test_integer_endpoints_keep_fractional_polymatroids(self):
        ranks = RankBounds(3, range(8), {0: 0, 1: 1, 2: 1, 4: 1, 7: 2})
        # The average of U(1,3) and U(3,3) has rank 3/2 on each pair.
        witness = [Fraction(min(s.bit_count(), 1) + s.bit_count(), 2)
                   for s in range(8)]
        for s, value in enumerate(witness):
            lo, hi = ranks.rank_bounds(s)
            self.assertLessEqual(lo, value)
            self.assertLessEqual(value, hi)
        self.assertEqual(ranks.rank_bounds(mask([1, 2])), (1, 2))

    def test_inconsistent_submodular_data_are_rejected(self):
        with self.assertRaisesRegex(AssertionError, "Inconsistent rank bounds"):
            from_rank_table([0, 1, 1, 3])

    def test_uncertain_contracted_rank_is_rejected(self):
        ranks = RankBounds(3, range(8), {0: 0, 1: 1, 2: 1, 4: 1, 7: 2})
        with self.assertRaisesRegex(AssertionError, "contracted set must have a fixed rank"):
            projection(ranks, mask([1, 2]), [3])

    def test_upper_endpoints_must_supply_a_feasible_witness(self):
        # These constraints are feasible (take 1||2, 3||4, and independent 5),
        # but the coordinatewise upper bounds are not jointly feasible.
        fixed = {0: 0, 31: 3, mask([1, 2, 3]): 2, mask([1, 2, 4]): 2}
        fixed.update({1 << i: 1 for i in range(5)})
        with self.assertRaisesRegex(AssertionError, "exact feasible witness"):
            RankBounds(5, range(32), fixed)


class CompleteVerificationTests(unittest.TestCase):
    def test_six_extensions_reproduce_all_minor_ranks_and_display_incidences(self):
        with TemporaryDirectory() as directory, redirect_stdout(StringIO()):
            output = Path(directory)
            main(output)
            # The CSV reader verifies every rank against the recovered bases.
            actual = load_matroid(output / "projection.csv")
            expected = load_matroid(ROOT / "data/final_matroid.json")
            self.assertEqual(actual.ground, expected.ground)
            self.assertEqual(actual.bases, expected.bases)
            stages = json.loads((output / "extensions.json").read_text())
            original = load_matroid(ROOT / "data/dual_matroid.json")
            saved = json.loads((ROOT / "data/im_steps.json").read_text())
            self.assertEqual(compact_steps(stages, original), saved)

    def test_optimized_python_cannot_bypass_verification_on_import(self):
        result = subprocess.run([sys.executable, "-O", "-c", "import im_projection"],
                                cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("assertions are part of the verification", result.stderr)


if __name__ == "__main__":
    unittest.main()
