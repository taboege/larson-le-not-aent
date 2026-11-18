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
from im_projection import (CONTRACT, IM_LINES, RETAIN, ClosureSystem, Inequality,
                           RankBounds, RankOracle, UndeterminedRankError,
                           check_desargues, check_matroid, from_rank_table,
                           im_extension_oracle, main, subsets)


ROOT = Path(__file__).resolve().parents[1]


def bounds(points, fixed):
    return RankBounds(ClosureSystem(subsets(points)),
                      {frozenset(s): value for s, value in fixed.items()})


class ClosureSystemTests(unittest.TestCase):
    def test_implications_apply_to_every_superset(self):
        original = ClosureSystem(subsets((2, 5, 9)))
        refined = original.imply(iter([(frozenset({2}), frozenset({2, 5}))]))
        self.assertEqual(refined({2}), {2, 5})
        self.assertEqual(refined({2, 9}), {2, 5, 9})
        self.assertEqual(refined({9}), {9})
        self.assertEqual(original({2, 9}), {2, 9})

    def test_adjoined_point_is_spanned_by_the_declared_sets(self):
        original = ClosureSystem(subsets((2, 5, 9)))
        extended = original.adjoin(11, iter(((2, 5), (5, 9))))
        self.assertEqual(extended({2, 5}), {2, 5, 11})
        self.assertEqual(extended({5, 9}), {5, 9, 11})
        self.assertEqual(extended({2, 9}), {2, 9})
        self.assertEqual(extended({11}), {11})
        self.assertEqual(original.ground, {2, 5, 9})
        with self.assertRaises(AssertionError):
            original.adjoin(5, ())

    def test_non_intersection_closed_family_is_rejected(self):
        with self.assertRaises(AssertionError):
            ClosureSystem(((), (1, 2), (2, 3), (1, 2, 3)))


class InequalityTests(unittest.TestCase):
    def test_submodularity_tightens_each_of_its_four_variables(self):
        inequality = Inequality(('a', 'b'), ('c', 'd'))
        # a+b >= c+d, with all but one of the ranks fixed.
        expected = {'a': (3, 10), 'b': (2, 10), 'c': (0, 6), 'd': (0, 2)}
        for variable, interval in expected.items():
            with self.subTest(variable=variable):
                lower = dict(a=4, b=3, c=5, d=1)
                upper = lower.copy()
                lower[variable], upper[variable] = 0, 10
                inequality.tighten(lower, upper)
                self.assertEqual((lower[variable], upper[variable]), interval)
        self.assertTrue(inequality.holds(dict(a=4, b=3, c=5, d=1)))
        self.assertFalse(inequality.holds(dict(a=2, b=3, c=5, d=1)))

    def test_monotonicity_uses_the_same_rule_with_exact_fractions(self):
        inequality = Inequality(('larger',), ('smaller',))
        lower = dict(smaller=Fraction(1, 3), larger=0)
        upper = dict(smaller=1, larger=Fraction(2, 3))
        inequality.tighten(lower, upper)
        self.assertEqual(lower['larger'], Fraction(1, 3))
        self.assertEqual(upper['smaller'], Fraction(2, 3))

    def test_overlapping_rank_variables_are_rejected(self):
        with self.assertRaises(AssertionError):
            Inequality(('a', 'b'), ('b', 'c'))


class RankBoundsTests(unittest.TestCase):
    def test_exact_fractional_input(self):
        # A modular polymatroid with singleton ranks 1/3 and 2/3.
        table = dict(zip(subsets((1, 2)), [0, Fraction(1, 3), Fraction(2, 3), 1]))
        ranks = from_rank_table(table)
        self.assertEqual([ranks.rank_bounds(s) for s in table],
                         [(r, r) for r in table.values()])

    def test_floating_point_input_is_rejected(self):
        with self.assertRaisesRegex(AssertionError, "exact rational"):
            bounds((1,), {(): 0, (1,): 1.0})

    def test_submodularity_forces_ranks_and_parallel_closure(self):
        # Points 1 and 2 are parallel, and point 3 is independent of them.
        ranks = bounds((1, 2, 3), {(): 0, (1,): 1, (2,): 1, (1, 2): 1, (3,): 1, (1, 2, 3): 2})
        expected = [0, 1, 1, 1, 1, 2, 2, 2]
        self.assertEqual([ranks.rank_bounds(s) for s in subsets((1, 2, 3))],
                         [(r, r) for r in expected])
        compressed = ranks.compress()
        self.assertEqual(compressed.closure({1}), {1, 2})
        self.assertEqual(compressed.closure({1, 3}), {1, 2, 3})
        self.assertEqual([compressed.rank_bounds(s) for s in subsets((1, 2, 3))],
                         [(r, r) for r in expected])

    def test_integer_endpoints_keep_fractional_polymatroids(self):
        ranks = bounds((1, 2, 3), {(): 0, (1,): 1, (2,): 1, (3,): 1, (1, 2, 3): 2})
        # The average of U(1,3) and U(3,3) has rank 3/2 on each pair.
        witness = {s: Fraction(min(len(s), 1) + len(s), 2) for s in subsets((1, 2, 3))}
        for s, value in witness.items():
            lo, hi = ranks.rank_bounds(s)
            self.assertLessEqual(lo, value)
            self.assertLessEqual(value, hi)
        self.assertEqual(ranks.rank_bounds({1, 2}), (1, 2))
        self.assertEqual(ranks.closure({1, 2}), {1, 2})
        self.assertNotIn(frozenset({1, 2}), ranks.fixed)

    def test_inconsistent_submodular_data_are_rejected(self):
        with self.assertRaisesRegex(AssertionError, "Inconsistent rank bounds"):
            from_rank_table(dict(zip(subsets((1, 2)), [0, 1, 1, 3])))

    def test_rank_table_must_contain_every_subset(self):
        with self.assertRaises(AssertionError):
            from_rank_table({frozenset(): 0, frozenset({2, 9}): 2})

    def test_loops_are_in_the_closure_of_the_empty_set(self):
        table = dict(zip(subsets((2, 9)), [0, 0, 1, 1]))
        ranks = from_rank_table(table)
        self.assertEqual(ranks.closure(()), {2})
        self.assertEqual([ranks.rank_bounds(s) for s in table],
                         [(r, r) for r in table.values()])

    def test_im_extension_requires_rank_four_for_the_three_line_union(self):
        # U(3,6): all line ranks and pairwise union ranks have the required
        # values, but the union of all three lines still has rank three.
        ranks = from_rank_table({s: min(len(s), 3) for s in subsets(range(1, 7))})
        with self.assertRaises(AssertionError):
            ranks.im_extend(7, ((1, 2), (3, 4), (5, 6)))

    def test_upper_endpoints_must_supply_a_feasible_witness(self):
        # These constraints are feasible (take 1||2, 3||4, and independent 5),
        # but the coordinatewise upper bounds are not jointly feasible.
        fixed = {(): 0, (1, 2, 3, 4, 5): 3, (1, 2, 3): 2, (1, 2, 4): 2}
        fixed.update({(p,): 1 for p in range(1, 6)})
        with self.assertRaisesRegex(AssertionError, "exact feasible witness"):
            bounds(range(1, 6), fixed)


class RankOracleTests(unittest.TestCase):
    def test_fractional_seed_and_minor_keep_the_original_labels(self):
        weights = {2: Fraction(1, 3), 7: Fraction(2, 3), 11: Fraction(4, 3)}
        rank = im_extension_oracle(weights, lambda s: sum(weights[p] for p in s), ())
        minor = rank.minor({7}, {2, 11})
        self.assertEqual(minor.ground, {2, 11})
        self.assertEqual(minor(()), 0)
        self.assertEqual(minor({2}), Fraction(1, 3))
        self.assertEqual(minor({11}), Fraction(4, 3))
        self.assertEqual(minor({2, 11}), Fraction(5, 3))
        with self.assertRaisesRegex(ValueError, 'Unknown elements'):
            minor({7})

    def test_minor_requires_disjoint_subsets_of_the_ground_set(self):
        rank = RankOracle({2, 7}, len)
        for contracted, retained in (({2}, {2}), ({9}, {7}), ({2}, {9})):
            with self.subTest(contracted=contracted, retained=retained):
                with self.assertRaisesRegex(ValueError, 'disjoint subsets'):
                    rank.minor(contracted, retained)

    def test_fractional_rank_is_allowed_by_the_oracle_but_is_not_a_matroid(self):
        rank = RankOracle({2}, lambda s: Fraction(len(s), 3))
        self.assertEqual(rank({2}), Fraction(1, 3))
        with self.assertRaises(AssertionError):
            check_matroid(rank)


class IMOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        seed = load_matroid(ROOT / 'data/dual_matroid.json')
        cls.ambient = im_extension_oracle(seed.ground, seed.rank, enumerate(IM_LINES, 11))
        cls.minor = cls.ambient.minor(CONTRACT, RETAIN)
        cls.expected = load_matroid(ROOT / 'data/final_matroid.json')

    def test_every_minor_rank_matches_the_saved_matroid(self):
        for s in subsets(RETAIN):
            self.assertEqual(self.minor(s), self.expected.rank(s), sorted(s))

    def test_uncertain_ambient_rank_and_contraction_raise_an_error(self):
        uncertain = {1, 3, 9, 14}
        with self.assertRaisesRegex(UndeterminedRankError, r'not uniquely determined: \[3, 4\]'):
            self.ambient(uncertain)
        with self.assertRaises(UndeterminedRankError):
            self.ambient.minor(uncertain, {2})

    def test_desargues_consumer_accepts_an_independent_basis_rank_oracle(self):
        with redirect_stdout(StringIO()):
            check_desargues(RankOracle(self.expected.ground, self.expected.rank))


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

    def test_standalone_script_reproduces_the_minor_without_other_project_files(self):
        with TemporaryDirectory() as directory:
            script = Path(directory) / 'im_projection.py'
            script.write_bytes((ROOT / 'im_projection.py').read_bytes())
            output = Path(directory) / 'output'
            result = subprocess.run([sys.executable, '-I', '-B', str(script), '--output', str(output)],
                                    cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(sorted(p.name for p in output.iterdir()), ['extensions.json', 'projection.csv'])
            actual = load_matroid(output / 'projection.csv')
            expected = load_matroid(ROOT / 'data/final_matroid.json')
            self.assertEqual(actual.ground, expected.ground)
            self.assertEqual(actual.bases, expected.bases)

    def test_optimized_python_cannot_bypass_verification_on_import(self):
        result = subprocess.run([sys.executable, "-O", "-c", "import im_projection"],
                                cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("assertions are part of the verification", result.stderr)


if __name__ == "__main__":
    unittest.main()
