"""Check the oracle boundary, minors, and agreement with the table verifier."""

from contextlib import redirect_stdout
from fractions import Fraction
from io import StringIO
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

from draw_matroid import load_matroid
from im_projection import CONTRACT, IM_LINES, RETAIN, main as table_main, subsets
from im_projection2 import (RankOracle, UndeterminedRankError, check_desargues, check_matroid,
                            im_extension_oracle, main)


ROOT = Path(__file__).resolve().parents[1]


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

    def test_optimized_python_cannot_bypass_verification(self):
        result = subprocess.run([sys.executable, '-O', '-c', 'import im_projection2'],
                                cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('assertions are part of the verification', result.stderr)


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

    def test_both_versions_export_identical_projection_tables(self):
        with TemporaryDirectory() as directory, redirect_stdout(StringIO()):
            first, second = Path(directory) / 'first', Path(directory) / 'second'
            table_main(first)
            main(second)
            self.assertEqual((first / 'projection.csv').read_bytes(),
                             (second / 'projection.csv').read_bytes())


if __name__ == '__main__':
    unittest.main()
