"""Check the planar theorem's hypotheses, including a degenerate counterexample."""

from contextlib import redirect_stdout
from io import StringIO
from itertools import combinations
import json
from pathlib import Path
import unittest

from im_projection import RETAIN, RankOracle, check_desargues


ROOT = Path(__file__).resolve().parents[1]


def field_rank(columns, prime=101):
    """Small exact Gaussian elimination over a prime field."""
    rows = [list(column) for column in columns]
    rank = 0
    for j in range(3):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][j] % prime), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][j] % prime, -1, prime)
        rows[rank] = [(x * inverse) % prime for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank:
                coefficient = rows[i][j]
                rows[i] = [(x - coefficient * y) % prime
                           for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


class PlanarDesarguesTests(unittest.TestCase):
    def test_final_matroid_passes_the_explicit_planar_hypotheses(self):
        data = json.loads((ROOT / "data/final_matroid.json").read_text())
        bases = list(map(frozenset, data["bases"]))
        rank = lambda subset: max(len(subset & basis) for basis in bases)
        with redirect_stdout(StringIO()):
            check_desargues(RankOracle(RETAIN, rank))

    def test_representable_degeneracy_is_not_a_desargues_obstruction(self):
        # The F_101 columns in the README, with the verifier's role labels.
        # The sides a2-a3 and b2-b3 coincide, allowing an arbitrary x1.
        columns = {
            4: (2, 0, 1),   # O
            2: (0, 1, 1),   # a1
            15: (0, 0, 1),  # a2
            12: (1, 0, 1),  # a3
            6: (2, 1, 2),   # b1
            16: (3, 0, 1),  # b2
            8: (4, 0, 1),   # b3
            5: (5, 0, 1),   # x1
            14: (2, 3, 5),  # x2
            1: (0, 3, 4),   # x3
        }
        rank = lambda subset: field_rank([columns[p] for p in subset])
        self.assertTrue(all(rank(pair) == 2 for pair in combinations(RETAIN, 2)))
        incidences = [(4, 2, 6), (4, 15, 16), (4, 12, 8),
                      (15, 12, 5), (16, 8, 5), (2, 12, 14), (6, 8, 14),
                      (2, 15, 1), (6, 16, 1)]
        self.assertTrue(all(rank(triple) == 2 for triple in incidences))
        self.assertEqual(rank([2, 15, 12]), 3)
        self.assertEqual(rank([6, 16, 8]), 3)
        self.assertEqual(rank([1, 5, 14]), 3)
        with self.assertRaisesRegex(AssertionError, "distinct corresponding side lines"):
            check_desargues(RankOracle(RETAIN, rank))


if __name__ == "__main__":
    unittest.main()
