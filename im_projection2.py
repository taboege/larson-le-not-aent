"""Verify the IM minor through a partial rank oracle.

Run: python im_projection2.py [--output im_projection2_output]
This comparison version shares the exact bounds engine and input data with
im_projection.py. It propagates bounds once; a rank query returns a value
only when both endpoints agree. Uncertain ambient ranks raise an error.
The minor and Desargues checks use only the rank oracle interface.
"""

import argparse
import csv
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from time import perf_counter

from im_projection import NONBASES, IM_LINES, CONTRACT, RETAIN, from_rank_table, subsets


class UndeterminedRankError(ValueError):
    """The inequalities do not force a unique rank for the queried set."""


class RankOracle:
    """A callable rank function on subsets of an explicit ground set."""

    def __init__(self, ground, evaluate):
        self.ground = frozenset(ground)
        self._evaluate = lru_cache(None)(evaluate)

    def __call__(self, points):
        s = frozenset(points)
        if not s <= self.ground:
            raise ValueError(f"Unknown elements in rank query: {sorted(s - self.ground)}.")
        return self._evaluate(s)

    def minor(self, contracted, retained):
        """Return the oracle S -> r(S union contracted) - r(contracted)."""
        contracted, retained = frozenset(contracted), frozenset(retained)
        if contracted & retained or not (contracted | retained) <= self.ground:
            raise ValueError("Contracted and retained sets must be disjoint subsets of the ground set.")
        offset = self(contracted)
        return RankOracle(retained, lambda s: self(s | contracted) - offset)


def im_extension_oracle(ground, seed_rank, extensions):
    """Propagate a seed rank oracle through (new point, three lines) steps."""
    bounds = from_rank_table({s: seed_rank(s) for s in subsets(ground)})
    for point, lines in extensions:
        bounds, _ = bounds.im_extend(point, lines)
        bounds = bounds.compress()

    def rank(s):
        lower, upper = bounds.rank_bounds(s)
        if lower != upper:
            raise UndeterminedRankError(
                f"Rank of {sorted(s)} is not uniquely determined: [{lower}, {upper}].")
        return lower

    return RankOracle(bounds.closure.ground, rank)


def check_matroid(rank):
    """Check all matroid rank axioms using only the oracle's public interface."""
    assert rank(()) == 0
    for s in subsets(rank.ground):
        value = rank(s)
        assert value.denominator == 1 and 0 <= value <= len(s)
        for p in rank.ground:
            assert rank(s | {p}) >= value
        for p, q in combinations(rank.ground, 2):
            assert rank(s | {p}) + rank(s | {q}) >= value + rank(s | {p, q})


def check_desargues(rank):
    """Check the complete minor and the planar theorem using rank queries only."""
    points = rank.ground
    assert points == frozenset(RETAIN)
    check_matroid(rank)
    assert rank(points) == 3
    r = lambda *ps: rank(ps)
    assert all(r(p) == 1 for p in points)
    assert all(r(p, q) == 2 for p, q in combinations(points, 2))
    o, a, b, x = 4, (2, 15, 12), (6, 16, 8), (5, 14, 1)
    assert r(*a) == r(*b) == 3
    assert all(r(o, a[i], b[i]) == 2 for i in range(3))
    assert all(r(o, a[i], a[j]) == 3 for i, j in combinations(range(3), 2)), \
        "Planar Desargues needs distinct corresponding side lines."
    for k, (i, j) in enumerate(((1, 2), (0, 2), (0, 1))):
        assert r(a[i], a[j], x[k]) == r(b[i], b[j], x[k]) == 2
    assert r(*x) == 3
    print(f"Desargues: O={o}, a={a}, b={b}, x={x}")
    print("The restriction is a simple rank-3 matroid; planar AE Desargues is violated.")


def main(output):
    start = perf_counter()
    nonbases = set(map(frozenset, NONBASES))
    bases = [frozenset(b) for b in combinations(range(1, 11), 4) if frozenset(b) not in nonbases]
    assert len(bases) == 190
    seed_rank = lambda s: max(len(s & b) for b in bases)
    ambient = im_extension_oracle(range(1, 11), seed_rank, enumerate(IM_LINES, 11))
    minor = ambient.minor(CONTRACT, RETAIN)
    check_desargues(minor)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "projection.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["local_mask", "points", "lower", "upper"])
        for local, s in enumerate(subsets(RETAIN)):
            value = minor(s)
            writer.writerow([local, " ".join(map(str, sorted(s))), str(value), str(value)])
    print(f"Wrote all 1024 exact minor ranks to {output}.")
    print(f"Exact arithmetic only; {perf_counter() - start:.2f} seconds.")


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("im_projection2_output"))
    main(parser.parse_args().output)


if __name__ == "__main__":
    cli()
