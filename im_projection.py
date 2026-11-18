"""Six IM extensions of Larson--Le's dual, verified by exact rank bounds.

Run: python im_projection.py [--output im_projection_output]
Only Python's standard library is needed.

Subsets are frozensets of element labels. A closure system identifies
h(S) with h(cl(S)); its closed sets need not be matroid flats.
This closure is defined on every subset using proved dependencies; ranks
of its closed sets are represented by intervals, not necessarily exact values.
For an intersection-closed family, monotonicity and submodularity on closed
sets imply both inequalities on all subsets: cl(A union B) equals
cl(cl(A) union cl(B)), and cl(A intersection B) is contained in
cl(A) intersection cl(B). Closure compression uses only proved equalities.

Bounds propagate by rearranging these inequalities, using exact arithmetic.
Integer interval endpoints do NOT require unknown ranks to be integers.
Propagation need not find optimal bounds, but proves every equality needed
here. The upper-bound vector is checked as an exact feasible witness at
every stage. No solver, tolerances, saved ranks, or certificates are used.
The partial rank oracle returns a value only when its bounds agree. Minor
construction and the final checks use only this oracle interface.
Run without Python's -O option: assertions are part of the verification.

The deductions use Shannon inequalities, not almost entropicity. The mathematical
ingredients outside this program are the AE three-line extension theorem,
minor closure, and nondegenerate planar Desargues (derived in the README
from Yashfe, arXiv:2509.06302, Lemma 16 and Theorem 18).
"""

import argparse
import csv
import json
from functools import lru_cache
from itertools import combinations
from numbers import Rational
from pathlib import Path
from time import perf_counter


if not __debug__:
    raise RuntimeError("Run without -O: assertions are part of the verification.")


def subsets(points):
    """All subsets, in binary order of the given labels for CSV local_mask."""
    result = [frozenset()]
    for point in points:
        result += [s | {point} for s in result]
    return result


class ClosureSystem:
    """Total closure from proved dependencies, not necessarily all dependencies.

    The family is intersection-closed and ordered by size, so its first
    superset of S is the unique least closed superset of S.
    """

    def __init__(self, closed_sets):
        family = frozenset(map(frozenset, closed_sets))
        self.ground = frozenset().union(*family)
        assert self.ground in family
        assert all(a & b in family for a, b in combinations(family, 2))
        self.closed = tuple(sorted(family, key=lambda s: (len(s), sorted(s))))
        self._closure = lru_cache(None)(lambda s: next(f for f in self.closed if s <= f))

    def __call__(self, points):
        s = frozenset(points)
        assert s <= self.ground
        return self._closure(s)

    def imply(self, implications):
        """Keep precisely the closed sets respecting every implication A -> B."""
        implications = tuple(implications)
        return ClosureSystem(s for s in self.closed
                             if all(not a <= s or b <= s for a, b in implications))

    def adjoin(self, point, spans):
        """Adjoin a point spanned by each specified old set."""
        spans = tuple(map(frozenset, spans))
        assert point not in self.ground and all(s <= self.ground for s in spans)
        return ClosureSystem([f | {point} for f in self.closed] +
                             [f for f in self.closed if not any(s <= f for s in spans)])


class Inequality:
    """sum h(positive) >= sum h(negative), with distinct rank variables."""

    def __init__(self, positive, negative):
        self.positive, self.negative = tuple(positive), tuple(negative)
        terms = self.positive + self.negative
        assert len(set(terms)) == len(terms)

    def tighten(self, lower, upper):
        left_max = sum(upper[s] for s in self.positive)
        right_min = sum(lower[s] for s in self.negative)
        for s in self.positive:
            lower[s] = max(lower[s], right_min - left_max + upper[s])
        for s in self.negative:
            upper[s] = min(upper[s], left_max - right_min + lower[s])

    def holds(self, rank):
        return sum(rank[s] for s in self.positive) >= sum(rank[s] for s in self.negative)


class RankBounds:
    """Sound rank intervals from monotonicity and submodularity.

    Inputs are exact rational ranks, including the empty and full sets.
    For these data the upper bounds also form a feasible rank function;
    stop if that witness check fails, even if the intervals are consistent.
    """

    def __init__(self, closure, fixed):
        self.closure = closure
        closed, ground = closure.closed, closure.ground
        assert all(isinstance(value, Rational) for value in fixed.values()), \
            "Ranks must be exact rational numbers."
        known = {}
        for s, value in fixed.items():
            f = closure(s)
            assert f not in known or known[f] == value
            known[f] = value
        assert known[closure(())] == 0
        lo = self.lower = {s: known.get(s, 0) for s in closed}
        hi = self.upper = {s: known.get(s, known[ground]) for s in closed}

        inclusions = dict.fromkeys((a, closure(a | {p})) for a in closed for p in sorted(ground - a))
        inequalities = [Inequality((b,), (a,)) for a, b in inclusions]
        inequalities += [Inequality((a, b), (closure(a | b), a & b))
                         for a, b in combinations(closed, 2) if not (a <= b or b <= a)]
        while True:
            previous = lo.copy(), hi.copy()
            for inequality in inequalities:
                inequality.tighten(lo, hi)
            assert all(lo[s] <= hi[s] for s in closed), "Inconsistent rank bounds."
            if previous == (lo, hi):
                break

        # Here the upper endpoints themselves witness feasibility exactly.
        assert all(inequality.holds(hi) for inequality in inequalities), \
            "The upper bounds must form an exact feasible witness."
        self.fixed = {s: lo[s] for s in closed if lo[s] == hi[s]}

    def rank_bounds(self, s):
        f = self.closure(s)
        return self.lower[f], self.upper[f]

    def compress(self):
        """Use forced equalities h(A)=h(B), A subset B, to compress closure.

        Such an equality gives B subset cl(A). Retaining only closed sets
        respecting these implications propagates them to EVERY superset.
        Nonconstant intervals continue to allow real ranks.
        """
        fixed = self.fixed
        implications = [(a, b) for a in fixed for b in fixed
                        if a < b and fixed[a] == fixed[b]]
        return RankBounds(self.closure.imply(implications), fixed)

    def im_extend(self, point, lines):
        """Add z of rank one spanned by each of three rank-two sets.

        Check the seven premises for every allowed polymatroid: each line has
        rank 2, each pair of lines has union rank 3, and the total rank is 4.
        The returned bounds include every polymatroid extension satisfying
        h(z)=1 and h(L+z)=h(L), for each input line L.
        """
        lines = tuple(map(frozenset, lines))
        assert len(lines) == 3
        premises = [(frozenset().union(*group), size + 1)
                    for size in (1, 2, 3) for group in combinations(lines, size)]
        for s, value in premises:
            assert self.rank_bounds(s) == (value, value), (sorted(s), self.rank_bounds(s))
        closure = self.closure.adjoin(point, lines)
        return RankBounds(closure, {**self.fixed, frozenset({point}): 1}), premises


def from_rank_table(rank):
    """Start with any rational polymatroid rank table, keyed by frozensets."""
    ground = frozenset().union(*rank)
    assert set(rank) == set(subsets(ground))
    closed = ClosureSystem(s for s in rank if all(rank[s | {p}] > rank[s] for p in ground - s))
    # Passing every original rank also detects inconsistent identifications.
    return RankBounds(closed, rank)


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
        """Return S -> r(S union contracted) - r(contracted), with fixed offset."""
        contracted, retained = frozenset(contracted), frozenset(retained)
        if contracted & retained or not (contracted | retained) <= self.ground:
            raise ValueError("Contracted and retained sets must be disjoint subsets of the ground set.")
        offset = self(contracted)
        return RankOracle(retained, lambda s: self(s | contracted) - offset)


def im_extension_oracle(ground, seed_rank, extensions, on_extension=None):
    """Propagate a seed rank oracle through (new point, three lines) steps.

    The optional observer receives each step's compressed bounds for reporting.
    """
    bounds = from_rank_table({s: seed_rank(s) for s in subsets(ground)})
    for point, lines in extensions:
        lines = tuple(map(frozenset, lines))
        bounds, premises = bounds.im_extend(point, lines)
        bounds = bounds.compress()
        if on_extension is not None:
            on_extension(point, lines, premises, bounds)

    def rank(s):
        lower, upper = bounds.rank_bounds(s)
        if lower != upper:
            raise UndeterminedRankError(
                f"Rank of {sorted(s)} is not uniquely determined: [{lower}, {upper}].")
        return lower

    return RankOracle(bounds.closure.ground, rank)


# The example-specific data starts here. N* has rank four; its other
# four-subsets are bases (Larson--Le, arXiv:2609.02834, equation (1)).
NONBASES = [
    (1,2,3,4), (1,2,5,6), (1,3,6,7), (2,3,5,7), (2,4,6,7),
    (2,5,8,9), (2,5,8,10), (2,5,9,10), (2,8,9,10),
    (3,4,5,6), (3,4,5,8), (3,4,6,8), (3,5,6,8), (4,5,6,8),
    (4,5,7,10), (4,6,9,10), (4,7,8,9), (5,6,7,9),
    (5,8,9,10), (6,7,8,10),
]
IM_LINES = [
    ((4,6),  (5,8),  (9,10)),  # add 11
    ((4,8),  (5,6),  (7,9)),   # add 12
    ((2,11), (4,7),  (5,10)),  # add 13
    ((4,5),  (6,8),  (7,10)),  # add 14
    ((1,2),  (3,4),  (5,6)),   # add 15
    ((1,6),  (2,5),  (3,7)),   # add 16
]
CONTRACT = (13,)
RETAIN = (1, 2, 4, 5, 6, 8, 12, 14, 15, 16)


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
    # Coincident corresponding sides would allow arbitrary intersections.
    # The explicit planar theorem in the README excludes that degeneracy.
    assert all(r(o, a[i], a[j]) == 3 for i, j in combinations(range(3), 2)), \
        "Planar Desargues needs distinct corresponding side lines."
    for k, (i, j) in enumerate(((1, 2), (0, 2), (0, 1))):
        assert r(a[i], a[j], x[k]) == r(b[i], b[j], x[k]) == 2
    assert r(*x) == 3
    print(f"Desargues: O={o}, a={a}, b={b}, x={x}")
    print("All nine collinear triples have rank 2, but r(x1,x2,x3)=3.")
    print("The restriction is a simple rank-3 matroid; planar AE Desargues is violated.")


def main(output):
    start = perf_counter()
    stages = []

    def record(point, lines, premises, bounds):
        intervals = [(s, *bounds.rank_bounds(s)) for s in bounds.closure.closed]
        uncertain = sum(lo != hi for _, lo, hi in intervals)
        print(f"Add {point}: {len(intervals)} closed sets after compression; {uncertain} undetermined ranks.", flush=True)
        stages.append({"new_point": point, "lines": [sorted(line) for line in lines],
                       "premises": [{"set": sorted(s), "rank": v} for s, v in premises],
                       "closed_set_bounds": [{"set": sorted(s), "lower": str(lo), "upper": str(hi)}
                                             for s, lo, hi in intervals]})

    nonbases = set(map(frozenset, NONBASES))
    bases = [frozenset(b) for b in combinations(range(1, 11), 4) if frozenset(b) not in nonbases]
    assert len(bases) == 190
    seed_rank = lambda s: max(len(s & b) for b in bases)
    ambient = im_extension_oracle(range(1, 11), seed_rank, enumerate(IM_LINES, 11), record)
    minor = ambient.minor(CONTRACT, RETAIN)
    check_desargues(minor)
    output.mkdir(parents=True, exist_ok=True)
    (output / "extensions.json").write_text(json.dumps(stages, indent=2) + "\n")
    with (output / "projection.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["local_mask", "points", "lower", "upper"])
        for local, s in enumerate(subsets(RETAIN)):
            value = minor(s)
            writer.writerow([local, " ".join(map(str, sorted(s))), str(value), str(value)])
    print(f"Wrote six extensions and all 1024 coordinate bounds to {output}.")
    print(f"Exact arithmetic only; {perf_counter() - start:.2f} seconds.")


def cli():
    """Entry point for the installed command and direct script execution."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("im_projection_output"))
    args = parser.parse_args()
    main(args.output)


if __name__ == "__main__":
    cli()
