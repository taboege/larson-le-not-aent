"""Six IM extensions of Larson--Le's dual, computed from scratch by LP.

Install: python -m pip install .
Run:     python im_projection.py [--output im_projection_output] [--write-lps]

The generic code uses labels 1,...,n and integer bitmasks for subsets. Unknown
ranks are REAL variables. A closure system identifies h(S) with h(cl(S));
its closed sets need not be the flats of a matroid. This is just compression
by known functional dependencies. Monotonicity and submodularity on the
closed sets imply those inequalities on all subsets: cl(A union B) equals
cl(cl(A) union cl(B)), while cl(A intersection B) is contained in
cl(A) intersection cl(B).

HiGHS finds floating-point dual solutions. minimum() turns these into exact
rational lower bounds using nonnegative row multipliers and column bounds.
Only EXACTLY matching lower/upper bounds are used to fix ranks or compress.
No old program, saved ranks, or certificate is read.
Rational reconstruction failures stop the computation; they cannot silently
certify a false equality. Run without Python's -O option (checks use assert).

The LPs impose Shannon inequalities, not almost entropicity. The mathematical
ingredients outside this program are the AE three-line extension theorem,
minor closure, and nondegenerate planar Desargues (derived in the README
from Yashfe, arXiv:2509.06302, Lemma 16 and Theorem 18).
"""

import argparse
import csv
import json
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from time import perf_counter

import highspy
import numpy as np


def mask(points):
    return sum(1 << (p - 1) for p in set(points))


def labels(s):
    return [i + 1 for i in range(s.bit_length()) if s >> i & 1]


def rational(x):
    # This guesses a useful rational, not an equality with the float.
    # The subsequent exact inequalities make any such guess safe.
    return Fraction(float(x)).limit_denominator(1_000_000)


def highs_ok(status):
    if status != highspy.HighsStatus.kOk:
        raise RuntimeError(f"HiGHS returned {status}")


class RankLP:
    """All polymatroids satisfying a closure system and fixed rational ranks.

    The rank of the full ground set must be fixed. Every LP variable is the
    rank of one closed set, bounded between zero and that full rank.
    """

    def __init__(self, n, closed_sets, fixed):
        self.n, self.full = n, (1 << n) - 1
        self.closed = sorted(set(closed_sets), key=lambda s: (s.bit_count(), s))
        closed = set(self.closed)
        assert self.full in closed and all(0 <= s <= self.full for s in closed)
        assert all(a & b in closed for a, b in combinations(closed, 2))

        @lru_cache(None)
        def closure(s):
            assert 0 <= s <= self.full
            return next(f for f in self.closed if f & s == s)

        self.closure = closure
        self.index = {s: i for i, s in enumerate(self.closed)}
        self.fixed = {}
        for s, value in fixed.items():
            f, value = closure(s), Fraction(value)
            assert f not in self.fixed or self.fixed[f] == value
            self.fixed[f] = value
        assert self.fixed[closure(0)] == 0
        total = self.fixed[self.full]
        self.lower = [self.fixed.get(s, Fraction(0)) for s in self.closed]
        self.upper = [self.fixed.get(s, total) for s in self.closed]
        assert all(0 <= lo <= hi <= total for lo, hi in zip(self.lower, self.upper))

        rows = set()

        def add(terms):
            row = {}
            for s, coefficient in terms:
                j = self.index[closure(s)]
                row[j] = row.get(j, 0) + coefficient
            row = tuple(sorted((j, c) for j, c in row.items() if c))
            if row:
                rows.add(row)

        # Each row means sum(coefficient * rank) >= 0.
        for a in self.closed:
            for i in range(n):
                add(((a | 1 << i, 1), (a, -1)))
        for a, b in combinations(self.closed, 2):
            if a & b != a and a & b != b:
                add(((a, 1), (b, 1), (a | b, -1), (a & b, -1)))
        self.rows = sorted(rows)
        starts, indices, coefficients = [0], [], []
        for row in self.rows:
            for j, c in row:
                indices.append(j)
                coefficients.append(c)
            starts.append(len(indices))

        self.solver = h = highspy.Highs()
        for option, value in (("output_flag", False), ("threads", 1),
                              ("solver", "simplex")):
            highs_ok(h.setOptionValue(option, value))
        count = len(self.closed)
        highs_ok(h.addVars(count, np.array(self.lower, dtype=float),
                          np.array(self.upper, dtype=float)))
        highs_ok(h.addRows(len(self.rows), np.zeros(len(self.rows)),
                          np.full(len(self.rows), highspy.kHighsInf), len(indices),
                          np.array(starts, dtype=np.int32),
                          np.array(indices, dtype=np.int32),
                          np.array(coefficients, dtype=float)))
        self.solves = 0
        self.run()
        # An exact feasible point rules out vacuous bounds from an empty LP.
        witness = list(map(rational, h.getSolution().col_value))
        assert all(lo <= x <= hi for lo, x, hi in zip(self.lower, witness, self.upper))
        assert all(sum(c * witness[j] for j, c in row) >= 0 for row in self.rows)

    def run(self):
        highs_ok(self.solver.run())
        self.solves += 1
        if self.solver.getModelStatus() != highspy.HighsModelStatus.kOptimal:
            raise RuntimeError(self.solver.modelStatusToString(self.solver.getModelStatus()))

    def minimum(self, objective):
        """An exact lower bound for a linear objective (column, coefficient).

        If A*x >= 0 and y >= 0, then c*x >= (c - A^T*y)*x. Bound each
        residual term using its column's lower or upper bound. Thus even
        imperfectly reconstructed dual multipliers give a VALID bound.
        """
        costs = np.zeros(len(self.closed))
        residual = [Fraction(0)] * len(self.closed)
        for j, c in objective:
            costs[j], residual[j] = float(c), c
        highs_ok(self.solver.changeColsCost(len(costs), np.arange(len(costs), dtype=np.int32), costs))
        self.run()
        for row, multiplier in zip(self.rows, self.solver.getSolution().row_dual):
            if multiplier > 0:
                y = rational(multiplier)
                for j, c in row:
                    residual[j] -= y * c
        return sum(d * (lo if d >= 0 else hi)
                   for d, lo, hi in zip(residual, self.lower, self.upper))

    @lru_cache(None)
    def column_bounds(self, objective):
        lo = self.minimum(objective)
        hi = -self.minimum(tuple((j, -c) for j, c in objective))
        assert lo <= hi
        return lo, hi

    def bounds(self, terms):
        """Bound sum(c*h(S) for S,c in terms), combining identical closures."""
        objective, constant = {}, Fraction(0)
        for s, coefficient in terms:
            f, c = self.closure(s), Fraction(coefficient)
            if f in self.fixed:
                constant += c * self.fixed[f]
            else:
                j = self.index[f]
                objective[j] = objective.get(j, Fraction(0)) + c
        key = tuple(sorted((j, c) for j, c in objective.items() if c))
        lo, hi = self.column_bounds(key) if key else (0, 0)
        return constant + lo, constant + hi

    def rank_bounds(self, s):
        return self.bounds(((s, 1),))

    def compress(self):
        """Discover forced ranks by LP, then use h(A)=h(B), A subset B.

        Such an equality gives B subset cl(A). Retaining only closed sets
        respecting these implications propagates them to EVERY superset.
        Nonconstant rank intervals remain unconstrained real variables.
        """
        bounds = {s: self.rank_bounds(s) for s in self.closed}
        fixed = {s: lo for s, (lo, hi) in bounds.items() if lo == hi}
        implications = [(a, b) for a in fixed for b in fixed
                        if a != b and a & b == a and fixed[a] == fixed[b]]
        closed = [s for s in self.closed
                  if all(s & a != a or s & b == b for a, b in implications)]
        return RankLP(self.n, closed, fixed)

    def im_extend(self, lines):
        """Add z of rank one spanned by each of three rank-two sets.

        Check the seven premises universally over this LP: each line has
        rank 2, each pair of lines has union rank 3, and the total rank is 4.
        The returned LP includes every polymatroid extension satisfying
        h(z)=1 and h(L+z)=h(L), for each input line L.
        """
        assert len(lines) == 3
        a, b, c = lines
        premises = [(a, 2), (b, 2), (c, 2), (a | b, 3), (a | c, 3),
                    (b | c, 3), (a | b | c, 4)]
        for s, value in premises:
            assert self.rank_bounds(s) == (value, value), (labels(s), self.rank_bounds(s))
        z = 1 << self.n
        # Old dependencies still hold, and any set containing a line spans z.
        closed = ({f | z for f in self.closed} |
                  {f for f in self.closed if not any(f & line == line for line in lines)})
        return RankLP(self.n + 1, closed, {**self.fixed, z: 1}), premises

    def write_lp(self, path):
        for s, j in self.index.items():
            name = "h_" + ("_".join(map(str, labels(s))) or "empty")
            highs_ok(self.solver.passColName(j, name))
        for i in range(len(self.rows)):
            highs_ok(self.solver.passRowName(i, f"Shannon_{i}"))
        count = len(self.closed)
        highs_ok(self.solver.changeColsCost(count, np.arange(count, dtype=np.int32), np.zeros(count)))
        highs_ok(self.solver.writeModel(str(path)))


def from_rank_table(rank):
    """Start with any rational polymatroid rank table, indexed by subset mask."""
    n = (len(rank) - 1).bit_length()
    assert len(rank) == 1 << n
    closed = [s for s in range(1 << n)
              if all(rank[s | 1 << i] > rank[s] for i in range(n) if not s >> i & 1)]
    # Passing every original rank also detects inconsistent identifications.
    return RankLP(n, closed, dict(enumerate(rank)))


def projection(lp, contracted, retained):
    """All coordinate bounds of (h / contracted)|retained, in binary order.

    Local bit i denotes retained[i]. Optimize the DIFFERENCE, so this also
    works when h(contracted) is not fixed. Never subtract unrelated extrema.
    """
    assert len(set(retained)) == len(retained) and not contracted & mask(retained)
    result = []
    for local in range(1 << len(retained)):
        s = mask(p for i, p in enumerate(retained) if local >> i & 1)
        lo, hi = lp.bounds(((s | contracted, 1), (contracted, -1)))
        result.append((s, lo, hi))
    return result


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


def check_desargues(projected):
    # First check that ALL 1024 coordinates are forced integral ranks.
    assert len(projected) == 1024
    assert all(lo == hi and lo.denominator == 1 for _, lo, hi in projected)
    rank = {s: int(lo) for s, lo, _ in projected}
    points = list(RETAIN)
    assert rank[0] == 0 and rank[mask(points)] == 3
    for s, r in rank.items():
        assert 0 <= r <= s.bit_count()
        for p in points:
            assert rank[s | mask([p])] >= r
        for p, q in combinations(points, 2):
            assert rank[s | mask([p])] + rank[s | mask([q])] >= r + rank[s | mask([p, q])]
    r = lambda *ps: rank[mask(ps)]
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
    print("All corresponding side lines are distinct: r(O,ai,aj)=3 for i != j.")
    print("All nine collinear triples have rank 2, but r(x1,x2,x3)=3.")
    print("The restriction is a simple rank-3 matroid; planar AE Desargues is violated.")


def main(output, write_lps):
    start = perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    nonbases = set(map(mask, NONBASES))
    bases = [mask(b) for b in combinations(range(1, 11), 4) if mask(b) not in nonbases]
    assert len(bases) == 190
    lp = from_rank_table([max((s & b).bit_count() for b in bases) for s in range(1024)])
    print(f"Initial rank-4 matroid: 10 points, {len(lp.closed)} closed sets.", flush=True)
    stages, solves = [], lp.solves
    for lines in IM_LINES:
        new_point = lp.n + 1
        extended, premises = lp.im_extend(tuple(map(mask, lines)))
        print(f"\nAdd {new_point}: h({new_point})=1; h(L+{new_point})=h(L) for L in {lines}.", flush=True)
        print("Premises: " + "; ".join(f"h({labels(s)})={v}" for s, v in premises), flush=True)
        if write_lps:
            extended.write_lp(output / f"extension_{new_point}.lp")
        lp = extended.compress()
        intervals = [(s, *extended.rank_bounds(s)) for s in lp.closed]
        uncertain = sum(lo != hi for _, lo, hi in intervals)
        print(f"LP: {len(extended.closed)} variables, {len(extended.rows)} Shannon rows; "
              f"compressed to {len(lp.closed)} closed sets, {uncertain} undetermined ranks.", flush=True)
        stages.append({"new_point": new_point, "lines": lines,
                       "premises": [{"set": labels(s), "rank": v} for s, v in premises],
                       "closed_set_bounds": [{"set": labels(s), "lower": str(lo), "upper": str(hi)}
                                             for s, lo, hi in intervals]})
        solves += extended.solves + lp.solves
    (output / "extensions.json").write_text(json.dumps(stages, indent=2) + "\n")
    before = lp.solves
    projected = projection(lp, mask(CONTRACT), RETAIN)
    solves += lp.solves - before
    with (output / "projection.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["local_mask", "points", "lower", "upper"])
        for local, (s, lo, hi) in enumerate(projected):
            writer.writerow([local, " ".join(map(str, labels(s))), str(lo), str(hi)])
    print(f"\nContract {CONTRACT}; restrict to {RETAIN}.")
    check_desargues(projected)
    print(f"Wrote six extensions and all 1024 coordinate bounds to {output}.")
    print(f"{solves} LP solves (including feasibility); {perf_counter() - start:.2f} seconds.")


def cli():
    """Entry point for the installed command and direct script execution."""
    if not __debug__:
        raise RuntimeError("Run without -O: assertions are part of the verification.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("im_projection_output"))
    parser.add_argument("--write-lps", action="store_true", help="also export the six LP models")
    args = parser.parse_args()
    main(args.output, args.write_lps)


if __name__ == "__main__":
    cli()
