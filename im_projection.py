"""Six IM extensions of Larson--Le's dual, verified by exact rank bounds.

Run: python im_projection.py [--output im_projection_output]
Only Python's standard library is needed.

Sets are integer bitmasks, with bit i-1 denoting element i. A closure system
identifies h(S) with h(cl(S)); its closed sets need not be matroid flats.
For an intersection-closed family, monotonicity and submodularity on closed
sets imply both inequalities on all subsets: cl(A union B) equals
cl(cl(A) union cl(B)), and cl(A intersection B) is contained in
cl(A) intersection cl(B). Closure compression uses only proved equalities.

Bounds propagate by rearranging these inequalities, using exact arithmetic.
Integer interval endpoints do NOT require unknown ranks to be integers.
Propagation need not find optimal bounds, but proves every equality needed
here. The upper-bound vector is checked as an exact feasible witness at
every stage. No solver, tolerances, saved ranks, or certificates are used.
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


def mask(points):
    return sum(1 << (p - 1) for p in set(points))


def labels(s):
    return [i + 1 for i in range(s.bit_length()) if s >> i & 1]


class RankBounds:
    """Sound rank intervals from monotonicity and submodularity.

    Inputs are exact rational ranks, including the empty and full sets.
    For these data the upper bounds also form a feasible rank function;
    stop if that witness check fails, even if the intervals are consistent.
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
        assert all(isinstance(value, Rational) for value in fixed.values()), \
            "Ranks must be exact rational numbers."
        known = {}
        for s, value in fixed.items():
            f = closure(s)
            assert f not in known or known[f] == value
            known[f] = value
        assert known[closure(0)] == 0
        lo = self.lower = {s: known.get(s, 0) for s in self.closed}
        hi = self.upper = {s: known.get(s, known[self.full]) for s in self.closed}

        inclusions = sorted({(a, closure(a | 1 << i))
                             for a in self.closed for i in range(n)})
        submodular = [(a, b, closure(a | b), a & b)
                      for a, b in combinations(self.closed, 2) if a & b not in (a, b)]
        while True:
            previous = lo.copy(), hi.copy()
            for a, b in inclusions:
                lo[b] = max(lo[b], lo[a])
                hi[a] = min(hi[a], hi[b])
            # h(a) + h(b) >= h(c) + h(d), with c=cl(a union b), d=a intersect b.
            for a, b, c, d in submodular:
                lo[a] = max(lo[a], lo[c] + lo[d] - hi[b])
                lo[b] = max(lo[b], lo[c] + lo[d] - hi[a])
                hi[c] = min(hi[c], hi[a] + hi[b] - lo[d])
                hi[d] = min(hi[d], hi[a] + hi[b] - lo[c])
            assert all(lo[s] <= hi[s] for s in self.closed), "Inconsistent rank bounds."
            if previous == (lo, hi):
                break

        # Here the upper endpoints themselves witness feasibility exactly.
        assert all(hi[a] <= hi[b] for a, b in inclusions)
        assert all(hi[a] + hi[b] >= hi[c] + hi[d] for a, b, c, d in submodular), \
            "The upper bounds must form an exact feasible witness."
        self.fixed = {s: lo[s] for s in self.closed if lo[s] == hi[s]}

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
                        if a != b and a & b == a and fixed[a] == fixed[b]]
        closed = [s for s in self.closed
                  if all(s & a != a or s & b == b for a, b in implications)]
        return RankBounds(self.n, closed, fixed)

    def im_extend(self, lines):
        """Add z of rank one spanned by each of three rank-two sets.

        Check the seven premises for every allowed polymatroid: each line has
        rank 2, each pair of lines has union rank 3, and the total rank is 4.
        The returned bounds include every polymatroid extension satisfying
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
        return RankBounds(self.n + 1, closed, {**self.fixed, z: 1}), premises


def from_rank_table(rank):
    """Start with any rational polymatroid rank table, indexed by subset mask."""
    n = (len(rank) - 1).bit_length()
    assert len(rank) == 1 << n
    closed = [s for s in range(1 << n)
              if all(rank[s | 1 << i] > rank[s] for i in range(n) if not s >> i & 1)]
    # Passing every original rank also detects inconsistent identifications.
    return RankBounds(n, closed, dict(enumerate(rank)))


def projection(ranks, contracted, retained):
    """Bounds of (h / contracted)|retained, when h(contracted) is fixed.

    Local bit i denotes retained[i]. Subtract only the proved fixed rank;
    subtracting separate interval endpoints would not give valid bounds.
    """
    assert len(set(retained)) == len(retained) and not contracted & mask(retained)
    lower, upper = ranks.rank_bounds(contracted)
    assert lower == upper, "The contracted set must have a fixed rank."
    result = []
    for local in range(1 << len(retained)):
        s = mask(p for i, p in enumerate(retained) if local >> i & 1)
        lo, hi = ranks.rank_bounds(s | contracted)
        result.append((s, lo - lower, hi - lower))
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


def main(output):
    start = perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    nonbases = set(map(mask, NONBASES))
    bases = [mask(b) for b in combinations(range(1, 11), 4) if mask(b) not in nonbases]
    assert len(bases) == 190
    ranks = from_rank_table([max((s & b).bit_count() for b in bases) for s in range(1024)])
    print(f"Initial rank-4 matroid: 10 points, {len(ranks.closed)} closed sets.", flush=True)
    stages = []
    for lines in IM_LINES:
        new_point = ranks.n + 1
        extended, premises = ranks.im_extend(tuple(map(mask, lines)))
        print(f"\nAdd {new_point}: h({new_point})=1; h(L+{new_point})=h(L) for L in {lines}.", flush=True)
        print("Premises: " + "; ".join(f"h({labels(s)})={v}" for s, v in premises), flush=True)
        ranks = extended.compress()
        intervals = [(s, *ranks.rank_bounds(s)) for s in ranks.closed]
        uncertain = sum(lo != hi for _, lo, hi in intervals)
        print(f"Closed sets: {len(extended.closed)} before compression, {len(ranks.closed)} after; "
              f"{uncertain} undetermined ranks.", flush=True)
        stages.append({"new_point": new_point, "lines": lines,
                       "premises": [{"set": labels(s), "rank": v} for s, v in premises],
                       "closed_set_bounds": [{"set": labels(s), "lower": str(lo), "upper": str(hi)}
                                             for s, lo, hi in intervals]})
    (output / "extensions.json").write_text(json.dumps(stages, indent=2) + "\n")
    projected = projection(ranks, mask(CONTRACT), RETAIN)
    with (output / "projection.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["local_mask", "points", "lower", "upper"])
        for local, (s, lo, hi) in enumerate(projected):
            writer.writerow([local, " ".join(map(str, labels(s))), str(lo), str(hi)])
    print(f"\nContract {CONTRACT}; restrict to {RETAIN}.")
    check_desargues(projected)
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
