"""Spatial incidence drawings of a rank-four matroid and IM extensions.

All mathematical data are separate from the layout. A saved extensions.json
from im_projection.py can be reduced to the small display data file with
--write-steps. The drawing never selects an integral polymatroid completion.

The layout may replace explicitly named planes by curved sheets. Every other
original dependence is checked against its rational spatial coordinates.
Every displayed line and plane is also checked, with the same exceptions.
SVG and the self-contained, rotatable HTML viewer need only the standard
library. Optional PNG export uses CairoSVG.
"""

import argparse
from fractions import Fraction
from html import escape
from itertools import combinations, cycle
import json
import math
from pathlib import Path

from draw_matroid import load_matroid, svg_text


INK = "#172033"
LINE_COLORS = ("#265eab", "#c27710", "#16836b")
PLANE_COLORS = ("#ada1d0", "#8cb9d0", "#c6b47e")
STEP_COLORS = ("#2463bb", "#bd7408", "#087f8c", "#b23c7d", "#7448b6", "#48852f")


def circuits(matroid):
    """Enumerate minimal dependent sets, whose sizes are at most rank + 1."""
    result = []
    for size in range(1, matroid.r + 2):
        for subset in combinations(matroid.ground, size):
            subset = frozenset(subset)
            if matroid.rank(subset) < size and all(
                    matroid.rank(subset - {e}) == size - 1 for e in subset):
                result.append(sorted(subset))
    return result


def affine_rank(points):
    """Exact rank of homogeneous coordinate rows (x,y,z,1)."""
    rows = [[*(Fraction(x) for x in p), Fraction(1)] for p in points]
    rank = 0
    for column in range(4):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        value = rows[rank][column]
        rows[rank] = [x / value for x in rows[rank]]
        for i in range(rank + 1, len(rows)):
            value = rows[i][column]
            rows[i] = [x - value * y for x, y in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def compact_steps(records, matroid):
    """Extract only the incidences used in the figure from verified rank bounds.

    The compressed closure system identifies each set with its least closed
    superset. Only exactly fixed ranks are accepted here. This reader is a
    display-data check, not a substitute for rerunning the proof.
    """
    result = []
    for record in records:
        point, lines = record["new_point"], record["lines"]
        closed = sorted(record["closed_set_bounds"], key=lambda row: len(row["set"]))

        def fixed_closure(points, rank):
            row = next(row for row in closed if set(points) <= set(row["set"]))
            if not Fraction(row["lower"]) == Fraction(row["upper"]) == rank:
                raise ValueError(f"The saved bounds do not force rank {rank} for {points}.")
            return row["set"]

        fixed_closure([point], 1)
        fixed_closure(sum(lines, []), 4)
        drawn_lines = [dict(generators=L, members=fixed_closure(L + [point], 2)) for L in lines]
        planes = []
        for i, j in combinations(range(3), 2):
            members = fixed_closure(lines[i] + lines[j] + [point], 3)
            original = sorted(set(members) & set(matroid.ground))
            planes.append(dict(lines=[i, j], members=members, original_flat=original))
        result.append(dict(new_point=point, lines=drawn_lines, planes=planes,
                           premises=record["premises"]))
    return {"description": "Selected forced incidences from the six IM extensions; not a complete ambient matroid.",
            "steps": result}


def vector(a, b, t=1):
    return [x + t * y for x, y in zip(a, b)]


def difference(a, b):
    return [x - y for x, y in zip(a, b)]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def unit(v):
    length = math.sqrt(dot(v, v))
    return [x / length for x in v]


def hull(points):
    """Indices of the convex hull of two-dimensional points, counterclockwise."""
    order = sorted(range(len(points)), key=lambda i: tuple(points[i]))

    def turn(i, j, k):
        a, b, c = points[i], points[j], points[k]
        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])

    def half(indices):
        out = []
        for i in indices:
            while len(out) > 1 and turn(out[-2], out[-1], i) <= 1e-10:
                out.pop()
            out.append(i)
        return out
    return half(order)[:-1] + half(reversed(order))[:-1]


def planar_patch(points):
    origin = points[0]
    x = unit(difference(points[1], origin))
    normal = next(cross(x, difference(p, origin)) for p in points[2:]
                  if dot(cross(x, difference(p, origin)), cross(x, difference(p, origin))) > 1e-12)
    y = unit(cross(unit(normal), x))
    local = [[dot(difference(p, origin), x), dot(difference(p, origin), y)] for p in points]
    boundary = [points[i] for i in hull(local)]
    center = [sum(p[j] for p in boundary) / len(boundary) for j in range(3)]
    return [vector(center, difference(p, center), 1.08) for p in boundary]


def curve_points(a, middle, b, count=48):
    """A quadratic arc with B(0)=a, B(1/2)=middle, B(1)=b."""
    control = [2*m - (x+y)/2 for x, m, y in zip(a, middle, b)]
    return [[(1-t)**2*x + 2*t*(1-t)*c + t*t*y for x, c, y in zip(a, control, b)]
            for t in (i/count for i in range(count+1))]


def plane_grid(boundary):
    """A faint, clipped Cartesian grid provides a plane's perspective cues."""
    origin = boundary[0]
    normal = unit(cross(difference(boundary[1], origin), difference(boundary[2], origin)))
    axis = min(range(3), key=lambda i: abs(normal[i]))
    e = [int(i == axis) for i in range(3)]
    x = unit(vector(e, normal, -dot(e, normal)))
    y = unit(cross(normal, x))
    uv = [[dot(difference(p, origin), x), dot(difference(p, origin), y)] for p in boundary]
    wires = []
    for axis in range(2):
        low, high = min(p[axis] for p in uv), max(p[axis] for p in uv)
        for fraction in (.2, .4, .6, .8):
            value, ends = low + fraction * (high - low), []
            for a, b in zip(uv, uv[1:] + uv[:1]):
                if (a[axis] - value) * (b[axis] - value) < 0:
                    ends.append(vector(a, difference(b, a), (value-a[axis])/(b[axis]-a[axis])))
            if len(ends) == 2:
                wires.append([vector(vector(origin, x, p[0]), y, p[1]) for p in ends])
    return wires


def curved_faces(perimeter, surface, subdivisions=10):
    """Triangulate a parameter polygon before applying its curved surface map."""
    center = [sum(p[i] for p in perimeter)/len(perimeter) for i in range(2)]
    result = []
    for a, b in zip(perimeter, perimeter[1:] + perimeter[:1]):
        def point(i, j):
            return surface(*vector(vector(center, difference(a, center), i/subdivisions),
                                   difference(b, center), j/subdivisions))
        for i in range(subdivisions):
            for j in range(subdivisions-i):
                result.append([point(i, j), point(i+1, j), point(i, j+1)])
                if i+j+2 <= subdivisions:
                    result.append([point(i+1, j), point(i+1, j+1), point(i, j+1)])
    return result


def depth_order(faces, toward_viewer):
    """Split intersecting sheets into a BSP tree and paint from back to front.

    Sorting whole transparent planes by their centers makes them look stacked.
    Splitting at their actual intersections allows the front sheet to change
    across an intersection. Curved facets are sorted inside each final cell.
    """
    splitter = next((f for f in faces if not f["curved"]), None)
    if splitter is None:
        return sorted(faces, key=lambda f: sum(dot(p, toward_viewer) for p in f["points"])/len(f["points"]))
    a, b, c = splitter["points"][:3]
    normal = unit(cross(difference(b, a), difference(c, a)))
    offset = dot(a, normal)
    below, on, above = [], [], []
    for face in faces:
        points = face["points"]
        values = [dot(p, normal)-offset for p in points]
        positive, negative = max(values) > 1e-8, min(values) < -1e-8
        if not positive and not negative:
            on.append(face)
        elif not positive:
            below.append(face)
        elif not negative:
            above.append(face)
        else:
            lo, hi = [], []
            for i, p in enumerate(points):
                q, v, w = points[(i+1) % len(points)], values[i], values[(i+1) % len(points)]
                if v <= 1e-8:
                    lo.append(p)
                if v >= -1e-8:
                    hi.append(p)
                if (v > 1e-8 and w < -1e-8) or (v < -1e-8 and w > 1e-8):
                    meet = vector(p, difference(q, p), v/(v-w))
                    lo.append(meet)
                    hi.append(meet)
            below.append(dict(face, points=lo))
            above.append(dict(face, points=hi))
    if dot(normal, toward_viewer) < 0:
        below, above = above, below
    return depth_order(below, toward_viewer) + on + depth_order(above, toward_viewer)


class SpatialDrawing:
    def __init__(self, matroid, steps, layout):
        if matroid.r != 4:
            raise ValueError("The spatial drawing needs a rank-four matroid.")
        self.matroid, self.steps, self.layout = matroid, steps["steps"], layout
        self.exact = {int(e): tuple(map(Fraction, p)) for e, p in layout["positions"].items()}
        self.points = {e: list(map(float, p)) for e, p in self.exact.items()}
        self.planes = [sorted(f) for f in matroid.flats(3) if len(f) >= 4]
        self.curved = {frozenset(p["flat"]): p for p in layout.get("curved_planes", [])}
        self.curved_lines = {frozenset(p["curve"][::2]): p["curve"] for p in self.curved.values()}
        self.check()

    def check(self):
        """Check the mathematical meaning and explicit limits of the spatial model."""
        original = set(self.matroid.ground)
        expected = original | {s["new_point"] for s in self.steps}
        if set(self.exact) != expected or any(len(p) != 3 for p in self.exact.values()):
            raise ValueError("Give exactly three rational coordinates for every original and added point.")
        if any(f not in map(frozenset, self.planes) for f in self.curved):
            raise ValueError("A curved sheet must name an actual original rank-three flat.")
        for flat, sheet in self.curved.items():
            straight, curve = sheet["straight"], sheet["curve"]
            if (len(straight) != 2 or len(curve) != 3 or len(set(curve)) != 3
                    or set(straight) | set(curve[::2]) != flat
                    or set(straight) & set(curve[::2]) or not set(curve) <= expected
                    or affine_rank([self.exact[e] for e in straight + [curve[1]]]) != 2):
                raise ValueError("A curved sheet needs two disjoint spanning pairs and a common middle point.")
        missed = set()
        for size in range(1, 5):
            for subset in combinations(self.matroid.ground, size):
                actual = affine_rank([self.exact[e] for e in subset])
                wanted = self.matroid.rank(frozenset(subset))
                if actual != wanted:
                    if size != 4 or actual != 4 or wanted != 3 or not any(set(subset) <= f for f in self.curved):
                        raise ValueError(f"Unannotated spatial incidence error at {subset}.")
                    missed.add(frozenset(subset))
        self.relaxed_nonbases = sorted(map(sorted, missed))
        known = original.copy()
        for step in self.steps:
            z = step["new_point"]
            if z in known or len(step["lines"]) != 3 or len(step["planes"]) != 3:
                raise ValueError("An IM step must add one fresh point on three lines.")
            generators = [set(L["generators"]) for L in step["lines"]]
            if any(len(L) != 2 or not L <= known for L in generators):
                raise ValueError("Line generators must be pairs of previously available elements.")
            premises = [(L, 2) for L in generators]
            premises += [(a | b, 3) for a, b in combinations(generators, 2)]
            premises += [(set().union(*generators), 4)]
            if {(frozenset(p["set"]), p["rank"]) for p in step["premises"]} != {
                    (frozenset(s), r) for s, r in premises}:
                raise ValueError("An IM step must record precisely the seven rank premises.")
            known.add(z)
            for line in step["lines"]:
                members, pair = set(line["members"]), frozenset(line["generators"])
                if not members <= known or not pair | {z} <= members:
                    raise ValueError("Invalid membership in an IM line.")
                if pair not in self.curved_lines and affine_rank([self.exact[e] for e in members]) != 2:
                    raise ValueError(f"The displayed straight line {sorted(members)} is not straight.")
                if pair in self.curved_lines and members != set(self.curved_lines[pair]):
                    raise ValueError("A curved line must list its two endpoints and specified middle point.")
            for plane, (i, j) in zip(step["planes"], combinations(range(3), 2)):
                f, members = frozenset(plane["original_flat"]), set(plane["members"])
                if (plane["lines"] != [i, j] or f not in map(frozenset, self.planes)
                        or not members <= known or members & original != f
                        or not generators[i] | generators[j] | {z} <= members):
                    raise ValueError("Invalid membership in an IM plane.")
                if f not in self.curved and affine_rank([self.exact[e] for e in members]) != 3:
                    raise ValueError(f"The unannotated plane {sorted(members)} is not flat in space.")
                if f in self.curved and not members <= f | {self.curved[f]["curve"][1]}:
                    raise ValueError("The curved sheet does not cover all the listed elements.")

    def sheet(self, flat, members=None, color="#a7bacb", name="", mesh=False):
        curved = self.curved.get(frozenset(flat))
        wires, faces = [], []
        if curved:
            a, middle, b = [self.points[e] for e in curved["curve"]]
            c, d = [self.points[e] for e in curved["straight"]]
            direction = difference(d, c)
            u = [dot(difference(p, middle), direction) / dot(direction, direction) for p in (c, d)]
            parameters = [[u[0], 0], [u[1], 0], [0, -1], [0, 1]]
            perimeter = [parameters[i] for i in hull(parameters)]

            def surface(u, v):
                linear = [(y-x)/2 for x, y in zip(a, b)]
                quadratic = [(x+y)/2-m for x, y, m in zip(a, b, middle)]
                return [m + u*w + v*x + v*v*y for m, w, x, y in zip(middle, direction, linear, quadratic)]

            boundary = []
            for first, second in zip(perimeter, perimeter[1:] + perimeter[:1]):
                boundary.extend(surface(*vector(first, difference(second, first), i/16)) for i in range(16))
            # Two interior curves make the departure from a flat sheet visible.
            for fraction in (0.35, 0.7):
                outer = max(u, key=abs)
                wires.append([surface(fraction*outer, (1-fraction)*(-1+2*i/24)) for i in range(25)])
            if mesh:
                faces = curved_faces(perimeter, surface)
        else:
            boundary = planar_patch([self.points[e] for e in (members or flat)])
            if mesh:
                faces = [boundary]
        result = dict(boundary=boundary, wires=wires, color=color, curved=bool(curved),
                      flat=flat, name=name or "{" + ",".join(map(str, flat)) + "}")
        if mesh:
            result["faces"] = faces
        return result

    def line(self, record, color, previous=False):
        curve = self.curved_lines.get(frozenset(record["generators"]))
        if curve:
            points = curve_points(*(self.points[e] for e in curve))
        else:
            points = [self.points[e] for e in record["members"]]
            axis = max(range(3), key=lambda j: max(p[j] for p in points)-min(p[j] for p in points))
            points.sort(key=lambda p: p[axis])
            a, b = points[0], points[-1]
            points = [vector(a, difference(a, b), .06), vector(b, difference(b, a), .06)]
        return dict(points=points, color=color, previous=previous,
                    members=record["members"], generators=record["generators"])

    def scene(self, number):
        if number == 0:
            five_colors = cycle(("#d5bd8e", "#9cbbdf"))
            sheets = [self.sheet(f, color=(next(five_colors) if len(f) == 5 else
                                         "#c0a4d5" if frozenset(f) in self.curved else "#b7c4d1"))
                      for f in self.planes]
            default = [i for i, f in enumerate(self.planes) if len(f) == 5 or frozenset(f) in self.curved]
            covered = set().union(*(set(self.planes[i]) for i in default))
            # Show at least one sheet through every original element.
            for i, flat in enumerate(self.planes):
                if not set(flat) <= covered:
                    default.append(i)
                    covered.update(flat)
            visible, lines, title = set(self.matroid.ground), [], "The starting matroid D = N*"
            active, z = set(), None
        else:
            step = self.steps[number-1]
            z = step["new_point"]
            visible = set(self.matroid.ground) | {s["new_point"] for s in self.steps[:number]}
            active = set().union(*(set(L["members"]) for L in step["lines"]))
            sheets = [self.sheet(p["original_flat"], p["members"], PLANE_COLORS[i],
                                  f"L{p['lines'][0]+1} + L{p['lines'][1]+1}")
                      for i, p in enumerate(step["planes"])]
            default = list(range(3))
            lines = [self.line(L, color) for L, color in zip(step["lines"], LINE_COLORS)]
            previous = [self.line(L, "#b7c0cb", previous=True)
                        for earlier in self.steps[:number-1] for L in earlier["lines"]]
            lines = previous + lines
            title = f"{number}. Add point {z}"
        points = [dict(label=e, position=self.points[e], new=e == z,
                       original=e in self.matroid.ground,
                       active=number == 0 or e in active)
                  for e in sorted(visible)]
        return dict(number=number, title=title, point=z, sheets=sheets,
                    default_sheets=default, lines=lines, points=points)

    def overview(self):
        """All original circuit planes and all six colored IM constructions."""
        sheets = []
        for index, flat in enumerate(self.planes):
            members = set(flat)
            for step in self.steps:
                for plane in step["planes"]:
                    if plane["original_flat"] == flat:
                        members.update(plane["members"])
            sheet = self.sheet(flat, sorted(members), mesh=True)
            # Three coordinate planes give the viewer a spatial frame; the
            # remaining planes stay lighter so that all IM lines are legible.
            primary = len(flat) == 5 or flat == [4, 6, 9, 10]
            sheet.update(letter=chr(65+index), primary=primary)
            sheet["color"] = ("#ae8abb" if sheet["curved"] else
                              "#c5a06a" if flat == [2, 5, 8, 9, 10] else
                              "#7fabbc" if flat == [3, 4, 5, 6, 8] else
                              "#90a9bf" if primary else "#98a6b7")
            sheet["opacity"] = .29 if primary else .20 if sheet["curved"] else .055
            if primary:
                sheet["wires"] = plane_grid(sheet["boundary"])
            sheets.append(sheet)
        # A generating pair used twice is still one line. In particular the
        # full line 5,6,12,15 receives both the colors for steps 2 and 5.
        merged = {}
        for step, color in zip(self.steps, STEP_COLORS):
            for record in step["lines"]:
                pair = tuple(record["generators"])
                line = merged.setdefault(pair, dict(generators=list(pair), members=set(), colors=[], points_added=[]))
                line["members"].update(record["members"])
                line["colors"].append(color)
                line["points_added"].append(step["new_point"])
        lines = []
        for record in merged.values():
            record["members"] = sorted(record["members"])
            lines.append(dict(self.line(record, record["colors"][0]),
                              colors=record["colors"], points_added=record["points_added"]))
        return sheets, lines


class Camera:
    def __init__(self, points, view, box):
        yaw, pitch = map(math.radians, (view["yaw"], view["pitch"]))
        self.axes = ([math.cos(yaw), math.sin(yaw), 0],
                     [-math.sin(yaw)*math.sin(pitch), math.cos(yaw)*math.sin(pitch), math.cos(pitch)],
                     [math.sin(yaw)*math.cos(pitch), -math.cos(yaw)*math.cos(pitch), math.sin(pitch)])
        projected = [[dot(p, a) for a in self.axes] for p in points]
        low = [min(p[i] for p in projected) for i in range(2)]
        high = [max(p[i] for p in projected) for i in range(2)]
        self.middle = [(a+b)/2 for a, b in zip(low, high)]
        x, y, width, height = box
        self.scale = min(width/(high[0]-low[0]), height/(high[1]-low[1]))
        self.center = (x+width/2, y+height/2)

    def project(self, point):
        x, y, z = [dot(point, a) for a in self.axes]
        return self.center[0]+self.scale*(x-self.middle[0]), self.center[1]-self.scale*(y-self.middle[1]), z


def path(points, camera, close=False):
    p = [camera.project(point) for point in points]
    return "M " + " L ".join(f"{x:.3f},{y:.3f}" for x, y, _ in p) + (" Z" if close else "")


def svg_start(width, height, title):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img">', f'<title>{escape(title)}</title>',
            '<desc>Spatial incidence schematic. Shaded sheets are rank-three flats. '
            'The explicitly curved sheet is not an ordinary Euclidean plane.</desc>',
            '<g font-family="DejaVu Sans, sans-serif">']


def segment_distance(point, a, b):
    direction = difference(b[:2], a[:2])
    delta = difference(point[:2], a[:2])
    length2 = dot(direction, direction)
    t = max(0, min(1, dot(delta, direction)/length2)) if length2 else 0
    v = difference(delta, [t*x for x in direction])
    return math.hypot(*v)


def crossing_gaps(projected, marked_points):
    """At an unmarked crossing, cut a short gap in the spatially farther line."""
    gaps = [[] for _ in projected]

    def det(a, b):
        return a[0]*b[1]-a[1]*b[0]

    for i, j in combinations(range(len(projected)), 2):
        for a, b in zip(projected[i], projected[i][1:]):
            for c, d in zip(projected[j], projected[j][1:]):
                u, v, delta = difference(b, a), difference(d, c), difference(c, a)
                denominator = det(u, v)
                if abs(denominator) < 1e-9:
                    continue
                t, s = det(delta, v)/denominator, det(delta, u)/denominator
                if not (0 <= t <= 1 and 0 <= s <= 1):
                    continue
                p, q = vector(a, u, t), vector(c, v, s)
                if abs(p[2]-q[2]) < 1e-7 or any(math.dist(p[:2], m[:2]) < 18 for m in marked_points):
                    continue
                farther = i if p[2] < q[2] else j
                angle = abs(denominator)/(math.hypot(*u[:2])*math.hypot(*v[:2]))
                gaps[farther].append((p[0], p[1], min(11, 4/max(angle, .1))))
    return gaps


def gapped_path(projected, holes):
    """Trim projected segments by circular gaps, using portable SVG subpaths.

    Explicit subpaths give identical crossings in SVG and PNG; luminance
    masks are interpreted differently by different SVG renderers.
    """
    commands, previous = [], None
    for a, b in zip(projected, projected[1:]):
        a, b = a[:2], b[:2]
        direction = difference(b, a)
        length2 = dot(direction, direction)
        if length2 < 1e-18:
            continue
        intervals = [(0, 1)]
        for x, y, radius in holes:
            delta = difference(a, [x, y])
            linear = dot(delta, direction)
            discriminant = linear*linear - length2*(dot(delta, delta)-radius*radius)
            if discriminant <= 0:
                continue
            root = math.sqrt(discriminant)
            low, high = (-linear-root)/length2, (-linear+root)/length2
            remaining = []
            for start, stop in intervals:
                if high <= start or low >= stop:
                    remaining.append((start, stop))
                else:
                    if start < low:
                        remaining.append((start, low))
                    if high < stop:
                        remaining.append((high, stop))
            intervals = remaining
        for start, stop in intervals:
            first, last = vector(a, direction, start), vector(a, direction, stop)
            if previous is None or math.dist(first, previous) > 1e-6:
                commands.append(f'M {first[0]:.3f},{first[1]:.3f}')
            commands.append(f'L {last[0]:.3f},{last[1]:.3f}')
            previous = last
    return " ".join(commands)


def plane_tags(sheets, lines, points, camera):
    """Place A–L just outside sheet edges, away from point labels and lines."""
    occupied = [(x+12, y-10) for x, y, _ in points]
    segments = [(a, b) for line in lines for a, b in zip(line, line[1:])]
    tags = []
    for sheet in sheets:
        boundary = [camera.project(p)[:2] for p in sheet["boundary"]]
        center = [sum(p[i] for p in boundary)/len(boundary) for i in range(2)]
        candidates = []
        stride = max(1, len(boundary)//12)
        for a, b in zip(boundary[::stride], boundary[::stride][1:] + boundary[::stride][:1]):
            for fraction in (.2, .5, .8):
                p = vector(a, difference(b, a), fraction)
                direction = unit(difference(p, center))
                p = vector(p, direction, 12)
                clearance = min([math.dist(p, q)-34 for q in occupied] +
                                [segment_distance(p, a, b)-13 for a, b in segments])
                candidates.append((clearance, p))
        _, p = max(candidates, key=lambda item: item[0])
        occupied.append(p)
        tags.append((sheet, p))
    return tags


def render_overview(drawing, sheets, lines, camera):
    out = []
    faces = [dict(points=face, color=s["color"], opacity=s["opacity"], curved=s["curved"])
             for s in sheets for face in s["faces"]]
    for face in depth_order(faces, camera.axes[2]):
        out.append(f'<path d="{path(face["points"], camera, True)}" fill="{face["color"]}" '
                   f'fill-opacity="{face["opacity"]}" stroke="none"/>')
    for sheet in sheets:
        members = " ".join(map(str, sheet["flat"]))
        out.append(f'<path data-plane="{members}" data-letter="{sheet["letter"]}" '
                   f'd="{path(sheet["boundary"], camera, True)}" fill="none" '
                   f'stroke="{sheet["color"]}" stroke-opacity=".7" stroke-width="1.1">'
                   f'<title>Plane {sheet["letter"]}: {sheet["name"]}</title></path>')
        for wire in sheet["wires"]:
            out.append(f'<path d="{path(wire, camera)}" fill="none" stroke="{sheet["color"]}" '
                       'stroke-opacity=".45" stroke-width=".7"/>')
    projected = [[camera.project(p) for p in line["points"]] for line in lines]
    points = [camera.project(p) for p in drawing.points.values()]
    gaps = crossing_gaps(projected, points)
    paths = [gapped_path(p, holes) for p, holes in zip(projected, gaps)]
    # All white under-strokes precede all colors, so they cannot erase a line
    # that should be in front. The gaps carry the depth at each crossing.
    for data in paths:
        out.append(f'<path d="{data}" fill="none" stroke="white" '
                   'stroke-width="5.5" stroke-opacity=".85"/>')
    for line, data in zip(lines, paths):
        members = " ".join(map(str, line["members"]))
        added = " ".join(map(str, line["points_added"]))
        pair = " ".join(map(str, line["generators"]))
        out.append(f'<g data-line="{members}" data-generators="{pair}" data-im-points="{added}">')
        for index, color in enumerate(line["colors"]):
            dash = ' stroke-dasharray="9 9"' if index else ''
            out.append(f'<path d="{data}" fill="none" stroke="{color}" stroke-width="3.2"{dash}/>')
        out.append('</g>')
    for sheet, (x, y) in plane_tags(sheets, projected, points, camera):
        out.append(f'<g data-plane-label="{sheet["letter"]}"><rect x="{x-11:.3f}" y="{y-12:.3f}" '
                   f'width="22" height="24" rx="4" fill="white" fill-opacity=".94" '
                   f'stroke="{sheet["color"]}" stroke-opacity=".65"/>')
        out.append(svg_text(x-5.7, y+5.6, sheet["letter"], 16, "#586b7c", halo=False) + '</g>')
    for e, p in sorted(drawing.points.items()):
        x, y, _ = camera.project(p)
        original = e in drawing.matroid.ground
        color = INK if original else STEP_COLORS[next(i for i, s in enumerate(drawing.steps) if s["new_point"] == e)]
        if not original:
            out.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="12" fill="white"/>')
        out.append(f'<circle data-point="{e}" cx="{x:.3f}" cy="{y:.3f}" r="{5.5 if original else 7.8}" '
                   f'fill="{color if original else "white"}" stroke="{color}" stroke-width="2.5"/>')
        dx, dy = drawing.layout.get("label_offsets", {}).get(str(e), [13, -13])
        out.append(svg_text(x+dx, y+dy, e, 22 if original else 24, color))
    return "\n".join(out)


def original_svg(drawing):
    """One publication figure: original planes, all circuits, and all IM points."""
    all_circuits = circuits(drawing.matroid)
    four = [c for c in all_circuits if len(c) == 4]
    five = [c for c in all_circuits if len(c) == 5]
    covered = {tuple(c) for flat in drawing.planes for c in combinations(flat, 4)}
    if len(four)+len(five) != len(all_circuits) or covered != set(map(tuple, four)):
        raise ValueError("The compact circuit description requires a rank-four paving matroid.")
    sheets, lines = drawing.overview()
    visible = list(drawing.points.values()) + [p for sheet in sheets for p in sheet["boundary"]]
    camera = Camera(visible, drawing.layout["view"], (65, 140, 1140, 930))
    out = svg_start(1680, 1320, "Larson–Le’s dual: all circuit planes and six IM extensions")
    out += ['<rect width="1680" height="1320" fill="white"/>',
            '<metadata id="matroid-circuits">' + escape(json.dumps(all_circuits)) + '</metadata>',
            svg_text(45, 52, "D = N* with its six Ingleton–Main extensions", 32, halo=False),
            svg_text(45, 85, "Rank 4  ·  10 original elements  ·  6 added points  ·  all 12 nontrivial planes", 18, "#596b7d", halo=False),
            render_overview(drawing, sheets, lines, camera),
            svg_text(1280, 145, "IM points and generating pairs", 19, halo=False)]
    for i, (step, color) in enumerate(zip(drawing.steps, STEP_COLORS)):
        y = 190 + 64*i
        out.append(f'<circle cx="1290" cy="{y-7}" r="7" fill="white" stroke="{color}" stroke-width="2.5"/>')
        out.append(svg_text(1310, y, step["new_point"], 24, color, halo=False))
        pairs = ["{" + ",".join(map(str, line["generators"])) + "}" for line in step["lines"]]
        out.append(svg_text(1310, y+25, "   ".join(pairs), 17, color, halo=False))
    out += [svg_text(1280, 579, "The line 5–6 carries both 12 and 15.", 16, "#596b7d", halo=False),
            svg_text(1280, 608, "Within each color: pairwise rank 3;", 16, "#596b7d", halo=False),
            svg_text(1280, 631, "the three lines together have rank 4.", 16, "#596b7d", halo=False),
            svg_text(1280, 688, "A–L: all 12 nontrivial planes", 19, halo=False)]
    for i, sheet in enumerate(sheets):
        y = 724+28*i
        out.append(f'<rect x="1280" y="{y-16}" width="23" height="22" rx="3" '
                   f'fill="{sheet["color"]}" fill-opacity=".22"/>')
        out.append(svg_text(1286, y, sheet["letter"], 16, "#586b7c", halo=False))
        out.append(svg_text(1318, y, sheet["name"], 17, halo=False))
        if sheet["curved"]:
            out.append(svg_text(1492, y, "curved", 15, "#866496", halo=False))
    out += ['<path d="M 45,1130 H 1635" stroke="#d6dee7"/>',
            svg_text(45, 1167, f"Circuits of D:  {len(four)} of size 4;  {len(five)} of size 5.", 20, halo=False),
            svg_text(45, 1199, "4-element circuits: all four-element subsets of a plane A–L.  5-element circuits: all five-element subsets containing none of these 4-circuits.", 17, halo=False),
            svg_text(45, 1237, "Filled points 1–10 are original; open colored points 11–16 are IM additions. Unmarked crossings are not elements.", 16, "#596b7d", halo=False),
            svg_text(45, 1266, "Only L = {6,7,8,10} and the line 7–14–10 are curved. The curved sheet encodes an abstract plane, not Euclidean coplanarity.", 16, "#596b7d", halo=False),
            svg_text(45, 1295, "The extensions are shown only through their selected forced incidences; this diagram does not choose integral matroid completions.", 16, "#596b7d", halo=False),
            '</g></svg>']
    return "\n".join(out) + "\n"


def interactive_html(drawing):
    data = dict(scenes=[drawing.scene(i) for i in range(7)], points=list(drawing.points.values()),
                view=drawing.layout["view"], steps=drawing.steps)
    encoded = json.dumps(data, separators=(",", ":")).replace("<", "\\u003c")
    return VIEWER.replace("/*SCENE_DATA*/", encoded)


VIEWER = r'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Larson–Le: six IM extensions</title>
<style>
:root{color-scheme:light;font-family:system-ui,sans-serif;color:#172033;background:#f7f8fb}
body{margin:0;padding:24px;max-width:1500px;margin-inline:auto}h1{font-size:26px;margin:0 0 8px}
p{line-height:1.55}button{font:inherit;color:inherit;background:white;border:1px solid #cad2dd;border-radius:7px;padding:9px 13px;cursor:pointer}
button[aria-pressed=true]{background:#243d61;color:white;border-color:#243d61}button:focus-visible{outline:3px solid #b53b69;outline-offset:2px}
nav{display:flex;gap:8px;flex-wrap:wrap;margin:20px 0}.layout{display:grid;grid-template-columns:minmax(400px,1fr) 300px;gap:24px}
.plot{background:white;border:1px solid #dce2eb;border-radius:12px;overflow:hidden}svg{display:block;width:100%;touch-action:none;cursor:grab;user-select:none}
aside{font-size:15px}h2{font-size:20px;margin-top:0}.plane-list{display:flex;flex-wrap:wrap;gap:6px}.plane-list button{font-size:13px;padding:7px}
.line-label{font-weight:600;margin:10px 0}.equation{font-family:Georgia,serif;font-size:18px;line-height:1.9}.subtle{color:#596779;font-size:14px}.controls{display:flex;gap:14px;align-items:center;flex-wrap:wrap;margin-top:14px}
.note{border-left:3px solid #ac94c6;padding-left:12px}a{color:#265eab}@media(max-width:850px){body{padding:14px}.layout{grid-template-columns:1fr}.plot svg{min-height:400px}aside{display:grid;grid-template-columns:1fr 1fr;gap:16px}}
</style>
<h1>Larson–Le’s dual and its six IM extensions</h1>
<p class="subtle">Drag to rotate the spatial schematic. Select an extension, then inspect its three pairwise planes.</p>
<nav id="steps" aria-label="Extension step"></nav>
<div class="layout"><div><div class="plot"><svg id="view" viewBox="0 0 1000 780" role="img" aria-label="Rotatable spatial incidence diagram"></svg></div>
<div class="controls"><button id="reset">Reset view</button><button id="download">Save this view as SVG</button><label><input id="earlier" type="checkbox"> Earlier extension lines</label></div></div>
<aside><section><h2 id="title"></h2><div id="lines"></div><div id="ranks" class="equation"></div></section>
<section><h3>Shaded planes</h3><div id="planes" class="plane-list"></div><p id="members" class="subtle"></p></section>
<section class="note"><p>The sheet {6,7,8,10} is curved. At step 4, the line through 7, 14 and 10 is curved too.</p><p class="subtle">These are abstract incidences, not a straight-line realization. Only the selected forced incidences of the extensions are shown; their other ranks may be fractional.</p></section>
<p class="subtle">Filled dots: original elements.<br>Open dots: earlier additions.<br>Magenta: the newly added point.<br>Unmarked crossings are not elements.</p>
</aside></div>
<script>
"use strict";
const data=/*SCENE_DATA*/;
const svg=document.getElementById("view"), colors=["#265eab","#c27710","#16836b"];
let current=1, selected=null, yaw=data.view.yaw*Math.PI/180, pitch=data.view.pitch*Math.PI/180;
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0), esc=s=>String(s).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll('"',"&quot;");
function camera(){
 const axes=[[Math.cos(yaw),Math.sin(yaw),0],[-Math.sin(yaw)*Math.sin(pitch),Math.cos(yaw)*Math.sin(pitch),Math.cos(pitch)],[Math.sin(yaw)*Math.cos(pitch),-Math.cos(yaw)*Math.cos(pitch),Math.sin(pitch)]];
 const points=data.points.map(p=>axes.map(a=>dot(p,a))), low=[0,1].map(i=>Math.min(...points.map(p=>p[i]))), high=[0,1].map(i=>Math.max(...points.map(p=>p[i])));
 const scale=Math.min(860/(high[0]-low[0]),620/(high[1]-low[1]));
 return p=>{const q=axes.map(a=>dot(p,a));return [500+scale*(q[0]-(low[0]+high[0])/2),390-scale*(q[1]-(low[1]+high[1])/2),q[2]]};
}
function render(){
 const scene=data.scenes[current], project=camera(), path=(points,closed=false)=>"M "+points.map(p=>project(p).slice(0,2).map(x=>x.toFixed(3)).join(",")).join(" L ")+(closed?" Z":"");
 const indices=selected===null?scene.default_sheets:[selected];
 const sheets=indices.map(i=>scene.sheets[i]).sort((a,b)=>a.boundary.reduce((s,p)=>s+project(p)[2],0)/a.boundary.length-b.boundary.reduce((s,p)=>s+project(p)[2],0)/b.boundary.length);
 let out=`<title>${esc(scene.title)}</title><g font-family="DejaVu Sans,sans-serif">`;
 for(const sheet of sheets){out+=`<path d="${path(sheet.boundary,true)}" fill="${sheet.color}" fill-opacity=".23" stroke="${sheet.color}" stroke-opacity=".7" stroke-width="1.4"/>`;for(const wire of sheet.wires)out+=`<path d="${path(wire)}" fill="none" stroke="${sheet.color}" stroke-opacity=".6"/>`;}
 for(const line of scene.lines){if(line.previous&&!document.getElementById("earlier").checked)continue;out+=`<path d="${path(line.points)}" fill="none" stroke="${line.color}" stroke-width="${line.previous?1.2:3.5}" stroke-opacity="${line.previous?.55:1}" stroke-linecap="round"/>`;}
 for(const p of [...scene.points].sort((a,b)=>a.new-b.new)){const [x,y]=project(p.position),color=p.new?"#b53b69":p.active?"#172033":"#8793a4";if(p.new)out+=`<circle cx="${x}" cy="${y}" r="16" fill="${color}" fill-opacity=".14"/>`;out+=`<circle cx="${x}" cy="${y}" r="${p.new?8:5.5}" fill="${p.original||p.new?color:'white'}" stroke="${color}" stroke-width="1.5"/><text x="${x+12}" y="${y-11}" font-size="${p.new?24:20}" fill="${color}" stroke="white" stroke-width="5" paint-order="stroke">${p.label}</text>`;}
 svg.innerHTML=out+"</g>";
}
function choosePlane(i){selected=i;document.querySelectorAll("#planes button").forEach((b,j)=>b.setAttribute("aria-pressed",String(j===(i===null?0:i+1))));const s=data.scenes[current];document.getElementById("members").textContent=i===null?(current?"Each sheet contains two of the highlighted lines. The three lines together have rank 4.":"The overview shows the two five-point planes, the curved sheet, and a plane through point 1."):"Original plane: {"+s.sheets[i].flat.join(", ")+"}";render();}
function chooseStep(i){current=i;selected=null;const scene=data.scenes[i];document.getElementById("title").textContent=scene.title;document.querySelectorAll("#steps button").forEach((b,j)=>b.setAttribute("aria-pressed",String(i===j)));
 document.getElementById("lines").innerHTML=i?data.steps[i-1].lines.map((l,j)=>`<div class="line-label" style="color:${colors[j]}">L${j+1} = cl{${l.generators.join(", ")}}</div>`).join(""):"<p>Rank 4 · 10 elements · 190 bases.<br>All 12 nontrivial planes can be selected below.</p>";
 document.getElementById("ranks").innerHTML=i?`h(Lᵢ) = 2<br>h(Lᵢ ∪ Lⱼ) = 3 &nbsp; (i ≠ j)<br>h(L₁ ∪ L₂ ∪ L₃) = 4<br><span style="color:#b53b69">${scene.point} ∈ cl(L₁) ∩ cl(L₂) ∩ cl(L₃)<br>h({${scene.point}}) = 1</span>`:"Every original line has two points. All remaining planes have three points.";
 const buttons=document.getElementById("planes");buttons.replaceChildren();[i?"All three":"Overview",...scene.sheets.map(s=>s.name)].forEach((label,j)=>{const b=document.createElement("button");b.textContent=label;b.onclick=()=>choosePlane(j?j-1:null);buttons.append(b)});choosePlane(null);
}
["Original",...data.steps.map((s,i)=>`${i+1}: add ${s.new_point}`)].forEach((label,i)=>{const b=document.createElement("button");b.textContent=label;b.onclick=()=>chooseStep(i);document.getElementById("steps").append(b)});
let drag=null;svg.addEventListener("pointerdown",e=>{drag=[e.clientX,e.clientY];svg.setPointerCapture(e.pointerId)});svg.addEventListener("pointermove",e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*.007;pitch=Math.max(-1.35,Math.min(1.35,pitch+(e.clientY-drag[1])*.007));drag=[e.clientX,e.clientY];render()});svg.addEventListener("pointerup",()=>drag=null);svg.addEventListener("pointercancel",()=>drag=null);
document.getElementById("reset").onclick=()=>{yaw=data.view.yaw*Math.PI/180;pitch=data.view.pitch*Math.PI/180;render()};document.getElementById("earlier").onchange=render;
document.addEventListener("keydown",e=>{if(e.key==="ArrowRight"||e.key==="ArrowLeft"){e.preventDefault();chooseStep((current+(e.key==="ArrowRight"?1:6))%7)}});
document.getElementById("download").onclick=()=>{const copy=svg.cloneNode(true);copy.setAttribute("xmlns","http://www.w3.org/2000/svg");copy.setAttribute("width","1000");copy.setAttribute("height","780");const url=URL.createObjectURL(new Blob([copy.outerHTML],{type:"image/svg+xml"}));const a=document.createElement("a");a.href=url;a.download=`rank4-step-${current}.svg`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};
chooseStep(1);
</script></html>
'''


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matroid", type=Path, help="complete original basis JSON")
    parser.add_argument("steps", type=Path, help="extensions.json from the verifier, or compact display data")
    parser.add_argument("--layout", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=Path("figures"))
    parser.add_argument("--write-steps", type=Path, help="write compact display data extracted from extensions.json")
    parser.add_argument("--png", action="store_true")
    parser.add_argument("--scale", type=float, default=2)
    args = parser.parse_args()
    matroid = load_matroid(args.matroid)
    steps = json.loads(args.steps.read_text())
    if isinstance(steps, list):
        steps = compact_steps(steps, matroid)
    drawing = SpatialDrawing(matroid, steps, json.loads(args.layout.read_text()))
    args.output.mkdir(parents=True, exist_ok=True)
    for name, content in (("rank4_original.svg", original_svg(drawing)),
                          ("rank4_interactive.html", interactive_html(drawing))):
        (args.output / name).write_text(content, encoding="utf-8")
        if args.png and name.endswith(".svg"):
            import cairosvg
            if not math.isfinite(args.scale) or args.scale <= 0:
                parser.error("--scale must be positive and finite")
            cairosvg.svg2png(bytestring=content.encode(), write_to=str(args.output / name.replace(".svg", ".png")), scale=args.scale)
    if args.write_steps:
        args.write_steps.parent.mkdir(parents=True, exist_ok=True)
        args.write_steps.write_text(json.dumps(steps, indent=2) + "\n")
    print(f"Checked all original ranks through size four, {3*len(drawing.steps)} IM lines, and {3*len(drawing.steps)} pairwise planes.")
    print(f"Explicitly curved original nonbases: {drawing.relaxed_nonbases}")
    print(f"Wrote the combined spatial drawing and the offline viewer to {args.output}.")


if __name__ == "__main__":
    cli()
