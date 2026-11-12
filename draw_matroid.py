"""Draw a small matroid of rank at most four as an SVG incidence diagram.

Input: JSON with integer ground_set labels and a complete list of bases, or
the verifier's projection.csv. Only the Python standard library is needed.
Nontrivial rank-two flats are drawn as curves. For rank four, nontrivial
rank-three flats get separate plane panels, avoiding ambiguous 3D projection.
Positions and Bezier controls may be supplied in a separate layout JSON.
Flats are computed from the bases. Optional dashed guides are annotations,
not incidences. SVG output has a transparent background; PNG needs CairoSVG.
"""

import argparse
import csv
import json
import math
from fractions import Fraction
from functools import lru_cache
from html import escape
from itertools import combinations
from pathlib import Path


class Matroid:
    def __init__(self, ground_set, bases, name="Matroid"):
        self.ground = tuple(ground_set)
        if any(type(e) is not int for e in self.ground):
            raise ValueError("Ground-set labels must be integers.")
        if len(set(self.ground)) != len(self.ground):
            raise ValueError("Ground-set labels must be distinct.")
        self.name = str(name)
        rows = [tuple(b) for b in bases]
        if not rows or any(len(set(b)) != len(b) or any(type(e) is not int for e in b) for b in rows):
            raise ValueError("Supply a nonempty basis family without repeated elements.")
        self.bases = frozenset(map(frozenset, rows))
        self.r = len(rows[0])
        if self.r > 4:
            raise ValueError("This renderer supports ranks zero through four.")
        if len(self.bases) != len(rows) or any(
                len(b) != self.r or not b <= set(self.ground) for b in self.bases):
            raise ValueError("Bases must be distinct, equicardinal subsets of the ground set.")
        # Basis exchange is checked before interpreting the input as a matroid.
        for a in self.bases:
            for b in self.bases:
                for e in a - b:
                    if not any((a - {e}) | {f} in self.bases for f in b - a):
                        raise ValueError(f"Basis exchange fails for {sorted(a)}, {sorted(b)}, {e}.")

    @lru_cache(None)
    def rank(self, subset):
        return max(len(subset & b) for b in self.bases)

    def flats(self, rank):
        """Every rank-k flat is the closure of some independent k-subset."""
        result = set()
        for row in combinations(self.ground, rank):
            s = frozenset(row)
            if self.rank(s) == rank:
                result.add(frozenset(e for e in self.ground if self.rank(s | {e}) == rank))
        return sorted(result, key=lambda f: tuple(e for e in self.ground if e in f))

    def as_dict(self):
        return {"name": self.name, "ground_set": list(self.ground),
                "bases": sorted(map(sorted, self.bases))}


def load_matroid(path):
    if path.suffix.lower() != ".csv":
        data = json.loads(path.read_text(encoding="utf-8"))
        return Matroid(data["ground_set"], data["bases"], data.get("name", path.stem))
    ranks = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            points = list(map(int, row["points"].split()))
            s = frozenset(points)
            lo, hi = Fraction(row["lower"]), Fraction(row["upper"])
            if lo != hi or lo.denominator != 1:
                raise ValueError("Projection ranks must be exact, equal integer bounds.")
            if len(s) != len(points) or s in ranks:
                raise ValueError("Repeated elements or subset in the projection table.")
            ranks[s] = int(lo)
    ground = sorted(set().union(*ranks))
    if len(ranks) != 1 << len(ground):
        raise ValueError("The projection must contain every subset exactly once.")
    rank = ranks[frozenset(ground)]
    bases = [sorted(s) for s, value in ranks.items() if len(s) == value == rank]
    matroid = Matroid(ground, bases, path.stem)
    if any(matroid.rank(s) != value for s, value in ranks.items()):
        raise ValueError("The projection table is not the rank function of its proposed bases.")
    return matroid


def xy(value):
    x, y = map(float, value)
    if not math.isfinite(x) or not math.isfinite(y):
        raise ValueError("Layout coordinates must be finite.")
    return x, y


def circle_positions(points, center, radius):
    if len(points) == 1:
        return {points[0]: center}
    return {p: (center[0] + radius * math.cos(2 * math.pi * i / len(points) - math.pi / 2),
                center[1] + radius * math.sin(2 * math.pi * i / len(points) - math.pi / 2))
            for i, p in enumerate(points)}


def line_path(order, positions, controls=None):
    """Join the listed points, with optional cubic Bezier controls per segment."""
    controls = controls if controls is not None else [None] * (len(order) - 1)
    if len(controls) != len(order) - 1:
        raise ValueError("A curve needs one control entry per segment.")
    x, y = positions[order[0]]
    commands = [f"M {x:g},{y:g}"]
    for point, control in zip(order[1:], controls):
        x, y = positions[point]
        if control is None:
            commands.append(f"L {x:g},{y:g}")
        else:
            first, second = map(xy, control)
            commands.append(f"C {first[0]:g},{first[1]:g} {second[0]:g},{second[1]:g} {x:g},{y:g}")
    return " ".join(commands)


def svg_text(x, y, value, size=17, color="#172033", halo=True):
    attributes = f'x="{x:g}" y="{y:g}" font-size="{size}"'
    content = escape(str(value))
    # A white halo keeps text readable over light or dark backgrounds.
    # Two layers work in SVG renderers that do not implement paint-order.
    outline = (f'<text {attributes} fill="none" stroke="white" stroke-width="5" '
               f'stroke-linejoin="round">{content}</text>') if halo else ""
    return outline + f'<text {attributes} fill="{escape(color)}">{content}</text>'


def draw_svg(matroid, layout=None):
    layout = layout or {}
    width, main_height = xy(layout.get("size", [900, 800]))
    if min(width, main_height) < 400:
        raise ValueError("The drawing canvas must be at least 400 by 400.")
    loops = frozenset(e for e in matroid.ground if matroid.rank(frozenset([e])) == 0)
    classes = [f - loops for f in matroid.flats(1)]
    representatives = [next(e for e in matroid.ground if e in f) for f in classes]
    point_labels = {p: " = ".join(map(str, sorted(f))) for p, f in zip(representatives, classes)}
    project = lambda f: tuple(p for p in representatives if p in f)
    lines = [project(f) for f in matroid.flats(2) if len(project(f)) >= 3]
    planes = [project(f) for f in matroid.flats(3)
              if matroid.r == 4 and len(project(f)) >= 4]
    positions = circle_positions(representatives, (width / 2, main_height / 2 - 30),
                                 min(width, main_height) * 0.32)
    if "positions" in layout:
        positions = {int(p): xy(v) for p, v in layout["positions"].items()}
        if set(positions) != set(representatives):
            raise ValueError("Positions must specify each parallel-class representative exactly once.")
    labels = {**point_labels, **{int(p): str(v) for p, v in layout.get("labels", {}).items()}}
    offsets = {int(p): xy(v) for p, v in layout.get("label_offsets", {}).items()}
    colors = {int(p): str(v) for p, v in layout.get("point_colors", {}).items()}
    if (set(labels) | set(offsets) | set(colors)) - set(representatives):
        raise ValueError("Unknown point in the layout labels, offsets, or colors.")
    styles = {}
    for style in layout.get("lines", []):
        key = frozenset(style["flat"])
        if key not in map(frozenset, lines) or key in styles:
            raise ValueError("A line style must name a distinct nontrivial rank-two flat.")
        styles[key] = style
    panel_height = 360
    total_height = main_height + panel_height * math.ceil(len(planes) / 2)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:g}" height="{total_height:g}" '
           f'viewBox="0 0 {width:g} {total_height:g}" role="img">',
           f'<title>{escape(matroid.name)}</title>',
           '<desc>Matroid incidence diagram. Only marked points are elements. '
           'Curves mark rank-two flats; rank-four plane panels mark rank-three flats.</desc>',
           '<g font-family="DejaVu Sans, sans-serif">',
           svg_text(30, 38, matroid.name, 25),
           svg_text(30, 67, f"Rank {matroid.r}; {len(matroid.ground)} elements; {len(matroid.bases)} bases", 17)]

    def panel(points, coords, visible_lines, custom=False):
        if custom:
            for triangle in layout.get("triangles", []):
                vertices = triangle["points"]
                if (len(vertices) != 3 or len(set(vertices)) != 3
                        or not set(vertices) <= set(points)
                        or matroid.rank(frozenset(vertices)) != 3):
                    raise ValueError("A highlighted triangle must be an independent triple.")
                polygon = " ".join(f"{coords[p][0]:g},{coords[p][1]:g}" for p in vertices)
                color = escape(triangle["color"])
                members = " ".join(map(str, vertices))
                svg.append(f'<polygon data-triangle="{members}" points="{polygon}" '
                           f'fill="{color}" fill-opacity="0.10"/>')
        for line in visible_lines:
            style = styles.get(frozenset(line), {}) if custom else {}
            order = style.get("order", sorted(line, key=lambda p: coords[p]))
            if len(order) != len(line) or set(order) != set(line):
                raise ValueError("A curve's order must contain exactly the points of its flat.")
            path = line_path(order, coords, style.get("controls"))
            color = escape(style.get("color", "#64748b"))
            line_width = float(style.get("width", 2.6))
            if not math.isfinite(line_width) or line_width <= 0:
                raise ValueError("Line widths must be positive and finite.")
            members = " ".join(map(str, sorted(line)))
            svg.append(f'<path data-flat="{members}" d="{path}" fill="none" '
                       f'stroke="{color}" stroke-width="{line_width:g}" stroke-linejoin="round" '
                       f'stroke-linecap="round"><title>Rank-two flat: {members}</title></path>')
        if custom:
            # Outlines distinguish the triangle sides from their extensions.
            for triangle in layout.get("triangles", []):
                if triangle.get("outline"):
                    vertices = triangle["points"]
                    polygon = " ".join(f"{coords[p][0]:g},{coords[p][1]:g}" for p in vertices)
                    color = escape(triangle["color"])
                    svg.append(f'<polygon points="{polygon}" fill="none" stroke="{color}" '
                               'stroke-width="3.5" stroke-linejoin="round"/>')
            for guide in layout.get("guides", []):
                order = guide["points"]
                if len(order) < 2 or len(set(order)) != len(order) or not set(order) <= set(points):
                    raise ValueError("A dashed guide needs distinct known points.")
                path = line_path(order, coords, guide.get("controls"))
                if "extend" in guide:
                    if len(order) != 2 or guide.get("controls"):
                        raise ValueError("Only a straight two-point guide can be extended.")
                    before, after = xy(guide["extend"])
                    if min(before, after) < 0:
                        raise ValueError("Guide extension lengths must be nonnegative.")
                    a, b = (coords[p] for p in order)
                    distance = math.dist(a, b)
                    if distance == 0:
                        raise ValueError("A straight guide needs distinct positions.")
                    direction = tuple((b[i] - a[i]) / distance for i in range(2))
                    start = tuple(a[i] - before * direction[i] for i in range(2))
                    end = tuple(b[i] + after * direction[i] for i in range(2))
                    path = f"M {start[0]:g},{start[1]:g} L {end[0]:g},{end[1]:g}"
                color = escape(guide.get("color", "#be3455"))
                members = " ".join(map(str, order))
                svg.append(f'<path data-guide="{members}" d="{path}" fill="none" '
                           f'stroke="{color}" stroke-width="2.8" stroke-dasharray="9 7" '
                           f'stroke-linecap="round"><title>{escape(guide.get("label", "Guide, not a matroid line"))}</title></path>')
        for p in points:
            x, y = coords[p]
            dx, dy = offsets.get(p, (11, -11)) if custom else (9, -9)
            color = escape(colors.get(p, "#172033")) if custom else "#172033"
            svg.append(f'<circle data-point="{p}" cx="{x:g}" cy="{y:g}" r="5.5" '
                       f'fill="{color}" stroke="white" stroke-width="1.6"/>')
            svg.append(svg_text(x + dx, y + dy, labels[p] if custom else point_labels[p],
                                18 if custom else 14, halo=True))

    panel(representatives, positions, lines, custom=True)
    for annotation in layout.get("annotations", []):
        x, y = xy(annotation["position"])
        svg.append(svg_text(x, y, annotation["text"], annotation.get("size", 17),
                            annotation.get("color", "#172033"), halo=True))
    if loops:
        svg.append(svg_text(30, main_height - 110, "Loops: " + ", ".join(map(str, sorted(loops)))))
    legend = layout.get("legend", [])
    for i, item in enumerate(legend):
        x, y = 30 + i * (width - 60) / max(1, len(legend)), main_height - 60
        dash = ' stroke-dasharray="6 4"' if item.get("dashed") else ""
        svg.append(f'<path d="M {x:g},{y-5:g} h 22" stroke="{escape(item["color"])}" stroke-width="3"{dash}/>')
        svg.append(svg_text(x + 30, y, item["label"], 14))
    caption = ("Each panel below is a rank-3 flat; three-point planes are implicit." if planes else
               "Curves mark rank-2 flats. Two-point lines are implicit. Unmarked crossings are not elements.")
    if layout.get("guides"):
        caption = "Solid curves mark rank-2 flats; dashed guides are annotations, not matroid lines."
    caption = layout.get("caption", caption)
    svg.append(svg_text(30, main_height - 25, caption, 13))
    # Rank four: show every plane with >=4 geometric points in its own panel.
    # Restricting each panel to that flat makes membership unambiguous.
    for i, plane in enumerate(planes):
        left, top = (i % 2) * width / 2, main_height + (i // 2) * panel_height
        members = " ".join(map(str, sorted(plane)))
        svg.append(f'<g data-plane="{members}"><rect x="{left+10:g}" y="{top:g}" '
                   f'width="{width/2-20:g}" height="{panel_height-15}" rx="14" '
                   'fill="#eff6ff" stroke="#bfdbfe"/>')
        svg.append(svg_text(left + 25, top + 27, "Plane {" + ", ".join(map(str, plane)) + "}", 14))
        coords = circle_positions(plane, (left + width / 4, top + panel_height / 2 + 10),
                                  min(width / 4 - 65, panel_height / 2 - 70))
        panel(plane, coords, [line for line in lines if set(line) <= set(plane)])
        svg.append('</g>')
    svg.extend(['</g>', '</svg>'])
    return "\n".join(svg) + "\n"


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matroid", type=Path, help="basis JSON or projection.csv")
    parser.add_argument("--layout", type=Path, help="optional positions and curve controls")
    parser.add_argument("--output", type=Path, default=Path("matroid.svg"))
    parser.add_argument("--write-json", type=Path, help="also export the computed complete basis family")
    parser.add_argument("--title", help="name used in the drawing and exported JSON")
    parser.add_argument("--png", type=Path, help="also render a transparent PNG (requires CairoSVG)")
    parser.add_argument("--scale", type=float, default=3, help="PNG pixels per SVG unit (default: 3)")
    args = parser.parse_args()
    matroid = load_matroid(args.matroid)
    if args.title:
        matroid.name = args.title
    layout = json.loads(args.layout.read_text(encoding="utf-8")) if args.layout else None
    drawing = draw_svg(matroid, layout)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(drawing, encoding="utf-8")
    if args.write_json:
        args.write_json.parent.mkdir(parents=True, exist_ok=True)
        # One basis per line keeps the complete list easy to audit.
        data = matroid.as_dict()
        basis_text = ",\n".join("    " + json.dumps(b) for b in data["bases"])
        encoded = ('{\n  "name": ' + json.dumps(data["name"], ensure_ascii=False)
                   + ',\n  "ground_set": ' + json.dumps(data["ground_set"])
                   + ',\n  "bases": [\n' + basis_text + '\n  ]\n}\n')
        args.write_json.write_text(encoded, encoding="utf-8")
    if args.png:
        if not math.isfinite(args.scale) or args.scale <= 0:
            parser.error("--scale must be positive and finite")
        try:
            import cairosvg
        except ImportError:
            parser.error("PNG export needs CairoSVG: python -m pip install '.[drawing]'")
        args.png.parent.mkdir(parents=True, exist_ok=True)
        cairosvg.svg2png(bytestring=drawing.encode("utf-8"), write_to=str(args.png), scale=args.scale)
    print(f"Rank {matroid.r}; {len(matroid.ground)} elements; {len(matroid.bases)} bases. Wrote {args.output}.")


if __name__ == "__main__":
    cli()
