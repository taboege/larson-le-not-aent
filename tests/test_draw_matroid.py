"""Checks of the mathematical input and the incidences encoded by the drawing."""

import csv
from itertools import combinations
import json
from pathlib import Path
import struct
import tempfile
import unittest
import xml.etree.ElementTree as ET

from draw_matroid import Matroid, draw_svg, load_matroid


ROOT = Path(__file__).resolve().parents[1]
SVG = {"s": "http://www.w3.org/2000/svg"}


class DrawMatroidTests(unittest.TestCase):
    def test_basis_exchange_rejects_nonmatroid(self):
        with self.assertRaisesRegex(ValueError, "Basis exchange"):
            Matroid([1, 2, 3, 4], [[1, 2], [3, 4]])

    def test_loops_and_parallel_elements(self):
        matroid = Matroid([1, 2, 3, 4], [[1, 3], [2, 3]])
        root = ET.fromstring(draw_svg(matroid))
        self.assertEqual({p.attrib["data-point"] for p in root.findall(".//s:circle", SVG)}, {"1", "3"})
        text = " ".join(root.itertext())
        self.assertIn("1 = 2", text)
        self.assertIn("Loops: 4", text)
        self.assertNotIn("data-flat", ET.tostring(root, encoding="unicode"))

    def test_rank_zero_and_one(self):
        for bases, expected_points in (([[]], 0), ([[1], [2]], 1)):
            root = ET.fromstring(draw_svg(Matroid([1, 2], bases)))
            self.assertEqual(len(root.findall(".//s:circle", SVG)), expected_points)

    def test_rank_two_draws_one_four_point_line(self):
        root = ET.fromstring(draw_svg(Matroid([1, 2, 3, 4], list(combinations(range(1, 5), 2)))))
        self.assertEqual([p.attrib["data-flat"] for p in root.findall(".//s:path[@data-flat]", SVG)],
                         ["1 2 3 4"])

    def test_fano_has_seven_three_point_lines(self):
        lines = {frozenset(s) for s in [(1, 2, 3), (1, 4, 5), (1, 6, 7),
                 (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 5, 6)]}
        matroid = Matroid(range(1, 8), [s for s in combinations(range(1, 8), 3) if frozenset(s) not in lines])
        self.assertEqual(set(matroid.flats(2)), lines)

    def test_rank_four_plane_panel(self):
        # U(3,4) plus a coloop: the only nontrivial plane is {1,2,3,4}.
        matroid = Matroid(range(1, 6), [list(s) + [5] for s in combinations(range(1, 5), 3)])
        root = ET.fromstring(draw_svg(matroid))
        self.assertEqual([p.attrib["data-plane"] for p in root.findall(".//s:g[@data-plane]", SVG)],
                         ["1 2 3 4"])

    def test_projection_requires_a_complete_exact_rank_table(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "projection.csv"

            def write_rows(change=None):
                with path.open("w", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(["points", "lower", "upper"])
                    for mask in range(8):
                        points = [i + 1 for i in range(3) if mask >> i & 1]
                        lo = hi = min(2, len(points))
                        if mask == 3 and change:
                            lo, hi = change
                        writer.writerow([" ".join(map(str, points)), lo, hi])

            write_rows()
            self.assertEqual(len(load_matroid(path).bases), 3)
            write_rows(("3/2", 2))
            with self.assertRaisesRegex(ValueError, "exact, equal integer"):
                load_matroid(path)
            write_rows((0, 0))
            with self.assertRaisesRegex(ValueError, "not the rank function"):
                load_matroid(path)

    def test_final_diagram_records_all_lines_and_only_a_dashed_missing_axis(self):
        matroid = load_matroid(ROOT / "data/final_matroid.json")
        layout = json.loads((ROOT / "figures/final_layout.json").read_text())
        drawing = draw_svg(matroid, layout)
        root = ET.fromstring(drawing)
        lines = [frozenset(map(int, p.attrib["data-flat"].split()))
                 for p in root.findall(".//s:path[@data-flat]", SVG)]
        self.assertEqual(len(matroid.bases), 104)
        self.assertEqual(len(lines), 10)
        self.assertEqual(set(lines), {f for f in matroid.flats(2) if len(f) >= 3})
        self.assertEqual(sorted(map(len, lines)), [3] * 8 + [4] * 2)
        guides = root.findall(".//s:path[@data-guide]", SVG)
        self.assertEqual(len(guides), 1)
        axis = frozenset(map(int, guides[0].attrib["data-guide"].split()))
        self.assertEqual(axis, frozenset([5, 14]))
        self.assertEqual(matroid.rank(axis), 2)
        self.assertEqual(matroid.rank(axis | {1}), 3)
        # The two guide points share a horizontal axis; the third is visibly off it.
        positions = layout["positions"]
        self.assertEqual(positions["5"][1], positions["14"][1])
        self.assertGreater(positions["1"][1] - positions["5"][1], 50)
        self.assertIn("stroke-dasharray", guides[0].attrib)
        triangles = root.findall(".//s:polygon[@data-triangle]", SVG)
        self.assertEqual({frozenset(map(int, p.attrib["data-triangle"].split())) for p in triangles},
                         {frozenset([2, 15, 12]), frozenset([6, 16, 8])})
        self.assertFalse(root.findall(".//s:rect", SVG))  # no opaque background
        self.assertEqual(drawing, (ROOT / "figures/final_matroid.svg").read_text())
        bad_layout = {"lines": [{"flat": [1, 5, 14]}]}
        with self.assertRaisesRegex(ValueError, "rank-two flat"):
            draw_svg(matroid, bad_layout)

    def test_png_is_high_resolution_rgba(self):
        data = (ROOT / "figures/final_matroid.png").read_bytes()
        self.assertEqual(data[:8], b"\x89PNG\r\n\x1a\n")
        width, height, depth, color, _, _, _ = struct.unpack(">IIBBBBB", data[16:29])
        self.assertEqual((width, height, depth, color), (3540, 2940, 8, 6))


if __name__ == "__main__":
    unittest.main()
