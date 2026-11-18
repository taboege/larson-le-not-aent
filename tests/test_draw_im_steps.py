"""Check the source data, geometric exceptions, and displayed IM incidences."""

import ast
import copy
from itertools import combinations
import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET

from draw_matroid import load_matroid
from draw_im_steps import (SpatialDrawing, affine_rank, original_svg, interactive_html,
                           crossing_gaps, gapped_path, depth_order)


ROOT = Path(__file__).resolve().parents[1]
NS = {"s": "http://www.w3.org/2000/svg"}


class DrawIMStepsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.matroid = load_matroid(ROOT / "data/dual_matroid.json")
        cls.steps = json.loads((ROOT / "data/im_steps.json").read_text())
        cls.layout = json.loads((ROOT / "figures/im_spatial_layout.json").read_text())
        cls.drawing = SpatialDrawing(cls.matroid, cls.steps, cls.layout)

    def test_display_inputs_match_the_verifiers_inputs(self):
        # Read only the literal inputs, without running the verification.
        module = ast.parse((ROOT / "im_projection.py").read_text())
        inputs = {node.targets[0].id: ast.literal_eval(node.value) for node in module.body
                  if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)
                  and node.targets[0].id in ("NONBASES", "IM_LINES")}
        nonbases = {frozenset(s) for s in inputs["NONBASES"]}
        self.assertEqual(self.matroid.bases,
                         {frozenset(s) for s in combinations(range(1, 11), 4)} - nonbases)
        self.assertEqual([[list(L) for L in step] for step in inputs["IM_LINES"]],
                         [[line["generators"] for line in step["lines"]] for step in self.steps["steps"]])
        self.assertEqual(len(self.matroid.flats(2)), 45)
        self.assertEqual(len(self.drawing.planes), 12)

    def test_only_the_explicit_exception_fails_spatial_coplanarity(self):
        self.assertEqual(self.drawing.relaxed_nonbases, [[6, 7, 8, 10]])
        points = self.drawing.exact
        self.assertEqual(self.matroid.rank(frozenset([6, 7, 8, 10])), 3)
        self.assertEqual(affine_rank([points[e] for e in (6, 7, 8, 10)]), 4)
        nonstraight = []
        for step in self.steps["steps"]:
            for line in step["lines"]:
                if affine_rank([points[e] for e in line["members"]]) != 2:
                    nonstraight.append(line["members"])
        self.assertEqual(nonstraight, [[7, 10, 14]])

    def test_unannotated_coordinate_error_is_rejected(self):
        layout = copy.deepcopy(self.layout)
        layout["positions"]["4"][2] = "1/100"
        with self.assertRaisesRegex(ValueError, "Unannotated spatial incidence"):
            SpatialDrawing(self.matroid, self.steps, layout)

    def test_incorrect_rank_premise_is_rejected(self):
        steps = copy.deepcopy(self.steps)
        steps["steps"][0]["premises"][-1]["rank"] = 3
        with self.assertRaisesRegex(ValueError, "seven rank premises"):
            SpatialDrawing(self.matroid, steps, self.layout)

    def test_curved_sheet_must_cover_its_declared_plane(self):
        layout = copy.deepcopy(self.layout)
        layout["curved_planes"][0]["straight"] = [6, 9]
        with self.assertRaisesRegex(ValueError, "two disjoint spanning pairs"):
            SpatialDrawing(self.matroid, self.steps, layout)

    def test_combined_drawing_shows_every_im_triple_in_its_points_color(self):
        root = ET.fromstring(original_svg(self.drawing))
        points = {int(p.attrib['data-point']): p for p in root.findall('.//s:circle[@data-point]', NS)}
        self.assertEqual(set(points), set(range(1, 17)))
        self.assertEqual(len({points[e].attrib['stroke'] for e in range(11, 17)}), 6)
        lines = root.findall('.//s:g[@data-line]', NS)
        self.assertEqual(len(lines), 17)  # 18 uses; the pair 5,6 occurs twice.
        for step in self.steps['steps']:
            z = step['new_point']
            selected = [line for line in lines if str(z) in line.attrib['data-im-points'].split()]
            self.assertEqual({tuple(map(int, line.attrib['data-generators'].split())) for line in selected},
                             {tuple(line['generators']) for line in step['lines']})
            for line in selected:
                self.assertIn(str(z), line.attrib['data-line'].split())
                self.assertIn(points[z].attrib['stroke'],
                              [p.attrib['stroke'] for p in line.findall('s:path', NS)])
        shared = next(line for line in lines if line.attrib['data-generators'] == '5 6')
        self.assertEqual(shared.attrib['data-line'], '5 6 12 15')
        self.assertEqual(shared.attrib['data-im-points'], '12 15')

    def test_all_original_circuits_are_accounted_for(self):
        root = ET.fromstring(original_svg(self.drawing))
        shown = {frozenset(map(int, p.attrib['data-plane'].split()))
                 for p in root.findall('.//s:path[@data-plane]', NS)}
        self.assertEqual(shown, {frozenset(s) for s in self.drawing.planes})
        nonbases = {frozenset(c) for c in combinations(self.matroid.ground, 4)} - self.matroid.bases
        self.assertEqual({frozenset(c) for f in shown for c in combinations(f, 4)}, nonbases)
        five = {frozenset(c) for c in combinations(self.matroid.ground, 5)
                if not any(n <= set(c) for n in nonbases)}
        self.assertEqual((len(nonbases), len(five)), (20, 140))
        stored = json.loads(root.find('.//s:metadata[@id="matroid-circuits"]', NS).text)
        self.assertEqual(set(map(frozenset, stored)), nonbases | five)

    def test_crossing_gaps_follow_spatial_depth_and_preserve_marked_meets(self):
        lines = [[(-100, 0, 0), (100, 0, 0)], [(0, -100, 1), (0, 100, 1)]]
        gaps = crossing_gaps(lines, [])
        self.assertEqual(list(map(len, gaps)), [1, 0])
        self.assertEqual(gapped_path(lines[0], gaps[0]).count('M '), 2)
        self.assertEqual(gapped_path(lines[1], gaps[1]).count('M '), 1)
        self.assertEqual(crossing_gaps(lines, [(0, 0, 0)]), [[], []])
        self.assertEqual(list(map(len, crossing_gaps(lines[::-1], []))), [0, 1])

    def test_intersecting_planes_change_depth_order_at_their_intersection(self):
        level = dict(name='level', curved=False, points=[[-1,-1,0], [1,-1,0], [1,1,0], [-1,1,0]])
        tilted = dict(name='tilted', curved=False, points=[[-1,-1,-1], [1,-1,1], [1,1,1], [-1,1,-1]])
        ordered = depth_order([level, tilted], [0, 0, 1])
        self.assertEqual([f['name'] for f in ordered], ['tilted', 'level', 'tilted'])
        self.assertTrue(all(p[0] <= 0 for p in ordered[0]['points']))
        self.assertTrue(all(p[0] >= 0 for p in ordered[-1]['points']))

    def test_checked_in_artifacts_are_reproducible(self):
        for name, content in [('rank4_original.svg', original_svg(self.drawing)),
                              ('rank4_interactive.html', interactive_html(self.drawing))]:
            self.assertEqual((ROOT / 'figures' / name).read_text(), content)


if __name__ == '__main__':
    unittest.main()
