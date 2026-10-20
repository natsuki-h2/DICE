"""Tests for source-population generation and matching."""

import unittest
from collections import namedtuple

from bliss_alma.mockdata import uniform_grid_positions
from bliss_alma.recovery import line_completeness, match_catalogue_to_positions


MockSource = namedtuple("MockSource", ["peak_index"])


class PopulationTests(unittest.TestCase):
    def test_uniform_grid_positions(self):
        positions = uniform_grid_positions(
            center=(50, 50), spacing_pixels=10, radius_pixels=15
        )
        self.assertEqual(len(positions), 9)
        self.assertIn((50, 50), positions)
        self.assertNotIn((30, 30), positions)

    def test_unique_greedy_matching(self):
        inputs = [(10, 10), (20, 20)]
        catalogue = [MockSource((11, 10)), MockSource((19, 20))]
        matches = match_catalogue_to_positions(inputs, catalogue, 2.0)
        self.assertEqual(len(matches), 2)
        self.assertEqual(matches[0].detection_index, 0)
        self.assertEqual(matches[1].detection_index, 1)

    def test_line_completeness_matches_in_space_frequency_and_spw(self):
        injected = [
            {"spw": 0, "x": 10, "y": 10, "channel": 50},
            {"spw": 0, "x": 20, "y": 20, "channel": 80},
            {"spw": 1, "x": 10, "y": 10, "channel": 50},
        ]
        detections = [
            {"spw": 0, "x": 10.5, "y": 10, "channel": 51},
            {"spw": 0, "x": 20, "y": 20, "channel": 95},
            {"spw": 0, "x": 10, "y": 10, "channel": 50.5},
        ]
        result = line_completeness(
            injected, detections,
            max_spatial_distance_pixels=2.0,
            max_channel_distance=3.0,
        )
        self.assertEqual(result.n_injected, 3)
        self.assertEqual(result.n_recovered, 1)
        self.assertAlmostEqual(result.completeness, 1 / 3)
        self.assertEqual(result.recovered, (True, False, False))

    def test_line_completeness_does_not_reuse_one_detection(self):
        injected = [
            {"x": 10, "y": 10, "channel": 50},
            {"x": 11, "y": 10, "channel": 50},
        ]
        detections = [{"x": 10.5, "y": 10, "channel": 50}]
        result = line_completeness(
            injected, detections,
            max_spatial_distance_pixels=2.0,
            max_channel_distance=1.0,
        )
        self.assertEqual(result.n_recovered, 1)


if __name__ == "__main__":
    unittest.main()
