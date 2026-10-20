"""Tests for configurable SciPy source extraction."""

import unittest

import numpy as np

from bliss_alma.extraction import estimate_rms, extract_sources


class ExtractionTests(unittest.TestCase):
    def test_positive_seed_and_flood_detection(self):
        data = np.zeros((9, 9), dtype=float)
        data[4, 4] = 6.0
        data[4, 3] = 4.0
        data[4, 5] = 4.0
        data[3, 4] = 4.0
        data[5, 4] = 4.0

        result = extract_sources(
            data,
            rms=1.0,
            seed_sigma=5.0,
            flood_sigma=3.0,
            min_pixels=5,
        )

        self.assertEqual(len(result.catalogue), 1)
        self.assertEqual(result.catalogue[0].peak_index, (4, 4))
        self.assertEqual(result.catalogue[0].peak_value, 6.0)
        self.assertEqual(result.catalogue[0].n_pixels, 5)

    def test_negative_detection_is_symmetric(self):
        data = np.zeros((9, 9), dtype=float)
        data[4, 4] = -6.0
        data[4, 3] = -4.0
        data[4, 5] = -4.0
        data[3, 4] = -4.0
        data[5, 4] = -4.0

        result = extract_sources(
            data,
            rms=1.0,
            polarity="negative",
            min_pixels=5,
        )

        self.assertEqual(len(result.catalogue), 1)
        self.assertEqual(result.catalogue[0].peak_index, (4, 4))
        self.assertEqual(result.catalogue[0].peak_value, 6.0)

    def test_rms_methods_are_selectable(self):
        data = np.array([-1.0, 0.0, 1.0, 100.0])
        mad = estimate_rms(data, method="mad")
        std = estimate_rms(data, method="std")
        rms = estimate_rms(data, method="rms")
        self.assertLess(mad, std)
        self.assertLess(std, rms)

    def test_absolute_thresholds_override_sigma(self):
        data = np.zeros((5, 5), dtype=float)
        data[2, 2] = 2.0
        result = extract_sources(
            data,
            rms=100.0,
            seed_threshold=1.5,
            flood_threshold=1.0,
            min_pixels=1,
        )
        self.assertEqual(len(result.catalogue), 1)

    def test_detection_mask_and_border_rejection(self):
        data = np.zeros((7, 7), dtype=float)
        data[0, 0:2] = 6.0
        data[3, 3] = 6.0

        masked = np.zeros(data.shape, dtype=bool)
        masked[3, 3] = True
        result = extract_sources(
            data,
            rms=1.0,
            flood_sigma=5.0,
            min_pixels=1,
            detection_mask=masked,
            exclude_border=True,
        )
        self.assertEqual(len(result.catalogue), 0)


if __name__ == "__main__":
    unittest.main()
