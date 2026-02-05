"""
test_element_yields
==========

This file contains unit tests for the yield interpolation methods
which are contained in element_yields
"""

import numpy as np
from numpy.testing import assert_approx_equal

import arsenal_gear


def assert_array_approx_equal(actual_arr, desired_arr, significant=7):
    for actual, desired in zip(actual_arr, desired_arr):
        assert_approx_equal(actual, desired, significant=significant)


class TestLimongiChieffi2018:
    """Tests for yields from Limongi & Chieffi (2018)"""

    lc18 = arsenal_gear.element_yields.LimongiChieffi2018()

    def test_monotinically_increasing_params(self):
        """Ensure monotonically increasing parameter grid"""
        for param in self.lc18.wind.params:
            assert np.all(
                param[1:] >= param[:-1]
            ), f"Parameter {param} for wind in {self.lc18.name} is not stricktly increasing."
        for param in self.lc18.ccsn.params:
            assert np.all(
                param[1:] >= param[:-1]
            ), f"Parameter {param} for ccsn in {self.lc18.name} is not stricktly increasing."

    def test_tables(self):
        """Make sure that online tables did not change by checking random yields."""

        elements = ["H", "C", "O", "Fe"]
        models = [
            (150, 3.236e-5, 13.0),  # 013d150
            (0, 1.345e-2, 15.0),  # 015a000
            (300, 3.236e-5, 20.0),  # 020d300
            (150, 3.236e-3, 25.0),
        ]  # 025b150

        ylds_tot = [
            [4.8816e0, 2.7393e-1, 1.4593e0, 7.2746e-2],
            [6.8845e0, 2.2174e-1, 5.8573e-1, 9.4084e-2],
            [6.7269e0, 6.4658e-1, 2.9605e0, 7.6375e-2],
            [7.8504e0, 7.0529e-1, 4.9064e0, 7.6535e-2],
        ]

        for model, yld_tot in zip(models, ylds_tot):
            ccsn_yld = self.lc18.ccsn.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            wind_yld = self.lc18.wind.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            loaded_yld = np.array([yc + yw for yc, yw in zip(ccsn_yld, wind_yld)])

            assert_array_approx_equal(loaded_yld.flatten(), yld_tot, significant=5)


class TestNuGrid:

    nugrid = arsenal_gear.element_yields.NuGrid()
    p16 = arsenal_gear.element_yields.Pignatari2016()
    r18 = arsenal_gear.element_yields.Ritter2018()
    bat = arsenal_gear.element_yields.Battino20192021()

    def test_monotinically_increasing_params(self):
        """Ensure monotonically increasing parameter grid"""
        for yld_set in [self.p16, self.r18, self.nugrid]:
            for param in yld_set.agb.params:
                assert np.all(
                    param[1:] >= param[:-1]
                ), f"Parameter {param} for agb in {yld_set.name} is not stricktly increasing."
            for param in yld_set.wind.params:
                assert np.all(
                    param[1:] >= param[:-1]
                ), f"Parameter {param} for wind in {yld_set.name} is not stricktly increasing."
            for param in yld_set.ccsn.params:
                assert np.all(
                    param[1:] >= param[:-1]
                ), f"Parameter {param} for ccsn in {yld_set.name} is not stricktly increasing."

        for param in self.bat.agb.params:
            assert np.all(
                param[1:] >= param[:-1]
            ), f"Parameter {param} for agb in {self.bat.name} is not stricktly increasing."

    def test_p16_tables(self):
        """Make sure that online tables did not change by checking random yields."""
        elements = ["H", "C", "O", "Fe"]
        models_lo = [(0.01, 2.0), (0.02, 3.0)]
        models_hi = [(0.01, 15.0), (0.02, 25.0)]

        ylds_lo = [
            [9.551e-1, 1.853e-2, 1.389e-2, 9.961e-4],
            [1.542e0, 3.955e-2, 3.603e-2, 3.377e-3],
        ]
        ylds_hi = [
            [1.920e0, 3.995e-3, 1.256e-2, 1.906e-3],
            [7.194e0, 2.787e-2, 8.769e-2, 1.619e-2],
        ]
        ylds_sn = [
            [4.569e0, 1.546e-1, 2.069e-1, 1.980e-1],
            [1.963e0, 4.365e-1, 7.289e-1, 9.694e-3],
        ]

        for model, yld in zip(models_lo, ylds_lo):
            agb_yld = self.p16.agb.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            assert_array_approx_equal(agb_yld.flatten(), yld, significant=4)

        for model, yld in zip(models_hi, ylds_hi):
            wind_yld = self.p16.wind.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            assert_array_approx_equal(wind_yld.flatten(), yld, significant=4)

        for model, yld in zip(models_hi, ylds_sn):
            wind_yld = self.p16.ccsn.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            assert_array_approx_equal(wind_yld.flatten(), yld, significant=4)

    def test_r18_tables(self):
        """Make sure that online tables did not change by checking random yields."""
        elements = ["H", "C", "O", "Fe"]
        models_lo = [(0.01, 2.0), (0.02, 3.0)]
        models_hi = [(0.01, 15.0), (0.02, 25.0)]

        ylds_lo = [
            [9.551e-1, 1.852e-2, 1.389e-2, 9.961e-4],
            [1.543e0, 3.956e-2, 3.604e-2, 3.379e-3],
        ]
        ylds_hi = [
            [6.590e0, 1.600e-1, 8.904e-1, 7.633e-2],
            [9.483e0, 2.324e-1, 9.370e-1, 2.499e-2],
        ]
        ylds_wi = [
            [1.169e0, 2.508e-3, 7.522e-3, 1.177e-3],
            [7.642e0, 2.894e-2, 9.239e-2, 1.752e-2],
        ]

        for model, yld in zip(models_lo, ylds_lo):
            agb_yld = self.r18.agb.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            assert_array_approx_equal(agb_yld.flatten(), yld, significant=3)

        for model, yld in zip(models_hi, ylds_wi):
            wind_yld = self.r18.wind.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            assert_array_approx_equal(wind_yld.flatten(), yld, significant=3)

        for model, yld_tot in zip(models_hi, ylds_hi):
            ccsn_yld = self.r18.ccsn.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            wind_yld = self.r18.wind.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            loaded_yld = np.array(
                [yc + yw for yc, yw in zip(ccsn_yld.flatten(), wind_yld.flatten())]
            )
            assert_array_approx_equal(loaded_yld, yld_tot, significant=3)

    def test_bat_tables(self):
        """Make sure that online tables did not change by checking random yields."""
        elements = ["H", "C", "O"]
        models = [(0.01, 2.0), (0.001, 3.0)]

        ylds_lo = [
            [7.334349323233e-1, 9.950594153e-3, 8.95887761e-3],
            [1.2695745259991285, 3.9231768329611e-3, 1.630396175869021e-3],
        ]

        for model, yld in zip(models, ylds_lo):
            agb_yld = self.bat.agb.get_yld(
                elements=elements,
                params=model,
                interpolate="nearest",
                extrapolate=False,
            )
            assert_array_approx_equal(agb_yld.flatten(), yld, significant=5)
