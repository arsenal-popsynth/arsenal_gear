"""
imf
==========

This submodule contains all the code required to sample from an IMF.
"""

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from pyerf import erf as _scalar_erf
from pyerf import erfinv as _scalar_erfinv
from scipy.stats import rv_continuous

erf = np.vectorize(_scalar_erf)
erfinv = np.vectorize(_scalar_erfinv)

__all__ = ["IMF", "Salpeter", "MillerScalo"]


class IMF(rv_continuous):
    """
    This class is the superclass of all Initial Mass Functions (IMFs)

    :param min_mass: Least massive star in the IMF
    :type min_mass: astropy mass unit
    :param max_mass: Most massive star in the IMF
    :type max_mass: astropy mass unit
    :param name: Name for the scipy.stats rv_continuous instance
    :type name: str
    :param seed: Random seed for sampling
    :type seed: None, int, numpy.random.Generator, or numpy.random.RandomState
    """

    def __init__(
        self, min_mass: Quantity["mass"], max_mass: Quantity["mass"], name="", seed=None
    ):
        self.min_mass: float = min_mass.to(u.Msun).value
        self.max_mass: float = max_mass.to(u.Msun).value
        if (self.min_mass <= 0) or (self.max_mass <= 0):
            raise ValueError(
                f"(min_mass,max_mass) ({min_mass}, {max_mass}) must both be positive. Stars cannot have zero or negative mass."
            )
        if self.min_mass >= self.max_mass:
            raise ValueError(
                f"min_mass ({min_mass}) must be less than max_mass ({max_mass})."
            )
        super().__init__(a=self.min_mass, b=self.max_mass, name=name, seed=seed)

    def sample_mass(self, mtot: Quantity["mass"]) -> Quantity["mass"]:
        """
        Draw a sample from the IMF with target total mass

        :param mtot: Targer total mass of the sample
        :type mtot: Quantity["mass"]
        :return: List of masses of stars
        :rtype: Quantity["mass"]
        """
        N_samp = round((mtot / self.mean()).value)
        return self.sample(N_samp)

    def sample(self, N: int) -> Quantity["mass"]:
        """
        Draw a sample from the IMF with a specific number of stars

        :param N: Number of stars to draw
        :type N: int
        :return: List of masses of stars
        :rtype: Quantity["mass"]
        """
        return self.rvs(size=N) * u.Msun


class Salpeter(IMF):
    """
    A simple, classic Salpeter 1955 (slope 2.35) IMF.

    :param min_mass: Least massive star in the IMF
    :type min_mass: astropy mass unit
    :param max_mass: Most massive star in the IMF
    :type max_mass: astropy mass unit
    :param alpha: The IMF slope (note that the slope applied is -alpha, so this
                  should be positive)
    :type alpha: float
    :param seed: Random seed for sampling
    :type seed: None, int, numpy.random.Generator, or numpy.random.RandomState
    """

    def __init__(
        self,
        min_mass: Quantity["mass"] = 0.08 * u.Msun,
        max_mass: Quantity["mass"] = 100 * u.Msun,
        alpha: float = 2.35,
        seed=None,
    ):
        self.alpha = alpha
        self.name = "Salpeter"
        assert alpha >= 0
        super().__init__(min_mass, max_mass, self.name, seed=seed)

    def _pdf(self, x: np.float64, *args) -> np.float64:
        upper = np.power(self.max_mass, 1 - self.alpha)
        lower = np.power(self.min_mass, 1 - self.alpha)
        norm = (upper - lower) / (1 - self.alpha)
        return np.power(x, -self.alpha) / norm

    def _ppf(self, q: np.float64, *args) -> np.float64:
        upper = np.power(self.max_mass, 1 - self.alpha)
        lower = np.power(self.min_mass, 1 - self.alpha)
        return (q * (upper - lower) + lower) ** (1.0 / (1 - self.alpha))


class PiecewisePowerLaw(IMF):
    def __init__(self, min_mass, max_mass, masses, betas, name="", seed=None):
        """
        Generic class for N-part piecewise power law IMFs.
        :param masses: List/Array of transition masses (e.g., [0.5, 1.0])
        :param betas: List/Array of slopes (e.g., [-1.3, -2.3, -2.7])
        """
        super().__init__(min_mass, max_mass, name=name, seed=seed)
        # define the slopes
        self.betas = np.array(betas)

        # create the boundaries
        m_pts = [m.to(u.Msun).value if hasattr(m, "unit") else m for m in masses]
        m_pts = list(np.sort(m_pts))
        if len(m_pts) != len(self.betas) - 1:
            raise ValueError(
                f"Number of transition masses (masses={masses}) must be one less than the number of slopes (betas={betas})."
            )
        if min_mass.to(u.Msun).value > m_pts[0] or max_mass.to(u.Msun).value < m_pts[-1]:
            raise ValueError(
                f"User-specified mass range ({min_mass}, {max_mass}) must encompass all transition masses ({masses})."
            )
        m_pts = np.array([min_mass.to(u.Msun).value] + m_pts + [max_mass.to(u.Msun).value])

        # calculate stellar continuity constants
        self.alphas = np.ones(len(self.betas))
        for i in range(1, len(self.betas)):
            self.alphas[i] = self.alphas[i - 1] * (
                m_pts[i] ** (self.betas[i - 1] - self.betas[i])
            )

        self.weights = []
        # calculate weights (or areas); weight is 0 if segment falls outside user's mass range
        # determines self.m_pts which creates the true user-defined boundaries
        for i in range(len(self.betas)):
            b_plus_1 = self.betas[i] + 1
            w = (self.alphas[i] / b_plus_1) * (
                m_pts[i + 1]**b_plus_1 - m_pts[i]**b_plus_1
            )
            self.weights.append(w)

        self.m_pts = m_pts
        self.weights = np.array(self.weights)
        self.total_area = np.sum(self.weights)
        self.cum_weights = np.cumsum(self.weights) / self.total_area

    def _pdf(self, x, *args):
        # check which segment x falls into based on m_pts
        conditions = [
            (x >= self.m_pts[i]) & (x < self.m_pts[i + 1])
            for i in range(len(self.betas))
        ]
        conditions[-1] = (x >= self.m_pts[-2]) & (x <= self.m_pts[-1])

        results = [
            self.alphas[i] * np.power(x, self.betas[i]) for i in range(len(self.betas))
        ]
        return np.select(conditions, results, default=0.0) / self.total_area

    def _ppf(self, q, *args):
        """Generalized inverse CDF (sampling) logic."""
        is_scalar = np.isscalar(q)
        q = np.atleast_1d(q)

        conditions = []
        prev_q = 0.0
        for i in range(len(self.cum_weights)):
            conditions.append((q >= prev_q) & (q <= self.cum_weights[i]))
            prev_q = self.cum_weights[i]

        results = []
        for i in range(len(self.betas)):
            b_plus_1 = self.betas[i] + 1
            w_prior = self.weights[:i].sum() if i > 0 else 0.0
            res = (
                (q * self.total_area - w_prior) * b_plus_1 / self.alphas[i]
                + self.m_pts[i] ** b_plus_1
            ) ** (1.0 / b_plus_1)
            results.append(res)

        final = np.select(conditions, results, default=np.nan)
        return final[0] if is_scalar else final


class Kroupa(PiecewisePowerLaw):
    """
    Kroupa (2001) / Kroupa (1993) IMF implementation

    :param min_mass: Least massive star in the IMF
    :type min_mass: astropy mass unit
    :param max_mass: Most massive star in the IMF
    :type max_mass: astropy mass unit
    :param seed: Random seed for sampling
    :type seed: None, int, numpy.random.Generator, or numpy.random.RandomState

    """

    def __init__(
        self,
        min_mass: Quantity["mass"] = 0.08 * u.Msun,
        max_mass: Quantity["mass"] = 100.0 * u.Msun,
        seed=None,
    ):
        super().__init__(
            min_mass=min_mass,
            max_mass=max_mass,
            masses=[0.5, 1.0],
            betas=[-1.3, -2.3, -2.7],
            name="Kroupa",
            seed=seed,
        )


class MillerScalo(PiecewisePowerLaw):
    """
    Miller & Scalo (1979) IMF implementation.

    :param min_mass: Least massive star in the IMF
    :type min_mass: astropy mass unit
    :param max_mass: Most massive star in the IMF
    :type max_mass: astropy mass unit
    :param seed: Random seed for sampling
    :type seed: None, int, numpy.random.Generator, or numpy.random.RandomState

    """

    def __init__(
        self,
        min_mass: Quantity["mass"] = 0.1 * u.Msun,
        max_mass: Quantity["mass"] = 100.0 * u.Msun,
        seed=None,
    ):
        super().__init__(
            min_mass=min_mass,
            max_mass=max_mass,
            masses=[1.0, 10.0],
            betas=[-0.4, -1.5, -2.3],
            name="MillerScalo",
            seed=seed,
        )


class Chabrier(IMF):
    def __init__(
        self,
        min_mass: u.Quantity = 0.08 * u.Msun,
        max_mass: u.Quantity = 100.0 * u.Msun,
        beta: float = -2.3,
        seed=None,
    ):
        self.name = "Chabrier"
        self.A = 0.158
        self.sigma = 0.69
        self.mc = 0.079
        val_at_1 = (self.A / np.log(10)) * np.exp(
            -((np.log10(1.0) - np.log10(self.mc)) ** 2) / (2 * self.sigma**2)
        )
        self.A_high = val_at_1
        self.beta = beta

        self.min_mass = min_mass.to(u.Msun).value
        self.max_mass = max_mass.to(u.Msun).value
        self.m_transition = 1.0

        self.E_min = erf(
            (np.log10(self.min_mass) - np.log10(self.mc)) / (np.sqrt(2) * self.sigma)
        )
        self.E_trans = erf(
            (np.log10(self.m_transition) - np.log10(self.mc))
            / (np.sqrt(2) * self.sigma)
        )

        self.area_low = self._integrate_low()
        self.area_high = self._integrate_high()
        self.total_area = self.area_low + self.area_high

        super().__init__(min_mass, max_mass, self.name, seed=seed)

    def _integrate_low(self):
        if self.min_mass >= self.m_transition:
            return 0.0
        actual_min = max(self.min_mass, 0.01)  # 0.01 as a safety floor
        E_min_actual = erf(
            (np.log10(actual_min) - np.log10(self.mc)) / (np.sqrt(2) * self.sigma)
        )
        factor = self.A * np.sqrt(np.pi / 2) * self.sigma
        return factor * (self.E_trans - E_min_actual)

    def _integrate_high(self):
        if self.max_mass <= self.m_transition:
            return 0.0
        actual_start = max(self.min_mass, self.m_transition)
        actual_end = self.max_mass
        p = self.beta + 1
        return (self.A_high / p) * (actual_end**p - actual_start**p)

    def _pdf(self, x: np.ndarray, *args) -> np.ndarray:
        is_scalar = np.isscalar(x)
        x = np.atleast_1d(x)
        condlist = [
            (x >= self.min_mass) & (x < self.m_transition),
            (x >= self.m_transition) & (x <= self.max_mass),
        ]

        res_low = (self.A / (x * np.log(10))) * np.exp(
            -((np.log10(x) - np.log10(self.mc)) ** 2) / (2 * self.sigma**2)
        )
        res_high = self.A_high * np.power(x, self.beta)

        result = np.select(condlist, [res_low, res_high], default=0.0) / self.total_area
        return result[0] if is_scalar else result

    def _ppf(self, q: np.ndarray, *args) -> np.ndarray:
        """inverse CDF for sampling using np.select."""
        q = np.atleast_1d(q)
        q_break = self.area_low / self.total_area

        # Segment 1: Low mass log-normal
        q_scaled_low = q / q_break
        target_erf = q_scaled_low * (self.E_trans - self.E_min) + self.E_min
        res_low = 10 ** (
            np.sqrt(2) * self.sigma * erfinv(target_erf) + np.log10(self.mc)
        )

        # Segment 2: High mass power law
        q_scaled_high = (q - q_break) / (1.0 - q_break)
        p = self.beta + 1
        res_high = (q_scaled_high * (self.max_mass**p - 1.0**p) + 1.0**p) ** (1.0 / p)

        condlist = [q < q_break, (q >= q_break) & (q <= 1.0)]
        return np.select(condlist, [res_low, res_high], default=np.nan)
